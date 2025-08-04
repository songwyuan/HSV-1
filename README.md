# HSV-1

## 目录准备
```
#创建主病毒目录结构（如果尚不存在）

mkdir -p /mnt/alamo01/users/yuansongwei7/download_dna/HSV-1/U87MG

#为两个GSE项目创建标准化目录

for gse in GSE235568; do
    mkdir -p /mnt/alamo01/users/yuansongwei7/download_dna/HSV-1//U87MG/${gse}/{scripts,raw_sra,fastq_files,fastqc_results,cleaned_data,alignment_results,analysis_results,logs}
done
```

##  download_GSE237079(GSE236646)
```
#!/bin/bash
set -euo pipefail
#设置环境变量
export PATH="/mnt/alamo01/users/chenyun730/micromamba/envs/sra-tools/bin:$PATH"
BASE_DIR="/mnt/alamo01/users/yuansongwei7/download_dna/HSV-1/GSE236646"
mkdir -p "${BASE_DIR}"/{raw_sra,fastq_files,fastqc_results,logs}

#样本列表（按实验条件分组）

MOI_0001=(
    SRR25181599 SRR25181600 SRR25181601  # Day3
    SRR25181588 SRR25181589 SRR25181590  # Day5
    SRR25181560 SRR25181561 SRR25181562  # Day7
)

MOI_001=(
    SRR25181585 SRR25181586 SRR25181587  # Day3
    SRR25181564 SRR25181594 SRR25181595  # Day5
    SRR25181557 SRR25181558 SRR25181559  # Day7
)

UNINFECTED=(
    SRR25181571 SRR25181580 SRR25181581  # Day3对照
    SRR25181568 SRR25181569 SRR25181570  # Day5对照
    SRR25181574 SRR25181575 SRR25181576  # Day7对照
)

ALL_SRR=("${MOI_0001[@]}" "${MOI_001[@]}" "${UNINFECTED[@]}")

#主处理函数
process_sample() {
    local srr="$1"
    local log_file="${BASE_DIR}/logs/${srr}.log"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START ${srr}" | tee -a "${log_file}"

    # 1. 下载SRA（带重试机制）
    if [[ ! -f "${BASE_DIR}/raw_sra/${srr}.sra" ]]; then
        for attempt in {1..3}; do
            echo "Download attempt ${attempt}" | tee -a "${log_file}"
            if wget -O "${BASE_DIR}/raw_sra/${srr}.sra" \
               "https://sra-pub-run-odp.s3.amazonaws.com/sra/${srr}/${srr}" \
               2>> "${log_file}"; then
                break
            else
                rm -f "${BASE_DIR}/raw_sra/${srr}.sra"
                sleep $((attempt * 10))
            fi
        done
    fi

    # 2. 转换FASTQ
    if [[ ! -f "${BASE_DIR}/fastq_files/${srr}_1.fastq.gz" ]]; then
        echo "Converting to FASTQ" | tee -a "${log_file}"
        fastq-dump --split-files "${BASE_DIR}/raw_sra/${srr}.sra" \
            --outdir "${BASE_DIR}/fastq_files" \
            --gzip \
            --skip-technical \
            --dumpbase \
            >> "${log_file}" 2>&1

        # 验证输出
        if [[ ! -f "${BASE_DIR}/fastq_files/${srr}_1.fastq.gz" ]]; then
            echo "ERROR: FASTQ conversion failed" | tee -a "${log_file}"
            exit 1
        fi
    fi

    # 3. 质控
    if [[ ! -f "${BASE_DIR}/fastqc_results/${srr}_1_fastqc.html" ]]; then
        echo "Running FastQC" | tee -a "${log_file}"
        fastqc -t 4 \
            "${BASE_DIR}/fastq_files/${srr}_1.fastq.gz" \
            "${BASE_DIR}/fastq_files/${srr}_2.fastq.gz" \
            -o "${BASE_DIR}/fastqc_results" \
            >> "${log_file}" 2>&1
    fi

    # 清理原始数据（保留注释，按需启用）
    # rm -f "${BASE_DIR}/raw_sra/${srr}.sra"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] FINISHED ${srr}" | tee -a "${log_file}"
}

#主程序（并行控制）
MAX_JOBS=8
for srr in "${ALL_SRR[@]}"; do
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 10
    done
    process_sample "$srr" &
done

wait  # 等待所有后台任务完成

#生成汇总报告
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generating MultiQC report"
multiqc "${BASE_DIR}/fastqc_results" -o "${BASE_DIR}/fastqc_results" \
    > "${BASE_DIR}/logs/multiqc.log" 2>&1

#最终校验
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Validation"
ls -lh "${BASE_DIR}"/fastq_files/*.gz | wc -l | \
    tee -a "${BASE_DIR}/logs/summary.log"
echo "Total samples processed: ${#ALL_SRR[@]}" \
    | tee -a "${BASE_DIR}/logs/summary.log"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ALL DONE"
```


## pinepline.sh(paired_single)
#paired
```
#!/bin/bash

########################################
# RNA-seq 自动流程（适用于多目录 PE 数据）
# 作者: ChatGPT + Yuansongwei7
########################################

# -------- 参数设置 --------
RAW_DIR="/mnt/alamo01/users/yuansongwei7/download_dna/HSV-1"  # 原始数据根目录（递归）
WORK_DIR="/mnt/alamo01/users/yuansongwei7/download_dna/HSV-1_pipeline_output"
CLEAN_DIR="$WORK_DIR/clean_reads"
ALIGN_DIR="$WORK_DIR/alignment"
QC_DIR="$WORK_DIR/qc"
combine_symbol_counts="/mnt/alamo01/projects/Group_Wang/script/function/combine_symbol_counts_STAR.R"

Species="Homo sapiens"
REF="/mnt/alamo01/users/yuansongwei7/genome_index/hg38/star_index"
GTF="/mnt/alamo01/users/yuansongwei7/genome_index/hg38/gencode.v43.annotation.gtf"

THREADS=64
featurecount_thread=32
FASTP_THREADS=64
use_R="/mnt/alamo01/users/chenyun730/micromamba/envs/R441/bin/Rscript"

# -------- 创建目录 --------
mkdir -p $WORK_DIR $QC_DIR $CLEAN_DIR $ALIGN_DIR

# -------- Step 1: 查找所有样本 --------
echo "🔍 正在查找所有双端 FASTQ 文件..."
find $RAW_DIR -type f -name "*_1.fastq.gz" | while read F1; do
    F2=${F1/_1.fastq.gz/_2.fastq.gz}
    if [ ! -f "$F2" ]; then
        echo "⚠️ 缺失配对文件: $F2，跳过！"
        continue
    fi

    SAMPLE=$(basename "$F1")
    SAMPLE=${SAMPLE%%_1.fastq.gz}

    echo "🎯 处理样本: $SAMPLE"

    # ---- Step 1: Fastp 清洗 ----
    CLEAN_R1="$CLEAN_DIR/${SAMPLE}_clean_R1.fq.gz"
    CLEAN_R2="$CLEAN_DIR/${SAMPLE}_clean_R2.fq.gz"
    DONE_FLAG="$CLEAN_DIR/${SAMPLE}.done"

    if [ -f "$DONE_FLAG" ]; then
        echo "⏩ [$SAMPLE] 已清洗，跳过"
    else
        echo "🧼 正在清洗 [$SAMPLE] ..."
        fastp -i "$F1" -I "$F2" -o "$CLEAN_R1" -O "$CLEAN_R2" \
            --detect_adapter_for_pe -q 25 -u 20 -e 20 -r -W 5 -M 30 \
            --length_required 50 -h "$CLEAN_DIR/${SAMPLE}.html" \
            -j "$CLEAN_DIR/${SAMPLE}.json" -w $FASTP_THREADS > "$CLEAN_DIR/${SAMPLE}_fastp.log" 2>&1

        if [ $? -eq 0 ]; then
            touch "$DONE_FLAG"
            echo "✅ [$SAMPLE] fastp 完成"
        else
            echo "❌ [$SAMPLE] fastp 失败" | tee "$CLEAN_DIR/${SAMPLE}.error"
            continue
        fi
    fi

    # ---- Step 2: 比对 ----
    ALIGN_SAMPLE_DIR="$ALIGN_DIR/$SAMPLE"
    mkdir -p "$ALIGN_SAMPLE_DIR"
    cd "$ALIGN_SAMPLE_DIR"

    ALIGN_DONE="$ALIGN_SAMPLE_DIR/${SAMPLE}.done"
    if [ -f "$ALIGN_DONE" ]; then
        echo "⏩ [$SAMPLE] 已比对，跳过"
        continue
    fi

    echo "🧬 [$SAMPLE] STAR 比对中..."
    STAR --runThreadN $THREADS \
        --genomeDir $REF \
        --readFilesIn "$CLEAN_R1" "$CLEAN_R2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${SAMPLE}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM NM MD \
        --outFilterMultimapNmax 10 --outFilterMismatchNmax 5 \
        --outFilterScoreMinOverLread 0.8 \
        --alignIntronMin 20 --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
        --sjdbScore 1 --twopassMode Basic \
        --twopass1readsN -1 --quantMode GeneCounts > star.log 2>&1

    if [ $? -eq 0 ]; then
        echo "✅ [$SAMPLE] STAR 比对完成"
        # ---- Step 3: featureCounts 定量 ----
        featureCounts -T $featurecount_thread -a $GTF \
            -o "${SAMPLE}_gene_counts.txt" "${SAMPLE}_Aligned.sortedByCoord.out.bam" > featurecounts.log 2>&1
        echo "✅ [$SAMPLE] featureCounts 完成"
        touch "$ALIGN_DONE"
    else
        echo "❌ [$SAMPLE] STAR 比对失败" | tee "${SAMPLE}.error"
        continue
    fi

    cd - > /dev/null
done

# -------- Step 4: 整合表达矩阵 --------
echo "📊 汇总所有基因表达矩阵..."
cd "$ALIGN_DIR"
$use_R $combine_symbol_counts "$ALIGN_DIR" "$Species" "$GTF"
echo "✅ 表达矩阵整合完成"

echo "🎉 全部分析流程结束！"

(rnaseq) yuansongwei7@mgt01:/mnt/alamo01/users/yuansongwei7/download_dna/HSV-1/RAW264.7/GSE204893
$

```




