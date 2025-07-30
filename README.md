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


## 




