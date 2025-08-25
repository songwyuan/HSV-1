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
    | tee -a "${BASE_DIR}/logs/summary.log"入眠-

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ALL DONE"
```


## pinepline.sh(paired_all_paired.sh)
```
$ cat RNAseq_pipeline_PE.sh
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

    ALIGN_DONE="$ALIGN_SAMPLE_DIR/${SAMPLE}.align.done"
    FEATURE_DONE="$ALIGN_SAMPLE_DIR/${SAMPLE}.feature.done"

    # ---- Step 2: 比对 ----
    if [ ! -f "$ALIGN_DONE" ]; then
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
            touch "$ALIGN_DONE"
        else
            echo "❌ [$SAMPLE] STAR 比对失败" | tee "${SAMPLE}.error"
            continue
        fi
    else
        echo "⏩ [$SAMPLE] 已比对，跳过比对步骤"
    fi

    # ---- Step 3: featureCounts 定量 ----
    if [ ! -f "$FEATURE_DONE" ]; then
        echo "📏 [$SAMPLE] featureCounts 定量中..."
        featureCounts -p -B -C -T $featurecount_thread -a $GTF \
            -o "${SAMPLE}_gene_counts.txt" "${SAMPLE}_Aligned.sortedByCoord.out.bam" > featurecounts.log 2>&1
        if [ $? -eq 0 ]; then
            echo "✅ [$SAMPLE] featureCounts 完成"
            touch "$FEATURE_DONE"
        else
            echo "❌ [$SAMPLE] featureCounts 失败" | tee "${SAMPLE}.feature.error"
        fi
    else
        echo "⏩ [$SAMPLE] 已完成 featureCounts，跳过"
    fi


    cd - > /dev/null
done

# -------- Step 4: 整合表达矩阵 --------
echo "📊 汇总所有基因表达矩阵..."
cd "$ALIGN_DIR"
$use_R $combine_symbol_counts "$ALIGN_DIR" "$Species" "$GTF"
echo "✅ 表达矩阵整合完成"

echo "🎉 全部分析流程结束！"


```
#单端数据分析
```
#!/bin/bash

########################################
# 自动化 RNA-seq 流程（单端版）
# FastQC → fastp → STAR + featureCounts → gene count matrix
# 作者: ChatGPT + Yuansongwei7（优化版）
########################################

# -------- 参数设置 --------
RAW_DIR="/mnt/alamo01/users/yuansongwei7/download_dna/HSV-1"
WORK_DIR="/mnt/alamo01/users/yuansongwei7/download_dna/HSV-1_pipeline_output_single"
CLEAN_DIR="$WORK_DIR/clean_reads"
ALIGN_DIR="$WORK_DIR/alignment"
QC_DIR="$WORK_DIR/qc"

Species="Homo sapiens"
REF="/mnt/alamo01/users/yuansongwei7/genome_index/hg38/star_index"
GTF="/mnt/alamo01/users/yuansongwei7/genome_index/hg38/gencode.v43.annotation.gtf"
combine_symbol_counts="/mnt/alamo01/users/yuansongwei7/scripts/combine_symbol_counts_STAR.R"
use_R="/mnt/alamo01/users/chenyun730/micromamba/envs/R441/bin/Rscript"

THREADS=64
FASTP_THREADS=64
featurecount_thread=32

mkdir -p $WORK_DIR $QC_DIR $CLEAN_DIR $ALIGN_DIR

########################################
# Step 1: FastQC + MultiQC
########################################
echo "📊 Step 1: Running FastQC..."
QC_REPORT="$QC_DIR/raw/multiqc_report.html"
if [ -f "$QC_REPORT" ]; then
    echo "✅ QA报告已存在，跳过FastQC与MultiQC：$QC_REPORT"
else
    mkdir -p $QC_DIR/raw
    fastqc -o $QC_DIR/raw -t $THREADS $RAW_DIR/*.fastq.gz
    multiqc -o $QC_DIR/raw $QC_DIR/raw
    echo "✅ FastQC与MultiQC完成！"
fi

########################################
# Step 2: fastp 清洗（单端）
########################################
echo "🧼 Step 2: Cleaning reads with fastp..."
for fq in $RAW_DIR/*.fastq.gz; do
    SAMPLE=$(basename "$fq" .fastq.gz)
    DONE_FLAG="$CLEAN_DIR/${SAMPLE}.done"

    # 跳过双端数据（*_1.fastq.gz / *_2.fastq.gz）
    if [[ "$fq" =~ _[12]\.fastq\.gz$ ]]; then
        echo "⚠️ 检测到双端文件 [$fq]，跳过"
        continue
    fi

    if [ -f "$DONE_FLAG" ]; then
        echo "⏩ [$SAMPLE] 已完成 fastp，跳过"
        continue
    fi

    echo "🧹 [$SAMPLE] fastp 清洗..."
    if fastp -i "$fq" -o "$CLEAN_DIR/${SAMPLE}_clean.fq.gz" \
        -q 25 -u 20 -e 20 -r -W 5 -M 30 \
        --length_required 50 \
        -h "$CLEAN_DIR/${SAMPLE}.html" \
        -j "$CLEAN_DIR/${SAMPLE}.json" \
        -w $FASTP_THREADS > "$CLEAN_DIR/${SAMPLE}_fastp.log" 2>&1; then
        touch "$DONE_FLAG"
        echo "✅ [$SAMPLE] fastp 完成"
    else
        echo "❌ [$SAMPLE] fastp 失败"
    fi
done

# fastp QC 合并
QC_REPORT2="$CLEAN_DIR/multiqc_report.html"
if [ ! -f "$QC_REPORT2" ]; then
    multiqc $CLEAN_DIR -o $CLEAN_DIR
fi

########################################
# Step 3: STAR + featureCounts
########################################
echo "🧬 Step 3: STAR + featureCounts..."
for fq in $CLEAN_DIR/*_clean.fq.gz; do
    SAMPLE=$(basename "$fq" _clean.fq.gz)
    ALIGN_SAMPLE_DIR=$ALIGN_DIR/$SAMPLE
    DONE_FLAG="$ALIGN_SAMPLE_DIR/${SAMPLE}.done"

    if [ -f "$DONE_FLAG" ]; then
        echo "⏩ [$SAMPLE] 已完成比对，跳过"
        continue
    fi

    mkdir -p $ALIGN_SAMPLE_DIR
    cd $ALIGN_SAMPLE_DIR

    echo "🧬 [$SAMPLE] STAR 比对..."
    if STAR --runThreadN $THREADS \
        --genomeDir $REF \
        --readFilesIn "$fq" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${SAMPLE}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM NM MD \
        --outFilterMultimapNmax 10 --outFilterMismatchNmax 5 \
        --outFilterScoreMinOverLread 0.8 \
        --alignIntronMin 20 --alignIntronMax 1000000 \
        --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
        --sjdbScore 1 --twopassMode Basic \
        --twopass1readsN -1 --quantMode GeneCounts > star.log 2>&1; then
        echo "✅ [$SAMPLE] STAR 完成"

        echo "🧮 [$SAMPLE] featureCounts..."
        featureCounts -T $featurecount_thread -a $GTF \
            -o "${SAMPLE}_gene_counts.txt" "${SAMPLE}_Aligned.sortedByCoord.out.bam" > featurecounts.log 2>&1
        echo "✅ [$SAMPLE] featureCounts 完成"

        touch "$DONE_FLAG"
    else
        echo "❌ [$SAMPLE] STAR 失败"
    fi
    cd - > /dev/null
done

########################################
# Step 4: 表达矩阵整合
########################################
echo "🟢 Step 4: 合并基因表达矩阵..."
cd $ALIGN_DIR
$use_R $combine_symbol_counts "$ALIGN_DIR" "$Species" "$GTF"
echo "✅ gene symbol 聚合完成！"

########################################
echo "🎉 所有样本流程完成！"

```

#combine_symbol_counts.R

```
#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(biomaRt))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(txdbmaker))

# 获取命令行参数（去掉 Rscript 路径本身）
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("❌ 请输入 alignment与物种信息 两个参数")
}

align_path <- args[1]
species = args[2]
gtf_path = args[3]

if(species == "Homo sapiens") {
  ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  symbol = "hgnc_symbol"
}else if(species == "Mus musculus"){
  ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  symbol = "mgi_symbol"
}

# 获取所有featureCounts输出文件的路径
files <- list.files(pattern = "gene_counts.txt$",path = align_path,recursive = T)

# 初始化一个空的数据框，用于存储基因计数矩阵
gene_counts <- data.frame()

# 遍历每个文件，读取数据并合并
for (file in files) {
  # 读取每个featureCounts的输出文件
  count_data <- read.table(file, header = TRUE, sep = "\t", comment.char = "#",check.names = F)
  
  # 提取Gene ID和计数列，确保计数列名称为"Count"
  count_data_subset <- count_data[, c(1,7)]
  
  # 将Gene ID列作为行名
  rownames(count_data_subset) <- count_data_subset$Geneid
  
  # 删除Geneid列，因为它已经作为行名
  count_data_subset$Geneid <- NULL
  
  # 获取样本名称作为列名
  sample_name <- strsplit(file,"\\/")[[1]][1]  # 从文件名中提取样本名称
  
  # 合并数据
  if (ncol(gene_counts) == 0) {
    # 如果是第一次合并，直接将数据放入gene_counts
    gene_counts <- count_data_subset
    colnames(gene_counts) <- sample_name
  } else {
    # 否则，将当前样本的计数数据添加到gene_counts矩阵中
    gene_counts <- cbind(gene_counts, count_data_subset)
    colnames(gene_counts)[ncol(gene_counts)] <- sample_name
  }
}

# 保存合并后的基因表达矩阵
write.csv(gene_counts, "transcript_count_matrix.csv", quote = FALSE, row.names = TRUE)

# 使用biomaRt查询Gene Symbol
gene_symbols <- getBM(attributes = c('ensembl_gene_id', symbol),
                      filters = 'ensembl_gene_id',
                      values = rownames(gene_counts),
                      mart = ensembl)
colnames(gene_symbols)[2] <- "symbol"

# 将Gene Symbol合并到基因计数矩阵中
gene_counts$GeneSymbol <- gene_symbols[,2][match(rownames(gene_counts), gene_symbols$ensembl_gene_id)]

# 聚合 & 清理
gene_counts <- gene_counts[-which(gene_counts$GeneSymbol == ""),]
gene_counts <- na.omit(gene_counts)
gene_matrix <- gene_counts %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), sum)) %>%
  as.data.frame()

rownames(gene_matrix) <- gene_matrix$GeneSymbol
gene_matrix <- gene_matrix[, -1]

# 输出
write.csv(gene_matrix, "gene_count_matrix_symbol_merged.csv", quote = FALSE)
cat("✅ Symbol Count聚合完成，输出文件：", paste(align_path,"gene_count_matrix_symbol_merged.csv",sep = "/"), "\n")

####

# 1. 建立TxDb对象，GTF文件路径根据你实际位置修改
txdb <- makeTxDbFromGFF(gtf_path, format="gtf")

# 2. 提取基因的外显子区
exonsByGene <- exonsBy(txdb, by="gene")

# 3. 计算每个基因外显子的宽度总和（合并重叠）
gene_lengths <- sum(width(reduce(exonsByGene)))

# 4. 转成数据框，gene_id 和 length
gene_length_df <- data.frame(
  ensembl = names(gene_lengths),
  length = as.numeric(gene_lengths)
)

gene_length_df <- left_join(gene_length_df, gene_symbols, by = c("ensembl" = "ensembl_gene_id")) %>%
  filter(symbol != "") %>%  # 去除没有 symbol 的
  group_by(symbol) %>%
  summarise(symbol_length = mean(length))

symbol_lengths <- gene_length_df$symbol_length
names(symbol_lengths) <- gene_length_df$symbol

use_gene <- intersect(gene_length_df$symbol,rownames(gene_matrix))

rpk <- sweep(gene_matrix[use_gene,], 1, symbol_lengths[use_gene] / 1000, "/")

tpm <- sweep(rpk, 2, colSums(rpk), "/") * 1e6

# 输出
write.csv(tpm, "gene_tpm_matrix_symbol_merged.csv", quote = FALSE)
cat("✅ Symbol TPM聚合完成，输出文件：", paste(align_path,"gene_tpm_matrix_symbol_merged.csv",sep = "/"), "\n")

```

# extractjson
```
#!/bin/env Rscript
# 修订点：
# 1) 路径改到你的 pipeline 输出：PAIRED 与 SINGLE 分别指向 HSV-1_pipeline_output / HSV-1_pipeline_output_single
# 2) 修复 unknown/na 判断；修复 TPM 列名匹配变量错误；增加列缺失告警
# 3) 允许 LibraryLayout 大小写差异（PAIRED / SINGLE / SINGLE-END）

Sys.setenv(TMPDIR = "/mnt/alamo01/projects/Group_Wang/sampledata/tmp_r")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("❌ Usage: Rscript generate_metadata_batch_from_table.R <input_csv_path>")
}

suppressPackageStartupMessages({
  library(jsonlite)
  library(readr)
  library(stringr)
  library(dplyr)
  library(fs)
})

input_csv <- args[1]
df <- read_csv(input_csv, show_col_types = FALSE)

# ---- 你的 pipeline 输出根目录（按 LibraryLayout 选择） ----
PAIRED_ALIGN_DIR <- "/mnt/alamo01/users/yuansongwei7/download_dna/HSV-1_pipeline_output/alignment"
SINGLE_ALIGN_DIR <- "/mnt/alamo01/users/yuansongwei7/download_dna/HSV-1_pipeline_output_single/alignment"

# 小工具：标准化 layout
norm_layout <- function(x){
  lx <- toupper(str_trim(as.character(x)))
  if (lx %in% c("PAIRED", "PE")) return("PAIRED")
  if (lx %in% c("SINGLE", "SE", "SINGLE-END", "SINGLE END")) return("SINGLE")
  return(lx)
}

for (i in seq_len(nrow(df))) {
  row <- df[i, ]

  number        <- sprintf("%03d", row$Number)
  accession     <- row$Accession
  virus         <- row$Virus
  cell_line     <- row$CellLine
  raw_time      <- row$`Infection Time`
  moi_str       <- row$VirusDose
  library_layout<- norm_layout(row$LibraryLayout)
  species       <- row$Species
  strain        <- row$VirusStrain
  safe_strain   <- fs::path_sanitize(strain, replacement = "_")

  # ---- 规范化 infection_time / MOI ----
  raw_time_lc <- tolower(str_trim(as.character(raw_time)))
  infection_time <- if (str_detect(raw_time_lc, "unknown|unkown|na|n/a")) "unknown" else toupper(str_trim(as.character(raw_time)))

  moi_lc <- tolower(str_trim(as.character(moi_str)))
  MOI <- if (str_detect(moi_lc, "unknown|unkown|na|n/a")) "unknown" else str_extract(as.character(moi_str), "[0-9.]+")
  if (is.na(MOI) || MOI == "") MOI <- "unknown"

  # ---- 输出 JSON 目录 ----
  time_moi <- paste0(infection_time, "_MOI", MOI)
  json_out_dir <- file.path("/mnt/alamo01/projects/Group_Wang/sampledata", virus, cell_line, paste0(time_moi, "_", accession, "_", safe_strain))
  dir_create(json_out_dir)

  # ---- 解析分组 ID ----
  split_ids <- function(x){
    if (is.null(x) || is.na(x)) return(character(0))
    unlist(str_split(x, "\\s*[,\\n]\\s*"))
  }
  mock_gsm      <- split_ids(row$Mock)
  infected_gsm  <- split_ids(row$Infected)
  SRR_mock      <- split_ids(row$SRR_mock)
  SRR_infected  <- split_ids(row$SRR_infected)

  # ---- 元数据 JSON ----
  version <- paste0("R-", R.version$major, ".", R.version$minor)
  today <- Sys.Date()
  metadata <- list(
    experiment = list(
      virus = virus,
      type = row$Type,
      cell_line = cell_line,
      strain = strain,
      infection_time = infection_time,
      MOI = MOI,
      LibraryLayout = library_layout,
      replicates = row$Replicates,
      control_condition = "Mock",
      Mock = mock_gsm,
      Infected = infected_gsm,
      SRR_mock = SRR_mock,
      SRR_infected = SRR_infected
    ),
    data_source = list(
      gse = accession,
      comparison = paste("Mock vs Infected ", row$Virus, "(", row$VirusStrain, ")", sep = ""),
      organism = species,
      reference = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", accession)
    ),
    analysis = list(
      tool = "DESeq2",
      version = version,
      date = as.character(today),
      parameters = list(
        pvalue_cutoff = 0.05,
        log2fc_cutoff = 1,
        normalization = "rlog"
      )
    )
  )
  json_path <- file.path(json_out_dir, "meta_analysis_info.json")
  write_json(metadata, json_path, pretty = TRUE, auto_unbox = TRUE)
  message("✅ JSON written: ", json_path)

  # ---- 选择矩阵文件：按 LibraryLayout 走你当前 pipeline 的输出 ----
  align_dir <- if (library_layout == "PAIRED") PAIRED_ALIGN_DIR else SINGLE_ALIGN_DIR
  matrix_path <- file.path(align_dir, "gene_count_matrix_symbol_merged.csv")
  tpm_path    <- file.path(align_dir, "gene_tpm_matrix_symbol_merged.csv")

  if (!file.exists(matrix_path)) {
    message("❌ Matrix not found at: ", matrix_path, "  (layout=", library_layout, ")  跳过该行")
    next
  }
  if (!file.exists(tpm_path)) {
    message("❌ TPM matrix not found at: ", tpm_path, "  (layout=", library_layout, ")  跳过该行")
    next
  }

  # ---- 选择使用 GSM 还是 SRR ----
  use_gsm <- length(SRR_mock) > 0 && all(tolower(SRR_mock) == "merged")
  if (use_gsm) {
    mock_ids     <- mock_gsm
    infected_ids <- infected_gsm
  } else {
    mock_ids     <- SRR_mock
    infected_ids <- SRR_infected
  }
  # 兜底：都为空时提示并跳过
  if (length(c(mock_ids, infected_ids)) == 0) {
    message("❌ No sample IDs (GSM/SRR) provided for accession: ", accession, "  跳过该行")
    next
  }

  # ---- 读取与挑列：Counts ----
  count_matrix <- read_csv(matrix_path, show_col_types = FALSE)
  # 第一列改名为 Gene
  if (ncol(count_matrix) > 0) {
    names(count_matrix)[1] <- "Gene"
  }
  colnames_to_extract <- unique(c(mock_ids, infected_ids))
  found_cols_counts <- intersect(colnames(count_matrix), colnames_to_extract)
  if (length(found_cols_counts) == 0) {
    message("⚠️ None of requested columns found in COUNT matrix for ", accession, "  (IDs: ", paste(colnames_to_extract, collapse=", "), ")")
  }
  expression_matrix <- count_matrix %>% select(any_of(c("Gene", found_cols_counts)))
  expression_out_path <- file.path(json_out_dir, "count_matrix.csv")
  write_csv(expression_matrix, expression_out_path)

  # ---- 读取与挑列：TPM（修正：用 tpm_matrix 的列名匹配）----
  tpm_matrix <- read_csv(tpm_path, show_col_types = FALSE)
  if (ncol(tpm_matrix) > 0) {
    names(tpm_matrix)[1] <- "Gene"
  }
  found_cols_tpm <- intersect(colnames(tpm_matrix), colnames_to_extract)
  if (length(found_cols_tpm) == 0) {
    message("⚠️ None of requested columns found in TPM matrix for ", accession)
  }
  expression_tpm_matrix <- tpm_matrix %>% select(any_of(c("Gene", found_cols_tpm)))
  tpm_out_path <- file.path(json_out_dir, "tpm_matrix.csv")
  write_csv(expression_tpm_matrix, tpm_out_path)

  # ---- 组信息 JSON ----
  rld_metadata <- list(Mock = mock_ids, Infected = infected_ids)
  rld_json_path <- file.path(json_out_dir, "meta_group_info.json")
  write_json(rld_metadata, rld_json_path, pretty = TRUE, auto_unbox = TRUE)

  message("📊 Saved: ", expression_out_path, " | ", tpm_out_path, " | ", rld_json_path)
}

```



