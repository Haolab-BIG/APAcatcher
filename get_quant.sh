#!/bin/bash

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <index_path> <fastq_directory> <output_directory>"
    exit 1
fi

# 从命令行参数中获取索引路径、fastq 文件目录和输出目录
index="$1"
fastq_dir="$2"
output_dir="$3"
threads=12

# 检查输出目录是否存在，不存在则创建
if [[ ! -d "$output_dir" ]]; then
    echo "Output directory does not exist. Creating directory: $output_dir"
    mkdir -p "$output_dir"
fi

extensions=("fastq.gz" "fq.gz" "fastq" "fq")

for ext in "${extensions[@]}"; do
    for file1 in "$fastq_dir"/*_1*."$ext"; do
        # 检查是否有匹配的文件
        if [[ ! -e "$file1" ]]; then
            continue
        fi

        # 获取样本名，去除 _1 和扩展名
        sample_name=$(basename "$file1")
        sample_name="${sample_name%%_1*}"
	echo $sample_name
        # 使用通配符匹配 _2 和相同的扩展名
        file2=$(find "$fastq_dir" -maxdepth 1 -type f -name "${sample_name}_2*.$ext")
        
	if [[ -z "$file2" ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Warning: Pair file for $file1 not found!"
            continue
        fi

        echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting processing for sample: $sample_name"

        # 运行 salmon quant 命令，并将输出结果放入指定的输出目录
        salmon quant -i "$index" -l IU -1 "$file1" -2 "$file2" \
        -o "${output_dir}/${sample_name}" -p "$threads" --validateMappings --useEM

        if [[ $? -eq 0 ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Completed processing for sample: $sample_name"
        else
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Error occurred while processing sample: $sample_name"
        fi
    done
done

echo "$(date '+%Y-%m-%d %H:%M:%S') - All samples processed."

