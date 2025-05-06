#!/bin/bash

# ====================================================
# 脚本名称: combine_bed_files.sh
# 功能描述: 遍历指定文件夹中的所有 BED 文件，合并、排序并处理，
#          最终生成 final_site.bed 文件。
# 使用方法: ./combine_bed_files.sh /path/to/bed/files /path/to/output
# ====================================================

# 设置脚本在遇到错误时立即退出
set -e

# 检查是否提供了输入文件夹和输出路径
if [ $# -ne 2 ]; then
    echo "用法: $0 /path/to/bed/files /path/to/output"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# 确保输出目录存在
mkdir -p "$OUTPUT_DIR"

# 进入输出目录
cd "$OUTPUT_DIR" || { echo "无法进入目录 $OUTPUT_DIR"; exit 1; }

# 创建临时文件用于合并所有 BED 文件
temp_file=$(mktemp)

# 遍历输入目录中的所有 .bed 文件并合并
for bed_file in "$INPUT_DIR"/*.bed; do
    if [[ -f "$bed_file" ]]; then
        cat "$bed_file" >> "$temp_file"
    else
        echo "未找到 BED 文件在 $INPUT_DIR 中。"
        exit 1
    fi
done

# 合并、排序、合并重叠区域并处理最后的输出
sort -k1,1 -k2,2n "$temp_file" | \
bedtools merge -i stdin -d 70 -s -c 4,5,6 -o distinct,first,distinct | \
awk 'BEGIN {FS=OFS="\t"} {
  if ($6 == "+") $2 = $3;
  else if ($6 == "-") $3 = $2;
  # 清理列5中的多余数据，仅保留第一个有效值
  split($5, arr, ","); $5 = arr[1]; 
  print
}' > combind_site.bed
echo "生成文件: combind_site.bed"

# 清理临时文件
rm "$temp_file"

echo "所有处理完成。"

