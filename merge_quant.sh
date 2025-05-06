#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Usage: $0 <folder1_list.txt> <base_directory> <output.txt>" 
    exit 1
fi

folder_list_file="$1"
base_dir="$2"
output_file="$3"

if [ ! -f "$folder_list_file" ]; then
    echo "File $folder_list_file does not exist."
    exit 1
fi

if [ ! -d "$base_dir" ]; then
    echo "Directory $base_dir does not exist."
    exit 1
fi

first_folder=$(head -n 1 "$folder_list_file" | tr -d '\r\n')  # 去掉回车符和换行符
if [ -d "$base_dir/$first_folder" ]; then
    awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3}' "$base_dir/$first_folder/quant.sf" > merge.txt
else
    echo "Error: First directory $base_dir/$first_folder does not exist."
    exit 1
fi

# 遍历 folder_list_file 文件中的文件夹路径
while IFS= read -r folder_name; do
    folder_name=$(echo "$folder_name" | tr -d '\r\n')  # 去掉回车符和换行符
    folder_path="$base_dir/$folder_name"

    if [ -d "$folder_path" ]; then
        column_name="${folder_name}_TPM"
        
        # 提取第四列并创建文件
        awk -v col="$column_name" 'BEGIN {FS=OFS="\t"} NR==1 {print col} NR>1 {print $4}' "$folder_path/quant.sf" > "${folder_name}_col4.txt"
        
        # 合并数据
        paste merge.txt "${folder_name}_col4.txt" > merge_temp.txt
        mv merge_temp.txt merge.txt
        
        # 清理中间文件
        rm "${folder_name}_col4.txt"
    else
        echo "Warning: Directory $folder_path does not exist. Skipping."
    fi
done < "$folder_list_file"

mv merge.txt "$output_file"
echo "Merging complete. Output written to $output_file."

