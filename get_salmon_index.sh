#!/bin/bash

# 设置脚本在遇到错误时立即退出
set -e

# 检查是否提供了足够的参数
if [ "$#" -lt 6 ]; then
    echo "Usage: $0 <final_site_bed> <refseq_utr_bed> <refseq_last_bed> <hg38_fa> <output_fa> <output_index>"
    exit 1
fi

# 输入参数
final_site_bed=$1          # 用户提供的 final_site.bed 文件
refseq_utr_bed=$2           # 用户提供的 RefSeq_UTR_final.bed 文件
refseq_last_bed=$3          # 用户提供的 refseq_last_exon.bed 文件
hg38_fa=$4                 # 用户提供的 UCSC_hg38.fa 文件
output_fa=$5               # 用户提供的 3UTRisoforms_sequences.fa 输出文件
output_index=$6            # 用户提供的 Salmon 索引输出目录

# 设置 Salmon 使用的线程数（根据需要调整）
THREADS=8

# 为 final_site_bed 添加递增的唯一 ID 并生成 final_site_with_id.bed
echo "Generating final_site_with_id.bed with unique incremental IDs..."
awk 'BEGIN {OFS="\t"} 
     { 
         $4 = NR;  # 将第四列替换为行号（1, 2, 3, ...）
         print 
     }' "$final_site_bed" > final_site_with_id.bed

# 生成 final_site_7col.bed
echo "Generating final_site_7col.bed..."
awk -F'\t' 'BEGIN {OFS="\t"} {
    if ($6 == "+") $7 = $3; 
    else if ($6 == "-") $7 = $2; 
    print $0
}' final_site_with_id.bed > final_site_7col.bed

# 生成 RefSeq_UTR_final_7col.bed
echo "Generating RefSeq_UTR_final_7col.bed..."
awk -F'\t' 'BEGIN {OFS="\t"} {
    if ($6 == "+") $7 = $2; 
    else if ($6 == "-") $7 = $3; 
    print $0
}' "$refseq_utr_bed" > RefSeq_UTR_final_7col.bed

# 生成 3UTRisoforms_temp.bed
echo "Generating 3UTRisoforms_temp.bed..."
awk 'BEGIN {FS=OFS="\t"} 
    NR==FNR { ref[$4] = $7; next } 
    $5 in ref { $8 = ref[$5] } 
    { print $0}' RefSeq_UTR_final_7col.bed final_site_7col.bed | \
awk 'BEGIN {FS=OFS="\t"} {
    if ($7 < $8) { 
        $2 = $7; 
        $3 = $8; 
    } else { 
        $2 = $8; 
        $3 = $7; 
    } 
    $7 = ""; 
    $8 = ""; 
    sub(/\t+$/, ""); 
    print 
}' | \
awk 'BEGIN {FS=OFS="\t"} {
    $4 = $5; 
    $7 = $5; 
    $5 = ($2 > $3) ? ($2 - $3) : ($3 - $2); 
    print 
}' > 3UTRisoforms_temp.bed

# 过滤出长度 ≥ 100 的条目，并生成 3UTRisoforms.bed
echo "Filtering 3UTRisoforms_temp.bed for lengths >= 100 and generating 3UTRisoforms.bed..."
awk 'BEGIN {FS=OFS="\t"} {
    id = NR      # 生成唯一的递增 ID
    $4 = id      # 更新第4列为唯一 ID
    $7 = ""      # 清空第7列
    NF = 6       # 确保只有6列
    if ($5 >= 100) 
        print
}' 3UTRisoforms_temp.bed > 3UTRisoforms.bed

# 删除临时文件
rm 3UTRisoforms_temp.bed

# 根据唯一 ID 匹配 final_site_bed 和 3UTRisoforms.bed
echo "Matching 3UTRisoforms.bed with final_site_with_id.bed based on unique ID..."
awk 'BEGIN {FS=OFS="\t"}
    NR==FNR { id[$4] = $0; next }
    $4 in id { 
        split(id[$4], a, "\t") 
        # 根据需要组合列，这里假设保留 final_site 的前7列，加上 3UTRisoforms 的前6列
        print a[1], a[2], a[3], a[4], a[5], a[6] 
    }' final_site_7col.bed 3UTRisoforms.bed > matched_final_site.bed
rm 3UTRisoforms.bed
rm final_site_7col.bed
rm RefSeq_UTR_final_7col.bed
# 使用 refseq_last_bed 生成 Last_Exon_3UTRisoforms.bed
echo "Generating Last_Exon_3UTRisoforms.bed using refseq_last_bed..."
echo "Generating final_site_7col.bed..."
awk -F'\t' 'BEGIN {OFS="\t"} {
    if ($6 == "+") $7 = $3; 
    else if ($6 == "-") $7 = $2; 
    print $0
}' matched_final_site.bed > final_site_7col.bed

# 检查输入数据，确保没有非法字符或空值
echo "Cleaning RefSeq_UTR_final_7col.bed..."
awk -F'\t' 'BEGIN {OFS="\t"} {
    # 检查是否有空列，或者列数不足
    if (NF < 6 || $6 == "" || $2 == "" || $3 == "") next;
    # 根据正负链调整第7列
    if ($6 == "+") $7 = $2;
    else if ($6 == "-") $7 = $3;
    print $0;
}' "$refseq_last_bed" > RefSeq_UTR_final_7col.bed

# 生成 3UTRisoforms_temp.bed，并处理边界问题
echo "Generating Last 3UTRisoforms_temp.bed..."
awk 'BEGIN {FS=OFS="\t"}
    NR==FNR {ref[$4] = $7; next}
    $5 in ref {
        $8 = ref[$5];
    }
    { print $0 }
' RefSeq_UTR_final_7col.bed final_site_7col.bed | \
awk 'BEGIN {FS=OFS="\t"} {
    # 检查是否有空值或非法数据
    if ($7 == "" || $8 == "") next;
    # 处理边界值
    if ($7 < $8) {
        $2 = $7;
        $3 = $8;
    } else {
        $2 = $8;
        $3 = $7;
    }
    $7 = "";
    $8 = "";
    sub(/\t+$/, "");  # 清理多余制表符
    print;
}' | \
awk 'BEGIN {FS=OFS="\t"} {
    # 确保第四列唯一，避免重复
    if ($2 == "" || $3 == "") next;  # 跳过无效行
    $4 = $5;
    $7 = $5;  # 第七列与第五列相同
    $5 = ($2 > $3) ? ($2 - $3) : ($3 - $2);  # 计算区间长度
    print;
}' > 3UTRisoforms_temp.bed

# 修正拼接和列裁剪问题
echo "Fixing 3UTRisoforms.bed format..."
awk 'BEGIN {FS=OFS="\t"} {
    if (NF < 6) next;  # 跳过列数不足的行
    $4 = $7":"$6;  # 拼接列
    $7 = "";       # 清空第七列
    NF = 6;        # 强制裁剪到6列
    print;
}' 3UTRisoforms_temp.bed > 3UTRisoforms.bed

sort -k1,1 -k2,2n 3UTRisoforms.bed > 3UTRisoforms_sorted.bed
mv 3UTRisoforms_sorted.bed 3UTRisoforms.bed

# 提取 3'UTR 序列并生成输出 FASTA 文件
echo "Extracting sequences to $output_fa..."
bedtools getfasta -fi "$hg38_fa" -fo "$output_fa" -bed 3UTRisoforms.bed -name
# 检查 Salmon 是否在 PATH 中
# 构建 Salmon 索引
echo "Building Salmon index at $output_index..."
salmon index -t "$output_fa" -k 31 -i "$output_index" -p "$THREADS"

# 清理中间文件（可选）
rm final_site_with_id.bed
rm final_site_7col.bed
rm RefSeq_UTR_final_7col.bed
rm matched_final_site.bed
rm 3UTRisoforms_temp.bed
echo "Pipeline complete."
