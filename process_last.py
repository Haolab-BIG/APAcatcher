import pandas as pd
import argparse

def process_group(group):
    """处理每个分组，按照链的不同进行过滤"""
    # 处理正链（+）
    if group["strand"].iloc[0] == "+":
        # 查找最大值和第二大值的差异
        max_value = group["end"].max()
        second_max_value = group[group["end"] != max_value]["end"].max()
        if max_value - second_max_value <= 100:
            # 保留最大值，删除第二大的
            group = group[group["end"] != second_max_value]

    # 处理负链（-）
    elif group["strand"].iloc[0] == "-":
        # 查找最小值和第二小值的差异
        min_value = group["start"].min()
        second_min_value = group[group["start"] != min_value]["start"].min()
        if second_min_value - min_value <= 100:
            # 保留最小值，删除第二小的
            group = group[group["start"] != second_min_value]

    return group

def process_bed_file(input_file, output_file):
    """读取BED文件，处理后保存输出文件"""
    # 读取BED文件
    df = pd.read_csv(input_file, sep="\t", header=None, names=["chr", "start", "end", "score", "Transcript", "strand"])

    # 按照 Transcript 和 start 排序
    df_sorted = df.sort_values(by=["Transcript", "start"])

    # 按照 Transcript 分组并应用处理函数
    df_processed = df_sorted.groupby("Transcript").apply(process_group).reset_index(drop=True)

    # 保存处理后的数据到输出文件
    df_processed_sorted = df_processed.sort_values(by=["chr","start","end"])
    df_processed_sorted.to_csv(output_file, sep="\t", header=False, index=False)

    print(f"处理完成，结果已保存到 {output_file}")

def main():
    # 使用 argparse 解析命令行参数
    parser = argparse.ArgumentParser(description="处理 BED 文件，删除符合条件的行")
    parser.add_argument("--input_file", type=str, help="输入 BED 文件的路径")
    parser.add_argument("--output_file", type=str, help="输出处理后的文件路径")
    
    args = parser.parse_args()

    # 调用处理函数
    process_bed_file(args.input_file, args.output_file)

if __name__ == "__main__":
    main()

