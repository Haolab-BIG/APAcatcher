#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: [PengCheng]
Date: 2024-11-08

描述:
    该脚本读取用户提供的任意数量的样本组文件和合并后的数据文件，计算TPM均值、usage均值，并导出多个CSV文件供进一步分析使用。
"""

import argparse
import pandas as pd
import numpy as np
import os
import logging
import sys


def setup_logging(log_file='script.log'):
    """配置日志记录"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='3\' UTR分析专业化脚本')
    parser.add_argument('--group_files', nargs='+', required=True, help='多个样本组文件路径')
    parser.add_argument('--merge_file', required=True, help='合并后的数据文件路径')
    parser.add_argument('--output_dir', default='output', help='输出CSV文件的目录')
    parser.add_argument("--length", default=-1,type=int ,help="filter length")
    return parser.parse_args()


def extract_group_name(file_path):
    """根据文件路径提取组名（不包含扩展名）"""
    return os.path.splitext(os.path.basename(file_path))[0]


def filter_groups_based_on_max_length(group, length):
    grouped = group.groupby('Transcript')
    print(length)
    # 遍历每个分组
    for group_name, group_data in grouped:
        max_length = group_data['Length'].max()
        if max_length < length:
            group = group[group['Transcript'] != group_name]
    logging.info("处理Transcript列，删除短3UTR")
    return group


def load_group_samples(group_file):
    """从文件中加载组样本名"""
    try:
        with open(group_file, 'r') as file:
            samples = file.read().splitlines()
            logging.info(f"加载组样本: {group_file}，样本数: {len(samples)}")
            return samples
    except Exception as e:
        logging.error(f"无法读取文件 {group_file}: {e}")
        sys.exit(1)


def add_group_TPM_mean_column(df, group_samples, group_name):
    """计算并添加每组的TPM均值列"""
    selected_columns = [col for col in df.columns if
                        any(sample in col for sample in group_samples) and col.endswith('_TPM')]
    if not selected_columns:
        logging.warning(f"没有找到符合条件的TPM列，组名: {group_name}")
    df[f"{group_name}_TPM_mean"] = df[selected_columns].mean(axis=1)


def add_group_usage_mean_column(df, group_samples, group_name):
    """计算并添加每组的usage均值列"""
    selected_columns = [col for col in df.columns if
                        any(sample in col for sample in group_samples) and col.endswith('_usage')]
    if not selected_columns:
        logging.warning(f"没有找到符合条件的usage列，组名: {group_name}")
    df[f"{group_name}_usage_mean"] = df[selected_columns].mean(axis=1)


def process_transcript_column(df, column_name):
    """处理Transcript列，并去重"""
    df = df.copy()
    df = df.reset_index(drop=True)
    df = df.drop_duplicates(subset=[column_name]).reset_index(drop=True)
    return df


def calculate_usage(group):
    """计算每个样本的usage"""
    tpm_columns = [col for col in group.columns if col.endswith('_TPM')]
    usage_columns = [col.replace('_TPM', '_usage') for col in tpm_columns]

    if len(group) == 1:
        for u_col in usage_columns:
            group[u_col] = 1.0
    else:
        tpm_sum = group[tpm_columns].sum()
        for t_col, u_col in zip(tpm_columns, usage_columns):
            if tpm_sum[t_col] == 0:
                group[u_col] = 0
            else:
                group[u_col] = group[t_col] / tpm_sum[t_col]
    return group


def calculate_other(group):
    usage_columns = [col for col in group.columns if col.endswith('_usage')]
    average_UTRlength_columns = [col.replace('_usage', '_averageLength') for col in usage_columns]
    index_UTR_columns = [col.replace('_usage', '_indexUTR') for col in usage_columns]
    PDUI_columns = [col.replace('_usage', '_PDUI') for col in usage_columns]
    PPUI_columns = [col.replace('_usage', '_PPUI') for col in usage_columns]
    # 用来保存计算结果的字典
    results = {}

    if len(group) == 1:
        for l_col, i_col, p_col,pp_col in zip(average_UTRlength_columns, index_UTR_columns, PDUI_columns,PPUI_columns):
            results[l_col] = group["Length"]
            results[i_col] = 1.0
            results[p_col] = 1.0
            results[pp_col] = 0.0
    else:
        # 当组内有多行数据时，按样本计算
        for u_col, l_col, i_col, p_col,pp_col in zip(usage_columns, average_UTRlength_columns, index_UTR_columns,
                                              PDUI_columns,PPUI_columns):
            sum_result = (group[u_col] * group["Length"]).sum()
            results[l_col] = [sum_result] * len(group)
            results[i_col] = [sum_result / group["Length"].max()] * len(group)
            # 找到 Length 最大值所在的索引，并取对应的 u_col 值
            max_length_index = group['Length'].idxmax()
            min_length_index = group["Length"].idxmin()
            p_result = group.loc[max_length_index, u_col]
            pp_result = group.loc[min_length_index,u_col]
            results[p_col] = [p_result] * len(group)
            results[pp_col] = [pp_result] * len(group)
    results_df = pd.DataFrame(results, index=group.index)
    group = pd.concat([group, results_df], axis=1)
    return group


def get_csv(df, output_dir):
    """处理并导出各类CSV文件"""
    base_col = ["Name", "Length", 'Transcript', 'start', 'end', 'strand']

    usage_col = base_col + [col for col in df.columns if col.endswith('_usage')]
    tpm_col = base_col + [col for col in df.columns if col.endswith('_TPM')]
    averageLength_col = base_col + [col for col in df.columns if col.endswith('_averageLength')]
    PDUI_col = base_col + [col for col in df.columns if col.endswith('_PDUI')]
    index_col = base_col + [col for col in df.columns if col.endswith('_indexUTR')]
    PPUI_col = base_col + [col for col in df.columns if col.endswith('_PPUI')]
    os.makedirs(output_dir, exist_ok=True)

    df_usage = df[usage_col]
    df_tpm = df[tpm_col]
    df_averageLength = process_transcript_column(df[averageLength_col], 'Transcript')
    df_PDUi = process_transcript_column(df[PDUI_col], 'Transcript')
    df_PPUI = process_transcript_column(df[PPUI_col],"Transcript")
    df_indexUTR = process_transcript_column(df[index_col], 'Transcript')
    #df_PDUi = df[PDUI_col]
    #df_PPUI = df[PPUI_col]
    #df_indexUTR = df[index_col]
    df_indexUTR = df_indexUTR.replace(0, np.nan)
    df_PDUi = df_PDUi.replace(0, np.nan)
    df_PPUI = df_PPUI.replace(0,np.nan)

    df_usage.to_csv(os.path.join(output_dir, "3UTR_usage.txt"), sep="\t", index=False)
    df_tpm.to_csv(os.path.join(output_dir, "TPM.txt"), sep="\t", index=False)
    df_averageLength.to_csv(os.path.join(output_dir, "3UTR_averageLength.txt"), sep="\t", index=False)
    df_PDUi.to_csv(os.path.join(output_dir, "PDUI.txt"), sep="\t", index=False)
    df_indexUTR.to_csv(os.path.join(output_dir, "3UTR_index.txt"), sep="\t", index=False)
    df_PPUI.to_csv(os.path.join(output_dir, "PPUI.txt"), sep="\t", index=False)
    logging.info(f"所有CSV文件已保存至 {output_dir}")


def main():
    """主函数"""
    setup_logging()
    args = parse_arguments()

    logging.info("脚本开始运行")

    # 读取数据
    try:
        data_raw = pd.read_csv(args.merge_file, sep='\t')
        logging.info(f"成功读取合并文件: {args.merge_file}, 行数: {data_raw.shape[0]}, 列数: {data_raw.shape[1]}")
    except Exception as e:
        logging.error(f"无法读取合并文件 {args.merge_file}: {e}")
        sys.exit(1)
    data_raw['tmp'] = data_raw['Name'].str.replace("::", ":")
    data_split = data_raw['tmp'].str.split(':', expand=True)
    data_split.columns = ['Transcript', 'strand', 'chr', 'position']
    data_split[['start', 'end']] = data_split['position'].str.split('-', expand=True)
    data_raw = pd.concat([data_raw, data_split[['Transcript', 'strand', 'chr', 'start', 'end']]], axis=1)
    logging.info("完成Name列的分割和提取信息")
    if args.length >0:
        data_raw = filter_groups_based_on_max_length(data_raw, args.length)
        print(data_raw.shape[0])
        data_raw = data_raw.reset_index(drop=True)
    else:
        data_raw = data_raw.reset_index(drop=True)
    # 加载所有样本组
    groups = {}
    for group_file in args.group_files:
        group_name = extract_group_name(group_file)
        if group_name in groups:
            logging.error(f"组名重复: {group_name}，请确保每个组文件具有唯一的文件名")
            sys.exit(1)
        groups[group_name] = load_group_samples(group_file)

    # 计算并添加每组的TPM均值列
    for group_name, samples in groups.items():
        add_group_TPM_mean_column(data_raw, samples, group_name)

    # 构建筛选条件：所有组的TPM均值均小于5
    filter_condition = " & ".join([f"(data_raw['{group_name}_TPM_mean'] < 5)" for group_name in groups.keys()])
    df_filtered = data_raw[~eval(filter_condition)].reset_index(drop=True)
    logging.info(f"筛选后数据行数: {df_filtered.shape[0]}")

    # 分割Name列，提取Transcript等信息

    # 根据Transcript进行分组并计算usage
    grouped = df_filtered.groupby('Transcript')
    df_usage = grouped.apply(calculate_usage).reset_index(drop=True)
    logging.info("完成usage的计算")

    # 计算每组的usage均值
    for group_name, samples in groups.items():
        add_group_usage_mean_column(df_usage, samples, group_name)

    filter_condition_usage = " & ".join(
        [f"(df_usage['{group_name}_usage_mean'] < 0.05)" for group_name in groups.keys()])
    df_usage_filtered = df_usage[~eval(filter_condition_usage)].reset_index(drop=True)
    logging.info(f"筛选后usage数据行数: {df_usage_filtered.shape[0]}")

    try:
        df_usage_filtered["Length"] = df_usage_filtered["Length"].astype(int)
        logging.info("将Length列转换为整数类型")
    except Exception as e:
        logging.error(f"转换Length列类型失败: {e}")
        sys.exit(1)

    # 计算其他指标
    grouped_final = df_usage_filtered.groupby('Transcript')
    final_result = grouped_final.apply(calculate_other).reset_index(drop=True)
    logging.info("计算其他指标完成")

    # 导出CSV文件
    get_csv(final_result, args.output_dir)

    logging.info("脚本运行完成")


if __name__ == "__main__":
    main()

