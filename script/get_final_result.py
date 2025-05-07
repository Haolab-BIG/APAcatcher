"""
Author: [ChengPeng]
"""
import argparse
import pandas as pd
import numpy as np
import os
import logging
import sys


def setup_logging(log_file='script.log'):
    """Configure logging settings"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Professional 3\' UTR analysis script')
    parser.add_argument('--group_files', nargs='+', required=True, help='Paths to multiple sample group files')
    parser.add_argument('--merge_file', required=True, help='Path to merged data file')
    parser.add_argument('--output_dir', default='output', help='Directory to save output CSV files')
    parser.add_argument("--length", default=-1, type=int, help="Filter by minimum length")
    return parser.parse_args()


def extract_group_name(file_path):
    """Extract group name from file path (without extension)"""
    return os.path.splitext(os.path.basename(file_path))[0]


def filter_groups_based_on_max_length(group, length):
    """Filter groups by maximum transcript length"""
    grouped = group.groupby('Transcript')
    print(length)
    # Iterate through each group
    for group_name, group_data in grouped:
        max_length = group_data['Length'].max()
        if max_length < length:
            group = group[group['Transcript'] != group_name]
    logging.info("Processing Transcript column, removing short 3' UTRs")
    return group


def load_group_samples(group_file):
    """Load sample names from group file"""
    try:
        with open(group_file, 'r') as file:
            samples = file.read().splitlines()
            logging.info(f"Loaded group samples: {group_file}, sample count: {len(samples)}")
            return samples
    except Exception as e:
        logging.error(f"Failed to read file {group_file}: {e}")
        sys.exit(1)


def add_group_TPM_mean_column(df, group_samples, group_name):
    """Calculate and add TPM mean column for each group"""
    selected_columns = [col for col in df.columns if
                        any(sample in col for sample in group_samples) and col.endswith('_TPM')]
    if not selected_columns:
        logging.warning(f"No matching TPM columns found, group name: {group_name}")
    df[f"{group_name}_TPM_mean"] = df[selected_columns].mean(axis=1)


def add_group_usage_mean_column(df, group_samples, group_name):
    """Calculate and add usage mean column for each group"""
    selected_columns = [col for col in df.columns if
                        any(sample in col for sample in group_samples) and col.endswith('_usage')]
    if not selected_columns:
        logging.warning(f"No matching usage columns found, group name: {group_name}")
    df[f"{group_name}_usage_mean"] = df[selected_columns].mean(axis=1)


def process_transcript_column(df, column_name):
    """Process and deduplicate Transcript column"""
    df = df.copy()
    df = df.reset_index(drop=True)
    df = df.drop_duplicates(subset=[column_name]).reset_index(drop=True)
    return df


def calculate_usage(group):
    """Calculate usage for each sample"""
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
    """Calculate additional metrics for each group"""
    usage_columns = [col for col in group.columns if col.endswith('_usage')]
    average_UTRlength_columns = [col.replace('_usage', '_averageLength') for col in usage_columns]
    index_UTR_columns = [col.replace('_usage', '_indexUTR') for col in usage_columns]
    PDUI_columns = [col.replace('_usage', '_PDUI') for col in usage_columns]
    PPUI_columns = [col.replace('_usage', '_PPUI') for col in usage_columns]
    # Dictionary to store calculation results
    results = {}

    if len(group) == 1:
        for l_col, i_col, p_col, pp_col in zip(average_UTRlength_columns, index_UTR_columns, PDUI_columns, PPUI_columns):
            results[l_col] = group["Length"]
            results[i_col] = 1.0
            results[p_col] = 1.0
            results[pp_col] = 0.0
    else:
        # Calculate when multiple rows exist in group
        for u_col, l_col, i_col, p_col, pp_col in zip(usage_columns, average_UTRlength_columns, index_UTR_columns,
                                                      PDUI_columns, PPUI_columns):
            sum_result = (group[u_col] * group["Length"]).sum()
            results[l_col] = [sum_result] * len(group)
            results[i_col] = [sum_result / group["Length"].max()] * len(group)
            # Find index of max/min length and get corresponding usage value
            max_length_index = group['Length'].idxmax()
            min_length_index = group["Length"].idxmin()
            p_result = group.loc[max_length_index, u_col]
            pp_result = group.loc[min_length_index, u_col]
            results[p_col] = [p_result] * len(group)
            results[pp_col] = [pp_result] * len(group)
    results_df = pd.DataFrame(results, index=group.index)
    group = pd.concat([group, results_df], axis=1)
    return group


def get_csv(df, output_dir):
    """Process and export various CSV files"""
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
    df_PPUI = process_transcript_column(df[PPUI_col], "Transcript")
    df_indexUTR = process_transcript_column(df[index_col], 'Transcript')

    df_indexUTR = df_indexUTR.replace(0, np.nan)
    df_PDUi = df_PDUi.replace(0, np.nan)
    df_PPUI = df_PPUI.replace(0, np.nan)

    df_usage.to_csv(os.path.join(output_dir, "3UTR_usage.txt"), sep="\t", index=False)
    df_tpm.to_csv(os.path.join(output_dir, "TPM.txt"), sep="\t", index=False)
    df_averageLength.to_csv(os.path.join(output_dir, "3UTR_averageLength.txt"), sep="\t", index=False)
    df_PDUi.to_csv(os.path.join(output_dir, "PDUI.txt"), sep="\t", index=False)
    df_indexUTR.to_csv(os.path.join(output_dir, "3UTR_index.txt"), sep="\t", index=False)
    df_PPUI.to_csv(os.path.join(output_dir, "PPUI.txt"), sep="\t", index=False)
    logging.info(f"All CSV files saved to {output_dir}")


def main():
    """Main function"""
    setup_logging()
    args = parse_arguments()

    logging.info("Script execution started")

    # Load data
    try:
        data_raw = pd.read_csv(args.merge_file, sep='\t')
        logging.info(f"Successfully loaded merged file: {args.merge_file}, rows: {data_raw.shape[0]}, columns: {data_raw.shape[1]}")
    except Exception as e:
        logging.error(f"Failed to load merged file {args.merge_file}: {e}")
        sys.exit(1)
    data_raw['tmp'] = data_raw['Name'].str.replace("::", ":")
    data_split = data_raw['tmp'].str.split(':', expand=True)
    data_split.columns = ['Transcript', 'strand', 'chr', 'position']
    data_split[['start', 'end']] = data_split['position'].str.split('-', expand=True)
    data_raw = pd.concat([data_raw, data_split[['Transcript', 'strand', 'chr', 'start', 'end']]], axis=1)
    logging.info("Completed splitting and extracting information from Name column")
    if args.length > 0:
        data_raw = filter_groups_based_on_max_length(data_raw, args.length)
        print(data_raw.shape[0])
        data_raw = data_raw.reset_index(drop=True)
    else:
        data_raw = data_raw.reset_index(drop=True)
    # Load all sample groups
    groups = {}
    for group_file in args.group_files:
        group_name = extract_group_name(group_file)
        if group_name in groups:
            logging.error(f"Duplicate group name: {group_name}, ensure unique filenames for each group file")
            sys.exit(1)
        groups[group_name] = load_group_samples(group_file)

    # Calculate and add TPM mean columns for each group
    for group_name, samples in groups.items():
        add_group_TPM_mean_column(data_raw, samples, group_name)

    # Build filter condition: all group TPM means are less than 5
    filter_condition = " & ".join([f"(data_raw['{group_name}_TPM_mean'] < 5)" for group_name in groups.keys()])
    df_filtered = data_raw[~eval(filter_condition)].reset_index(drop=True)
    logging.info(f"Filtered data row count: {df_filtered.shape[0]}")

    # Split Name column and extract transcript information

    # Group by Transcript and calculate usage
    grouped = df_filtered.groupby('Transcript')
    df_usage = grouped.apply(calculate_usage).reset_index(drop=True)
    logging.info("Completed usage calculation")

    # Calculate group usage means
    for group_name, samples in groups.items():
        add_group_usage_mean_column(df_usage, samples, group_name)

    filter_condition_usage = " & ".join(
        [f"(df_usage['{group_name}_usage_mean'] < 0.05)" for group_name in groups.keys()])
    df_usage_filtered = df_usage[~eval(filter_condition_usage)].reset_index(drop=True)
    logging.info(f"Filtered usage data row count: {df_usage_filtered.shape[0]}")

    try:
        df_usage_filtered["Length"] = df_usage_filtered["Length"].astype(int)
        logging.info("Converted Length column to integer type")
    except Exception as e:
        logging.error(f"Failed to convert Length column: {e}")
        sys.exit(1)

    # Calculate additional metrics
    grouped_final = df_usage_filtered.groupby('Transcript')
    final_result = grouped_final.apply(calculate_other).reset_index(drop=True)
    logging.info("Completed additional metric calculations")

    # Export CSV files
    get_csv(final_result, args.output_dir)

    logging.info("Script execution completed")


if __name__ == "__main__":
    main()
