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
    parser.add_argument('--bed_no_lastexon', default=None,
                        help='BED file without lastexon; used to correct Length by removing lastexon portion')
    return parser.parse_args()


def extract_group_name(file_path):
    """Extract group name from file path (without extension)"""
    return os.path.splitext(os.path.basename(file_path))[0]


def filter_groups_based_on_max_length(group, length):
    """Filter groups by maximum transcript length"""
    grouped = group.groupby('Transcript')
    print(length)
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
    results = {}

    if len(group) == 1:
        for l_col, i_col, p_col, pp_col in zip(average_UTRlength_columns, index_UTR_columns, PDUI_columns, PPUI_columns):
            results[l_col] = group["Length"]
            results[i_col] = 1.0
            results[p_col] = 1.0
            results[pp_col] = 0.0
    else:
        for u_col, l_col, i_col, p_col, pp_col in zip(usage_columns, average_UTRlength_columns, index_UTR_columns,
                                                      PDUI_columns, PPUI_columns):
            sum_result = (group[u_col] * group["Length"]).sum()
            results[l_col] = [sum_result] * len(group)
            results[i_col] = [sum_result / group["Length"].max()] * len(group)
            max_length_index = group['Length'].idxmax()
            min_length_index = group["Length"].idxmin()
            p_result = group.loc[max_length_index, u_col]
            pp_result = group.loc[min_length_index, u_col]
            results[p_col] = [p_result] * len(group)
            results[pp_col] = [pp_result] * len(group)
    results_df = pd.DataFrame(results, index=group.index)
    group = pd.concat([group, results_df], axis=1)
    return group


def calculate_cpm(df):
    """
    Calculate CPM for each sample.
    Columns ending with '_TPM' are actually NumReads.
    Adds corresponding '_CPM' columns.
    """
    read_cols = [col for col in df.columns if col.endswith('_TPM')]
    for col in read_cols:
        total_reads = df[col].sum()
        cpm_col = col.replace('_TPM', '_CPM')
        if total_reads > 0:
            df[cpm_col] = df[col] / total_reads * 1e6
        else:
            df[cpm_col] = 0.0
        logging.info(f"  {col}: total reads = {total_reads:.0f}, CPM calculated → {cpm_col}")
    return df


def load_bed_no_lastexon(bed_path):
    bed_dict = {}
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            chrom, start, end, name, score, strand = (
                parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4], parts[5]
            )
            key = (name, strand)
            if key not in bed_dict:
                bed_dict[key] = []
            bed_dict[key].append({'chrom': chrom, 'start': start, 'end': end})
    logging.info(f"Loaded BED (no lastexon): {len(bed_dict)} transcript-strand entries from {bed_path}")
    return bed_dict


def fix_length_remove_lastexon(df, bed_dict):
    """
    For each row, correct the Length by removing the lastexon portion using the
    no-lastexon BED file:
      - Positive strand (+): replace start with BED start (lastexon is upstream/left)
      - Negative strand (-): replace end with BED end (lastexon is upstream/right)
    Length is recalculated as new_end - new_start.
    The 'start' and 'end' columns in df are also updated accordingly.
    Rows with no BED match are kept unchanged (with a warning).
    """
    fixed = 0
    not_found = set()

    for idx, row in df.iterrows():
        transcript = row['Transcript']
        strand = row['strand']
        orig_start = int(row['start'])
        orig_end = int(row['end'])

        key = (transcript, strand)
        if key not in bed_dict:
            not_found.add(transcript)
            continue

        # Pick the matching BED entry (by coordinate overlap if multiple)
        candidates = bed_dict[key]
        entry = None
        if len(candidates) == 1:
            entry = candidates[0]
        else:
            for c in candidates:
                if c['start'] <= orig_end and c['end'] >= orig_start:
                    entry = c
                    break
        if entry is None:
            not_found.add(transcript)
            continue

        if strand == '+':
            new_start = entry['start']
            new_end = orig_end
        else:  # '-'
            new_start = orig_start
            new_end = entry['end']

        new_length = new_end - new_start
        if new_length <= 0:
            logging.warning(
                f"Non-positive length after lastexon removal for {transcript} "
                f"({strand}): {new_start}-{new_end}, keeping original"
            )
            continue

        df.at[idx, 'start'] = str(new_start)
        df.at[idx, 'end'] = str(new_end)
        df.at[idx, 'Length'] = new_length
        fixed += 1

    if not_found:
        logging.warning(
            f"No BED match for {len(not_found)} transcript(s), keeping original lengths. "
            f"Examples: {list(not_found)[:5]}"
        )
    logging.info(f"Lastexon length correction: {fixed} / {len(df)} rows updated")
    return df


def get_csv(df, output_dir):
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
    setup_logging()
    args = parse_arguments()

    logging.info("Script execution started")


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
    data_raw = data_raw.drop(columns=['tmp'])
    logging.info("Completed splitting and extracting information from Name column")

    if args.bed_no_lastexon:
        bed_dict = load_bed_no_lastexon(args.bed_no_lastexon)
        data_raw = fix_length_remove_lastexon(data_raw, bed_dict)
    else:
        logging.info("--bed_no_lastexon not provided, skipping lastexon length correction")


    if args.length > 0:
        data_raw = filter_groups_based_on_max_length(data_raw, args.length)
        print(data_raw.shape[0])
        data_raw = data_raw.reset_index(drop=True)
    else:
        data_raw = data_raw.reset_index(drop=True)

    groups = {}
    for group_file in args.group_files:
        group_name = extract_group_name(group_file)
        if group_name in groups:
            logging.error(f"Duplicate group name: {group_name}, ensure unique filenames for each group file")
            sys.exit(1)
        groups[group_name] = load_group_samples(group_file)

    data_raw = calculate_cpm(data_raw)

    cpm_cols = [col for col in data_raw.columns if col.endswith('_CPM')]
    low_cpm_mask = data_raw[cpm_cols].mean(axis=1) < 1 
    df_filtered = data_raw[~low_cpm_mask].reset_index(drop=True)
    logging.info(f"Filtered data row count: {df_filtered.shape[0]}")
    logging.info(f"Removed {low_cpm_mask.sum()} rows where mean CPM across all samples CPM < 1")

    transcript_order = df_filtered['Transcript'].copy()
    grouped = df_filtered.groupby('Transcript', group_keys=False)
    df_usage = grouped.apply(calculate_usage).reset_index(drop=True)
    if 'Transcript' not in df_usage.columns:
        df_usage.insert(df_usage.columns.get_loc('strand'), 'Transcript', transcript_order.values)
    logging.info("Completed usage calculation")

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

    transcript_order_final = df_usage_filtered['Transcript'].copy()
    grouped_final = df_usage_filtered.groupby('Transcript', group_keys=False)
    final_result = grouped_final.apply(calculate_other).reset_index(drop=True)
    if 'Transcript' not in final_result.columns:
        final_result.insert(final_result.columns.get_loc('strand'), 'Transcript', transcript_order_final.values)
    logging.info("Completed additional metric calculations")

    get_csv(final_result, args.output_dir)

    logging.info("Script execution completed")


if __name__ == "__main__":
    main()
