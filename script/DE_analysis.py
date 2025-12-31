import pandas as pd
from scipy.stats import wilcoxon, chi2_contingency, ranksums
import os
import sys
import argparse
import logging
from statsmodels.stats.multitest import multipletests
import time


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Automated APA Analysis Script: Compare TPM and Index UTR between two groups."
    )
    parser.add_argument(
        "-t", "--tpm_file",
        required=True,
        help="Path to the TPM file (tab-separated)."
    )
    parser.add_argument(
        "-i", "--index_file",
        required=True,
        help="Path to the Index UTR file (tab-separated)."
    )
    parser.add_argument(
        "-a", "--group_a",
        required=True,
        help="Path to the Group A sample list file (txt)."
    )
    parser.add_argument(
        "-b", "--group_b",
        required=True,
        help="Path to the Group B sample list file (txt)."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output results file (tab-separated)."
    )
    parser.add_argument(
        "-m", "--model",
        required=True,
        type=str,
        help="Index file or PUDI file."
    )
    return parser.parse_args()


def extract_group_name(file_path):
    """Extract group name from the file name (without extension)."""
    return os.path.splitext(os.path.basename(file_path))[0]


def load_group_samples(group_file):
    """Load group samples from a txt file."""
    try:
        with open(group_file, 'r') as file:
            samples = file.read().splitlines()
            logging.info(f"Loaded {len(samples)} samples from {group_file}")
            return samples
    except Exception as e:
        logging.error(f"Error reading group file {group_file}: {e}")
        sys.exit(1)


def add_group_TPM_sum_column(df, group_samples, group_name):
    """Calculate and add the group's TPM sum column."""
    # Find columns that match the group's samples and end with '_TPM'
    selected_columns = [col for col in df.columns if
                        any(sample in col for sample in group_samples) and col.endswith('_TPM')]
    if not selected_columns:
        logging.warning(f"No matching TPM columns found for group '{group_name}'.")
        df[f"{group_name}_TPM_sum"] = 0
    else:
        df[f"{group_name}_TPM_sum"] = (df[selected_columns].sum(axis=1) + 1).astype('int32')
        logging.info(f"Added TPM_sum column for group '{group_name}' using columns: {selected_columns}")


def add_group_TPM_mean_column(df, group_samples, group_name):
    """Calculate and add the group's TPM  mean column."""
    # Find columns that match the group's samples and end with 'TPM'
    selected_columns = [col for col in df.columns if
                        any(sample in col for sample in group_samples) and col.endswith('_TPM')]
    if not selected_columns:
        logging.warning(f"No matching TPM columns found for group '{group_name}'.")
        df[f"{group_name}_TPM_mean"] = 0
    else:
        df[f"{group_name}_TPM_mean"] = df[selected_columns].mean(axis=1)
        logging.info(f"Added TPM_mean column for group '{group_name}' using columns: {selected_columns}")


def add_group_index_mean_column(df, group_samples, group_name):
    """Calculate and add the group's Index UTR mean column."""
    # Find columns that match the group's samples and end with '_indexUTR'
    selected_columns = [col for col in df.columns if
                        any(sample in col for sample in group_samples) and col.endswith('_indexUTR')]
    if not selected_columns:
        logging.warning(f"No matching Index UTR columns found for group '{group_name}'.")
        df[f"{group_name}_indexUTR_mean"] = 0
    else:
        df[f"{group_name}_indexUTR_mean"] = df[selected_columns].mean(axis=1)
        logging.info(f"Added indexUTR_mean column for group '{group_name}' using columns: {selected_columns}")


def add_group_pdui_mean_column(df, group_samples, group_name):
    """Calculate and add the group's PDUI UTR mean column."""
    # Find columns that match the group's samples and end with '_PDUI'
    selected_columns = [col for col in df.columns if
                        any(sample in col for sample in group_samples) and col.endswith('_PDUI')]
    if not selected_columns:
        logging.warning(f"No matching Index UTR columns found for group '{group_name}'.")
        df[f"{group_name}_PDUI_mean"] = 0
    else:
        df[f"{group_name}_PDUI_mean"] = df[selected_columns].mean(axis=1)
        logging.info(f"Added PDUI_mean column for group '{group_name}' using columns: {selected_columns}")


def get_group_TPM_column(df, samples, group_name):
    """
    Find columns in the DataFrame that correspond to the samples for a given group.

    :param df: DataFrame containing columns with sample data
    :param samples: List of sample names for the group
    :param group_name: Name of the group (used for logging)
    :return: List of column names matching the group's sample names
    """
    # Find columns in the DataFrame that correspond to the group's samples
    index_columns = [f"{sample}_TPM" for sample in samples]

    # Debugging print to confirm column names
    print(f"Checking for columns in group '{group_name}': {index_columns}")

    # Verify columns exist
    missing_columns = [col for col in index_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in DataFrame for group '{group_name}': {missing_columns}")

    return index_columns


def get_group_index_column(df, samples, group_name):
    """
    Find columns in the DataFrame that correspond to the samples for a given group.

    :param df: DataFrame containing columns with sample data
    :param samples: List of sample names for the group
    :param group_name: Name of the group (used for logging)
    :return: List of column names matching the group's sample names
    """
    # Find columns in the DataFrame that correspond to the group's samples
    index_columns = [f"{sample}_indexUTR" for sample in samples]

    # Debugging print to confirm column names
    print(f"Checking for columns in group '{group_name}': {index_columns}")

    # Verify columns exist
    missing_columns = [col for col in index_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in DataFrame for group '{group_name}': {missing_columns}")

    return index_columns


def get_group_pdui_column(df, samples, group_name):
    """
    Find columns in the DataFrame that correspond to the samples for a given group.

    :param df: DataFrame containing columns with sample data
    :param samples: List of sample names for the group
    :param group_name: Name of the group (used for logging)
    :return: List of column names matching the group's sample names
    """
    # Find columns in the DataFrame that correspond to the group's samples
    index_columns = [f"{sample}_PDUI" for sample in samples]

    # Debugging print to confirm column names
    print(f"Checking for columns in group '{group_name}': {index_columns}")

    # Verify columns exist
    missing_columns = [col for col in index_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in DataFrame for group '{group_name}': {missing_columns}")

    return index_columns


def chi_square_test(group, A_col, B_col):
    """
    Perform Chi-square test between two categorical variables.
    Note: This function assumes A_col and B_col are categorical counts.
    """
    # 如果 group 只有一行数据，跳过卡方检验
    if group.shape[0] == 1:
        return pd.Series({"chi2_p_value": 1.0})

    contingency_table = [group[A_col].values, group[B_col].values]

    chi2, p, _, _ = chi2_contingency(contingency_table, correction=True)

    return pd.Series({"chi2_p_value": p})


def wilcoxon_rank_sum_test(group, A_col, B_col):
    """
    Perform Wilcoxon rank-sum (Mann-Whitney U) test between two independent groups.
    """
    # Perform Mann-Whitney U test
    try:
        stat, p = ranksums(group[A_col].values.flatten(), group[B_col].values.flatten(),nan_policy='omit')
        return pd.Series({"wilcoxon_p_value": p})
    except Exception as e:
        logging.error(f"Error performing Wilcoxon test: {e}")
        return pd.Series({"wilcoxon_p_value": 1.0})


def main():
    setup_logging()
    args = parse_arguments()

    # Load TPM and Index files
    logging.info("Loading TPM file...")
    try:
        tpm_df = pd.read_csv(args.tpm_file, sep='\t')
        if 'Transcript' not in tpm_df.columns:
            logging.error("TPM file must contain a 'Transcript' column.")
            sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading TPM file: {e}")
        sys.exit(1)

    logging.info("Loading Index UTR file...")
    try:
        index_df = pd.read_csv(args.index_file, sep='\t')
        if 'Transcript' not in index_df.columns:
            logging.error("Index file must contain a 'Transcript' column.")
            sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading Index UTR file: {e}")
        sys.exit(1)

    # Load group samples
    groups = {}
    group_files = [args.group_a, args.group_b]
    for group_file in group_files:
        group_name = extract_group_name(group_file)
        if group_name in groups:
            logging.error(f"Duplicate group name '{group_name}'. Ensure each group file has a unique name.")
            sys.exit(1)
        groups[group_name] = load_group_samples(group_file)

    # Add TPM_sum and indexUTR_mean columns for each group
    for group_name, samples in groups.items():
        add_group_TPM_sum_column(tpm_df, samples, group_name)
        add_group_TPM_mean_column(tpm_df,samples, group_name)

    if args.model in ["pdui","PDUI"]:
        for group_name, samples in groups.items():
            add_group_pdui_mean_column(index_df, samples, group_name)
    else:
        for group_name, samples in groups.items():
            add_group_index_mean_column(index_df, samples, group_name)

    # Merge TPM and Index dataframes on 'Transcript'
    logging.info("Merging TPM and Index UTR data...")
    # Define comparison groups (Assuming two groups)
    if len(groups) != 2:
        logging.error("This script currently supports exactly two groups for comparison.")
        sys.exit(1)

    group_names = list(groups.keys())
    A_group, B_group = group_names[0], group_names[1]

    A_TPM_sum_col = f"{A_group}_TPM_sum"
    B_TPM_sum_col = f"{B_group}_TPM_sum"

    A_TPM_mean_col = f"{A_group}_TPM_mean"
    B_TPM_mean_col = f"{B_group}_TPM_mean"

    if args.model in ["pdui","PDUI"]:
        A_index_mean_col = f"{A_group}_PDUI_mean"
        B_index_mean_col = f"{B_group}_PDUI_mean"
        group_a_index_columns = get_group_pdui_column(index_df, groups[A_group], A_group)
        group_b_index_columns = get_group_pdui_column(index_df, groups[B_group], B_group)
        group_a_index_columns = list(map(str, group_a_index_columns))
        group_b_index_columns = list(map(str, group_b_index_columns))

    else:
        A_index_mean_col = f"{A_group}_indexUTR_mean"
        B_index_mean_col = f"{B_group}_indexUTR_mean"
        group_a_index_columns = get_group_index_column(index_df, groups[A_group], A_group)
        group_b_index_columns = get_group_index_column(index_df, groups[B_group], B_group)
        group_a_index_columns = list(map(str, group_a_index_columns))
        group_b_index_columns = list(map(str, group_b_index_columns))

    group_a_tpm_columns = get_group_TPM_column(tpm_df, groups[A_group], A_group)
    group_b_tpm_columns = get_group_TPM_column(tpm_df, groups[B_group], B_group)
    # Prepare for statistical testing
    logging.info("Performing statistical tests...")
    group_a_tpm_columns = list(map(str, group_a_tpm_columns))
    group_b_tpm_columns = list(map(str, group_b_tpm_columns))
    # Initialize results dataframe
    results_df = tpm_df[['Transcript', A_TPM_sum_col, B_TPM_sum_col] + group_a_tpm_columns + group_b_tpm_columns].copy()
    # Perform Wilcoxon rank-sum test for each transcript

    results_wilcoxon = index_df.groupby("Transcript").apply(
        lambda group: wilcoxon_rank_sum_test(group, group_a_index_columns, group_b_index_columns)).reset_index()

    results_chi2 = results_df.groupby("Transcript").apply(
        lambda group: chi_square_test(group, A_TPM_sum_col, B_TPM_sum_col)).reset_index()

    logging.info("Applying FDR correction to p-values...")
    # Remove transcripts where p-values are NaN
    valid_chi2_p = results_chi2['chi2_p_value'].dropna()
    _, corrected_chi2_p, _, _ = multipletests(valid_chi2_p, method='fdr_bh')
    results_chi2.loc[results_chi2["chi2_p_value"].notna(), "corrected_chi2_p_value"] = corrected_chi2_p

    # For Wilcoxon p-values
    valid_wilcoxon_p = results_wilcoxon['wilcoxon_p_value'].dropna()
    _, corrected_wilcoxon_p, _, _ = multipletests(valid_wilcoxon_p, method='fdr_bh')
    results_wilcoxon.loc[
        results_wilcoxon["wilcoxon_p_value"].notna(), "corrected_wilcoxon_p_value"] = corrected_wilcoxon_p
    #TMP_mean_delta_df = pd.DataFrame()
    # Calculate delta indexUTR mean
    logging.info("Calculating delta indexUTR mean between groups...")
    temp_1 = A_index_mean_col + '-' + B_index_mean_col
    #temp_2 = A_TPM_mean_col + '-' + B_TPM_mean_col
    results_chi2["wilcoxon_p_value"] = results_wilcoxon["wilcoxon_p_value"]
    results_chi2["corrected_wilcoxon_p_value"] = results_wilcoxon["corrected_wilcoxon_p_value"]
    results_chi2[temp_1] = index_df[A_index_mean_col] - index_df[B_index_mean_col]
    #TMP_mean_delta_df[temp_2] = tpm_df[A_TPM_mean_col]-tpm_df[B_TPM_mean_col]
    # Merge corrected p-values back to the main results

    # Save results to output file
    logging.info(f"Saving results to {args.output}...")
    results_chi2.to_csv(args.output, sep='\t', index=False)
    logging.info("Analysis complete.")


if __name__ == "__main__":
    main()

