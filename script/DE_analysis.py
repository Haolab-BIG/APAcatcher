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
        "-n", "--numreads_file",
        required=True,
        help="Path to the NumReads file (tab-separated). Used for chi-square test."
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
        help="Index file or PDUI file."
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



def add_group_NumReads_sum_column(df, group_samples, group_name):
    """Calculate and add the group's NumReads sum column for chi-square test."""
    target_cols = {f"{sample}_NumReads" for sample in group_samples}
    selected_columns = [col for col in df.columns if col in target_cols]

    if not selected_columns:
        logging.warning(f"No matching NumReads columns found for group '{group_name}'.")
        df[f"{group_name}_NumReads_sum"] = 1   # avoid zero counts in chi-square
    else:
        # Round to int and add 1 as pseudocount
        df[f"{group_name}_NumReads_sum"] = (
            df[selected_columns].sum(axis=1).round().astype('int32') + 1
        )
        logging.info(f"Added NumReads_sum column for group '{group_name}' using columns: {selected_columns}")



def add_group_TPM_mean_column(df, group_samples, group_name):
    """Calculate and add the group's TPM mean column."""
    target_cols = {f"{sample}_TPM" for sample in group_samples}
    selected_columns = [col for col in df.columns if col in target_cols]
    if not selected_columns:
        logging.warning(f"No matching TPM columns found for group '{group_name}'.")
        df[f"{group_name}_TPM_mean"] = 0
    else:
        df[f"{group_name}_TPM_mean"] = df[selected_columns].mean(axis=1)
        logging.info(f"Added TPM_mean column for group '{group_name}' using columns: {selected_columns}")


def add_group_index_mean_column(df, group_samples, group_name):
    """Calculate and add the group's Index UTR mean column."""
    target_cols = {f"{sample}_indexUTR" for sample in group_samples}
    selected_columns = [col for col in df.columns if col in target_cols]
    if not selected_columns:
        logging.warning(f"No matching Index UTR columns found for group '{group_name}'.")
        df[f"{group_name}_indexUTR_mean"] = 0
    else:
        df[f"{group_name}_indexUTR_mean"] = df[selected_columns].mean(axis=1)
        logging.info(f"Added indexUTR_mean column for group '{group_name}' using columns: {selected_columns}")


def add_group_pdui_mean_column(df, group_samples, group_name):
    """Calculate and add the group's PDUI mean column."""
    target_cols = {f"{sample}_PDUI" for sample in group_samples}
    selected_columns = [col for col in df.columns if col in target_cols]
    if not selected_columns:
        logging.warning(f"No matching PDUI columns found for group '{group_name}'.")
        df[f"{group_name}_PDUI_mean"] = 0
    else:
        df[f"{group_name}_PDUI_mean"] = df[selected_columns].mean(axis=1)
        logging.info(f"Added PDUI_mean column for group '{group_name}' using columns: {selected_columns}")


def get_group_TPM_column(df, samples, group_name):
    """Return TPM column names for a group, raising if any are missing."""
    index_columns = [f"{sample}_TPM" for sample in samples]
    print(f"Checking for columns in group '{group_name}': {index_columns}")
    missing_columns = [col for col in index_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in DataFrame for group '{group_name}': {missing_columns}")
    return index_columns


def get_group_index_column(df, samples, group_name):
    """Return indexUTR column names for a group, raising if any are missing."""
    index_columns = [f"{sample}_indexUTR" for sample in samples]
    print(f"Checking for columns in group '{group_name}': {index_columns}")
    missing_columns = [col for col in index_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in DataFrame for group '{group_name}': {missing_columns}")
    return index_columns


def get_group_pdui_column(df, samples, group_name):
    """Return PDUI column names for a group, raising if any are missing."""
    index_columns = [f"{sample}_PDUI" for sample in samples]
    print(f"Checking for columns in group '{group_name}': {index_columns}")
    missing_columns = [col for col in index_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in DataFrame for group '{group_name}': {missing_columns}")
    return index_columns



def get_group_NumReads_column(df, samples, group_name):
    """Return NumReads column names for a group, raising if any are missing."""
    numreads_columns = [f"{sample}_NumReads" for sample in samples]
    print(f"Checking for NumReads columns in group '{group_name}': {numreads_columns}")
    missing_columns = [col for col in numreads_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing NumReads columns in DataFrame for group '{group_name}': {missing_columns}")
    return numreads_columns



def chi_square_test(group, A_col, B_col):
    """
    Perform Chi-square test using raw count columns (NumReads sums).
    A_col / B_col should be integer count columns, NOT TPM.
    """
    if group.shape[0] == 1:
        return pd.Series({"chi2_p_value": 1.0})

    contingency_table = [group[A_col].values, group[B_col].values]
    chi2, p, _, _ = chi2_contingency(contingency_table, correction=True)
    return pd.Series({"chi2_p_value": p})


def wilcoxon_rank_sum_test(group, A_col, B_col):
    """
    Perform Wilcoxon rank-sum (Mann-Whitney U) test between two independent groups.
    Returns p=1.0 when:
      - Either group has no valid (non-NaN) values
      - All values across both groups are identical (zero-variance / all-tie),
        which would cause a divide-by-zero in the z-score calculation
    """
    import warnings
    import numpy as np

    a_vals = pd.to_numeric(group[A_col].values.flatten(), errors='coerce')
    b_vals = pd.to_numeric(group[B_col].values.flatten(), errors='coerce')

    # Drop NaNs
    a_valid = a_vals[~np.isnan(a_vals)]
    b_valid = b_vals[~np.isnan(b_vals)]

    # Need at least 1 value per group
    if len(a_valid) == 0 or len(b_valid) == 0:
        return pd.Series({"wilcoxon_p_value": 1.0})

    # All values identical across both groups → zero variance → skip
    combined = np.concatenate([a_valid, b_valid])
    if np.all(combined == combined[0]):
        return pd.Series({"wilcoxon_p_value": 1.0})

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            stat, p = ranksums(a_valid, b_valid, nan_policy='omit')
        return pd.Series({"wilcoxon_p_value": p})
    except RuntimeWarning:
        # Catches divide-by-zero from heavy ties not caught by the check above
        return pd.Series({"wilcoxon_p_value": 1.0})
    except Exception as e:
        logging.error(f"Error performing Wilcoxon test: {e}")
        return pd.Series({"wilcoxon_p_value": 1.0})


def main():
    setup_logging()
    args = parse_arguments()
    logging.info("Loading TPM file...")
    try:
        tpm_df = pd.read_csv(args.tpm_file, sep='\t')
        if 'Transcript' not in tpm_df.columns:
            logging.error("TPM file must contain a 'Transcript' column.")
            sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading TPM file: {e}")
        sys.exit(1)

    logging.info("Loading NumReads file...")
    try:
        numreads_df = pd.read_csv(args.numreads_file, sep='\t')
        if 'Transcript' not in numreads_df.columns:
            logging.error("NumReads file must contain a 'Transcript' column.")
            sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading NumReads file: {e}")
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
    for group_file in [args.group_a, args.group_b]:
        group_name = extract_group_name(group_file)
        if group_name in groups:
            logging.error(f"Duplicate group name '{group_name}'. Ensure each group file has a unique name.")
            sys.exit(1)
        groups[group_name] = load_group_samples(group_file)

    if len(groups) != 2:
        logging.error("This script currently supports exactly two groups for comparison.")
        sys.exit(1)

    group_names = list(groups.keys())
    A_group, B_group = group_names[0], group_names[1]

    for group_name, samples in groups.items():
        add_group_NumReads_sum_column(numreads_df, samples, group_name)

    # TPM mean columns (kept for reporting, not for chi-square)
    for group_name, samples in groups.items():
        add_group_TPM_mean_column(tpm_df, samples, group_name)

    # Index / PDUI mean columns
    if args.model in ["pdui", "PDUI"]:
        for group_name, samples in groups.items():
            add_group_pdui_mean_column(index_df, samples, group_name)
    else:
        for group_name, samples in groups.items():
            add_group_index_mean_column(index_df, samples, group_name)

    # ── Column name aliases ─────────
    A_NumReads_sum_col = f"{A_group}_NumReads_sum"   
    B_NumReads_sum_col = f"{B_group}_NumReads_sum"   
    A_TPM_mean_col     = f"{A_group}_TPM_mean"
    B_TPM_mean_col     = f"{B_group}_TPM_mean"

    if args.model in ["pdui", "PDUI"]:
        A_index_mean_col    = f"{A_group}_PDUI_mean"
        B_index_mean_col    = f"{B_group}_PDUI_mean"
        group_a_index_columns = list(map(str, get_group_pdui_column(index_df, groups[A_group], A_group)))
        group_b_index_columns = list(map(str, get_group_pdui_column(index_df, groups[B_group], B_group)))
    else:
        A_index_mean_col    = f"{A_group}_indexUTR_mean"
        B_index_mean_col    = f"{B_group}_indexUTR_mean"
        group_a_index_columns = list(map(str, get_group_index_column(index_df, groups[A_group], A_group)))
        group_b_index_columns = list(map(str, get_group_index_column(index_df, groups[B_group], B_group)))

    group_a_tpm_columns = list(map(str, get_group_TPM_column(tpm_df, groups[A_group], A_group)))
    group_b_tpm_columns = list(map(str, get_group_TPM_column(tpm_df, groups[B_group], B_group)))

    group_a_numreads_columns = list(map(str, get_group_NumReads_column(numreads_df, groups[A_group], A_group)))
    group_b_numreads_columns = list(map(str, get_group_NumReads_column(numreads_df, groups[B_group], B_group)))

    logging.info("Performing statistical tests...")

    chi2_base_df = numreads_df[
        ['Transcript', A_NumReads_sum_col, B_NumReads_sum_col]
    ].copy()

    results_chi2 = chi2_base_df.groupby("Transcript").apply(
        lambda group: chi_square_test(group, A_NumReads_sum_col, B_NumReads_sum_col)
    ).reset_index()


    # Wilcoxon rank-sum on index/PDUI values (unchanged)
    results_wilcoxon = index_df.groupby("Transcript").apply(
        lambda group: wilcoxon_rank_sum_test(group, group_a_index_columns, group_b_index_columns)
    ).reset_index()

    logging.info("Applying FDR correction to p-values...")

    valid_chi2_p = results_chi2['chi2_p_value'].dropna()
    _, corrected_chi2_p, _, _ = multipletests(valid_chi2_p, method='fdr_bh')
    results_chi2.loc[results_chi2["chi2_p_value"].notna(), "corrected_chi2_p_value"] = corrected_chi2_p

    valid_wilcoxon_p = results_wilcoxon['wilcoxon_p_value'].dropna()
    _, corrected_wilcoxon_p, _, _ = multipletests(valid_wilcoxon_p, method='fdr_bh')
    results_wilcoxon.loc[
        results_wilcoxon["wilcoxon_p_value"].notna(), "corrected_wilcoxon_p_value"
    ] = corrected_wilcoxon_p

    logging.info("Calculating delta index mean between groups...")
    delta_col = A_index_mean_col + '-' + B_index_mean_col

    results_chi2["wilcoxon_p_value"]           = results_wilcoxon["wilcoxon_p_value"]
    results_chi2["corrected_wilcoxon_p_value"] = results_wilcoxon["corrected_wilcoxon_p_value"]
    results_chi2[delta_col]                    = index_df[A_index_mean_col] - index_df[B_index_mean_col]

    # ── Append TPM mean columns to results for reference ─────────────────────
    # Use merge on Transcript to avoid length mismatch (tpm_df may have more
    # rows than results_chi2 because results_chi2 is already deduplicated by
    # groupby("Transcript"), whereas tpm_df retains one row per isoform).
    tpm_mean_ref = (
        tpm_df[['Transcript', A_TPM_mean_col, B_TPM_mean_col]]
        .drop_duplicates(subset='Transcript')
        .reset_index(drop=True)
    )
    results_chi2 = results_chi2.merge(tpm_mean_ref, on='Transcript', how='left')


    logging.info(f"Saving results to {args.output}...")
    results_chi2.to_csv(args.output, sep='\t', index=False)
    logging.info("Analysis complete.")



if __name__ == "__main__":
    main()

