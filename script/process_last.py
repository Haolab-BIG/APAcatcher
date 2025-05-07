# Author: [ChengPeng]
import pandas as pd
import argparse
from typing import Optional

def process_transcript_group(group: pd.DataFrame) -> pd.DataFrame:
    """Process transcript groups by strand-specific filtering.
    Args:
        group: DataFrame containing genomic coordinates for a single transcript
    
    Returns:
        Filtered DataFrame with removed entries based on proximity criteria
    """
    if group.empty or len(group) == 1:
        return group

    strand = group["strand"].iloc[0]
    
    if strand == "+":
        sorted_ends = group["end"].sort_values(ascending=False).unique()
        if len(sorted_ends) > 1 and (sorted_ends[0] - sorted_ends[1] <= 100):
            second_max_end = sorted_ends[1]
            return group[group["end"] != second_max_end]
            
    elif strand == "-":
        sorted_starts = group["start"].sort_values().unique()
        if len(sorted_starts) > 1 and (sorted_starts[1] - sorted_starts[0] <= 100):
            second_min_start = sorted_starts[1]
            return group[group["start"] != second_min_start]
    
    return group

def process_bed_file(input_path: str, output_path: str) -> None:
    """Process BED file through genomic coordinate filtering.
    
    Args:
        input_path: Path to input BED file
        output_path: Path for output processed BED file
    """
    # Load data with type validation
    column_names = [
        "chr", "start", "end", 
        "score", "transcript_id", "strand"
    ]
    
    df = pd.read_csv(
        input_path,
        sep="\t",
        header=None,
        names=column_names,
        dtype={
            "chr": "category",
            "start": "int32",
            "end": "int32",
            "strand": "category"
        }
    )

    # Process genomic coordinates
    processed_df = (
        df.sort_values(["transcript_id", "start"])
          .groupby("transcript_id", group_keys=False)
          .apply(process_transcript_group)
    )

    # Save sorted results maintaining BED format
    (
        processed_df.sort_values(["chr", "start", "end"])
        [column_names]  # Maintain original column order
        .to_csv(output_path, sep="\t", header=False, index=False)
    )

    print(f"Processing complete. Results saved to: {output_path}")

def main() -> None:
    """Command-line interface for BED file processing."""
    parser = argparse.ArgumentParser(
        description="BED File Processor - Remove proximal genomic features by strand-specific criteria",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "-i", "--input", 
        required=True,
        help="Input BED file path",
        metavar="INPUT.bed"
    )
    
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output BED file path",
        metavar="OUTPUT.bed"
    )

    args = parser.parse_args()
    
    try:
        process_bed_file(args.input, args.output)
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        raise SystemExit(1) from e

if __name__ == "__main__":
    main()
