#!/usr/bin/env bash

# ====================================================
# Script Name : combine_bed_files.sh
# Description : Iterate through all BED files in a specified directory,
#               merge, sort, process intervals, and output final_site.bed
# Usage       : combine.sh -i <input_dir> -o <output_dir> [-d <merge_distance>]
# Options     :
#   -i INPUT_DIR       Path to directory containing BED files (required)
#   -o OUTPUT_DIR      Path to directory for output files (required)
#   -d MERGE_DISTANCE  Maximum distance between intervals to merge (default: 70)
#   -h                 Show help message and exit
# ====================================================

set -euo pipefail

# Default parameters
MERGE_DISTANCE=70

# Print usage information
usage() {
  cat << EOF
Usage: $0 -i <input_dir> -o <output_dir> [-d <merge_distance>]

Options:
  -i INPUT_DIR       Path to directory containing BED files (required)
  -o OUTPUT_DIR      Path to directory for output files (required)
  -d MERGE_DISTANCE  Maximum separation (in bases) to merge intervals (default: ${MERGE_DISTANCE})
  -h                 Show this help message and exit
EOF
  exit 1
}

# Parse command-line options
while getopts ":i:o:d:h" opt; do
  case ${opt} in
    i ) INPUT_DIR=${OPTARG} ;;
    o ) OUTPUT_DIR=${OPTARG} ;;
    d ) MERGE_DISTANCE=${OPTARG} ;;
    h ) usage ;;
    \? ) echo "Invalid option: -${OPTARG}" >&2 ; usage ;;
    : ) echo "Option -${OPTARG} requires an argument." >&2 ; usage ;;
  esac
done

# Ensure required arguments are provided
if [[ -z "${INPUT_DIR:-}" || -z "${OUTPUT_DIR:-}" ]]; then
  echo "Error: Both -i and -o options are required." >&2
  usage
fi

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}" || { echo "Cannot change to directory ${OUTPUT_DIR}" >&2; exit 1; }

# Temporary file for merging BED contents
temp_file=$(mktemp)
trap 'rm -f "$temp_file"' EXIT

# Concatenate all BED files
shopt -s nullglob
for bed_file in "${INPUT_DIR}"/*.bed; do
  if [[ -r "${bed_file}" ]]; then
    cat "${bed_file}" >> "${temp_file}"
  else
    echo "Warning: Cannot read ${bed_file}, skipping." >&2
  fi
done
shopt -u nullglob

# Merge, sort, merge nearby intervals, and post-process
sort -k1,1 -k2,2n "${temp_file}" | \
bedtools merge -i stdin -d ${MERGE_DISTANCE} -s -c 4,5,6 -o distinct,first,distinct | \
awk 'BEGIN {FS=OFS="\t"} {
  # For '+' strand, set end = start; for '-' strand, set start = end
  if ($6 == "+") $3 = $2 + 1;
  else if ($6 == "-") $2 = $3 - 1;
  # Retain only the first comma-delimited value in column 5
  split($5, arr, ","); $5 = arr[1];
  print
}' > combined_site.bed

echo "Generated file: combined_site.bed"
echo "Processing complete."
