#!/bin/bash

# ====================================================
# SCRIPT: combine_bed_files.sh
# DESCRIPTION: Traverse BED files in specified directory, merge, sort,
#              and process to generate final_site.bed
# USAGE: ./combine_bed_files.sh -i <input_dir> -o <output_dir> [OPTIONS]
# ====================================================

set -e

# Initialize variables with default values
INPUT_DIR=""
OUTPUT_DIR=""
MERGE_DISTANCE=70
SORT_MODE="asc"
SCRIPT_NAME=$(basename "$0")

# Display help information
usage() {
    cat <<EOF
${SCRIPT_NAME} - Merge and process BED files

USAGE:
    ${SCRIPT_NAME} -i <INPUT_DIR> -o <OUTPUT_DIR> [OPTIONS]

OPTIONS:
    -i, --input-dir      Directory containing input BED files (required)
    -o, --output-dir     Output directory for processed files (required)
    -d, --merge-dist     Maximum distance between features to merge (default: 70)
    -s, --sort-mode      Sorting mode: 'asc' for ascending, 'desc' for descending (default: asc)
    -h, --help           Show this help message

EXAMPLES:
    ${SCRIPT_NAME} -i ./bed_files -o ./results
    ${SCRIPT_NAME} -i /data/bed -o /output -d 50 -s desc
EOF
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input-dir)
            INPUT_DIR="$2"
            shift
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift
            ;;
        -d|--merge-dist)
            MERGE_DISTANCE="$2"
            shift
            ;;
        -s|--sort-mode)
            SORT_MODE="$2"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: Unknown parameter: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
    shift
done

# Validate required parameters
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required arguments" >&2
    usage >&2
    exit 1
fi

# Validate sort mode
if [[ "$SORT_MODE" != "asc" && "$SORT_MODE" != "desc" ]]; then
    echo "Error: Invalid sort mode. Use 'asc' or 'desc'" >&2
    exit 1
fi

# Configure sort arguments
if [[ "$SORT_MODE" == "asc" ]]; then
    SORT_ARGS="-k1,1 -k2,2n"
else
    SORT_ARGS="-k1,1 -k2,2nr"
fi

# Create output directory if not exists
mkdir -p "$OUTPUT_DIR" || { echo "Error: Failed to create output directory" >&2; exit 1; }

# Change to output directory
cd "$OUTPUT_DIR" || { echo "Error: Cannot access output directory" >&2; exit 1; }

# Create temporary file
TMP_FILE=$(mktemp) || { echo "Error: Failed to create temporary file" >&2; exit 1; }

# Merge BED files
FILE_COUNT=0
for BED_FILE in "${INPUT_DIR}"/*.bed; do
    if [[ -f "$BED_FILE" ]]; then
        cat "$BED_FILE" >> "$TMP_FILE"
        ((FILE_COUNT++))
    fi
done

if [[ $FILE_COUNT -eq 0 ]]; then
    echo "Error: No BED files found in input directory" >&2
    exit 1
fi

# Process BED data
echo "Processing ${FILE_COUNT} BED files..."
sort $SORT_ARGS "$TMP_FILE" | \
bedtools merge -i stdin -d "$MERGE_DISTANCE" -s -c 4,5,6 -o distinct,first,distinct | \
awk 'BEGIN {FS=OFS="\t"} {
    if ($6 == "+") $2 = $3;
    else if ($6 == "-") $3 = $2;
    split($5, arr, ","); $5 = arr[1];
    print
}' > combined_sites.bed

# Cleanup temporary file
rm "$TMP_FILE"

echo "Successfully generated: ${OUTPUT_DIR}/combined_sites.bed"
