#!/bin/bash

# ====================================================
# Author: [ChengPeng]
# SCRIPT: merge_quant.sh
# DESCRIPTION: Merge TPM columns from multiple Salmon quant.sf files
# USAGE: ./merge_quant.sh [OPTIONS] -l SAMPLE_LIST -b BASE_DIR -o OUTPUT_FILE
# ====================================================

set -euo pipefail

# Configuration
TEMP_FILES=()
SCRIPT_NAME=$(basename "$0")

# Initialize parameters
SAMPLE_LIST=""
BASE_DIR=""
OUTPUT_FILE=""
KEEP_TEMP=false

# Cleanup handler
cleanup() {
    if [[ "$KEEP_TEMP" == false ]]; then
        for file in "${TEMP_FILES[@]}"; do
            [[ -f "$file" ]] && rm -f "$file"
        done
    fi
}

trap cleanup EXIT

# Display help information
show_usage() {
    cat <<EOF
${SCRIPT_NAME} - Salmon Quant Merge Utility

USAGE:
    ${SCRIPT_NAME} [OPTIONS] -l SAMPLE_LIST -b BASE_DIR -o OUTPUT_FILE

MANDATORY ARGUMENTS:
    -l, --sample-list    File containing list of sample directories
    -b, --base-dir       Base directory containing sample folders
    -o, --output         Merged output file path

OPTIONS:
    -k, --keep-temp      Keep intermediate temporary files
    -h, --help           Show this help message
    -v, --version        Display version information

EXAMPLES:
    Basic usage:
    ${SCRIPT_NAME} -l sample_list.txt -b /data/quant -o merged_tpm.txt

    Keep temporary files:
    ${SCRIPT_NAME} -l samples.txt -b /path/to/quants -o matrix.txt --keep-temp
EOF
}

# Parse command-line arguments
parse_arguments() {
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            -l|--sample-list)
                SAMPLE_LIST="$2"
                shift
                ;;
            -b|--base-dir)
                BASE_DIR="$2"
                shift
                ;;
            -o|--output)
                OUTPUT_FILE="$2"
                shift
                ;;
            -k|--keep-temp)
                KEEP_TEMP=true
                ;;
            -h|--help)
                show_usage
                exit 0
                ;;
            *)
                echo "Error: Unknown parameter: $1" >&2
                show_usage >&2
                exit 1
                ;;
        esac
        shift
    done
}

# Validate input parameters
validate_inputs() {
    [[ -z "$SAMPLE_LIST" ]] && { echo "Error: Sample list required" >&2; exit 1; }
    [[ -z "$BASE_DIR" ]] && { echo "Error: Base directory required" >&2; exit 1; }
    [[ -z "$OUTPUT_FILE" ]] && { echo "Error: Output file required" >&2; exit 1; }

    [[ -f "$SAMPLE_LIST" ]] || { echo "Error: Sample list not found: $SAMPLE_LIST" >&2; exit 1; }
    [[ -d "$BASE_DIR" ]] || { echo "Error: Base directory not found: $BASE_DIR" >&2; exit 1; }
}

# Logging system
log() {
    local level=$1
    local message=$2
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${level^^}] $message" >&2
}

# Core processing functions
initialize_matrix() {
    local first_sample=$(head -n 1 "$SAMPLE_LIST" | tr -d '\r\n')
    local first_quant="${BASE_DIR}/${first_sample}/quant.sf"

    [[ -d "${BASE_DIR}/${first_sample}" ]] || { 
        log ERROR "Initial sample directory missing: ${BASE_DIR}/${first_sample}"
        exit 1
    }

    [[ -f "$first_quant" ]] || {
        log ERROR "Quant file missing: $first_quant"
        exit 1
    }

    awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3}' "$first_quant" > merge.txt
    TEMP_FILES+=("merge.txt")
}

merge_samples() {
    while IFS= read -r sample; do
        sample=$(echo "$sample" | tr -d '\r\n')
        local sample_dir="${BASE_DIR}/${sample}"
        local quant_file="${sample_dir}/quant.sf"
        local tmp_column="${sample}_tpm.txt"

        [[ -d "$sample_dir" ]] || {
            log WARNING "Sample directory missing: $sample_dir"
            continue
        }

        [[ -f "$quant_file" ]] || {
            log WARNING "Quant file missing: $quant_file"
            continue
        }

        log INFO "Processing sample: $sample"
        
        # Extract TPM column with header
        awk -v col="${sample}_TPM" '
            BEGIN {FS=OFS="\t"} 
            NR == 1 {print col} 
            NR > 1 {print $4}
        ' "$quant_file" > "$tmp_column"
        
        # Validate column dimensions
        local main_lines=$(wc -l < merge.txt)
        local column_lines=$(wc -l < "$tmp_column")
        
        if [[ "$main_lines" -ne "$column_lines" ]]; then
            log ERROR "Dimension mismatch: $sample (main: $main_lines vs column: $column_lines)"
            rm "$tmp_column"
            continue
        fi

        paste merge.txt "$tmp_column" > merge_temp.txt
        mv merge_temp.txt merge.txt
        rm "$tmp_column"
        
    done < "$SAMPLE_LIST"
}

# Main workflow
main() {
    parse_arguments "$@"
    validate_inputs
    
    log INFO "Starting TPM matrix merge"
    log INFO "System Configuration:"
    log INFO "  Sample List: $SAMPLE_LIST"
    log INFO "  Base Directory: $BASE_DIR"
    log INFO "  Output File: $OUTPUT_FILE"

    initialize_matrix
    merge_samples

    mv merge.txt "$OUTPUT_FILE"
    log SUCCESS "Merge completed successfully. Output: $OUTPUT_FILE"
}

main "$@"
