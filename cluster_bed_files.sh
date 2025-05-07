#!/bin/bash

# ====================================================
# SCRIPT: cluster_bed_files.sh
# DESCRIPTION: Batch process BED files through sorting, merging,
#              filtering, and generate final_site.bed
# USAGE: ./cluster_bed_files.sh [OPTIONS] -i INPUT_DIR -o OUTPUT_DIR
# ====================================================

set -eo pipefail

# Configuration
DEFAULT_MERGE_DISTANCE=70
MIN_REPLICATE_COUNT=2  
TEMP_FILES=()
SCRIPT_NAME=$(basename "$0")

# Initialize parameters
INPUT_DIR=""
OUTPUT_DIR=""
KEEP_TEMP=false

# Cleanup handler
cleanup() {
    if [[ "$KEEP_TEMP" == false ]]; then
        for file in "${TEMP_FILES[@]}"; do
            if [[ -f "$file" ]]; then
                rm -f "$file"
            fi
        done
    fi
}

trap cleanup EXIT

# Display help information
show_usage() {
    cat <<EOF
${SCRIPT_NAME} - BED Processing Pipeline

USAGE:
    ${SCRIPT_NAME} [OPTIONS] -i INPUT_DIR -o OUTPUT_DIR

MANDATORY ARGUMENTS:
    -i, --input-dir      Directory containing input BED files
    -o, --output-dir     Output directory for results

OPTIONS:
    -d, --merge-dist     Merging distance (default: ${DEFAULT_MERGE_DISTANCE})
    -c, --min-count      Minimum replicate count (default: ${MIN_REPLICATE_COUNT})
    -k, --keep-temp      Keep intermediate files
    -h, --help           Show this help message


EXAMPLES:
    Basic usage:
    ${SCRIPT_NAME} -i ./bed_files -o ./results

    Custom parameters:
    ${SCRIPT_NAME} -i /data/bed -o /output \\
                   -d 100 -c 3 --keep-temp
EOF
}

# Parse command-line arguments
parse_arguments() {
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
                DEFAULT_MERGE_DISTANCE="$2"
                shift
                ;;
            -c|--min-count)
                MIN_REPLICATE_COUNT="$2"
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
    if ! [[ "$MIN_REPLICATE_COUNT" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: Invalid min-count value: must be positive integer" >&2
        exit 1
    fi

    if ! [[ "$DEFAULT_MERGE_DISTANCE" =~ ^[0-9]+$ ]]; then
        echo "Error: Invalid merge distance: must be non-negative integer" >&2
        exit 1
    fi
    
    if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
        echo "Error: Missing required arguments" >&2
        show_usage >&2
        exit 1
    fi

    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "Error: Input directory not found: $INPUT_DIR" >&2
        exit 1
    fi

    mkdir -p "$OUTPUT_DIR" || {
        echo "Error: Failed to create output directory" >&2
        exit 1
    }
}

# Logging system
log() {
    local level=$1
    local message=$2
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${level^^}] $message"
}

# Core processing functions
process_individual_beds() {
    local -a processed_files
    log INFO "Processing individual BED files"

    for bed_file in "${INPUT_DIR}"/*.bed; do
        [[ -f "$bed_file" ]] || continue
        [[ "$bed_file" == *_keeplast.bed ]] && continue

        local base_name=$(basename "${bed_file%.*}")
        local output_file="${OUTPUT_DIR}/${base_name}_selfcluster_keeplast.bed"

        log DEBUG "Processing: $bed_file"
        
        {
            sort -k1,1 -k2,2n "$bed_file" | \
            bedtools merge -i stdin -d "$DEFAULT_MERGE_DISTANCE" -s -c 4,5,6 -o distinct | \
            awk 'BEGIN {FS=OFS="\t"} {
                if ($6 == "+") $2 = $3
                else if ($6 == "-") $3 = $2
                print
            }'
        } > "$output_file"

        processed_files+=("$output_file")
        TEMP_FILES+=("$output_file")
    done

    if [[ ${#processed_files[@]} -eq 0 ]]; then
        log ERROR "No valid BED files found in input directory"
        exit 1
    fi
}

merge_group_a() {
    log INFO "Generating groupA.bed"
    local group_a="${OUTPUT_DIR}/groupA.bed"
    cat "${OUTPUT_DIR}"/*_selfcluster_keeplast.bed > "$group_a"
    TEMP_FILES+=("$group_a")
}

filter_group_a() {
    log INFO "Filtering groupA.bed"
    local input="${OUTPUT_DIR}/groupA.bed"
    local output="${OUTPUT_DIR}/groupA_nonunique_keeplast.bed"

    bedtools sort -i "$input" | \
    bedtools merge -i stdin -d "$DEFAULT_MERGE_DISTANCE" -s -c 4,5,6 -o count,distinct,distinct | \
    awk -v min="$MIN_REPLICATE_COUNT" 'BEGIN {FS=OFS="\t"} $4 >= min' | \
    awk 'BEGIN {FS=OFS="\t"} {
        if ($6 == "+") $2 = $3
        else if ($6 == "-") $3 = $2
        print
    }' > "$output"

    TEMP_FILES+=("$output")
}

generate_final_output() {
    log INFO "Creating final output"
    local base_name=$(basename "$INPUT_DIR")
    local final_output="${OUTPUT_DIR}/${base_name}_final_site.bed"
    
    cp "${OUTPUT_DIR}/groupA_nonunique_keeplast.bed" "$final_output"
    log SUCCESS "Final output created: $final_output"
}

# Main workflow
main() {
    parse_arguments "$@"
    validate_inputs
    log INFO "Starting BED processing pipeline"
    log INFO "System Configuration:"
    log INFO "  Input Directory: $INPUT_DIR"
    log INFO "  Output Directory: $OUTPUT_DIR"
    log INFO "  Merge Distance: $DEFAULT_MERGE_DISTANCE"

    process_individual_beds
    merge_group_a
    filter_group_a
    generate_final_output

    log INFO "Pipeline completed successfully"
}

main "$@"
