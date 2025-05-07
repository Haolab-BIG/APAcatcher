#!/bin/bash

# ====================================================
# SCRIPT: salmon_quantification.sh
# DESCRIPTION: Perform RNA-seq quantification using Salmon
# USAGE: ./salmon_quantification.sh -i INDEX -f FASTQ_DIR -o OUTPUT_DIR [OPTIONS]
# ====================================================

set -euo pipefail

# Initialize configuration variables
INDEX_PATH=""
FASTQ_DIR=""
OUTPUT_DIR=""
THREADS=12
OVERWRITE=false
LOG_FILE="salmon_processing.log"
EXTENSIONS=("fastq.gz" "fq.gz" "fastq" "fq")
SCRIPT_NAME=$(basename "$0")

# Display help information
usage() {
    cat <<EOF
${SCRIPT_NAME} - Salmon quantification pipeline

USAGE:
    ${SCRIPT_NAME} -i INDEX -f FASTQ_DIR -o OUTPUT_DIR [OPTIONS]

REQUIRED PARAMETERS:
    -i, --index          Path to Salmon index
    -f, --fastq-dir      Directory containing FASTQ files
    -o, --output-dir     Output directory for quantification results

OPTIONS:
    -t, --threads        Number of threads to use (default: 12)
    -e, --extensions     Comma-separated file extensions to process
                        (default: fastq.gz,fq.gz,fastq,fq)
    -l, --log-file       Path to log file (default: salmon_processing.log)
    --overwrite          Overwrite existing output directories
    -h, --help           Show this help message

EXAMPLES:
    Basic usage:
    ${SCRIPT_NAME} -i /data/index -f /data/fastq -o /results

    Custom extensions and threads:
    ${SCRIPT_NAME} -i /data/index -f /data/fastq -o /results \\
                   -t 16 -e "fastq.gz,fq"
EOF
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--index)
            INDEX_PATH="$2"
            shift
            ;;
        -f|--fastq-dir)
            FASTQ_DIR="$2"
            shift
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift
            ;;
        -t|--threads)
            THREADS="$2"
            shift
            ;;
        -e|--extensions)
            IFS=',' read -ra EXTENSIONS <<< "$2"
            shift
            ;;
        -l|--log-file)
            LOG_FILE="$2"
            shift
            ;;
        --overwrite)
            OVERWRITE=true
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
check_path() {
    if [[ -z "$1" ]]; then
        echo "Error: Missing required parameter: $2" >&2
        usage >&2
        exit 1
    fi
    if [[ ! -e "$1" ]]; then
        echo "Error: Path does not exist: $1" >&2
        exit 1
    fi
}

check_path "$INDEX_PATH" "index"
check_path "$FASTQ_DIR" "fastq-dir"

# Initialize output directory and logging
mkdir -p "$OUTPUT_DIR" || {
    echo "Error: Failed to create output directory" >&2
    exit 1
}

exec > >(tee -a "$LOG_FILE") 2>&1

log() {
    local level=$1
    local message=$2
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${level^^}] $message"
}

# File pair validation functions
validate_file_pair() {
    local ext=$1
    local count=0
    
    while IFS= read -r -d $'\0' file; do
        ((count++))
        sample_name=$(basename "$file")
        sample_name="${sample_name%%_1*}"
        
        local pair_pattern="${sample_name}_2${file#*_1}"
        local pair_files=()
        
        while IFS= read -r -d $'\0' pf; do
            pair_files+=("$pf")
        done < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "$pair_pattern" -print0)

        case ${#pair_files[@]} in
            0)
                log WARNING "Missing pair for $file"
                ;;
            1)
                process_sample "$file" "${pair_files[0]}" "$sample_name"
                ;;
            *)
                log ERROR "Multiple pair candidates found for $file"
                ;;
        esac
    done < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "*_1*.$ext" -print0)
}

# Main processing function
process_sample() {
    local r1=$1
    local r2=$2
    local sample=$3
    local out_dir="${OUTPUT_DIR}/${sample}"

    if [[ -d "$out_dir" && "$OVERWRITE" == false ]]; then
        log WARNING "Skipping existing sample: $sample (use --overwrite to force)"
        return
    fi

    log INFO "Processing sample: $sample"
    log DEBUG "R1: $r1 | R2: $r2"

    mkdir -p "$out_dir"
    
    if salmon quant -i "$INDEX_PATH" -l IU \
        -1 "$r1" -2 "$r2" \
        -o "$out_dir" -p "$THREADS" \
        --validateMappings --useEM
    then
        log INFO "Completed successfully: $sample"
    else
        log ERROR "Failed processing: $sample"
        rm -rf "$out_dir"
    fi
}

# Main execution flow
log INFO "Starting Salmon quantification pipeline"
log INFO "System configuration:"
log INFO "  - CPU threads: $THREADS"
log INFO "  - Input directory: $FASTQ_DIR"
log INFO "  - Output directory: $OUTPUT_DIR"

for ext in "${EXTENSIONS[@]}"; do
    log INFO "Processing files with extension: $ext"
    validate_file_pair "$ext"
done

log INFO "Processing completed successfully"
