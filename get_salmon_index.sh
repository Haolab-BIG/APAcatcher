#!/bin/bash

# ====================================================
# SCRIPT: utr_isoform_processing.sh
# DESCRIPTION: Process 3'UTR isoforms and build Salmon index
# USAGE: ./utr_isoform_processing.sh [OPTIONS] REQUIRED_ARGS
# ====================================================

set -euo pipefail

# Configuration
THREADS=8
KMER_SIZE=31
TEMP_FILES=()

# Initialize parameters
FINAL_SITE_BED=""
REFSEQ_UTR_BED=""
REFSEQ_LAST_BED=""
HG38_FA=""
OUTPUT_FA=""
OUTPUT_INDEX=""
KEEP_TEMP=false

SCRIPT_NAME=$(basename "$0")
VERSION="1.0.0"

# Cleanup function
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
show_help() {
    cat <<EOF
${SCRIPT_NAME} v${VERSION} - 3'UTR Isoform Processing Pipeline

USAGE:
    ${SCRIPT_NAME} [OPTIONS] -f FINAL_BED -u UTR_BED -l LAST_EXON_BED -g GENOME_FA -o OUTPUT_FA -x INDEX_DIR

REQUIRED PARAMETERS:
    -f, --final-bed      Final site BED file
    -u, --utr-bed        RefSeq UTR BED file
    -l, --last-exon      RefSeq last exon BED file
    -g, --genome         hg38 reference genome FASTA
    -o, --output-fa      Output FASTA file
    -x, --index-dir      Salmon index directory

OPTIONS:
    -t, --threads        Processing threads (default: ${THREADS})
    -k, --kmer           Salmon k-mer size (default: ${KMER_SIZE})
    --keep-temp          Keep temporary files
    -h, --help           Show this help message
    -v, --version        Display version

EXAMPLES:
    Basic usage:
    ${SCRIPT_NAME} -f final_site.bed -u RefSeq_UTR.bed \\
                   -l refseq_last_exon.bed -g hg38.fa \\
                   -o isoforms.fa -x salmon_index

    Custom parameters:
    ${SCRIPT_NAME} -f input.bed -u utr.bed -l last_exon.bed \\
                   -g genome.fa -o output.fa -x index_dir \\
                   -t 16 -k 27 --keep-temp
EOF
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f|--final-bed)
            FINAL_SITE_BED="$2"
            shift
            ;;
        -u|--utr-bed)
            REFSEQ_UTR_BED="$2"
            shift
            ;;
        -l|--last-exon)
            REFSEQ_LAST_BED="$2"
            shift
            ;;
        -g|--genome)
            HG38_FA="$2"
            shift
            ;;
        -o|--output-fa)
            OUTPUT_FA="$2"
            shift
            ;;
        -x|--index-dir)
            OUTPUT_INDEX="$2"
            shift
            ;;
        -t|--threads)
            THREADS="$2"
            shift
            ;;
        -k|--kmer)
            KMER_SIZE="$2"
            shift
            ;;
        --keep-temp)
            KEEP_TEMP=true
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        -v|--version)
            echo "${SCRIPT_NAME} v${VERSION}"
            exit 0
            ;;
        *)
            echo "Error: Unknown parameter: $1" >&2
            show_help >&2
            exit 1
            ;;
    esac
    shift
done

# Validate required parameters
validate_file() {
    if [[ ! -f "$1" ]]; then
        echo "Error: Missing required file: $1" >&2
        exit 1
    fi
}

check_params() {
    local -a required=(
        "$FINAL_SITE_BED" 
        "$REFSEQ_UTR_BED" 
        "$REFSEQ_LAST_BED" 
        "$HG38_FA" 
        "$OUTPUT_FA" 
        "$OUTPUT_INDEX"
    )
    
    for param in "${required[@]}"; do
        if [[ -z "$param" ]]; then
            echo "Error: Missing required parameters" >&2
            show_help >&2
            exit 1
        fi
    done

    validate_file "$FINAL_SITE_BED"
    validate_file "$REFSEQ_UTR_BED"
    validate_file "$REFSEQ_LAST_BED"
    validate_file "$HG38_FA"
}

check_params

# Logging functions
log() {
    local level=$1
    local message=$2
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${level^^}] $message"
}

add_temp_file() {
    TEMP_FILES+=("$1")
}

# Core processing functions
process_final_site() {
    log INFO "Generating final_site_with_id.bed"
    awk 'BEGIN {OFS="\t"} { $4 = NR; print }' "$FINAL_SITE_BED" > final_site_with_id.bed
    add_temp_file "final_site_with_id.bed"

    log INFO "Creating 7-column final site BED"
    awk -F'\t' 'BEGIN {OFS="\t"} {
        if ($6 == "+") $7 = $3; 
        else if ($6 == "-") $7 = $2; 
        print
    }' final_site_with_id.bed > final_site_7col.bed
    add_temp_file "final_site_7col.bed"
}

process_utr_data() {
    log INFO "Processing UTR data"
    awk -F'\t' 'BEGIN {OFS="\t"} {
        if ($6 == "+") $7 = $2;
        else if ($6 == "-") $7 = $3;
        print
    }' "$REFSEQ_UTR_BED" > RefSeq_UTR_final_7col.bed
    add_temp_file "RefSeq_UTR_final_7col.bed"
}

generate_isoforms() {
    log INFO "Generating 3UTR isoforms"
    awk 'BEGIN {FS=OFS="\t"} NR==FNR { ref[$4] = $7; next } $5 in ref { $8 = ref[$5] } { print }' \
        RefSeq_UTR_final_7col.bed final_site_7col.bed | \
    awk 'BEGIN {FS=OFS="\t"} {
        if ($7 < $8) { $2 = $7; $3 = $8 } 
        else { $2 = $8; $3 = $7 }
        $7 = $8 = ""; sub(/\t+$/, "")
        print
    }' | \
    awk 'BEGIN {FS=OFS="\t"} {
        $4 = $5; $7 = $5
        $5 = ($2 > $3) ? ($2 - $3) : ($3 - $2)
        print
    }' > 3UTRisoforms_temp.bed
    add_temp_file "3UTRisoforms_temp.bed"

    log INFO "Filtering isoforms"
    awk 'BEGIN {FS=OFS="\t"} {
        if ($5 >= 100) {
            $4 = NR; $7 = ""; NF = 6
            print
        }
    }' 3UTRisoforms_temp.bed > 3UTRisoforms.bed
    add_temp_file "3UTRisoforms.bed"
}

build_index() {
    log INFO "Extracting sequences"
    bedtools getfasta -fi "$HG38_FA" -fo "$OUTPUT_FA" -bed 3UTRisoforms.bed -name

    log INFO "Building Salmon index"
    salmon index -t "$OUTPUT_FA" -k "$KMER_SIZE" -i "$OUTPUT_INDEX" -p "$THREADS"
}

# Main execution flow
main() {
    log INFO "Starting 3'UTR processing pipeline"
    log INFO "System configuration:"
    log INFO "  - Threads: $THREADS"
    log INFO "  - k-mer size: $KMER_SIZE"
    log INFO "  - Genome: $HG38_FA"

    process_final_site
    process_utr_data
    generate_isoforms
    build_index

    log INFO "Pipeline completed successfully"
}

main
