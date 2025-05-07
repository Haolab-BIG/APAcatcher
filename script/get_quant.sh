#!/bin/bash
# Author: [ChengPeng]
# Set script to exit on error
set -e

# Parse command line arguments with option flags
usage() {
    echo "Usage: $0 -i <index_path> -d <fastq_directory> -o <output_directory>"
    exit 1
}

# Initialize variables
index=""
fastq_dir=""
output_dir=""
threads=12  # Default number of threads

# Parse command line options
while getopts ":i:d:o:" opt; do
    case $opt in
        i) index="$OPTARG" ;;
        d) fastq_dir="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        \?) echo "Invalid option -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Validate required parameters
if [[ -z "$index" || -z "$fastq_dir" || -z "$output_dir" ]]; then
    echo "Missing required parameters"
    usage
fi

# Check if output directory exists, create if not
if [[ ! -d "$output_dir" ]]; then
    echo "Output directory does not exist. Creating directory: $output_dir"
    mkdir -p "$output_dir"
fi

# Supported FASTQ file extensions
extensions=("fastq.gz" "fq.gz" "fastq" "fq")

# Process paired-end FASTQ files
for ext in "${extensions[@]}"; do
    for file1 in "$fastq_dir"/*_1*."$ext"; do
        # Skip if file does not exist
        if [[ ! -e "$file1" ]]; then
            continue
        fi

        # Extract sample name by removing _1 and extension
        sample_name=$(basename "$file1")
        sample_name="${sample_name%%_1*}"
        echo "$sample_name"

        # Find corresponding R2 file using wildcard match
        file2=$(find "$fastq_dir" -maxdepth 1 -type f -name "${sample_name}_2*.$ext")
        
        if [[ -z "$file2" ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Warning: Pair file for $file1 not found!"
            continue
        fi

        echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting processing for sample: $sample_name"

        # Run Salmon quantification
        salmon quant -i "$index" -l IU -1 "$file1" -2 "$file2" \
        -o "${output_dir}/${sample_name}" -p "$threads" --validateMappings --useEM

        if [[ $? -eq 0 ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Completed processing for sample: $sample_name"
        else
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Error occurred while processing sample: $sample_name"
        fi
    done
done

echo "$(date '+%Y-%m-%d %H:%M:%S') - All samples processed."
