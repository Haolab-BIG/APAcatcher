#!/bin/bash

# ====================================================
# Script Name: cluster_bed_files.sh
# Description : Batch process multiple .bed files: sort, merge, filter,
#               generate a consolidated final_site.bed file, and clean up intermediate files.
# Usage       : ./cluster_bed_files.sh <BED_FILES_DIR> <OUTPUT_DIR>
# Example     : ./cluster_bed_files.sh /data/bed_files /data/results
# ====================================================

# Exit immediately on error
set -e

# Check for required input arguments
if [ $# -lt 2 ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <BED_FILES_DIRECTORY> <OUTPUT_DIRECTORY>"
    echo "       <BED_FILES_DIRECTORY>: Path to folder containing .bed files"
    echo "       <OUTPUT_DIRECTORY>   : Path to folder where results will be saved"
    exit 1
fi

# Assign input parameters\ nBED_DIR="$1"
OUTPUT_DIR="$2"

# Change to BED files directory
cd "$BED_DIR" || { echo "Error: Cannot change to directory '$BED_DIR'"; exit 1; }

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "Starting batch processing of .bed files in: $BED_DIR"

# Array to hold intermediate "keeplast" file paths
keeplast_files=()

echo "Step 1: Processing individual .bed files..."
for bed_file in *.bed; do
    # Skip already-processed files
    if [[ "$bed_file" == *_selfcluster_keeplast.bed ]]; then
        continue
    fi

    # Verify file exists
    if [[ ! -f "$bed_file" ]]; then
        echo "Warning: No .bed files found in directory."
        break
    fi

    # Derive sample prefix from filename (e.g., "sample1.bed" -> "sample1")
    sample_prefix="${bed_file%.*}"
    echo "Processing: $bed_file"

    # Define output filename
    keeplast_file="$OUTPUT_DIR/${sample_prefix}_selfcluster_keeplast.bed"

    # Sort, merge, and retain last base position per strand
    sort -k1,1 -k2,2n "$bed_file" | \
    bedtools merge -i stdin -d 70 -s -c 4,5,6 -o distinct | \
    awk 'BEGIN {FS=OFS="\t"} { if ($6 == "+") $2=$3; else if ($6 == "-") $3=$2; print }' > "$keeplast_file"

    echo "Generated: $keeplast_file"
    keeplast_files+=("$keeplast_file")
    echo "----------------------------------------"
done

# Ensure at least one intermediate file was created
if [ ${#keeplast_files[@]} -eq 0 ]; then
    echo "Error: No intermediate files generated. Exiting."
    exit 1
fi

echo "Step 2: Consolidating all keeplast files into groupA.bed..."
cat "${keeplast_files[@]}" > "$OUTPUT_DIR/groupA.bed"
echo "Generated: $OUTPUT_DIR/groupA.bed"

echo "Step 3: Sorting, merging, and filtering non-unique sites..."
bedtools sort -i "$OUTPUT_DIR/groupA.bed" | \
bedtools merge -i stdin -d 70 -s -c 4,5,6 -o count,distinct,distinct | \
awk '\$4 > 1' | \
awk 'BEGIN {FS=OFS="\t"} { if ($6 == "+") $2=$3; else if ($6 == "-") $3=$2; print }' > "$OUTPUT_DIR/groupA_nonunique_keeplast.bed"
echo "Generated: $OUTPUT_DIR/groupA_nonunique_keeplast.bed"

echo "Step 4: Creating final output file..."
base_name=$(basename "$BED_DIR")
final_output="$OUTPUT_DIR/${base_name}_final_site.bed"
cp "$OUTPUT_DIR/groupA_nonunique_keeplast.bed" "$final_output"
echo "Generated: $final_output"

echo "Step 5: Cleaning up intermediate files..."
rm -f "$OUTPUT_DIR/groupA.bed" "$OUTPUT_DIR/groupA_nonunique_keeplast.bed"
for file in "${keeplast_files[@]}"; do
    [[ -f "$file" ]] && rm -f "$file"
done

echo "Processing complete. Final BED file available at: $final_output"

