#!/bin/bash

# ====================================================
# Author: [ChengPeng]
# SCRIPT: get_salmon_index.sh
# DESCRIPTION: Process 3'UTR isoforms and build Salmon index
# USAGE: ./get_salmon_index.sh [OPTIONS] REQUIRED_ARGS
# ====================================================

set -e

# Parse command line arguments with option flags

usage() {
    echo "Usage: $0 -f <final_site_bed> -r <refseq_utr_bed> -l <refseq_last_bed> -g <hg38_fa> -o <output_fa> -i <output_index>"
    exit 1
}

while getopts ":f:r:l:g:o:i:" opt; do
    case $opt in
        f) final_site_bed="$OPTARG" ;;
        r) refseq_utr_bed="$OPTARG" ;;
        l) refseq_last_bed="$OPTARG" ;;
        g) hg38_fa="$OPTARG" ;;
        o) output_fa="$OPTARG" ;;
        i) output_index="$OPTARG" ;;
        \?) echo "Invalid option -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Validate required parameters
if [ -z "$final_site_bed" ] || [ -z "$refseq_utr_bed" ] || [ -z "$refseq_last_bed" ] || [ -z "$hg38_fa" ] || [ -z "$output_fa" ] || [ -z "$output_index" ]; then
    echo "Missing required parameters"
    usage
fi

# Generate final_site_with_id.bed with incremental IDs
echo "Generating final_site_with_id.bed with unique incremental IDs..."
awk 'BEGIN {OFS="\t"} 
     { 
         $4 = NR;  # Replace fourth column with row number (1, 2, 3, ...)
         print 
     }' "$final_site_bed" > final_site_with_id.bed

# Generate final_site_7col.bed
echo "Generating final_site_7col.bed..."
awk -F'\t' 'BEGIN {OFS="\t"} {
    if ($6 == "+") $7 = $3; 
    else if ($6 == "-") $7 = $2; 
    print $0
}' final_site_with_id.bed > final_site_7col.bed

# Generate RefSeq_UTR_final_7col.bed
echo "Generating RefSeq_UTR_final_7col.bed..."
awk -F'\t' 'BEGIN {OFS="\t"} {
    if ($6 == "+") $7 = $2; 
    else if ($6 == "-") $7 = $3; 
    print $0
}' "$refseq_utr_bed" > RefSeq_UTR_final_7col.bed

# Generate 3UTRisoforms_temp.bed
echo "Generating 3UTRisoforms_temp.bed..."
awk 'BEGIN {FS=OFS="\t"} 
    NR==FNR { ref[$4] = $7; next } 
    $5 in ref { $8 = ref[$5] } 
    { print $0}' RefSeq_UTR_final_7col.bed final_site_7col.bed | \
awk 'BEGIN {FS=OFS="\t"} {
    if ($7 < $8) { 
        $2 = $7; 
        $3 = $8; 
    } else { 
        $2 = $8; 
        $3 = $7; 
    } 
    $7 = ""; 
    $8 = ""; 
    sub(/\t+$/, ""); 
    print 
}' | \
awk 'BEGIN {FS=OFS="\t"} {
    $4 = $5; 
    $7 = $5; 
    $5 = ($2 > $3) ? ($2 - $3) : ($3 - $2); 
    print 
}' > 3UTRisoforms_temp.bed

# Filter entries with length ≥ 100 and generate 3UTRisoforms.bed
echo "Filtering 3UTRisoforms_temp.bed for lengths >= 100 and generating 3UTRisoforms.bed..."
awk 'BEGIN {FS=OFS="\t"} {
    id = NR      # Generate unique incremental ID
    $4 = id      # Update fourth column to unique ID
    $7 = ""      # Clear seventh column
    NF = 6       # Ensure only six columns
    if ($5 >= 100) 
        print
}' 3UTRisoforms_temp.bed > 3UTRisoforms.bed

# Clean up temporary files
rm 3UTRisoforms_temp.bed

# Match 3UTRisoforms.bed with final_site_bed using unique IDs
echo "Matching 3UTRisoforms.bed with final_site_with_id.bed based on unique ID..."
awk 'BEGIN {FS=OFS="\t"}
    NR==FNR { id[$4] = $0; next }
    $4 in id { 
        split(id[$4], a, "\t") 
        print a[1], a[2], a[3], a[4], a[5], a[6] 
    }' final_site_7col.bed 3UTRisoforms.bed > matched_final_site.bed
rm 3UTRisoforms.bed
rm final_site_7col.bed
rm RefSeq_UTR_final_7col.bed

# Generate Last_Exon_3UTRisoforms.bed using refseq_last_bed
echo "Generating final_site_7col.bed..."
awk -F'\t' 'BEGIN {OFS="\t"} {
    if ($6 == "+") $7 = $3; 
    else if ($6 == "-") $7 = $2; 
    print $0
}' matched_final_site.bed > final_site_7col.bed

# Clean RefSeq_UTR_final_7col.bed and check for valid data
echo "Cleaning RefSeq_UTR_final_7col.bed..."
awk -F'\t' 'BEGIN {OFS="\t"} {
    if (NF < 6 || $6 == "" || $2 == "" || $3 == "") next;
    if ($6 == "+") $7 = $2;
    else if ($6 == "-") $7 = $3;
    print $0;
}' "$refseq_last_bed" > RefSeq_UTR_final_7col.bed

# Generate 3UTRisoforms_temp.bed with boundary handling
echo "Generating Last 3UTRisoforms_temp.bed..."
awk 'BEGIN {FS=OFS="\t"}
    NR==FNR {ref[$4] = $7; next}
    $5 in ref {
        $8 = ref[$5];
    }
    { print $0 }
' RefSeq_UTR_final_7col.bed final_site_7col.bed | \
awk 'BEGIN {FS=OFS="\t"} {
    if ($7 == "" || $8 == "") next;
    if ($7 < $8) {
        $2 = $7;
        $3 = $8;
    } else {
        $2 = $8;
        $3 = $7;
    }
    $7 = "";
    $8 = "";
    sub(/\t+$/, ""); 
    print;
}' | \
awk 'BEGIN {FS=OFS="\t"} {
    if ($2 == "" || $3 == "") next; 
    $4 = $5;
    $7 = $5; 
    $5 = ($2 > $3) ? ($2 - $3) : ($3 - $2); 
    print;
}' > 3UTRisoforms_temp.bed

# Fix formatting issues in 3UTRisoforms.bed
echo "Fixing 3UTRisoforms.bed format..."
awk 'BEGIN {FS=OFS="\t"} {
    if (NF < 6) next; 
    $4 = $7":"$6; 
    $7 = "";       
    NF = 6;        
    print;
}' 3UTRisoforms_temp.bed > 3UTRisoforms.bed

sort -k1,1 -k2,2n 3UTRisoforms.bed > 3UTRisoforms_sorted.bed
mv 3UTRisoforms_sorted.bed 3UTRisoforms.bed

# Extract 3'UTR sequences and generate FASTA file
echo "Extracting sequences to $output_fa..."
bedtools getfasta -fi "$hg38_fa" -fo "$output_fa" -bed 3UTRisoforms.bed -name

# Build Salmon index
echo "Building Salmon index at $output_index..."
salmon index -t "$output_fa" -k 31 -i "$output_index" -p 8

# Clean up temporary files
rm final_site_with_id.bed
rm final_site_7col.bed
rm RefSeq_UTR_final_7col.bed
rm matched_final_site.bed
rm 3UTRisoforms_temp.bed
echo "Pipeline complete."
