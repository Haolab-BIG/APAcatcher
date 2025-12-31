#!/usr/bin/env bash

set -euo pipefail

INPUT_DIR=""
DATA_DIR=""
OUT_DIR=""

while getopts "i:d:o:" opt; do
  case $opt in
    i) INPUT_DIR="$OPTARG" ;;
    d) DATA_DIR="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    *)
      echo "Usage: $0 -i <index_dir> -d <data_dir> -o <output_dir>"
      exit 1
      ;;
  esac
done

if [[ -z "$INPUT_DIR" || -z "$DATA_DIR" || -z "$OUT_DIR" ]]; then
  echo "Error: -i -d -o are required."
  echo "Usage: $0 -i <index_dir> -d <data_dir> -o <output_dir>"
  exit 1
fi

mkdir -p "$OUT_DIR"

echo "===== Running quantification ====="
echo "Index directory: $INPUT_DIR"
echo "Data directory : $DATA_DIR"
echo "Output dir     : $OUT_DIR"

for fq in "$DATA_DIR"/*.fastq.gz "$DATA_DIR"/*.fq.gz "$DATA_DIR"/*.fastq "$DATA_DIR"/*.fq; do
    [[ -e "$fq" ]] || continue

    base=$(basename "$fq")
    sample="${base%%.*}"

    echo "Processing $sample ..."
    salmon quant \
        -i "$INPUT_DIR" \
        -l A \
        -r "$fq" \
        -o "$OUT_DIR/$sample" \
        -p 12 \
	--validateMappings \
	--useEM
done

echo "===== Done ====="
