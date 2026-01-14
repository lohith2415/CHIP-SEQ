#!/bin/bash
set -euo pipefail

############################################
# Detect paths (BASED ON YOUR LAYOUT)
############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Raw FASTQ location (your current folder)
RAW_CHIP_DIR="$(dirname "$SCRIPT_DIR")/chip"

# Target folder (same chip/, but organized)
CHIP_DIR="$RAW_CHIP_DIR"

############################################
# Sanity checks
############################################
if [[ ! -d "$RAW_CHIP_DIR" ]]; then
    echo "❌ chip directory not found:"
    echo "   $RAW_CHIP_DIR"
    exit 1
fi

FASTQ_FILES=$(find "$RAW_CHIP_DIR" -maxdepth 1 -type f -name "*.fastq.gz")

if [[ -z "$FASTQ_FILES" ]]; then
    echo "❌ No FASTQ files found in chip/"
    exit 1
fi

echo "🔍 Found FASTQ files:"
ls "$RAW_CHIP_DIR"/*.fastq.gz
echo

############################################
# Arrange FASTQ files
############################################
echo "📦 Arranging FASTQ files into chip/<condition>/ ..."

for fq in $FASTQ_FILES; do
    fname=$(basename "$fq")

    # Detect replicate
    rep=$(echo "$fname" | grep -oE "rep[0-9]+" || true)
    if [[ -z "$rep" ]]; then
        echo "⚠ Cannot detect replicate in $fname — skipping"
        continue
    fi

    # INPUT files
    if echo "$fname" | grep -qi "input"; then
        condition=$(echo "$fname" | sed -E 's/_?input.*//')
        out_dir="$CHIP_DIR/${condition}_input"
        out_file="${condition}_input_${rep}.fastq.gz"

    # CHIP files
    else
        condition=$(echo "$fname" | sed -E 's/_?(chip_)?rep[0-9]+.*//')
        out_dir="$CHIP_DIR/$condition"
        out_file="${condition}_${rep}.fastq.gz"
    fi

    mkdir -p "$out_dir"

    # Move instead of copy (cleaner)
    echo "➡ $fname → $out_dir/$out_file"
    mv "$fq" "$out_dir/$out_file"
done

echo
echo "🎉 FASTQ arrangement completed!"
echo
echo "📁 Final chip/ structure:"
find "$CHIP_DIR" -maxdepth 2 -type f -name "*.fastq.gz" | sort

