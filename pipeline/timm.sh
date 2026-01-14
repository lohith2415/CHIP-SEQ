#!/bin/bash
set -euo pipefail

############################################
# Detect paths
############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

CHIP_DIR="$PROJECT_ROOT/chip"
QC_DIR="$PROJECT_ROOT/qc"
TRIM_DIR="$PROJECT_ROOT/trimmed"

TRIMMOMATIC_JAR="$SCRIPT_DIR/Trimmomatic-0.39/trimmomatic-0.39.jar"

THREADS=${THREADS:-4}
MINLEN=30

mkdir -p "$QC_DIR" "$TRIM_DIR"

############################################
# Check Trimmomatic
############################################
if [[ ! -f "$TRIMMOMATIC_JAR" ]]; then
    echo "❌ Trimmomatic JAR not found:"
    echo "   $TRIMMOMATIC_JAR"
    exit 1
fi

############################################
# Find FASTQ files inside chip/
############################################
echo "🔍 Searching FASTQ files inside chip/ ..."
FASTQ_FILES=$(find "$CHIP_DIR" -type f -name "*.fastq.gz")

if [[ -z "$FASTQ_FILES" ]]; then
    echo "❌ No FASTQ files found in chip/"
    exit 1
fi

############################################
# FASTQC (run only if missing)
############################################
RUN_FASTQC=false
for fq in $FASTQ_FILES; do
    base=$(basename "$fq" .fastq.gz)
    if [[ ! -f "$QC_DIR/${base}_fastqc.html" ]]; then
        RUN_FASTQC=true
        break
    fi
done

if $RUN_FASTQC; then
    echo "📊 Running FastQC..."
    fastqc -t "$THREADS" $FASTQ_FILES -o "$QC_DIR"
    echo "✅ FastQC completed"
else
    echo "⏭ FastQC already exists — skipping"
fi

############################################
# Ask before trimming
############################################
echo
echo "👉 Press ENTER to start trimming (existing trimmed files will be skipped)"
echo "👉 Press Ctrl+C to cancel"
read

############################################
# TRIMMING (SKIP IF TRIM EXISTS)
############################################
echo "✂️ Starting trimming..."

for fq in $FASTQ_FILES; do
    # Path relative to chip/
    rel_path="${fq#$CHIP_DIR/}"

    # Mirror directory structure inside trimmed/
    out_dir="$TRIM_DIR/$(dirname "$rel_path")"
    mkdir -p "$out_dir"

    base=$(basename "$fq" .fastq.gz)
    out="$out_dir/${base}_trimmed.fastq.gz"

    # ✅ SKIP if trimmed file already exists
    if [[ -f "$out" ]]; then
        echo "⏭ Trimmed file exists — skipping:"
        echo "   $out"
        continue
    fi

    echo "🔹 Trimming: $rel_path"

    java -jar "$TRIMMOMATIC_JAR" SE \
        -threads "$THREADS" \
        -phred33 \
        "$fq" \
        "$out" \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:4:20 \
        MINLEN:$MINLEN
done

echo "🎉 Trimming completed (structure preserved, existing files skipped)"

