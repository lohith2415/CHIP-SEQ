#!/bin/bash
set -euo pipefail

############################################
# Detect paths
############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

BAM_DIR="$PROJECT_ROOT/bam"
CLEAN_BAM_DIR="$PROJECT_ROOT/bam_clean"

THREADS=${THREADS:-4}
MAPQ=30

mkdir -p "$CLEAN_BAM_DIR"

############################################
# Sanity check: required tools
############################################
for tool in samtools; do
    if ! command -v "$tool" &>/dev/null; then
        echo "❌ Required tool not found: $tool"
        echo "👉 Please install it before running this script"
        exit 1
    fi
done

############################################
# Find BAM files
############################################
echo "🔍 Searching BAM files..."
BAM_FILES=$(find "$BAM_DIR" -type f -name "*.sorted.bam")

if [[ -z "$BAM_FILES" ]]; then
    echo "❌ No BAM files found in $BAM_DIR"
    exit 1
fi

############################################
# Filter + Deduplicate
############################################
echo "🧬 Filtering & removing PCR duplicates..."

for bam in $BAM_FILES; do
    # Preserve directory structure
    rel_path="${bam#$BAM_DIR/}"
    out_dir="$CLEAN_BAM_DIR/$(dirname "$rel_path")"
    mkdir -p "$out_dir"

    base=$(basename "$bam" .sorted.bam)
    filt_bam="$out_dir/${base}.q${MAPQ}.bam"
    dedup_bam="$out_dir/${base}.dedup.bam"

    # STEP 1: MAPQ filtering
    if [[ ! -f "$filt_bam" ]]; then
        echo "🔹 Filtering MAPQ ≥ $MAPQ : $rel_path"
        samtools view -@ "$THREADS" -b -q "$MAPQ" "$bam" > "$filt_bam"
        samtools index "$filt_bam"
    else
        echo "⏭ MAPQ-filtered BAM exists — skipping"
    fi

    # STEP 2: Duplicate removal
    if [[ ! -f "$dedup_bam" ]]; then
        echo "🔹 Removing PCR duplicates..."

        samtools sort -@ "$THREADS" -n \
            -o "${filt_bam%.bam}.name_sorted.bam" "$filt_bam"

        samtools fixmate -m \
            "${filt_bam%.bam}.name_sorted.bam" \
            "${filt_bam%.bam}.fixmate.bam"

        samtools sort -@ "$THREADS" \
            -o "${filt_bam%.bam}.pos_sorted.bam" \
            "${filt_bam%.bam}.fixmate.bam"

        samtools markdup -r \
            "${filt_bam%.bam}.pos_sorted.bam" "$dedup_bam"

        samtools index "$dedup_bam"

        # Cleanup intermediate files
        rm -f "${filt_bam%.bam}".{name_sorted.bam,fixmate.bam,pos_sorted.bam}
    else
        echo "⏭ Deduplicated BAM exists — skipping"
    fi
done

echo "🎉 Filtering + duplicate removal completed!"
echo "📁 Clean BAM files are in: $CLEAN_BAM_DIR"

