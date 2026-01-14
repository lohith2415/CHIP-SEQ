#!/bin/bash
set -euo pipefail

############################################
# Detect paths
############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

############################################
# Auto-activate MACS2 environment
############################################
MACS2_ENV="$SCRIPT_DIR/tools/macs2_env"

if [[ -f "$MACS2_ENV/bin/activate" ]]; then
    source "$MACS2_ENV/bin/activate"
else
    echo "❌ MACS2 environment not found:"
    echo "   $MACS2_ENV"
    exit 1
fi

if ! command -v macs2 &>/dev/null; then
    echo "❌ macs2 not found even after activation"
    exit 1
fi

############################################
# Paths and parameters
############################################
BAM_CLEAN_DIR="$PROJECT_ROOT/bam_clean"
PEAK_DIR="$PROJECT_ROOT/peaks"

GENOME_SIZE=1.2e8        # Arabidopsis thaliana
QVALUE=0.01
EXTSIZE=200              # TF ChIP-seq best practice
THREADS=${THREADS:-4}

mkdir -p "$PEAK_DIR"

############################################
# Peak calling
############################################
echo "🧬 Starting MACS2 peak calling..."

for condition in cdf2 cdf2_pif4; do
    CHIP_DIR="$BAM_CLEAN_DIR/rep/$condition"
    INPUT_DIR="$BAM_CLEAN_DIR/input/$condition"

    if [[ ! -d "$CHIP_DIR" || ! -d "$INPUT_DIR" ]]; then
        echo "⏭ Skipping condition $condition (missing BAM directories)"
        continue
    fi

    for chip_bam in "$CHIP_DIR"/*.dedup.bam; do
        chip_base=$(basename "$chip_bam" .dedup.bam)

        # Extract replicate ID (rep1, rep2, rep3)
        rep_id=$(echo "$chip_base" | grep -o "rep[0-9]" || true)

        if [[ -z "$rep_id" ]]; then
            echo "⚠ Cannot determine replicate for $chip_base — skipping"
            continue
        fi

        # FLEXIBLE input BAM detection
        input_bam=$(ls "$INPUT_DIR"/*input*"$rep_id"*.dedup.bam 2>/dev/null | head -n1)

        if [[ -z "$input_bam" ]]; then
            echo "⚠ Input BAM not found for $chip_base — skipping"
            continue
        fi

        out_dir="$PEAK_DIR/$condition/$rep_id"
        mkdir -p "$out_dir"

        peak_name="${chip_base}_vs_input"

        if [[ -f "$out_dir/${peak_name}_peaks.narrowPeak" ]]; then
            echo "⏭ Peaks already exist — skipping:"
            echo "   $peak_name"
            continue
        fi

        echo "🔹 Calling peaks:"
        echo "   Condition : $condition"
        echo "   Replicate : $rep_id"
        echo "   ChIP      : $chip_bam"
        echo "   Input     : $input_bam"

        macs2 callpeak \
            -t "$chip_bam" \
            -c "$input_bam" \
            -f BAM \
            -g "$GENOME_SIZE" \
            -n "$peak_name" \
            --outdir "$out_dir" \
            -q "$QVALUE" \
            --nomodel \
            --extsize "$EXTSIZE" \
            --keep-dup all
    done
done

echo "🎉 MACS2 peak calling completed successfully!"

