#!/bin/bash
set -euo pipefail

############################################
# Detect paths
############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

GENOME_DATA_DIR="$PROJECT_ROOT/genome_data"
INDEX_DIR="$PROJECT_ROOT/genome_index"
INDEX_PREFIX="$INDEX_DIR/TAIR10"

TRIM_DIR="$PROJECT_ROOT/trimmed"
BAM_DIR="$PROJECT_ROOT/bam"

THREADS=${THREADS:-4}

mkdir -p "$INDEX_DIR" "$BAM_DIR"

############################################
# Detect genome (fa or fa.gz)
############################################
GENOME_FA=$(find "$GENOME_DATA_DIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) | head -n1)
GENOME_GZ=$(find "$GENOME_DATA_DIR" -maxdepth 1 -type f \( -name "*.fa.gz" -o -name "*.fasta.gz" \) | head -n1)

if [[ -z "$GENOME_FA" && -z "$GENOME_GZ" ]]; then
    echo "❌ No genome FASTA found in genome_data/"
    exit 1
fi

############################################
# Unzip genome if needed (KEEP original .gz)
############################################
if [[ -z "$GENOME_FA" && -n "$GENOME_GZ" ]]; then
    GENOME_FA="${GENOME_GZ%.gz}"
    if [[ ! -f "$GENOME_FA" ]]; then
        echo "🧬 Unzipping genome:"
        echo "   $GENOME_GZ"
        gunzip -c "$GENOME_GZ" > "$GENOME_FA"
        echo "✅ Genome unzipped to:"
        echo "   $GENOME_FA"
    else
        echo "⏭ Unzipped genome already exists — skipping unzip"
    fi
fi

############################################
# Build Bowtie2 index (only if missing)
############################################
if [[ ! -f "${INDEX_PREFIX}.1.bt2" ]]; then
    echo "🧬 Building Bowtie2 index..."
    echo "📌 Genome: $GENOME_FA"
    echo "📁 Index dir: $INDEX_DIR"

    bowtie2-build "$GENOME_FA" "$INDEX_PREFIX"
    echo "✅ Genome indexing completed"
else
    echo "⏭ Bowtie2 index already exists — skipping indexing"
fi

############################################
# Find trimmed FASTQ files
############################################
echo "🔍 Searching trimmed FASTQ files..."
FASTQ_FILES=$(find "$TRIM_DIR" -type f -name "*_trimmed.fastq.gz")

if [[ -z "$FASTQ_FILES" ]]; then
    echo "❌ No trimmed FASTQ files found"
    exit 1
fi

############################################
# Alignment (structure preserved)
############################################
echo "🧬 Starting alignment..."

for fq in $FASTQ_FILES; do
    rel_path="${fq#$TRIM_DIR/}"
    out_dir="$BAM_DIR/$(dirname "$rel_path")"
    mkdir -p "$out_dir"

    base=$(basename "$fq" _trimmed.fastq.gz)
    bam_out="$out_dir/${base}.sorted.bam"

    if [[ -f "$bam_out" ]]; then
        echo "⏭ BAM exists — skipping:"
        echo "   $bam_out"
        continue
    fi

    echo "🔹 Aligning: $rel_path"

    bowtie2 \
        -x "$INDEX_PREFIX" \
        -U "$fq" \
        -p "$THREADS" | \
    samtools view -bS - | \
    samtools sort -@ "$THREADS" -o "$bam_out"

    samtools index "$bam_out"
done

echo "🎉 Genome unzip + indexing + alignment completed!"

