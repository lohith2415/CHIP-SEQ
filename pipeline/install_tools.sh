#!/bin/bash
set -euo pipefail

############################################
# Paths
############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="$SCRIPT_DIR/tools"

FASTQC_DIR="$TOOLS_DIR/fastqc"
MACS2_ENV="$TOOLS_DIR/macs2_env"

# 👉 Trimmomatic OUTSIDE tools folder
TRIMMOMATIC_BASE="/opt/trimmomatic"
TRIMMOMATIC_DIR="$TRIMMOMATIC_BASE/Trimmomatic-0.39"
TRIMMOMATIC_JAR="$TRIMMOMATIC_DIR/trimmomatic-0.39.jar"

mkdir -p "$TOOLS_DIR"

echo "📁 Tools directory: $TOOLS_DIR"
echo "📁 Trimmomatic directory: $TRIMMOMATIC_DIR"
echo

############################################
# Helper: install apt package if missing
############################################
apt_install_if_missing () {
    local cmd=$1
    local pkg=$2

    if command -v "$cmd" &>/dev/null; then
        echo "⏭ $cmd already installed"
    else
        echo "📦 Installing $pkg..."
        sudo apt update
        sudo apt install -y "$pkg"
    fi
}

############################################
# SYSTEM TOOLS
############################################
echo "🔧 Checking system tools..."

apt_install_if_missing java openjdk-11-jre
apt_install_if_missing bowtie2 bowtie2
apt_install_if_missing samtools samtools
apt_install_if_missing wget wget
apt_install_if_missing unzip unzip
apt_install_if_missing python3.10 python3.10
apt_install_if_missing pip3 python3-pip

# Python venv support
if python3.10 -m venv --help &>/dev/null; then
    echo "⏭ python3.10-venv already available"
else
    echo "📦 Installing python3.10-venv..."
    sudo apt install -y python3.10-venv python3.10-dev
fi

############################################
# FASTQC (local install)
############################################
if [[ -x "$FASTQC_DIR/fastqc" ]]; then
    echo "⏭ FastQC already installed"
else
    echo "📦 Installing FastQC..."
    mkdir -p "$FASTQC_DIR"
    cd "$FASTQC_DIR"
    wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
    unzip -q fastqc_v0.12.1.zip
    chmod +x FastQC/fastqc
    ln -sf "$FASTQC_DIR/FastQC/fastqc" "$FASTQC_DIR/fastqc"
fi

############################################
# TRIMMOMATIC (GLOBAL install)
############################################
if [[ -f "$TRIMMOMATIC_JAR" ]]; then
    echo "⏭ Trimmomatic already installed at $TRIMMOMATIC_DIR"
else
    echo "📦 Installing Trimmomatic globally..."
    sudo mkdir -p "$TRIMMOMATIC_BASE"
    sudo chown "$USER":"$USER" "$TRIMMOMATIC_BASE"

    cd "$TRIMMOMATIC_BASE"
    wget -q http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
    unzip -q Trimmomatic-0.39.zip
fi

############################################
# MACS2 (Python venv)
############################################
if [[ -x "$MACS2_ENV/bin/macs2" ]]; then
    echo "⏭ MACS2 already installed"
else
    echo "📦 Installing MACS2 in Python venv..."
    python3.10 -m venv "$MACS2_ENV"
    source "$MACS2_ENV/bin/activate"
    pip install --upgrade pip
    pip install macs2
    deactivate
fi

############################################
# SUMMARY
############################################
echo
echo "🎉 ChIP-seq tool setup completed!"
echo
echo "🔹 FastQC:      $FASTQC_DIR/fastqc"
echo "🔹 Trimmomatic: java -jar $TRIMMOMATIC_JAR"
echo "🔹 Bowtie2:     $(command -v bowtie2)"
echo "🔹 SAMtools:   $(command -v samtools)"
echo "🔹 MACS2:      source $MACS2_ENV/bin/activate"

