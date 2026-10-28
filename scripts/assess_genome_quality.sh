#!/usr/bin/env bash

set -euo pipefail

# --------------------------------------------------
# Project: E. coli TEM beta-lactamase analysis
# Purpose: Assess genome assembly quality with QUAST
# --------------------------------------------------

# Move to the project root directory
cd "$(dirname "$0")/.." || exit 1

# Define directories
input_dir="data/genomes"
output_dir="results/quast"

# Check that QUAST is installed
if ! command -v quast.py >/dev/null 2>&1; then
    echo "ERROR: QUAST (quast.py) is not installed or not in PATH."
    exit 1
fi

# Check that genome directory exists
if [[ ! -d "$input_dir" ]]; then
    echo "ERROR: Genome directory not found: $input_dir"
    exit 1
fi

# Find genome assembly FASTA files
mapfile -d '' genome_files < <(
    find "$input_dir" -type f \
        \( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \) \
        -print0 | sort -z
)

# Ensure genome assemblies are available
if [[ ${#genome_files[@]} -eq 0 ]]; then
    echo "ERROR: No genome assembly FASTA files found."
    echo "Check $input_dir before running QUAST."
    exit 1
fi

# Create output directory
mkdir -p "$output_dir"

echo "========================================"
echo "E. coli genome assembly quality assessment"
echo "========================================"
echo "Genome assemblies found: ${#genome_files[@]}"
echo "Output directory: $output_dir"
echo

# Run QUAST on all selected genome assemblies
quast.py \
    "${genome_files[@]}" \
    -o "$output_dir"

echo
echo "========================================"
echo "QUAST analysis completed successfully."
echo "========================================"
echo "Reports saved in: $output_dir"
echo "Main report: $output_dir/report.html"
