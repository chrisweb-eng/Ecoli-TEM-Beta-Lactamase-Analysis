#!/usr/bin/env bash
set -euo pipefail

# E. coli TEM beta-lactamase comparative genomics
# Genome assembly quality assessment using QUAST

cd "$HOME/Ecoli_AMR_Project"

samples=(TEM-1 TEM-2 TEM-3 TEM-10 TEM-24 TEM-57)

output_dir="results/quast"
mkdir -p "$output_dir"

genomes=()
labels=()

echo "========================================"
echo "E. coli genome assembly quality assessment"
echo "========================================"

# Locate and validate the six genome assemblies
for sample in "${samples[@]}"; do
    echo "Checking $sample..."

    genome_dir="data/genomes/$sample/extracted/ncbi_dataset/data"

    if [[ ! -d "$genome_dir" ]]; then
        echo "ERROR: Missing genome directory for $sample"
        exit 1
    fi

    mapfile -t matches < <(
        find "$genome_dir" -type f -name "*_genomic.fna" | sort
    )

    if [[ ${#matches[@]} -ne 1 ]]; then
        echo "ERROR: Expected one genome FASTA for $sample"
        echo "Found: ${#matches[@]}"
        exit 1
    fi

    genome="${matches[0]}"

    if [[ ! -s "$genome" ]]; then
        echo "ERROR: Empty genome FASTA for $sample"
        exit 1
    fi

    genomes+=("$genome")
    labels+=("$sample")

    echo "Found: $genome"
done

echo
echo "Running QUAST on all six assemblies..."

# Create comma-separated labels
labels_csv=$(IFS=,; echo "${labels[*]}")

quast.py \
    "${genomes[@]}" \
    --labels "$labels_csv" \
    --threads 4 \
    --min-contig 500 \
    --output-dir "$output_dir"

echo
echo "========================================"
echo "QUAST analysis completed"
echo "========================================"

echo "Results directory: $output_dir"
echo "Combined report: $output_dir/report.tsv"
echo "HTML report: $output_dir/report.html"

if [[ -s "$output_dir/report.tsv" ]]; then
    echo
    echo "Assembly quality summary:"
    cat "$output_dir/report.tsv"
else
    echo "WARNING: QUAST report.tsv not found"
fi
