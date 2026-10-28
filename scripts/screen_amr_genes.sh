#!/usr/bin/env bash

set -uo pipefail

# =====================================================
# Project: Comparative Analysis of TEM Beta-Lactamases
# Purpose: Screen six E. coli genome assemblies
# Tool: NCBI AMRFinderPlus
# =====================================================

# Move to project root
cd "$(dirname "$0")/.." || exit 1

# Directories
genome_dir="data/genomes"
output_dir="results/amrfinder"
summary="metadata/tem_gene_verification.tsv"

mkdir -p "$output_dir" metadata

# Six selected TEM variants
variants=(
    "TEM-1"
    "TEM-2"
    "TEM-3"
    "TEM-10"
    "TEM-24"
    "TEM-57"
)

# Create summary table
printf "Sample\tExpected_gene\tDetected_TEM_genes\tStatus\n" \
    > "$summary"

# =====================================================
# Screen each genome
# =====================================================

for variant in "${variants[@]}"; do

    echo
    echo "========================================"
    echo "Screening $variant"
    echo "========================================"

    expected="bla${variant}"
    result="$output_dir/${variant}.tsv"

    # Locate downloaded genome FASTA
    fasta=$(find "$genome_dir/$variant/extracted" \
        -type f -name "*_genomic.fna" \
        -print -quit 2>/dev/null)

    if [[ -z "$fasta" || ! -s "$fasta" ]]; then

        echo "ERROR: Genome FASTA not found."

        printf "%s\t%s\t%s\t%s\n" \
            "$variant" "$expected" "NA" \
            "FASTA_MISSING" >> "$summary"

        continue
    fi

    echo "Genome: $fasta"
    echo "Running AMRFinderPlus..."

    # Run AMRFinderPlus on nucleotide assembly
    if ! amrfinder \
        -n "$fasta" \
        -o "$result"; then

        echo "ERROR: AMRFinderPlus failed."

        printf "%s\t%s\t%s\t%s\n" \
            "$variant" "$expected" "NA" \
            "AMRFINDER_FAILED" >> "$summary"

        continue
    fi

    # =================================================
    # Extract TEM genes from AMRFinderPlus results
    # Correct column: Element symbol
    # =================================================

    detected=$(awk -F '\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "Element symbol") {
                    col = i
                }
            }
            next
        }

        col > 0 && $col ~ /^blaTEM/ {
            print $col
        }
    ' "$result" | sort -u | paste -sd "," -)

    # =================================================
    # Compare detected genes with expected allele
    # =================================================

    if [[ -z "$detected" ]]; then

        detected="NONE"
        status="NOT_DETECTED"

    elif printf '%s\n' "$detected" |
         tr ',' '\n' |
         grep -Fxq "$expected"; then

        status="MATCH"

    else

        status="MISMATCH"

    fi

    # Save summary
    printf "%s\t%s\t%s\t%s\n" \
        "$variant" "$expected" "$detected" "$status" \
        >> "$summary"

    echo "Expected gene: $expected"
    echo "Detected genes: $detected"
    echo "Verification status: $status"

done

# =====================================================
# Display final results
# =====================================================

echo
echo "========================================"
echo "TEM GENE VERIFICATION SUMMARY"
echo "========================================"

column -t -s $'\t' "$summary"

echo
echo "Individual AMRFinderPlus results:"
echo "$output_dir"

echo "Summary file:"
echo "$summary"

echo
echo "Screening completed."
