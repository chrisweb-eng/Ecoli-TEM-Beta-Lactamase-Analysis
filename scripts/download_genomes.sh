#!/usr/bin/env bash

set -uo pipefail

# Move to the project directory
cd "$(dirname "$0")/.."

# Input and output directories
input="metadata/assembly_verification.tsv"
output_dir="data/genomes"
log="metadata/genome_download_log.tsv"

mkdir -p "$output_dir"

# Create download log
printf "TEM_variant\tAssembly\tStatus\n" > "$log"

# Check that verification results exist
if [[ ! -f "$input" ]]; then
    echo "ERROR: Assembly verification file not found."
    exit 1
fi

# Read verified assemblies
while IFS=$'\t' read -r variant accession status; do

    # Skip header
    [[ "$variant" == "TEM_variant" ]] && continue

    # Only process available assemblies
    [[ "$status" != "AVAILABLE" ]] && continue

    echo "=================================="
    echo "Processing $variant ($accession)"
    echo "=================================="

    genome_dir="$output_dir/$variant"
    mkdir -p "$genome_dir"

    zip_file="$genome_dir/genome.zip"
    extract_dir="$genome_dir/extracted"

    # Check whether a genome FASTA already exists
    if find "$extract_dir" -type f -name "*_genomic.fna" \
        -size +0c -print -quit 2>/dev/null | grep -q .; then

        echo "Genome already downloaded. Skipping."
        result="SKIPPED_EXISTING"

    else
        echo "Downloading genome..."

        rm -f "$zip_file"

        if datasets download genome accession "$accession" \
            --include genome \
            --filename "$zip_file"; then

            mkdir -p "$extract_dir"

            if unzip -oq "$zip_file" -d "$extract_dir" &&
               find "$extract_dir" -type f \
                   -name "*_genomic.fna" -size +0c \
                   -print -quit | grep -q .; then

                echo "Download and extraction successful."
                result="DOWNLOADED"

                rm -f "$zip_file"

            else
                echo "ERROR: Genome extraction or FASTA check failed."
                result="EXTRACTION_FAILED"
            fi

        else
            echo "ERROR: Download failed."
            result="DOWNLOAD_FAILED"
        fi
    fi

    printf "%s\t%s\t%s\n" \
        "$variant" "$accession" "$result" >> "$log"

done < "$input"

echo
echo "=================================="
echo "Genome download summary"
echo "=================================="

column -t -s $'\t' "$log"
