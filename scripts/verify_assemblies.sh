#!/usr/bin/env bash

set -uo pipefail

mkdir -p metadata

output="metadata/assembly_verification.tsv"

printf "TEM_variant\tAssembly\tStatus\n" > "$output"

assemblies=(
    "TEM-1:GCF_000010765.1"
    "TEM-2:GCA_032298435.1"
    "TEM-3:GCA_012578525.1"
    "TEM-10:GCA_026378375.1"
    "TEM-24:GCA_059152815.1"
    "TEM-57:GCA_061231245.1"
)

for entry in "${assemblies[@]}"; do

    variant="${entry%%:*}"
    accession="${entry#*:}"

    echo "Checking $variant ($accession)..."

    if datasets summary genome accession "$accession" \
        > /dev/null 2>&1; then

        status="AVAILABLE"

    else
        status="CHECK_FAILED"

    fi

    printf "%s\t%s\t%s\n" \
        "$variant" "$accession" "$status" >> "$output"

done

echo
echo "Verification results:"
column -t -s $'\t' "$output"
