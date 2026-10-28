#!/usr/bin/env bash
set -euo pipefail

GENE_DIR="results/tem_genes"
DB="reference/tem/TEM_blast_db"
OUT_DIR="results/blastx"

mkdir -p "$OUT_DIR"

SAMPLES=(TEM-1 TEM-2 TEM-3 TEM-10 TEM-24 TEM-57)

# Check required programs and database
command -v blastx >/dev/null || {
    echo "ERROR: BLASTX is not installed."
    exit 1
}

blastdbcmd -db "$DB" -dbtype prot -info >/dev/null || {
    echo "ERROR: Local TEM database not found."
    exit 1
}

for SAMPLE in "${SAMPLES[@]}"; do

    INPUT="$GENE_DIR/${SAMPLE}.fasta"
    OUTPUT="$OUT_DIR/${SAMPLE}_blastx.tsv"

    if [[ ! -s "$INPUT" ]]; then
        echo "ERROR: Missing file $INPUT"
        exit 1
    fi

    echo "================================="
    echo "Running local BLASTX for $SAMPLE"
    echo "================================="

    blastx \
        -query "$INPUT" \
        -db "$DB" \
        -evalue 1e-10 \
        -max_target_seqs 250 \
        -num_threads 2 \
        -outfmt "6 qseqid sacc stitle pident length qcovs evalue bitscore" \
        -out "$OUTPUT"

    echo "Completed: $SAMPLE"
    echo "Results: $OUTPUT"
    echo
done

echo "All six local BLASTX searches completed!"
