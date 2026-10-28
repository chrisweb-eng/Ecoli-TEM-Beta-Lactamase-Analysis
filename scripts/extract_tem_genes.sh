#!/usr/bin/env bash
set -euo pipefail

# Extract complete, exact blaTEM allele matches from
# six E. coli genome assemblies using AMRFinderPlus results.
#
# Requirements: Python 3 (standard library only)
#
# Outputs:
# results/tem_genes/TEM-X.fasta
# results/tem_genes/copies/TEM-X_copyN.fasta
# metadata/tem_extraction_summary.tsv
# metadata/tem_extraction_copies.tsv

cd "$(dirname "$0")/.."

OUT="results/tem_genes"
COPIES="$OUT/copies"
SUMMARY="metadata/tem_extraction_summary.tsv"
DETAILS="metadata/tem_extraction_copies.tsv"

mkdir -p "$OUT" "$COPIES" metadata

python - "$OUT" "$COPIES" "$SUMMARY" "$DETAILS" <<'PY'
import csv
import glob
import os
import sys
from collections import defaultdict

outdir, copiesdir, summary_path, details_path = sys.argv[1:]

samples = ["TEM-1", "TEM-2", "TEM-3",
           "TEM-10", "TEM-24", "TEM-57"]

def read_fasta(path):
    records = {}
    name = None
    chunks = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if name in records:
                        raise ValueError("Duplicate FASTA ID: " + name)
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if name is not None:
        if name in records:
            raise ValueError("Duplicate FASTA ID: " + name)
        records[name] = "".join(chunks).upper()

    return records

def reverse_complement(seq):
    table = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(table)[::-1]

def write_fasta(path, header, sequence):
    with open(path, "w") as handle:
        handle.write(">" + header + "\n")
        for i in range(0, len(sequence), 70):
            handle.write(sequence[i:i+70] + "\n")

summary_rows = []
detail_rows = []
failures = []

for sample in samples:
    expected = "bla" + sample
    report = f"results/amrfinder/{sample}.tsv"
    representative = os.path.join(outdir, sample + ".fasta")

    print("\n" + "=" * 55)
    print("Processing:", sample)
    print("Expected allele:", expected)

    if not os.path.isfile(report):
        print("ERROR: Missing AMRFinderPlus report:", report)
        failures.append(sample)
        continue

    genomes = sorted(glob.glob(
        f"data/genomes/{sample}/extracted/"
        "ncbi_dataset/data/*/*_genomic.fna"
    ))

    if len(genomes) != 1:
        print("ERROR: Expected one genome FASTA; found", len(genomes))
        failures.append(sample)
        continue

    genome_path = genomes[0]
    print("Genome:", genome_path)

    try:
        genome = read_fasta(genome_path)
        with open(report, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            hits = list(reader)
            columns = reader.fieldnames or []

        required = [
            "Contig id", "Start", "Stop", "Strand",
            "Element symbol", "Method",
            "% Coverage of reference",
            "% Identity to reference"
        ]
        missing = [c for c in required if c not in columns]
        if missing:
            raise ValueError("Missing report columns: " + str(missing))

        selected = []

        for hit in hits:
            if hit["Element symbol"] != expected:
                continue
            if hit["Method"] != "ALLELEX":
                continue

            coverage = float(hit["% Coverage of reference"])
            identity = float(hit["% Identity to reference"])

            if abs(coverage - 100) > 0.001:
                continue
            if abs(identity - 100) > 0.001:
                continue

            contig = hit["Contig id"]
            start = int(hit["Start"])
            stop = int(hit["Stop"])
            strand = hit["Strand"]

            if contig not in genome:
                raise ValueError("Contig not found: " + contig)

            if strand not in ("+", "-"):
                raise ValueError("Invalid strand: " + strand)

            left = min(start, stop)
            right = max(start, stop)

            contig_seq = genome[contig]

            if left < 1 or right > len(contig_seq):
                raise ValueError("Coordinates outside contig")

            # AMRFinderPlus reports 1-based inclusive coordinates.
            seq = contig_seq[left-1:right]

            if strand == "-":
                seq = reverse_complement(seq)

            if len(seq) % 3 != 0:
                raise ValueError(
                    "Extracted coding sequence not divisible by 3"
                )

            if any(base not in "ACGT" for base in seq):
                raise ValueError(
                    "Ambiguous nucleotide in extracted sequence"
                )

            selected.append({
                "contig": contig,
                "start": start,
                "stop": stop,
                "strand": strand,
                "coverage": coverage,
                "identity": identity,
                "sequence": seq
            })

        # Remove duplicate report rows for the same genomic locus.
        unique_loci = {}
        for hit in selected:
            key = (hit["contig"], hit["start"],
                   hit["stop"], hit["strand"])
            unique_loci[key] = hit

        selected = sorted(
            unique_loci.values(),
            key=lambda x: (x["contig"], x["start"], x["stop"])
        )

        if not selected:
            raise ValueError(
                "No complete exact " + expected + " allele found"
            )

        # Group copies by their actual nucleotide sequences.
        sequence_groups = defaultdict(list)
        for hit in selected:
            sequence_groups[hit["sequence"]].append(hit)

        print("Complete exact genomic copies:", len(selected))
        print("Distinct nucleotide sequences:", len(sequence_groups))

        # Save each genomic copy.
        for index, hit in enumerate(selected, start=1):
            copy_file = os.path.join(
                copiesdir, f"{sample}_copy{index}.fasta"
            )

            header = (
                f"{sample}_copy{index}"
                f"|{hit['contig']}"
                f"|{hit['start']}-{hit['stop']}"
                f"|strand={hit['strand']}"
            )

            write_fasta(copy_file, header, hit["sequence"])

            detail_rows.append([
                sample, expected, index,
                hit["contig"], hit["start"], hit["stop"],
                hit["strand"], len(hit["sequence"]),
                hit["coverage"], hit["identity"],
                copy_file
            ])

        # Preserve an existing representative if its sequence
        # matches one of the verified exact copies.
        if os.path.isfile(representative):
            old_records = read_fasta(representative)
            old_sequences = list(old_records.values())

            if (len(old_sequences) == 1 and
                    old_sequences[0] in sequence_groups):
                print("Existing representative verified; preserving:",
                      representative)
            else:
                raise ValueError(
                    "Existing representative differs from verified "
                    "copies. File preserved for manual inspection: "
                    + representative
                )
        elif len(sequence_groups) == 1:
            seq = next(iter(sequence_groups))
            write_fasta(representative, sample, seq)
            print("Representative saved:", representative)
        else:
            print(
                "WARNING: Copies differ in nucleotide sequence."
            )
            print(
                "No representative selected automatically."
            )

        summary_rows.append([
            sample, expected, len(selected),
            len(sequence_groups),
            len(selected[0]["sequence"])
            if len(sequence_groups) == 1 else "VARIABLE",
            "IDENTICAL" if len(sequence_groups) == 1
            else "DIFFERENT",
            "YES" if os.path.isfile(representative)
            else "REVIEW_REQUIRED"
        ])

    except (ValueError, OSError, KeyError) as exc:
        print("ERROR:", exc)
        failures.append(sample)
        continue

with open(summary_path, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "Sample", "Expected_allele", "Genomic_copies",
        "Distinct_nucleotide_sequences", "Length_bp",
        "Copy_comparison", "Representative_available"
    ])
    writer.writerows(summary_rows)

with open(details_path, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "Sample", "Allele", "Copy", "Contig",
        "Start", "Stop", "Strand", "Length_bp",
        "Reference_coverage_pct", "Reference_identity_pct",
        "FASTA_file"
    ])
    writer.writerows(detail_rows)

print("\n" + "=" * 55)
print("EXTRACTION SUMMARY")
print("=" * 55)

for row in summary_rows:
    print(
        f"{row[0]}: {row[2]} copies; "
        f"{row[3]} distinct sequences; "
        f"{row[5]}"
    )

print("\nSummary:", summary_path)
print("Copy details:", details_path)
print("FASTA directory:", outdir)

if failures:
    print("\nSamples requiring attention:", ", ".join(failures))
    sys.exit(1)

print("\nAll six samples processed successfully.")
PY
