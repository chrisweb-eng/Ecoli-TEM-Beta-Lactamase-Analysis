#!/usr/bin/env python3

import csv
from pathlib import Path
from itertools import combinations

# Input and output
input_file = Path("results/alignments/TEM_nucleotide_aligned.fasta")
output_dir = Path("results/variants")
output_dir.mkdir(parents=True, exist_ok=True)

# Expected sequence names
variants = [
    "TEM-1", "TEM-2", "TEM-3",
    "TEM-10", "TEM-24", "TEM-57"
]

# Read aligned FASTA
sequences = {}
header = None

with open(input_file) as f:
    for line in f:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            header = line[1:].split()[0]

            # TEM-1 has a longer FASTA header
            if header.startswith("TEM-1_blaTEM-1"):
                header = "TEM-1"

            if header in sequences:
                raise ValueError(f"Duplicate sequence: {header}")

            sequences[header] = ""

        else:
            sequences[header] += line.upper()

# Validate sequences
if set(sequences) != set(variants):
    raise ValueError(
        f"Unexpected sequence names: {list(sequences)}"
    )

lengths = {len(seq) for seq in sequences.values()}

if len(lengths) != 1:
    raise ValueError("Alignment lengths are inconsistent")

alignment_length = lengths.pop()

# Identify variable nucleotide positions
variable_rows = []
conserved_count = 0

for i in range(alignment_length):

    bases = [sequences[v][i] for v in variants]

    if any(base not in "ACGT" for base in bases):
        raise ValueError(
            f"Gap or ambiguous base at alignment position {i+1}"
        )

    if len(set(bases)) == 1:
        conserved_count += 1
    else:
        row = {"Position": i + 1}

        for variant in variants:
            row[variant] = sequences[variant][i]

        variable_rows.append(row)

# Save variable nucleotide positions
with open(
    output_dir / "TEM_variable_positions.csv",
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=["Position"] + variants
    )

    writer.writeheader()
    writer.writerows(variable_rows)

# Compare all 15 pairs
pairwise_rows = []

for v1, v2 in combinations(variants, 2):

    differences = sum(
        a != b
        for a, b in zip(sequences[v1], sequences[v2])
    )

    pairwise_rows.append({
        "Variant_1": v1,
        "Variant_2": v2,
        "Differences": differences,
        "Identity_Percent": round(
            100 * (alignment_length - differences)
            / alignment_length, 3
        )
    })

# Save pairwise comparisons
with open(
    output_dir / "TEM_pairwise_differences.csv",
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Variant_1", "Variant_2",
            "Differences", "Identity_Percent"
        ]
    )

    writer.writeheader()
    writer.writerows(pairwise_rows)

# Save summary
summary = [
    ["Alignment_Length", alignment_length],
    ["Number_of_Variants", len(variants)],
    ["Conserved_Positions", conserved_count],
    ["Variable_Positions", len(variable_rows)],
    ["Pairwise_Comparisons", len(pairwise_rows)]
]

with open(
    output_dir / "TEM_nucleotide_summary.csv",
    "w",
    newline=""
) as f:

    writer = csv.writer(f)
    writer.writerow(["Metric", "Value"])
    writer.writerows(summary)

# Print results
print("\nNUCLEOTIDE COMPARISON COMPLETED")
print("--------------------------------")
print(f"Variants analysed: {len(variants)}")
print(f"Alignment length: {alignment_length} bp")
print(f"Conserved positions: {conserved_count}")
print(f"Variable positions: {len(variable_rows)}")
print(f"Pairwise comparisons: {len(pairwise_rows)}")

print("\nPairwise nucleotide differences:")
for row in pairwise_rows:
    print(
        f"{row['Variant_1']} vs {row['Variant_2']}: "
        f"{row['Differences']} differences "
        f"({row['Identity_Percent']}% identity)"
    )

print("\nResults saved in results/variants/")
