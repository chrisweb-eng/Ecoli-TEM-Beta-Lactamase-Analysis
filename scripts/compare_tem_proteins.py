#!/usr/bin/env python3

import csv
from pathlib import Path
from itertools import combinations

input_file = Path(
    "results/alignments/TEM_protein_aligned.faa"
)

output_dir = Path("results/variants")
output_dir.mkdir(parents=True, exist_ok=True)

variants = [
    "TEM-1", "TEM-2", "TEM-3",
    "TEM-10", "TEM-24", "TEM-57"
]

sequences = {}
current = None

# Read aligned protein FASTA
with open(input_file) as f:
    for line in f:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            header = line[1:].split()[0]

            # Extract the variant name
            current = header.split("_")[0]

            if current not in variants:
                raise ValueError(
                    f"Unexpected sequence: {header}"
                )

            if current in sequences:
                raise ValueError(
                    f"Duplicate sequence: {current}"
                )

            sequences[current] = ""
        else:
            if current is None:
                raise ValueError("Sequence without header")
            sequences[current] += line.upper()

if set(sequences) != set(variants):
    raise ValueError(
        f"Expected six variants, found: {list(sequences)}"
    )

lengths = {len(s) for s in sequences.values()}

if len(lengths) != 1:
    raise ValueError("Different alignment lengths")

alignment_length = lengths.pop()

# Check amino acid characters
valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

for variant, sequence in sequences.items():
    unexpected = set(sequence) - valid_amino_acids
    if unexpected:
        raise ValueError(
            f"{variant} contains unexpected symbols: "
            f"{unexpected}"
        )

# Identify variable positions
variable_rows = []
conserved_count = 0

for i in range(alignment_length):

    residues = [
        sequences[variant][i]
        for variant in variants
    ]

    if len(set(residues)) == 1:
        conserved_count += 1
    else:
        row = {"Position": i + 1}

        for variant in variants:
            row[variant] = sequences[variant][i]

        variable_rows.append(row)

# Save variable amino acid positions
with open(
    output_dir / "TEM_amino_acid_variations.csv",
    "w",
    newline=""
) as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["Position"] + variants
    )
    writer.writeheader()
    writer.writerows(variable_rows)

# Pairwise protein differences
pairwise_rows = []

for v1, v2 in combinations(variants, 2):

    differences = sum(
        a != b
        for a, b in zip(sequences[v1], sequences[v2])
    )

    identity = (
        100 * (alignment_length - differences)
        / alignment_length
    )

    pairwise_rows.append({
        "Variant_1": v1,
        "Variant_2": v2,
        "Differences": differences,
        "Identity_Percent": round(identity, 3)
    })

with open(
    output_dir / "TEM_protein_pairwise_differences.csv",
    "w",
    newline=""
) as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Variant_1",
            "Variant_2",
            "Differences",
            "Identity_Percent"
        ]
    )
    writer.writeheader()
    writer.writerows(pairwise_rows)

# Save summary
with open(
    output_dir / "TEM_protein_summary.csv",
    "w",
    newline=""
) as f:
    writer = csv.writer(f)
    writer.writerow(["Metric", "Value"])
    writer.writerow(["Alignment_Length", alignment_length])
    writer.writerow(["Number_of_Variants", len(variants)])
    writer.writerow(["Conserved_Positions", conserved_count])
    writer.writerow(["Variable_Positions", len(variable_rows)])
    writer.writerow(["Pairwise_Comparisons", len(pairwise_rows)])

# Display results
print("\nPROTEIN COMPARISON COMPLETED")
print("--------------------------------")
print(f"Variants analysed: {len(variants)}")
print(f"Alignment length: {alignment_length}")
print(f"Conserved positions: {conserved_count}")
print(f"Variable positions: {len(variable_rows)}")

print("\nVariable amino acid positions:")
for row in variable_rows:
    print(row)

print("\nPairwise protein differences:")
for row in pairwise_rows:
    print(
        f"{row['Variant_1']} vs {row['Variant_2']}: "
        f"{row['Differences']} differences "
        f"({row['Identity_Percent']}% identity)"
    )

print("\nResults saved in results/variants/")
