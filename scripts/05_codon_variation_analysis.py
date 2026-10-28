from pathlib import Path
from itertools import combinations
import csv

# Project files
alignment = Path("results/alignments/TEM_nucleotide_aligned.fasta")
output_dir = Path("results/variants")
output_dir.mkdir(parents=True, exist_ok=True)

# Standard genetic code
genetic_code = {
    codon: amino_acid
    for amino_acid, codons in {
        "F": ["TTT", "TTC"],
        "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        "I": ["ATT", "ATC", "ATA"],
        "M": ["ATG"],
        "V": ["GTT", "GTC", "GTA", "GTG"],
        "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "P": ["CCT", "CCC", "CCA", "CCG"],
        "T": ["ACT", "ACC", "ACA", "ACG"],
        "A": ["GCT", "GCC", "GCA", "GCG"],
        "Y": ["TAT", "TAC"],
        "*": ["TAA", "TAG", "TGA"],
        "H": ["CAT", "CAC"],
        "Q": ["CAA", "CAG"],
        "N": ["AAT", "AAC"],
        "K": ["AAA", "AAG"],
        "D": ["GAT", "GAC"],
        "E": ["GAA", "GAG"],
        "C": ["TGT", "TGC"],
        "W": ["TGG"],
        "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "G": ["GGT", "GGC", "GGA", "GGG"],
    }.items()
    for codon in codons
}

# Read aligned FASTA sequences
sequences = {}
name = None

with alignment.open() as file:
    for line in file:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = ""
        else:
            sequences[name] += line.upper()

if len(sequences) != 6:
    raise ValueError("Expected exactly six TEM sequences.")

lengths = {len(seq) for seq in sequences.values()}
if len(lengths) != 1:
    raise ValueError("Sequences have different lengths.")

length = lengths.pop()
if length % 3 != 0:
    raise ValueError("Sequence length is not divisible by 3.")

if any(set(seq) - set("ACGT") for seq in sequences.values()):
    raise ValueError("Sequences contain gaps or ambiguous bases.")

names = list(sequences)
codon_count = length // 3

# Store every codon and translated amino acid
codons = {
    name: [
        sequences[name][i:i + 3]
        for i in range(0, length, 3)
    ]
    for name in names
}

amino_acids = {
    name: [genetic_code[codon] for codon in codons[name]]
    for name in names
}

# Identify variable codon positions across all six sequences
variable_positions = [
    pos
    for pos in range(codon_count)
    if len({codons[name][pos] for name in names}) > 1
]

# Output 1: Codons and amino acids at variable positions
codon_file = output_dir / "TEM_codon_variation.csv"

with codon_file.open("w", newline="") as file:
    writer = csv.writer(file)

    header = ["Codon_Position", "Variation_Type"]
    for name in names:
        header.extend([f"{name}_Codon", f"{name}_AA"])
    writer.writerow(header)

    for pos in variable_positions:
        distinct_aa = {amino_acids[name][pos] for name in names}
        variation_type = (
            "Synonymous_only"
            if len(distinct_aa) == 1
            else "Nonsynonymous_present"
        )

        row = [pos + 1, variation_type]
        for name in names:
            row.extend([codons[name][pos], amino_acids[name][pos]])
        writer.writerow(row)

# Output 2: All 15 pairwise comparisons
pairwise_file = output_dir / "TEM_pairwise_codon_classification.csv"

with pairwise_file.open("w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow([
        "Variant_1",
        "Variant_2",
        "Identical_Codons",
        "Synonymous_Codon_Differences",
        "Nonsynonymous_Codon_Differences",
        "Nucleotide_Differences",
        "Amino_Acid_Differences",
    ])

    for name1, name2 in combinations(names, 2):
        identical = 0
        synonymous = 0
        nonsynonymous = 0
        nucleotide_differences = 0
        amino_acid_differences = 0

        for pos in range(codon_count):
            c1 = codons[name1][pos]
            c2 = codons[name2][pos]
            a1 = amino_acids[name1][pos]
            a2 = amino_acids[name2][pos]

            nucleotide_differences += sum(
                x != y for x, y in zip(c1, c2)
            )

            if c1 == c2:
                identical += 1
            elif a1 == a2:
                synonymous += 1
            else:
                nonsynonymous += 1
                amino_acid_differences += 1

        writer.writerow([
            name1,
            name2,
            identical,
            synonymous,
            nonsynonymous,
            nucleotide_differences,
            amino_acid_differences,
        ])

print("Analysis completed.")
print(f"Sequences analysed: {len(names)}")
print(f"Codons per sequence: {codon_count}")
print(f"Variable codon positions: {len(variable_positions)}")
print(f"Pairwise comparisons: 15")
print(f"Saved: {codon_file}")
print(f"Saved: {pairwise_file}")
