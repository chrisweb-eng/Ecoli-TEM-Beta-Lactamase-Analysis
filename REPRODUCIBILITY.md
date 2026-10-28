# Reproducibility Guide

## Comparative Analysis of TEM Beta-Lactamase Gene Variants in Escherichia coli

This document describes how to reproduce the computational analyses performed in this project, including nucleotide sequence comparison, amino acid sequence comparison, codon-level variation analysis, and visualisation of sequence differences.

The study examines six TEM beta-lactamase variants: TEM-1, TEM-2, TEM-3, TEM-10, TEM-24, and TEM-57.

The repository contains the Bash, Python, and R scripts used in the analysis, together with reference sequences, selected intermediate files, summary tables, and visualisations.

Two levels of reproducibility are supported:

1. **Sequence-level analysis:** Reproduce the comparative analyses using the extracted TEM sequences and alignments provided in the repository.
2. **Genome-level workflow:** Repeat the upstream genome retrieval, quality assessment, antimicrobial resistance gene screening, and TEM sequence extraction using the documented NCBI assembly accessions.

The sequence-level workflow is the recommended starting point because it uses the existing repository files and does not require downloading complete bacterial genomes.

This project investigates sequence variation and does not independently establish antibiotic susceptibility, minimum inhibitory concentrations, or the functional effects of individual mutations.

## 1. Software Environment

This project uses Linux/Bash, Python, R and established bioinformatics software to investigate sequence variation among six selected TEM beta-lactamase variants in Escherichia coli.

The original analyses were developed using Ubuntu through Windows Subsystem for Linux (WSL).

### 1.1 Clone the repository

```bash
git clone https://github.com/chrisweb-eng/Ecoli-TEM-Beta-Lactamase-Analysis.git
cd Ecoli-TEM-Beta-Lactamase-Analysis
```

### 1.2 Set up the Conda environment

Install Miniconda or Anaconda, then create the environment using the supplied configuration:

```bash
conda env create -f environment.yml
conda activate ecoli_amr
```

The environment is intended to provide the principal command-line dependencies, including Python, NCBI Datasets CLI, QUAST, AMRFinderPlus, BLAST+, SeqKit and MAFFT.

The environment file must be committed to the repository before these setup commands can be used by other researchers.

### 1.3 Install R dependencies

The four visualisation scripts use R, including the `ggplot2` and `tidyr` packages.

If R is not included in the Conda environment, install it separately. Then run:

```r
install.packages(c("ggplot2", "tidyr"))
```

Other packages should be installed if required by the `library()` statements in the individual scripts.

### 1.4 Reference databases

AMRFinderPlus requires its appropriate reference database.

The BLASTX verification stage uses a local TEM protein reference database generated from:

`reference/tem/TEM_reference_proteins.fasta`

Create the local BLAST database using:

```bash
makeblastdb \
  -in reference/tem/TEM_reference_proteins.fasta \
  -dbtype prot \
  -out reference/tem/TEM_blast_db
```

Reference database versions can affect upstream identification results and should be recorded when repeating the analysis.

---

## 2. Repository Organisation

The repository contains the following principal directories:

```text
Ecoli-TEM-Beta-Lactamase-Analysis/
├── README.md
├── environment.yml
├── REPRODUCIBILITY.md
├── metadata/
├── reference/
│   └── tem/
├── scripts/
└── results/
    ├── alignments/
    ├── amrfinder/
    ├── blastx/
    ├── figures/
    ├── proteins/
    ├── quast/
    ├── tem_genes/
    └── variants/
```

The `metadata/` directory contains records documenting genome assembly selection, downloads, TEM gene verification and extraction.

The `scripts/` directory contains the Bash, Python and R scripts used throughout the analysis.

The `results/` directory contains selected intermediate files, comparative tables and final visualisations.

Complete genome assembly downloads are retained locally under `data/genomes/` and are not committed to the repository. They can be retrieved again using the recorded NCBI assembly accessions.

---

## 3. Input Data

Six representative TEM beta-lactamase coding sequences were selected from publicly available E. coli genome assemblies.

| TEM variant | NCBI assembly accession |
|-------------|-------------------------|
| TEM-1 | GCF_000010765.1 |
| TEM-2 | GCA_032298435.1 |
| TEM-3 | GCA_012578525.1 |
| TEM-10 | GCA_026378375.1 |
| TEM-24 | GCA_059152815.1 |
| TEM-57 | GCA_061231245.1 |

Assembly verification records are available in:

`metadata/assembly_verification.tsv`

One selected TEM coding sequence per assembly was included in the final comparative analysis.

Each selected nucleotide sequence contains 858 nucleotides, corresponding to a translated protein of 286 amino acids.

The extracted nucleotide sequences are available in:

`results/tem_genes/`

The corresponding translated protein sequences are available in:

`results/proteins/`

The published sequences and alignments provide the recommended starting point for reproducing the comparative analysis.

---

## 4. Upstream Genome Processing

This section documents the original genome-level workflow.

Researchers who wish to reproduce only the sequence comparisons can skip to Section 5.

### 4.1 Assembly verification

The selected NCBI genome assemblies were checked using:

```bash
bash scripts/verify_assemblies.sh
```

### 4.2 Genome retrieval

The assemblies were downloaded using:

```bash
bash scripts/download_genomes.sh
```

Download records are maintained in:

`metadata/genome_download_log.tsv`

### 4.3 Genome assembly quality assessment

QUAST was used to examine genome assembly quality.

```bash
bash scripts/assess_genome_quality.sh
```

Outputs are organised under:

`results/quast/`

### 4.4 Antimicrobial resistance gene screening

NCBI AMRFinderPlus was used to identify antimicrobial resistance-associated genes in the selected assemblies.

```bash
bash scripts/screen_amr_genes.sh
```

Results are organised under:

`results/amrfinder/`

### 4.5 TEM gene extraction

Selected TEM coding sequences were extracted using:

```bash
bash scripts/extract_tem_genes.sh
```

The extraction records include:

- `metadata/tem_gene_verification.tsv`
- `metadata/tem_extraction_summary.tsv`
- `metadata/tem_extraction_copies.tsv`

The final selected nucleotide sequences are stored under:

`results/tem_genes/`

### 4.6 BLASTX verification

The extracted nucleotide sequences were checked against the TEM reference protein collection using BLASTX.

After preparing the local BLAST database described in Section 1.4, the verification script can be invoked using:

```bash
bash scripts/run_blastx.sh
```

The results are organised under:

`results/blastx/`

### 4.7 Important upstream requirements

The upstream scripts may require downloaded assemblies, prepared reference databases and particular intermediate files.

Their expected input paths and processing requirements should be checked before execution in a fresh environment.

This upstream workflow is documented as a sequence of individual processing stages, not as an independently validated one-command pipeline.

---

## 5. Nucleotide and Protein Sequence Alignment

The extracted TEM coding sequences were aligned using MAFFT.

The translated TEM protein sequences were also aligned to support amino acid-level comparisons.

The relevant files are:

```text
results/alignments/TEM_combined.fasta
results/alignments/TEM_nucleotide_aligned.fasta
results/alignments/TEM_protein_aligned.faa
```

The final selected alignments contain six sequences, with no alignment gaps.

For direct reproduction of the comparative analyses, the existing aligned FASTA files should be used as input.

Researchers wishing to repeat the alignment stage should use the same selected nucleotide and protein sequences and record the MAFFT version and settings used.

---

## 6. Nucleotide Variation Analysis

The nucleotide comparison was performed using:

`scripts/compare_all_tem.py`

Run the script from the repository root:

```bash
python scripts/compare_all_tem.py
```

The corresponding result files include:

```text
results/variants/TEM_variable_positions.csv
results/variants/TEM_nucleotide_summary.csv
results/variants/TEM_pairwise_differences.csv
```

The analysis identified:

- 858 aligned nucleotide positions
- 846 conserved nucleotide positions
- 12 variable nucleotide positions
- 15 unique pairwise comparisons

Pairwise nucleotide differences ranged from 1 to 10 positions among the selected variants.

These counts describe sequence differences and variable alignment positions, not separate evolutionary mutation events.

---

## 7. Protein Variation Analysis

Amino acid sequence comparisons were performed using:

`scripts/compare_tem_proteins.py`

Run:

```bash
python scripts/compare_tem_proteins.py
```

The corresponding result files include:

```text
results/variants/TEM_amino_acid_variations.csv
results/variants/TEM_protein_summary.csv
results/variants/TEM_protein_pairwise_differences.csv
```

The analysis identified:

- 286 aligned amino acid positions
- 279 conserved amino acid positions
- 7 variable amino acid positions
- 15 unique pairwise comparisons

Pairwise protein differences ranged from 1 to 6 positions.

The seven variable amino acid positions were 37, 90, 102, 162, 235, 236 and 237.

These are positions in the selected aligned protein sequences and should not automatically be interpreted as conventional Ambler beta-lactamase residue numbers.

---

## 8. Codon-Level Variation Analysis

Codon-level analysis was performed using:

`scripts/05_codon_variation_analysis.py`

Run:

```bash
python scripts/05_codon_variation_analysis.py
```

The corresponding outputs include:

```text
results/variants/TEM_codon_variation.csv
results/variants/TEM_pairwise_codon_classification.csv
```

The analysis identified 12 variable codon positions across the six selected sequences.

Five positions showed synonymous-only variation:

6, 46, 132, 158 and 239.

Seven positions showed nonsynonymous variation:

37, 90, 102, 162, 235, 236 and 237.

For example, the codons TTT and TTC both encode phenylalanine. Their difference is therefore synonymous.

The classification describes variation across the selected sequences and does not establish the functional effect of individual substitutions.

---

## 9. Reproducing the Four Heatmaps

Four R scripts were developed to visualise the comparative results.

Run these commands from the repository root after the corresponding CSV input files are available.

### 9.1 Nucleotide variation heatmap

```bash
Rscript scripts/01_nucleotide_variation_heatmap.R
```

Expected output:

`results/figures/TEM_nucleotide_variation_heatmap.png`

### 9.2 Amino acid variation heatmap

```bash
Rscript scripts/02_amino_acid_variation_heatmap.R
```

Expected output:

`results/figures/TEM_amino_acid_variation_heatmap.png`

### 9.3 Pairwise nucleotide differences

```bash
Rscript scripts/03_pairwise_nucleotide_heatmap.R
```

Expected output:

`results/figures/TEM_pairwise_nucleotide_heatmap.png`

### 9.4 Pairwise protein differences

```bash
Rscript scripts/04_pairwise_protein_heatmap.R
```

Expected output:

`results/figures/TEM_pairwise_protein_heatmap.png`

Corresponding PDF figures are also available in the repository.

The heatmaps represent sequence variation and pairwise differences, not measured antibiotic resistance or enzyme activity.

---

## 10. Recommended Reproduction Procedure

The recommended procedure starts with the aligned nucleotide and protein sequences already supplied in the repository.

This avoids downloading complete bacterial genomes or repeating upstream sequence extraction when the objective is to reproduce the comparative results.

### Step 1: Prepare the software environment

Create and activate the Conda environment as described in Section 1.

Ensure R and the necessary packages are installed.

### Step 2: Confirm the input alignments

The following files should be present:

```text
results/alignments/TEM_nucleotide_aligned.fasta
results/alignments/TEM_protein_aligned.faa
```

### Step 3: Run the comparative analyses

```bash
python scripts/compare_all_tem.py
python scripts/compare_tem_proteins.py
python scripts/05_codon_variation_analysis.py
```

### Step 4: Regenerate the figures

```bash
Rscript scripts/01_nucleotide_variation_heatmap.R
Rscript scripts/02_amino_acid_variation_heatmap.R
Rscript scripts/03_pairwise_nucleotide_heatmap.R
Rscript scripts/04_pairwise_protein_heatmap.R
```

### Step 5: Compare the results

Compare the regenerated tables and figures with the reference outputs already included under:

```text
results/variants/
results/figures/
```

These commands document the intended execution order. A complete clean-environment reproduction test has not yet been performed.

---

## 11. Expected Results

Using the same selected sequences and analysis settings, the principal results should be consistent with the following:

| Measure | Expected result |
|---------|-----------------|
| TEM variants | 6 |
| Coding sequence length | 858 nucleotides |
| Protein sequence length | 286 amino acids |
| Conserved nucleotide positions | 846 |
| Variable nucleotide positions | 12 |
| Conserved amino acid positions | 279 |
| Variable amino acid positions | 7 |
| Synonymous-only variable codon positions | 5 |
| Nonsynonymous-variable codon positions | 7 |
| Unique pairwise comparisons | 15 |
| Final heatmaps | 4 |

Minor visual differences may occur between graphical devices or software versions without changing the underlying numerical results.

---

## 12. Limitations and Reproducibility Considerations

### 12.1 Representative sequences

The analysis includes one selected TEM coding sequence from each of six E. coli genome assemblies.

The findings should not be generalised to all TEM beta-lactamase variants or all E. coli isolates.

### 12.2 Reference database versions

AMRFinderPlus and BLAST results may depend on reference database versions and detection settings.

These should be recorded when repeating the upstream analysis.

### 12.3 External genome availability

The original assemblies are identified by their NCBI accessions.

Researchers should verify that the corresponding records remain available when retrieving the data.

### 12.4 Sequence numbering

The reported nucleotide, codon and amino acid positions correspond to the selected aligned sequences.

They should not automatically be treated as standardised beta-lactamase residue numbering.

### 12.5 Phenotypic interpretation

This project is a computational sequence-comparison study.

It does not include:

- Experimental antimicrobial susceptibility testing
- Minimum inhibitory concentration measurements
- Enzyme activity or kinetic assays
- Functional validation of individual substitutions

Consequently, the observed sequence differences alone cannot establish amoxicillin resistance levels or the clinical significance of individual TEM variants.

### 12.6 Reproducibility scope

The repository provides the scripts, selected input sequences, intermediate files, numerical outputs and figures needed to inspect and reproduce the main comparative analyses.

The Conda environment and documented commands support reproducibility, but the full workflow has not been independently validated from a fresh installation.

Exact software versions, external database versions and script-specific requirements should be recorded for more stringent reproducibility.

---

## 13. Data and Output Availability

The repository includes the following research materials:

| Directory | Contents |
|-----------|----------|
| `metadata/` | Assembly verification, download and TEM extraction records |
| `reference/tem/` | TEM protein reference collection |
| `scripts/` | Bash, Python and R analysis scripts |
| `results/alignments/` | Nucleotide and protein alignments |
| `results/tem_genes/` | Extracted TEM coding sequences |
| `results/proteins/` | Translated TEM proteins |
| `results/amrfinder/` | AMRFinderPlus results |
| `results/blastx/` | BLASTX verification results |
| `results/quast/` | Genome assembly quality reports |
| `results/variants/` | Comparative nucleotide, protein and codon tables |
| `results/figures/` | Four heatmaps in PNG and PDF formats |

The complete genome assembly downloads are excluded from the repository but can be retrieved using the accession numbers provided in Section 3.

---

## 14. Reuse and Citation

The scripts and analytical outputs are provided to support transparent inspection and reproducibility of this independent bioinformatics project.

Users should acknowledge the original NCBI data sources, the relevant bioinformatics tools and this repository when reusing the workflow.

Reuse permissions depend on any applicable repository licence and the terms of the underlying data sources.

Project repository:

https://github.com/chrisweb-eng/Ecoli-TEM-Beta-Lactamase-Analysis

---

## 15. Summary

This project provides a documented computational workflow for comparing nucleotide, amino acid and codon-level variation among six selected TEM beta-lactamase variants in Escherichia coli.

The recommended reproduction approach uses the published nucleotide and protein alignments to regenerate the comparative tables and four heatmaps.

An extended genome-level workflow is also documented for researchers wishing to repeat assembly retrieval, quality assessment, antimicrobial resistance gene screening, TEM sequence extraction and BLASTX verification.

The combination of documented input data, source scripts, intermediate results, comparative tables and visualisations supports transparency and reuse of the sequence-level findings.

