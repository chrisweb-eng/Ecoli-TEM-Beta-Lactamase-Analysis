# ============================================================
# Script 02: Amino Acid Variation Heatmap
# Project: Comparative Analysis of TEM Beta-Lactamase Variants
#
# Purpose:
# Visualise variable amino acid positions across six
# translated TEM beta-lactamase protein sequences.
#
# Input:
# results/variants/TEM_amino_acid_variations.csv
#
# Output:
# results/figures/TEM_amino_acid_variation_heatmap.png
# results/figures/TEM_amino_acid_variation_heatmap.pdf
#
# Run from project root:
# Rscript scripts/02_amino_acid_variation_heatmap.R
# ============================================================

# 1. Load required libraries
library(ggplot2)
library(tidyr)

# 2. Define input and output paths
input_file <- "results/variants/TEM_amino_acid_variations.csv"
output_dir <- "results/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Import amino acid variation data
protein_data <- read.csv(
  input_file,
  check.names = FALSE
)

# 4. Define TEM variant order
variant_order <- c(
  "TEM-1", "TEM-2", "TEM-3",
  "TEM-10", "TEM-24", "TEM-57"
)

# Check expected columns
stopifnot(
  all(c("Position", variant_order) %in%
        names(protein_data))
)

# 5. Convert wide-format data into long format
protein_long <- pivot_longer(
  protein_data,
  cols = all_of(variant_order),
  names_to = "Variant",
  values_to = "AminoAcid"
)

# 6. Set ordering for the heatmap
protein_long$Variant <- factor(
  protein_long$Variant,
  levels = rev(variant_order)
)

protein_long$Position <- factor(
  protein_long$Position,
  levels = sort(unique(protein_data$Position))
)

# 7. Define colours for amino acids
# Each amino acid receives a consistent colour
amino_acid_colours <- c(
  "A" = "#8FB8DE",
  "D" = "#E7A2A2",
  "E" = "#F1B47D",
  "G" = "#B5D8A6",
  "K" = "#D3B1E3",
  "Q" = "#F2D27B",
  "R" = "#9BCBC8",
  "S" = "#D5B5C6",
  "T" = "#A7B9DF"
)

# 8. Create amino acid variation heatmap
protein_plot <- ggplot(
  protein_long,
  aes(
    x = Position,
    y = Variant,
    fill = AminoAcid
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = AminoAcid),
    colour = "black",
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = amino_acid_colours,
    breaks = names(amino_acid_colours)
  ) +
  labs(
    title = "Amino Acid Variation Across Six TEM Beta-Lactamase Proteins",
    subtitle = "7 variable positions identified from a 286-aa alignment",
    x = "Protein alignment position (1-based)",
    y = "TEM variant",
    fill = "Amino acid"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      size = 11,
      colour = "grey35"
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 11
    ),
    axis.text.y = element_text(
      size = 11,
      face = "bold"
    ),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(15, 15, 15, 15)
  )

# 9. Save high-resolution PNG
ggsave(
  filename = file.path(
    output_dir,
    "TEM_amino_acid_variation_heatmap.png"
  ),
  plot = protein_plot,
  width = 11,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

# 10. Save vector PDF
ggsave(
  filename = file.path(
    output_dir,
    "TEM_amino_acid_variation_heatmap.pdf"
  ),
  plot = protein_plot,
  width = 11,
  height = 5,
  units = "in",
  bg = "white"
)

# 11. Print summary
cat("\nAmino acid variation heatmap completed.\n")
cat("TEM variants:", length(variant_order), "\n")
cat("Variable amino acid positions:", nrow(protein_data), "\n")
cat("Output directory:", output_dir, "\n")
