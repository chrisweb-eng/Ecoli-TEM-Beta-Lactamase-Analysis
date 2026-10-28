# ============================================================
# Script 01: Nucleotide Variation Heatmap
# Project: Comparative Analysis of TEM Beta-Lactamase Variants
#
# Purpose:
# Visualise variable nucleotide positions across six TEM genes.
#
# Input:
# results/variants/TEM_variable_positions.csv
#
# Output:
# results/figures/TEM_nucleotide_variation_heatmap.png
# results/figures/TEM_nucleotide_variation_heatmap.pdf
#
# Run from project root:
# Rscript scripts/01_nucleotide_variation_heatmap.R
# ============================================================

# 1. Load required libraries
library(ggplot2)
library(tidyr)

# 2. Define input and output paths
input_file <- "results/variants/TEM_variable_positions.csv"
output_dir <- "results/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Import nucleotide variation data
nucleotide_data <- read.csv(
  input_file,
  check.names = FALSE
)

# 4. Define TEM variant order
variant_order <- c(
  "TEM-1", "TEM-2", "TEM-3",
  "TEM-10", "TEM-24", "TEM-57"
)

# Check that all expected columns exist
stopifnot(
  all(c("Position", variant_order) %in%
        names(nucleotide_data))
)

# 5. Convert wide-format data into long format
nucleotide_long <- pivot_longer(
  nucleotide_data,
  cols = all_of(variant_order),
  names_to = "Variant",
  values_to = "Nucleotide"
)

# 6. Set ordering for the heatmap
nucleotide_long$Variant <- factor(
  nucleotide_long$Variant,
  levels = rev(variant_order)
)

nucleotide_long$Position <- factor(
  nucleotide_long$Position,
  levels = sort(unique(nucleotide_data$Position))
)

# 7. Define nucleotide colours
nucleotide_colours <- c(
  "A" = "#70B7A6",
  "T" = "#F1C46C",
  "G" = "#85A9D8",
  "C" = "#D99CB5"
)

# 8. Create nucleotide variation heatmap
nucleotide_plot <- ggplot(
  nucleotide_long,
  aes(
    x = Position,
    y = Variant,
    fill = Nucleotide
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = Nucleotide),
    colour = "black",
    size = 4.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = nucleotide_colours,
    breaks = c("A", "T", "G", "C")
  ) +
  labs(
    title = "Nucleotide Variation Across Six TEM Beta-Lactamase Genes",
    subtitle = "12 variable positions identified from an 858-bp alignment",
    x = "Nucleotide alignment position (1-based)",
    y = "TEM variant",
    fill = "Nucleotide"
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
      size = 10
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
    "TEM_nucleotide_variation_heatmap.png"
  ),
  plot = nucleotide_plot,
  width = 12,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

# 10. Save vector PDF
ggsave(
  filename = file.path(
    output_dir,
    "TEM_nucleotide_variation_heatmap.pdf"
  ),
  plot = nucleotide_plot,
  width = 12,
  height = 5,
  units = "in",
  bg = "white"
)

# 11. Print summary
cat("\nNucleotide variation heatmap completed.\n")
cat("TEM variants:", length(variant_order), "\n")
cat("Variable nucleotide positions:", nrow(nucleotide_data), "\n")
cat("Output directory:", output_dir, "\n")
