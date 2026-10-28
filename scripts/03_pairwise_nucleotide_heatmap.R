# ============================================================
# Script 03: Pairwise Nucleotide Difference Heatmap
# Project: Comparative Analysis of TEM Beta-Lactamase Variants
#
# Purpose:
# Compare nucleotide differences among six TEM genes
# using a symmetrical 6 x 6 heatmap.
#
# Input:
# results/variants/TEM_pairwise_differences.csv
#
# Output:
# results/figures/TEM_pairwise_nucleotide_heatmap.png
# results/figures/TEM_pairwise_nucleotide_heatmap.pdf
#
# Run:
# Rscript scripts/03_pairwise_nucleotide_heatmap.R
# ============================================================

# 1. Load required library
library(ggplot2)

# 2. Define file paths
input_file <- "results/variants/TEM_pairwise_differences.csv"
output_dir <- "results/figures"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 3. Define TEM variant order
variant_order <- c(
  "TEM-1",
  "TEM-2",
  "TEM-3",
  "TEM-10",
  "TEM-24",
  "TEM-57"
)

# 4. Import pairwise nucleotide comparison data
pairwise_data <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Original CSV columns:
# Variant_1, Variant_2, Differences, Identity_Percent

# Select only columns needed for the heatmap
pairwise_data <- pairwise_data[, c(
  "Variant_1",
  "Variant_2",
  "Differences"
)]

# Rename columns for plotting
names(pairwise_data) <- c(
  "Variant1",
  "Variant2",
  "Differences"
)

# Convert differences to numeric values
pairwise_data$Differences <- as.numeric(
  pairwise_data$Differences
)

# 5. Validate imported data
if (
  anyNA(pairwise_data) ||
  !all(pairwise_data$Variant1 %in% variant_order) ||
  !all(pairwise_data$Variant2 %in% variant_order)
) {
  stop("Invalid variant names or missing difference values.")
}

# Six variants should produce 15 unique comparisons
if (nrow(pairwise_data) != 15) {
  stop("Expected 15 unique pairwise comparisons.")
}

# Check for duplicate unordered pairs
pair_keys <- vapply(
  seq_len(nrow(pairwise_data)),
  function(i) {
    paste(
      sort(c(
        pairwise_data$Variant1[i],
        pairwise_data$Variant2[i]
      )),
      collapse = "|"
    )
  },
  character(1)
)

if (
  anyDuplicated(pair_keys) ||
  any(pairwise_data$Variant1 == pairwise_data$Variant2)
) {
  stop("Duplicate pairs or unexpected self-comparisons.")
}

# Validate difference counts
if (
  any(pairwise_data$Differences < 0) ||
  any(pairwise_data$Differences > 858) ||
  any(pairwise_data$Differences !=
        floor(pairwise_data$Differences))
) {
  stop("Invalid nucleotide difference counts.")
}

# 6. Generate reverse pairwise comparisons
reverse_data <- data.frame(
  Variant1 = pairwise_data$Variant2,
  Variant2 = pairwise_data$Variant1,
  Differences = pairwise_data$Differences
)

# 7. Add self-comparisons
# Every sequence is identical to itself
diagonal_data <- data.frame(
  Variant1 = variant_order,
  Variant2 = variant_order,
  Differences = 0
)

# 8. Combine all comparisons into a 6 x 6 matrix
heatmap_data <- rbind(
  pairwise_data,
  reverse_data,
  diagonal_data
)

# Set the display order of TEM variants
heatmap_data$Variant1 <- factor(
  heatmap_data$Variant1,
  levels = variant_order
)

heatmap_data$Variant2 <- factor(
  heatmap_data$Variant2,
  levels = rev(variant_order)
)

# Confirm that the matrix contains 36 cells
stopifnot(nrow(heatmap_data) == 36)

# 9. Create pairwise nucleotide heatmap
nucleotide_pairwise_plot <- ggplot(
  heatmap_data,
  aes(
    x = Variant1,
    y = Variant2,
    fill = Differences
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 1
  ) +
  geom_text(
    aes(label = Differences),
    colour = "black",
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_gradient(
    low = "#F0F5FA",
    high = "#5084B8",
    limits = c(
      0,
      max(heatmap_data$Differences)
    )
  ) +
  coord_equal() +
  labs(
    title = "Pairwise Nucleotide Differences Among Six TEM Genes",
    subtitle = "Number of differing positions in an 858-bp alignment",
    x = "TEM variant",
    y = "TEM variant",
    fill = "Differences"
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
      size = 11
    ),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      15, 15, 15, 15
    )
  )

# 10. Save high-resolution PNG
ggsave(
  filename = file.path(
    output_dir,
    "TEM_pairwise_nucleotide_heatmap.png"
  ),
  plot = nucleotide_pairwise_plot,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

# 11. Save vector PDF
ggsave(
  filename = file.path(
    output_dir,
    "TEM_pairwise_nucleotide_heatmap.pdf"
  ),
  plot = nucleotide_pairwise_plot,
  width = 9,
  height = 7,
  units = "in",
  bg = "white"
)

# 12. Print completion summary
cat("\nPairwise nucleotide heatmap completed.\n")
cat("TEM variants:", length(variant_order), "\n")
cat("Unique pairwise comparisons:", nrow(pairwise_data), "\n")
cat("Heatmap cells:", nrow(heatmap_data), "\n")
cat(
  "Maximum nucleotide differences:",
  max(heatmap_data$Differences),
  "\n"
)
cat("Output directory:", output_dir, "\n")

