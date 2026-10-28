# ============================================================
# Script 04: Pairwise Amino Acid Difference Heatmap
# Project: Comparative Analysis of TEM Beta-Lactamase Variants
#
# Purpose:
# Visualise pairwise amino acid differences among six
# TEM beta-lactamase protein sequences.
#
# Input:
# results/variants/TEM_protein_pairwise_differences.csv
#
# Outputs:
# results/figures/TEM_pairwise_protein_heatmap.png
# results/figures/TEM_pairwise_protein_heatmap.pdf
#
# Run from project root:
# Rscript scripts/04_pairwise_protein_heatmap.R
# ============================================================

# 1. Load required library
library(ggplot2)

# 2. Define file paths
input_file <- "results/variants/TEM_protein_pairwise_differences.csv"
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

# 4. Import pairwise protein comparison data
pairwise_data <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Original CSV columns:
# Variant_1, Variant_2, Differences, Identity_Percent

# Check that required columns exist
required_columns <- c(
  "Variant_1",
  "Variant_2",
  "Differences"
)

if (!all(required_columns %in% names(pairwise_data))) {
  stop("Required columns are missing from the CSV.")
}

# Select only the columns needed for plotting
pairwise_data <- pairwise_data[, required_columns]

# Rename columns
names(pairwise_data) <- c(
  "Variant1",
  "Variant2",
  "Differences"
)

# Convert differences to numeric
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

# Six sequences should produce 15 unique comparisons
if (nrow(pairwise_data) != 15) {
  stop("Expected 15 unique pairwise comparisons.")
}

# Check for duplicate unordered comparisons
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

# Validate amino acid difference counts
if (
  any(pairwise_data$Differences < 0) ||
  any(pairwise_data$Differences > 286) ||
  any(pairwise_data$Differences !=
        floor(pairwise_data$Differences))
) {
  stop("Invalid amino acid difference counts.")
}

# 6. Generate reverse comparisons
reverse_data <- data.frame(
  Variant1 = pairwise_data$Variant2,
  Variant2 = pairwise_data$Variant1,
  Differences = pairwise_data$Differences
)

# 7. Add self-comparisons
# Every protein has zero differences with itself
diagonal_data <- data.frame(
  Variant1 = variant_order,
  Variant2 = variant_order,
  Differences = 0
)

# 8. Construct complete 6 x 6 comparison matrix
heatmap_data <- rbind(
  pairwise_data,
  reverse_data,
  diagonal_data
)

# Set plotting order
heatmap_data$Variant1 <- factor(
  heatmap_data$Variant1,
  levels = variant_order
)

heatmap_data$Variant2 <- factor(
  heatmap_data$Variant2,
  levels = rev(variant_order)
)

# Confirm 36 heatmap cells
stopifnot(nrow(heatmap_data) == 36)

# 9. Create pairwise protein difference heatmap
protein_pairwise_plot <- ggplot(
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
    title = "Pairwise Amino Acid Differences Among Six TEM Proteins",
    subtitle = "Number of differing positions in a 286-aa alignment",
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
    "TEM_pairwise_protein_heatmap.png"
  ),
  plot = protein_pairwise_plot,
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
    "TEM_pairwise_protein_heatmap.pdf"
  ),
  plot = protein_pairwise_plot,
  width = 9,
  height = 7,
  units = "in",
  bg = "white"
)

# 12. Print completion summary
cat("\nPairwise amino acid heatmap completed.\n")
cat("TEM variants:", length(variant_order), "\n")
cat("Unique pairwise comparisons:", nrow(pairwise_data), "\n")
cat("Heatmap cells:", nrow(heatmap_data), "\n")
cat(
  "Maximum amino acid differences:",
  max(heatmap_data$Differences),
  "\n"
)
cat("Output directory:", output_dir, "\n")
