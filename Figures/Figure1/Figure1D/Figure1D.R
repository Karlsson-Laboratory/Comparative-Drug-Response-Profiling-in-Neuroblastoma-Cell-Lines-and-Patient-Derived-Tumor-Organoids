# Figure 1D: DSS heatmap of the 30 highest-ranked drugs.

output_dir <- "Figures/Figure1/Figure1D/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(pheatmap)

data_dir <- "data"

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_primary_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Define the model display order.
cell_line_cols <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y",
                    "SK.N.FI", "SK.N.DZ", "IMR32", "AMC772T", "NB129", "NB067")

# Identify bone-marrow control columns by the DSS_BM prefix.
bm_cols <- grep("^DSS_BM", names(dss), value = TRUE)

# Rank drugs by the input diff column and keep the top 30.
top30 <- dss %>%
  arrange(desc(diff)) %>%
  slice_head(n = 30)

# Build the DSS matrix in model-then-control order.
heatmap_matrix <- top30 %>%
  select(Drug_Name, all_of(bm_cols), all_of(cell_line_cols)) %>%
  column_to_rownames("Drug_Name") %>%
  as.matrix()

# Cap displayed DSS values at 30.
heatmap_matrix[heatmap_matrix > 30] <- 30

# Plot the heatmap without clustering.
panel_plot <- pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(0, 30, length.out = 101),

  main = "Top 30 Drugs"
)

# Save PDF and PNG figures.
ggplot2::ggsave(file.path(output_dir, "Figure1D.pdf"), plot = panel_plot$gtable,
    width = 11, height = 8, units = "in")
ggplot2::ggsave(file.path(output_dir, "Figure1D.png"), plot = panel_plot$gtable,
    width = 11, height = 8, units = "in", dpi = 300, bg = "white")
