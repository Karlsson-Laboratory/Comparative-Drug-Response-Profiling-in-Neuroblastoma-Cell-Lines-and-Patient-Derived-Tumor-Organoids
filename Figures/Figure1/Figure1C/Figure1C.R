# Figure 1C: Clustered DSS heatmap for drugs active in at least one model.

output_dir <- "Figures/Figure1/Figure1C/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(pheatmap)

data_dir <- "data"

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_primary_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Define neuroblastoma and bone-marrow control columns.
cell_line_cols <- c(
  "SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y",
  "SK.N.FI", "SK.N.DZ", "AMC772T", "IMR32", "NB129", "NB067",

  "DSS_BM1", "DSS_BM2", "DSS_BM3", "DSS_BM4", "DSS_BM5"
)

# Keep drugs with DSS > 10 in at least one included model or control.
dss_filt <- dss %>%
  filter(if_any(all_of(cell_line_cols), ~ . > 10))

# Build the drug-by-model DSS matrix.
heatmap_matrix <- dss_filt %>%
  select(Drug_Name, all_of(cell_line_cols)) %>%
  column_to_rownames("Drug_Name") %>%
  as.matrix()

# Cap displayed DSS values at 30.
heatmap_matrix[heatmap_matrix > 30] <- 30

# Cluster drugs and models; use Pearson correlation distance for columns.
panel_plot <- pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),

  breaks = seq(0, 30, length.out = 101),
  show_rownames = FALSE,
  clustering_distance_cols = 'correlation',
  main = "All Drugs (DSS > 10 in any model/BM)"
)

# Save PDF and PNG figures.
ggplot2::ggsave(file.path(output_dir, "Figure1C.pdf"), plot = panel_plot$gtable,
    width = 11, height = 8, units = "in")
ggplot2::ggsave(file.path(output_dir, "Figure1C.png"), plot = panel_plot$gtable,
    width = 11, height = 8, units = "in", dpi = 300, bg = "white")
