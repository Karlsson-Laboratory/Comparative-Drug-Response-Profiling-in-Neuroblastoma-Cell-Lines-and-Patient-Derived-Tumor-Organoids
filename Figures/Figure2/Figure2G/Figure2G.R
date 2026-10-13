# Figure 2G: DSS heatmap ranked by 3D-minus-control response.

output_dir <- "Figures/Figure2/Figure2G/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(pheatmap)

data_dir <- "data"

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Define control, 2D and 3D model columns.
bm_cols    <- c("DSS_BM1","DSS_BM2","DSS_BM3","DSS_BM4","DSS_BM5")
models_2D  <- c("SK.N.SH","SK.N.BE2","SK.N.AS","SH.SY5Y","SK.N.FI","SK.N.DZ","AMC772T")
models_3D  <- c("IMR32","NB129","NB067","LU_NB_1","LU_NB_2","LU_NB_3")
col_order  <- c(bm_cols, models_2D, models_3D)

# Use the first semicolon-separated target as the primary target.
first_target <- function(x) sapply(strsplit(x, ";"), function(xx) trimws(xx[1]))

# Standardize the Sepantronium bromide name.
dss$Drug_Name <- gsub("Sepantronium_bromide", "Sepan_bromide", dss$Drug_Name)

# Rank drugs by mean 3D-minus-control DSS and keep the top 20 primary targets.
dss_rank <- dss %>%
  mutate(
    mean_3D       = rowMeans(select(., all_of(models_3D)), na.rm = TRUE),
    mean_BM       = rowMeans(select(., all_of(bm_cols)),   na.rm = TRUE),
    diff_3D_BM    = mean_3D - mean_BM,
    Primary_Target = first_target(Gene_Targeted)
  ) %>%
  filter(!is.na(Primary_Target), Primary_Target != "", Primary_Target != "Na")

top_unique <- dss_rank %>%
  group_by(Primary_Target) %>%
  arrange(desc(diff_3D_BM)) %>%
  slice(1) %>%               # Keep the highest-ranked drug for each primary target.
  ungroup() %>%
  arrange(desc(diff_3D_BM)) %>%
  slice_head(n = 20)

# Build the DSS matrix in control, 2D and 3D order.
top_drug_order <- top_unique$Drug_Name
heat_mat <- top_unique %>%
  select(Drug_Name, all_of(col_order)) %>%
  mutate(Drug_Name = factor(Drug_Name, levels = top_drug_order)) %>%
  arrange(Drug_Name) %>%
  column_to_rownames("Drug_Name") %>%
  as.matrix()

# Cap displayed DSS values at 30.
heat_mat[heat_mat > 30] <- 30

# Label rows as Drug [Target].
rownames(heat_mat) <- paste0(
  rownames(heat_mat), " [",
  top_unique$Primary_Target[match(rownames(heat_mat), top_unique$Drug_Name)],
  "]"
)

# Plot the heatmap without clustering.
panel_plot <- pheatmap(
  heat_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  na_col = "#CCCCCC",
  main = "Top 20 by (3D − BM), unique targets | DSS (BM + 2D + 3D)",
  fontsize_row = 10
)

# Save PDF and PNG figures.
ggplot2::ggsave(file.path(output_dir, "Figure2G.pdf"), plot = panel_plot$gtable,
    width = 11, height = 8, units = "in")
ggplot2::ggsave(file.path(output_dir, "Figure2G.png"), plot = panel_plot$gtable,
    width = 11, height = 8, units = "in", dpi = 300, bg = "white")
