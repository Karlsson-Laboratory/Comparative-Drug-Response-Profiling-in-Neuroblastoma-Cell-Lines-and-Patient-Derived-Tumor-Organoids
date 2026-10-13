# Figure 2B: Spearman correlations between model drug-response profiles.

output_dir <- "Figures/Figure2/Figure2B/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(pheatmap)

data_dir <- "data"

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

cell_line_cols <- c(
  "SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y",
  "SK.N.FI", "SK.N.DZ", "IMR32", "AMC772T", "NB129", "NB067", "LU_NB_1", "LU_NB_2", "LU_NB_3"
)
bm_cols <- grep("^DSS_BM", names(dss), value = TRUE)
model_cols <- c(cell_line_cols, bm_cols)

# Define MYCN, TP53 and model-type annotations.
mycn_status <- tibble(
  model = model_cols,
  MYCN = case_when(
    model %in% c("IMR32", "NB129", "NB067", "LU_NB_1", "LU_NB_2", "LU_NB_3", "SK.N.DZ", "SK.N.BE2") ~ "AMP",
    model %in% c("SK.N.SH", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "AMC772T") ~ "non-AMP",
    grepl("^DSS_BM", model) ~ "non-AMP",
    TRUE ~ NA_character_
  ),
  TP53 = case_when(
    model %in% c("SK.N.BE2", "SK.N.AS", "SK.N.FI", "SK.N.DZ") ~ "MUT",
    model %in% c("AMC772T", "NB129", "NB067", "LU_NB_1", "LU_NB_2", "LU_NB_3","SK.N.SH","SH.SY5Y","IMR32") ~ "WT",
    grepl("^DSS_BM", model) ~ "WT",
    TRUE ~ NA_character_
  ),'Model type' = case_when(
    model %in% c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ", "IMR32") ~ "cell line",
    model %in% c("AMC772T", "NB129", "NB067") ~ "nPDO",
    model %in% c("LU_NB_1", "LU_NB_2", "LU_NB_3") ~ "nPDXO",
    grepl("^DSS_BM", model) ~ "cell line",
    TRUE ~ NA_character_
  )
  )
mycn_anno <- mycn_status %>% column_to_rownames("model")

anno_colors <- list(
  MYCN = c("AMP" = "firebrick", "non-AMP" = "white"),
  'Model type' = c(
    "cell line" = "#800020",
    "nPDO"     = "#d95f02",
    "nPDXO"    = "#7570b3"
  ),
  TP53 = c(
    "MUT" = "#0052B1",
    "WT"  = "white"
  ))

# Keep drugs with DSS >= 10 in at least one included model.
dss_filtered <- dss %>%
  filter(apply(select(., all_of(model_cols)), 1, function(x) any(x >= 10, na.rm = TRUE)))

# Build the drug-by-model DSS matrix.
model_matrix <- dss_filtered %>%
  select(all_of(model_cols)) %>%
  as.matrix()

# Calculate Spearman correlations between model response profiles.
cor_matrix <- cor(model_matrix, method = "spearman", use = "pairwise.complete.obs")

# Plot the clustered correlation heatmap with model annotations.
panel_plot <- pheatmap(
  cor_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  clustering_distance_cols = "correlation",
  clustering_distance_rows = "correlation",
  main = "Spearman Correlation (Drugs with at least one DSS ≥ 10)",
  annotation_col = mycn_anno,
  annotation_row = mycn_anno,
  annotation_colors = anno_colors
)

# Save PDF and PNG figures.
ggplot2::ggsave(file.path(output_dir, "Figure2B.pdf"), plot = panel_plot$gtable,
    width = 11, height = 8, units = "in")
ggplot2::ggsave(file.path(output_dir, "Figure2B.png"), plot = panel_plot$gtable,
    width = 11, height = 8, units = "in", dpi = 300, bg = "white")
