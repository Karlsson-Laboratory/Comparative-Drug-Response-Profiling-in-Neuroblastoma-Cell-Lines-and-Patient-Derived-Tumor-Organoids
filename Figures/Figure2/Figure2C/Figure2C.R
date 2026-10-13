# Figure 2C: Within-topology pairwise model correlations.

output_dir <- "Figures/Figure2/Figure2C/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggpubr)

data_dir <- "data"

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

cell_line_cols <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ", "IMR32", "AMC772T",
    "NB129", "NB067", "LU_NB_1", "LU_NB_2", "LU_NB_3")
bm_cols <- grep("^DSS_BM", names(dss), value = TRUE)
model_cols <- c(cell_line_cols, bm_cols)

# Keep drugs with DSS >= 10 in at least one included model.
dss_filtered <- dss %>%
    filter(apply(select(., all_of(model_cols)), 1, function(x) any(x >= 10, na.rm = TRUE)))

# Build the drug-by-model DSS matrix.
model_matrix <- dss_filtered %>%
    select(all_of(model_cols)) %>%
    as.matrix()

# Calculate Spearman correlations between model response profiles.
cor_matrix <- cor(model_matrix, method = "spearman", use = "pairwise.complete.obs")

# Define model groups.
models_2D <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ", "AMC772T")
models_3D <- c("IMR32", "NB129", "NB067", "LU_NB_1", "LU_NB_2", "LU_NB_3")

get_pairwise_corr <- function(models, cor_mat) {
    # Extract each unique model pair once.
    combs <- combn(models, 2, simplify = FALSE)
    tibble(model1 = sapply(combs, `[`, 1), model2 = sapply(combs, `[`, 2), correlation = sapply(combs,
        function(pair) cor_mat[pair[1], pair[2]]))
}

corr_2D <- get_pairwise_corr(models_2D, cor_matrix) %>%
    mutate(group = "2D")
corr_3D <- get_pairwise_corr(models_3D, cor_matrix) %>%
    mutate(group = "3D")
corr_long <- bind_rows(corr_2D, corr_3D)

# Compare within-topology correlations between the 2D and 3D groups.

panel_plot <- ggplot(corr_long, aes(x = group, y = correlation, fill = group)) + geom_boxplot(width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.8) + stat_compare_means(comparisons = list(c("2D",
    "3D")), method = "wilcox.test", label = "p.signif") + scale_fill_manual(values = c(`2D` = "steelblue",
    `3D` = "darkorange")) + labs(x = "Group", y = "Spearman Correlation", title = "Pairwise Spearman Correlations: 2D vs 3D Models") +
    theme_classic(base_size = 15) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

print(panel_plot)

# Save PDF and PNG figures.
ggplot2::ggsave(file.path(output_dir, "Figure2C.pdf"), plot = panel_plot,
    width = 11, height = 8, units = "in")
ggplot2::ggsave(file.path(output_dir, "Figure2C.png"), plot = panel_plot,
    width = 11, height = 8, units = "in", dpi = 300, bg = "white")
