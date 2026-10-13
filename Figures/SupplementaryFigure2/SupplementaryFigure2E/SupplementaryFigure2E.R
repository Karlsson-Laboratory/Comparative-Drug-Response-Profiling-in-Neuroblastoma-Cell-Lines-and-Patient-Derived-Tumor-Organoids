# Supplementary Figure 2E: Per-model mean DSS by topology.

output_dir <- "Figures/SupplementaryFigure2/SupplementaryFigure2E/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(ggrepel)
library(ggpubr)

data_dir <- "data"

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Define model groups.
models_2D <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ", "AMC772T")
models_3D <- c("NB129", "NB067", "LU_NB_1", "LU_NB_2", "LU_NB_3", "IMR32")

# Calculate mean DSS across drugs for each model.
model_means <- dss %>%
    select(all_of(c(models_2D, models_3D))) %>%
    summarise_all(mean, na.rm = TRUE) %>%
    pivot_longer(cols = everything(), names_to = "Model", values_to = "Mean_DSS") %>%
    mutate(Topology = case_when(Model %in% models_2D ~ "2D", Model %in% models_3D ~ "3D"))

# Plot per-model means by topology.
panel_plot <- ggplot(model_means, aes(x = Topology, y = Mean_DSS, color = Topology)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) + geom_text_repel(aes(label = Model), size = 4, min.segment.length = 0) +
    scale_color_manual(values = c(`2D` = "steelblue", `3D` = "darkorange")) + labs(title = "Mean DSS Score by Topology (Each dot = one model)",
    y = "Mean DSS (across all drugs)", x = NULL) + theme_classic(base_size = 14) + ggpubr::stat_compare_means(method = "wilcox.test")

print(panel_plot)

# Save PDF and PNG figures.
ggplot2::ggsave(file.path(output_dir, "SupplementaryFigure2E.pdf"), plot = panel_plot,
    width = 11, height = 8, units = "in")
ggplot2::ggsave(file.path(output_dir, "SupplementaryFigure2E.png"), plot = panel_plot,
    width = 11, height = 8, units = "in", dpi = 300, bg = "white")
