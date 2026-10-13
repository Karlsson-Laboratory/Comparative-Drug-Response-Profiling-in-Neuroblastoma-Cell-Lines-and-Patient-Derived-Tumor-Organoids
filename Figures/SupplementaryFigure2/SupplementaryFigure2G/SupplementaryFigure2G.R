# Supplementary Figure 2G: Per-model mean negative AUC by topology.

library(tidyverse)
library(ggpubr)
library(ggrepel)

data_dir <- "data"
output_dir <- "Figures/SupplementaryFigure2/SupplementaryFigure2G/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load PMC drug-response AUC tables.
auc_tumoroids <- read.delim(file.path(data_dir, "pmc_tumoroid_drug_response_auc.tsv"),
    stringsAsFactors = FALSE, na.strings = c("NA", "na", "Na", "N/A"))
auc_cell_lines <- read.delim(file.path(data_dir, "pmc_cell_line_drug_response_auc.tsv"),
    stringsAsFactors = FALSE, na.strings = c("NA", "na", "Na", "N/A"))

# Use numeric model columns; assign IMR32 to 3D and AMC772 to 2D.
tumoroid_cols <- names(auc_tumoroids)[vapply(auc_tumoroids, is.numeric, logical(1))]
cell_line_cols <- names(auc_cell_lines)[vapply(auc_cell_lines, is.numeric, logical(1))]
models_3D <- c(setdiff(tumoroid_cols, "amc772"), "imr32")
models_2D <- c(setdiff(cell_line_cols, "imr32"), "amc772")

# Calculate each model’s mean AUC over its available screened drugs.
# Omit missing values without imputing or restricting to the FO5 overlap.
model_means <- bind_rows(
    auc_tumoroids %>% summarise(across(all_of(tumoroid_cols), ~ mean(.x, na.rm = TRUE))) %>%
        pivot_longer(everything(), names_to = "Model", values_to = "Mean_AUC"),
    auc_cell_lines %>% summarise(across(all_of(cell_line_cols), ~ mean(.x, na.rm = TRUE))) %>%
        pivot_longer(everything(), names_to = "Model", values_to = "Mean_AUC")
) %>% mutate(Mean_negative_AUC = -Mean_AUC,
    Topology = factor(ifelse(Model %in% models_2D, "2D", "3D"), levels = c("2D", "3D")))
comparison <- wilcox.test(Mean_negative_AUC ~ Topology, data = model_means)
print(comparison)

panel_plot <- ggplot(model_means, aes(x = Topology, y = Mean_negative_AUC, color = Topology)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    geom_text_repel(aes(label = Model), size = 3.5, min.segment.length = 0) +
    scale_color_manual(values = c("2D" = "steelblue", "3D" = "darkorange")) +
    stat_compare_means(method = "wilcox.test") +
    theme_classic(base_size = 14) + theme(legend.position = "none") + labs(x = NULL, y = "Mean -AUC")

print(panel_plot)
# Save PDF and PNG figures.
ggsave(file.path(output_dir, "SupplementaryFigure2G.pdf"), plot = panel_plot,
    width = 7, height = 8, units = "in")
ggsave(file.path(output_dir, "SupplementaryFigure2G.png"), plot = panel_plot,
    width = 7, height = 8, units = "in", dpi = 300, bg = "white")
