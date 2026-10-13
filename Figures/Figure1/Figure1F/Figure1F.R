# Figure 1F: Per-model target mean DSS in neuroblastoma and controls.

library(tidyverse)
library(ggpubr)

data_dir <- "data"
output_dir <- "Figures/Figure1/Figure1F/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Define neuroblastoma and bone-marrow control models.
Neuroblastoma <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ", "AMC772T", "IMR32", "NB129", "NB067")
models_CTRL <- c("DSS_BM1", "DSS_BM2", "DSS_BM3", "DSS_BM4", "DSS_BM5")

selected_targets <- c("AURKA", "AURKB", "PLK1", "AURKC", "XPO1", "IGF1R")

# Calculate one mean DSS per model and target across target-annotated drugs.

dss_long <- dss %>% separate_rows(Gene_Targeted, sep = ";") %>%
    mutate(Gene_Targeted = trimws(Gene_Targeted)) %>% filter(Gene_Targeted %in% selected_targets) %>%
    pivot_longer(cols = all_of(c(models_CTRL, Neuroblastoma)), names_to = "Model", values_to = "DSS") %>%
    mutate(Group = ifelse(Model %in% Neuroblastoma, "NB", "CTRL")) %>% filter(!is.na(DSS))
model_means <- dss_long %>% group_by(Model, Gene_Targeted, Group) %>%
    summarise(Mean_DSS = mean(DSS, na.rm = TRUE), .groups = "drop") %>%
    mutate(Group = factor(Group, levels = c("CTRL", "NB")), Gene_Targeted = factor(Gene_Targeted, levels = selected_targets))
print(model_means %>% count(Gene_Targeted, Group))
# Compare per-model target means with unpaired two-sided Wilcoxon tests.
comparisons <- list(c("CTRL", "NB"))

panel_plot <- ggplot(model_means, aes(x = Group, y = Mean_DSS, color = Group)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
    facet_wrap(~ Gene_Targeted, ncol = 6) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", hide.ns = FALSE, step.increase = 0.2) +
    scale_color_manual(values = c("CTRL" = "grey40", "NB" = "steelblue", "2D" = "steelblue", "3D" = "darkorange")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    theme_classic(base_size = 14) + theme(legend.position = "none") + labs(x = "Model", y = "Mean DSS")

print(panel_plot)
# Save PDF and PNG figures.
ggsave(file.path(output_dir, "Figure1F.pdf"), plot = panel_plot,
    width = 14, height = 5, units = "in")
ggsave(file.path(output_dir, "Figure1F.png"), plot = panel_plot,
    width = 14, height = 5, units = "in", dpi = 300, bg = "white")
