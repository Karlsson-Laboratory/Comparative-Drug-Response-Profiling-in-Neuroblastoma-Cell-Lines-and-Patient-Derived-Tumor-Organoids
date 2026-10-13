# Supplementary Figure 1E: Clinical-drug responses in neuroblastoma and controls.

library(tidyverse)
library(ggpubr)

data_dir <- "data"
output_dir <- "Figures/SupplementaryFigure1/SupplementaryFigure1E/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_primary_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Define neuroblastoma and bone-marrow control models.
Neuroblastoma <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ", "AMC772T", "IMR32", "NB129", "NB067")
models_CTRL <- c("DSS_BM1", "DSS_BM2", "DSS_BM3", "DSS_BM4", "DSS_BM5")

selected_drugs <- c("Carboplatin", "Cisplatin", "Doxorubicin", "Etoposide", "Topotecan", "Vincristine")

dss_long <- dss %>%
    filter(Drug_Name %in% selected_drugs) %>%
    pivot_longer(cols = all_of(c(models_CTRL, Neuroblastoma)), names_to = "Model", values_to = "DSS") %>%
    mutate(Group = factor(ifelse(Model %in% Neuroblastoma, "NB", "CTRL"), levels = c("CTRL", "NB")),
        Drug_Name = factor(Drug_Name, levels = selected_drugs))

# Use unpaired two-sided Wilcoxon comparisons with default tie handling.
comparisons <- list(c("CTRL", "NB"))
panel_plot <- ggplot(dss_long, aes(x = Group, y = DSS, color = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    facet_wrap(~ Drug_Name, ncol = 3) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", hide.ns = FALSE) +
    scale_color_manual(values = c("CTRL" = "grey40", "NB" = "steelblue"), name = "Group") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    theme_classic(base_size = 14) + labs(x = "Group", y = "DSS")

print(panel_plot)
# Save PDF and PNG figures.
ggsave(file.path(output_dir, "SupplementaryFigure1E.pdf"), plot = panel_plot,
    width = 9, height = 7, units = "in")
ggsave(file.path(output_dir, "SupplementaryFigure1E.png"), plot = panel_plot,
    width = 9, height = 7, units = "in", dpi = 300, bg = "white")
