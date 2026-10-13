# Supplementary Figure 2D: Drugs with higher responses in 2D than controls and 3D.

library(tidyverse)
library(ggpubr)

data_dir <- "data"
output_dir <- "Figures/SupplementaryFigure2/SupplementaryFigure2D/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Define model groups.
models_2D <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ", "AMC772T")
models_3D <- c("IMR32", "NB129", "NB067", "LU_NB_1", "LU_NB_2", "LU_NB_3")
models_CTRL <- c("DSS_BM1", "DSS_BM2", "DSS_BM3", "DSS_BM4", "DSS_BM5")

# Calculate per-drug group means and unpaired two-sided Wilcoxon p-values.

drug_statistics <- dss %>% rowwise() %>%
    mutate(mean_CTRL = mean(c_across(all_of(models_CTRL)), na.rm = TRUE),
        mean_2D = mean(c_across(all_of(models_2D)), na.rm = TRUE),
        mean_3D = mean(c_across(all_of(models_3D)), na.rm = TRUE),
        p_CTRL_2D = wilcox.test(c_across(all_of(models_CTRL)), c_across(all_of(models_2D)))$p.value,
        p_2D_3D = wilcox.test(c_across(all_of(models_2D)), c_across(all_of(models_3D)))$p.value) %>%
    ungroup()
# Keep drugs with higher mean 2D DSS than both CTRL and 3D, with both p-values < 0.05.
selected <- drug_statistics %>%
    filter(mean_2D > mean_CTRL, mean_2D > mean_3D, p_CTRL_2D < 0.05, p_2D_3D < 0.05)
print(selected %>% select(Drug_Name, mean_CTRL, mean_2D, mean_3D, p_CTRL_2D, p_2D_3D))

# Set the drug display order and append other qualifying compounds.
drug_order <- c("MST-312", "8-amino-adenosine", "Enzastaurin", "AZD8055", "Vistusertib", "PD0325901")
drug_order <- c(intersect(drug_order, selected$Drug_Name), setdiff(selected$Drug_Name, drug_order))
dss_long <- selected %>%
    pivot_longer(cols = all_of(c(models_CTRL, models_2D, models_3D)), names_to = "Model", values_to = "DSS") %>%
    mutate(Group = case_when(Model %in% models_CTRL ~ "CTRL", Model %in% models_2D ~ "2D", Model %in% models_3D ~ "3D"),
        Group = factor(Group, levels = c("CTRL", "2D", "3D")), Drug_Name = factor(Drug_Name, levels = drug_order))

# Compare CTRL versus 2D and 2D versus 3D without multiple-testing correction.
comparisons <- list(c("CTRL", "2D"), c("2D", "3D"))
panel_plot <- ggplot(dss_long, aes(x = Group, y = DSS, color = Group)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.7) +
    facet_wrap(~ Drug_Name, ncol = 3) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", hide.ns = FALSE) +
    scale_color_manual(values = c("CTRL" = "grey40", "2D" = "steelblue", "3D" = "darkorange")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
    theme_classic(base_size = 14) + theme(legend.position = "none") + labs(x = "Group", y = "DSS")

print(panel_plot)
# Save PDF and PNG figures.
ggsave(file.path(output_dir, "SupplementaryFigure2D.pdf"), plot = panel_plot,
    width = 9, height = 7, units = "in")
ggsave(file.path(output_dir, "SupplementaryFigure2D.png"), plot = panel_plot,
    width = 9, height = 7, units = "in", dpi = 300, bg = "white")
