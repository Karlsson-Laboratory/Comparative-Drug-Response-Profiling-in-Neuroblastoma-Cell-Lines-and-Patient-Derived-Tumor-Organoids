# Figure 2H: Drug-specific DSS comparisons across control, 2D and 3D models.

library(tidyverse)
library(ggpubr)

data_dir <- "data"
output_dir <- "Figures/Figure2/Figure2H/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load integrated drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"),
    stringsAsFactors = FALSE)

# Define topology groups: AMC772T is in 2D and IMR32 is in 3D.
models_2D <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y",
    "SK.N.FI", "SK.N.DZ", "AMC772T")
models_3D <- c("IMR32", "NB129", "NB067", "LU_NB_1", "LU_NB_2", "LU_NB_3")
models_CTRL <- c("DSS_BM1", "DSS_BM2", "DSS_BM3", "DSS_BM4", "DSS_BM5")

selected_drugs <- c("Carboplatin", "Cisplatin", "Doxorubicin",
    "Etoposide", "Topotecan", "Vincristine")

# Each point represents one drug-model DSS value.
dss_long <- dss %>%
    filter(Drug_Name %in% selected_drugs) %>%
    pivot_longer(cols = all_of(c(models_CTRL, models_2D, models_3D)),
        names_to = "Model", values_to = "DSS") %>%
    mutate(
        Group = case_when(
            Model %in% models_CTRL ~ "CTRL",
            Model %in% models_2D ~ "2D",
            Model %in% models_3D ~ "3D"
        ),
        Group = factor(Group, levels = c("CTRL", "2D", "3D")),
        Drug_Name = factor(Drug_Name, levels = selected_drugs)
    )

# Use unpaired two-sided Wilcoxon tests without multiple-testing correction.

# With ties, use the default asymptotic p-value and continuity correction.
comparisons <- list(c("CTRL", "2D"), c("CTRL", "3D"), c("2D", "3D"))
statistics <- dss_long %>%
    group_by(Drug_Name) %>%
    group_modify(~ tibble(
        comparison = c("CTRL vs 2D", "CTRL vs 3D", "2D vs 3D"),
        p_value = c(
            wilcox.test(.x$DSS[.x$Group == "CTRL"], .x$DSS[.x$Group == "2D"])$p.value,
            wilcox.test(.x$DSS[.x$Group == "CTRL"], .x$DSS[.x$Group == "3D"])$p.value,
            wilcox.test(.x$DSS[.x$Group == "2D"], .x$DSS[.x$Group == "3D"])$p.value
        )
    )) %>%
    ungroup()
print(statistics)

# Plot drug-specific boxplots and model points.
# Bracket heights and axis limits control the display.
panel_plot <- ggplot(dss_long, aes(x = Group, y = DSS, color = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    facet_wrap(~ Drug_Name, ncol = 3) +
    stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        label = "p.signif",
        hide.ns = FALSE,
        label.y = c(40, 55, 70),
        step.increase = 0
    ) +
    scale_color_manual(values = c("CTRL" = "grey40", "2D" = "steelblue", "3D" = "darkorange")) +
    scale_y_continuous(breaks = c(0, 25, 50, 75),
        expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(ylim = c(0, 85)) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") +
    labs(y = "DSS", x = "Topology")

print(panel_plot)

# Save PDF and PNG figures.
ggsave(file.path(output_dir, "Figure2H.pdf"), plot = panel_plot,
    width = 7.5, height = 7, units = "in")
ggsave(file.path(output_dir, "Figure2H.png"), plot = panel_plot,
    width = 7.5, height = 7, units = "in", dpi = 300, bg = "white")
