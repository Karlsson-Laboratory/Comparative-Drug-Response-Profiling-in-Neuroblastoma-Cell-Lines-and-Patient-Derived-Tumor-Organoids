# Supplementary Figure 1C: DSS comparisons for the highest-ranked drugs.

output_dir <- "Figures/SupplementaryFigure1/SupplementaryFigure1C/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(ggpubr)

data_dir <- "data"

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_primary_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Define neuroblastoma and control model groups.
Neuroblastoma <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ", "AMC772T", "IMR32",
    "NB129", "NB067")
models_CTRL <- c("DSS_BM1", "DSS_BM2", "DSS_BM3", "DSS_BM4", "DSS_BM5")

# Calculate the neuroblastoma-minus-control mean DSS for each drug.
dss <- dss %>%
    mutate(mean_NB = rowMeans(select(., all_of(Neuroblastoma)), na.rm = TRUE), mean_CTRL = rowMeans(select(.,
        all_of(models_CTRL)), na.rm = TRUE), diff_NB_CTRL = mean_NB - mean_CTRL)

# Keep the 12 drugs with the largest mean DSS differences.
top12_drugs <- dss %>%
    arrange(desc(diff_NB_CTRL)) %>%
    slice_head(n = 12) %>%
    pull(Drug_Name)

# Reshape DSS values into drug-model observations.
dss_long <- dss %>%
    filter(Drug_Name %in% top12_drugs) %>%
    pivot_longer(cols = all_of(c(models_CTRL, Neuroblastoma)), names_to = "Model", values_to = "DSS") %>%
    mutate(Group = ifelse(Model %in% Neuroblastoma, "NB", "CTRL"), Drug_Name = factor(Drug_Name, levels = top12_drugs))

# Plot one facet per drug and one point per model.
panel_plot <- ggboxplot(dss_long, x = "Group", y = "DSS", color = "Group", palette = c("grey40", "steelblue"), outlier.shape = NA,
    facet.by = "Drug_Name", add = "jitter", ylab = "DSS", xlab = NULL) + theme_classic(base_size = 14) +
    labs(title = "Top 12 Most Effective Drugs: NB vs CTRL (Each dot = model value)") + stat_compare_means(comparisons = list(c("CTRL",
    "NB")), method = "wilcox.test", label = "p.signif", hide.ns = FALSE) + scale_y_continuous(expand = expansion(mult = c(0.05,
    0.25)))

print(panel_plot)

# Save PDF and PNG figures.
ggplot2::ggsave(file.path(output_dir, "SupplementaryFigure1C.pdf"), plot = panel_plot,
    width = 11, height = 8, units = "in")
ggplot2::ggsave(file.path(output_dir, "SupplementaryFigure1C.png"), plot = panel_plot,
    width = 11, height = 8, units = "in", dpi = 300, bg = "white")
