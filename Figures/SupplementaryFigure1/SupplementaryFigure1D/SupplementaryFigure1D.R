# Supplementary Figure 1D: ALK inhibitor responses with ALK status annotations.

library(tidyverse)
library(ggpubr)

data_dir <- "data"
output_dir <- "Figures/SupplementaryFigure1/SupplementaryFigure1D/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_primary_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Define neuroblastoma and bone-marrow control models.
Neuroblastoma <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ", "AMC772T", "IMR32", "NB129", "NB067")
models_CTRL <- c("DSS_BM1", "DSS_BM2", "DSS_BM3", "DSS_BM4", "DSS_BM5")

selected_drugs <- c("Alectinib", "ASP3026", "Brigatinib", "Ceritinib", "Crizotinib")

dss_long <- dss %>%
    filter(Drug_Name %in% selected_drugs) %>%
    pivot_longer(cols = all_of(c(models_CTRL, Neuroblastoma)), names_to = "Model", values_to = "DSS") %>%
    mutate(Group = factor(ifelse(Model %in% Neuroblastoma, "NB", "CTRL"), levels = c("CTRL", "NB")),
        Drug_Name = factor(Drug_Name, levels = selected_drugs))

# Use unpaired two-sided Wilcoxon comparisons with default tie handling.
comparisons <- list(c("CTRL", "NB"))

# ALK snv: SK.N.SH, SH.SY5Y and NB129; amplification: IMR32.
# Both alteration types share the displayed category ALK mut.
alk_snv <- c("SK.N.SH", "SH.SY5Y", "NB129")
alk_amplified <- "IMR32"
alk_altered <- c(alk_snv, alk_amplified)
stopifnot(all(alk_altered %in% names(dss)))
dss_long <- dss_long %>% mutate(ALK_status = ifelse(Model %in% alk_altered,
    "ALK mut", "WT/unknown"))
print(dss_long %>% count(Drug_Name, ALK_status))

panel_plot <- ggplot(dss_long, aes(x = Group, y = DSS)) +
    geom_boxplot(aes(color = Group), outlier.shape = NA) +

    # Draw white markers first; dodge red markers by model to separate tied DSS values.
    geom_point(data = dss_long %>% filter(ALK_status == "WT/unknown"), aes(fill = ALK_status),
        position = position_jitter(width = 0.2, height = 0, seed = 123),
        shape = 21, color = "grey30", size = 2.5) +
    geom_point(data = dss_long %>% filter(ALK_status == "ALK mut"), aes(fill = ALK_status, group = Model),
        position = position_dodge(width = 0.6), shape = 21, color = "grey30", size = 2.5) +
    scale_fill_manual(values = c("ALK mut" = "firebrick", "WT/unknown" = "white"), name = "ALK status") +
    facet_wrap(~ Drug_Name, ncol = 3) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", hide.ns = FALSE) +
    scale_color_manual(values = c("CTRL" = "grey40", "NB" = "steelblue"), name = "Group") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    theme_classic(base_size = 14) + labs(x = "Group", y = "DSS")

print(panel_plot)
# Save PDF and PNG figures.
ggsave(file.path(output_dir, "SupplementaryFigure1D.pdf"), plot = panel_plot,
    width = 9, height = 7, units = "in")
ggsave(file.path(output_dir, "SupplementaryFigure1D.png"), plot = panel_plot,
    width = 9, height = 7, units = "in", dpi = 300, bg = "white")
