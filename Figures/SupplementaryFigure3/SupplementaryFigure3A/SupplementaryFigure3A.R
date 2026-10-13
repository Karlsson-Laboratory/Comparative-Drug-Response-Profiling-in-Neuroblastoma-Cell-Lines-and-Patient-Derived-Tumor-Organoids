# Supplementary Figure 3A: Mean 2D and 3D DSS with drug-category annotations.

library(tidyverse)
library(ggrepel)

data_dir <- "data"
output_dir <- "Figures/SupplementaryFigure3/SupplementaryFigure3A/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Define model groups.
models_2D <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ", "AMC772T")
models_3D <- c("IMR32", "NB129", "NB067", "LU_NB_1", "LU_NB_2", "LU_NB_3")
models_CTRL <- c("DSS_BM1", "DSS_BM2", "DSS_BM3", "DSS_BM4", "DSS_BM5")

effective_2d <- c("Dasatinib", "LY3009120", "Trametinib", "Cilengitide", "AZD8055",
    "Saracatinib", "Vistusertib", "Ulixertinib", "Binimetinib", "Sabutoclax")
effective_3d <- c("Vinblastine", "LY-2874455", "Etoposide", "Vincristine", "Sepantronium bromide",
    "GSK-1070916", "Danusertib", "AZD1152-HQPA", "Topotecan", "OSU-03012")
effective_both <- c("E7820", "Omipalisib", "Alisertib", "MK-8745", "Vinorelbine",
    "BI 2536", "Rigosertib", "BMS-754807", "GSK-461364")

# Calculate group means and assign the specified drug categories.
df <- dss %>% mutate(mean_2D = rowMeans(select(., all_of(models_2D)), na.rm = TRUE),
    mean_3D = rowMeans(select(., all_of(models_3D)), na.rm = TRUE),
    drug_clean = gsub("_", " ", Drug_Name),
    category = case_when(drug_clean %in% effective_2d ~ "Effective in 2D",
        drug_clean %in% effective_3d ~ "Effective in 3D", drug_clean %in% effective_both ~ "Effective in both", TRUE ~ "Other")) %>%
    filter(is.finite(mean_2D), is.finite(mean_3D))
comparison <- cor.test(df$mean_2D, df$mean_3D, method = "spearman")
print(comparison)

panel_plot <- ggplot(df, aes(x = mean_2D, y = mean_3D)) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    geom_point(aes(color = category, shape = category), size = 2, alpha = 0.8) +
    geom_text_repel(data = df %>% filter(category != "Other"), aes(label = drug_clean, color = category),
        size = 3, seed = 123, show.legend = FALSE, max.overlaps = Inf) +
    scale_color_manual(values = c("Effective in 2D" = "steelblue", "Effective in 3D" = "darkorange", "Effective in both" = "darkgreen", "Other" = "grey75")) +
    scale_shape_manual(values = c("Effective in 2D" = 16, "Effective in 3D" = 17, "Effective in both" = 15, "Other" = 1)) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
        label = sprintf("Spearman rho = %.2f, p = %.2g", comparison$estimate, comparison$p.value), fontface = "bold") +
    theme_classic(base_size = 14) + labs(x = "2D", y = "3D", title = "Mean DSS", color = "Category", shape = "Category")

print(panel_plot)
# Save PDF and PNG figures.
ggsave(file.path(output_dir, "SupplementaryFigure3A.pdf"), plot = panel_plot,
    width = 11, height = 8, units = "in")
ggsave(file.path(output_dir, "SupplementaryFigure3A.png"), plot = panel_plot,
    width = 11, height = 8, units = "in", dpi = 300, bg = "white")
