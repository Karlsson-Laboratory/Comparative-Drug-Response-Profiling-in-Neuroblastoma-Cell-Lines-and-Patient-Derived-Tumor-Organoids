# Supplementary Figure 4E: Paired pre-formed and non-pre-formed drug responses.

library(tidyverse)
library(ggrepel)

data_dir <- "data"
output_dir <- "Figures/SupplementaryFigure4/SupplementaryFigure4E/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load drug sensitivity scores.
dss_all <- read.delim(file.path(data_dir, "focused_28_drug_integrated_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

effective_2d <- c("Dasatinib", "LY3009120", "Trametinib", "Cilengitide", "AZD8055",
    "Saracatinib", "Vistusertib", "Ulixertinib", "Binimetinib", "Sabutoclax")
effective_3d <- c("Vinblastine", "LY-2874455", "Etoposide", "Vincristine", "Sepantronium bromide",
    "GSK-1070916", "Danusertib", "AZD1152-HQPA", "Topotecan", "OSU-03012")
effective_both <- c("E7820", "Omipalisib", "Alisertib", "MK-8745", "Vinorelbine",
    "BI 2536", "Rigosertib", "BMS-754807", "GSK-461364")

dss_all <- dss_all %>% mutate(
    drug_clean = gsub("_", " ", DRUG_NAME),
    category = case_when(drug_clean %in% effective_2d ~ "Effective in 2D",
        drug_clean %in% effective_3d ~ "Effective in 3D", drug_clean %in% effective_both ~ "Effective in both", TRUE ~ "Other"),
    category = factor(category, levels = c("Effective in 2D", "Effective in 3D", "Effective in both", "Other"))
)
colors <- c("Effective in 2D" = "steelblue", "Effective in 3D" = "darkorange", "Effective in both" = "darkgreen", "Other" = "grey60")
shapes <- c("Effective in 2D" = 16, "Effective in 3D" = 17, "Effective in both" = 15, "Other" = 1)

# Define paired model columns for complete-case Pearson correlations.
# Primary-screen columns are non-pre-formed references; _pre columns are pre-formed values.
pairs <- tibble(
    model = c("NB067_primary_reference", "NB129_primary_reference", "SK_N_DZ_primary_reference"),
    x_col = c("NB067", "NB129", "SK.N.DZ"),
    y_col = c("NB067_pre", "NB129_pre", "SK.N.DZ_pre")
)

# Use paired nonmissing DSS values; label absolute differences > 10.
for (i in seq_len(nrow(pairs))) {
    col_x <- pairs$x_col[i]
    col_y <- pairs$y_col[i]
    df <- dss_all %>% filter(!is.na(.data[[col_x]]), !is.na(.data[[col_y]])) %>%
        mutate(diff = .data[[col_y]] - .data[[col_x]], label_outlier = abs(diff) > 10)
    correlation <- cor.test(df[[col_x]], df[[col_y]], method = "pearson")
    print(tibble(model = pairs$model[i], n = nrow(df), R = unname(correlation$estimate), p = correlation$p.value))
    panel_plot <- ggplot(df, aes(x = .data[[col_x]], y = .data[[col_y]])) +
        geom_point(aes(color = category, shape = category), size = 2.5, alpha = 0.85) +
        geom_abline(slope = 1, intercept = 0, color = "black") +
        geom_abline(slope = 1, intercept = 10, linetype = "dotted", color = "red") +
        geom_abline(slope = 1, intercept = -10, linetype = "dotted", color = "blue") +
        geom_text_repel(data = df %>% filter(label_outlier), aes(label = DRUG_NAME, color = category),
            size = 3.5, fontface = "bold", seed = 123, show.legend = FALSE) +
        scale_color_manual(values = colors, name = "Category") + scale_shape_manual(values = shapes, name = "Category") +
        theme_classic(base_size = 14) + labs(x = "Non-pre-formed", y = "Pre-formed", title = pairs$model[i],
            subtitle = sprintf("Pearson R = %.2f, p = %.2g", correlation$estimate, correlation$p.value))
    print(panel_plot)
    # Save PDF and PNG figures for each model pair.
    filename <- paste0("SupplementaryFigure4E_", pairs$model[i], "_categories")
    ggsave(file.path(output_dir, paste0(filename, ".pdf")), plot = panel_plot, width = 8, height = 6, units = "in")
    ggsave(file.path(output_dir, paste0(filename, ".png")), plot = panel_plot, width = 8, height = 6, units = "in", dpi = 300, bg = "white")
}
