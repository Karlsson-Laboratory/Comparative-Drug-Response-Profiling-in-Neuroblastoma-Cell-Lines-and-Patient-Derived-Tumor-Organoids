# Supplementary Figure 4C: Drug-specific mean DSS in NBOM and RPMI model groups.

library(tidyverse)
library(ggrepel)

data_dir <- "data"
output_dir <- "Figures/SupplementaryFigure4/SupplementaryFigure4C/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"),
    stringsAsFactors = FALSE)

# Medium groups contain different models; means are not paired medium-switch measurements.
models_RPMI <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ")
models_NBOM <- c("AMC772T", "NB129", "NB067", "LU_NB_1", "LU_NB_2", "LU_NB_3")
# Define drug target labels and families for plotting.
drug_labels <- tribble(
    ~Drug_Name, ~Target_Label, ~Family,
    "BMS-754807", "IGF1R", "IGF1R",
    "NVP-AEW541", "IGF1R", "IGF1R",
    "Linsitinib", "IGF1R/INSR", "IGF1R",
    "LY-2874455", "FGFR", "FGFR",
    "AZD4547", "FGFR", "FGFR",
    "Infigratinib", "FGFR", "FGFR",
    "Gefitinib", "EGFR", "EGFR/ERBB",
    "Erlotinib", "EGFR", "EGFR/ERBB",
    "Afatinib", "ERBB", "EGFR/ERBB"
)

# Calculate mean DSS across the models in each medium group.
drug_means <- dss %>% filter(Drug_Name %in% drug_labels$Drug_Name) %>%
    mutate(NBOM = rowMeans(select(., all_of(models_NBOM)), na.rm = TRUE),
        RPMI = rowMeans(select(., all_of(models_RPMI)), na.rm = TRUE)) %>%
    select(Drug_Name, NBOM, RPMI) %>% left_join(drug_labels, by = "Drug_Name")
print(drug_means)
plot_data <- drug_means %>% pivot_longer(cols = c(NBOM, RPMI),
    names_to = "Medium", values_to = "Mean_DSS") %>%
    mutate(Medium = factor(Medium, levels = c("NBOM", "RPMI")),
        Label = paste0(Drug_Name, " [", Target_Label, "]"))

panel_plot <- ggplot(plot_data, aes(x = Medium, y = Mean_DSS, color = Family)) +
    geom_line(aes(group = Drug_Name), linewidth = 0.7) + geom_point(size = 3) +
    geom_text_repel(data = plot_data %>% filter(Medium == "NBOM"), aes(label = Label),
        size = 3.5, seed = 123, direction = "y", nudge_x = -0.05, max.overlaps = Inf,
        show.legend = FALSE) +
    scale_color_manual(values = c("IGF1R" = "steelblue", "FGFR" = "darkorange", "EGFR/ERBB" = "red")) +
    theme_bw(base_size = 14) + labs(x = NULL, y = "Mean DSS", color = NULL)
print(panel_plot)
# Save PDF and PNG figures.
ggsave(file.path(output_dir, "SupplementaryFigure4C.pdf"), plot = panel_plot,
    width = 9, height = 8, units = "in")
ggsave(file.path(output_dir, "SupplementaryFigure4C.png"), plot = panel_plot,
    width = 9, height = 8, units = "in", dpi = 300, bg = "white")
