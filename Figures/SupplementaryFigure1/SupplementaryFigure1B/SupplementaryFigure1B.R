# Supplementary Figure 1B: Number of drugs annotated to each target.

output_dir <- "Figures/SupplementaryFigure1/SupplementaryFigure1B/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(RColorBrewer)

data_dir <- "data"

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Count drugs per target and keep targets with more than five drugs.
targets <- dss %>%
    separate_rows(Gene_Targeted, sep = ";") %>%
    mutate(Gene_Targeted = trimws(Gene_Targeted)) %>%
    filter(!is.na(Gene_Targeted) & Gene_Targeted != "" & Gene_Targeted != "Na") %>%
    count(Gene_Targeted) %>%
    filter(n > 5) %>%
    arrange(n, Gene_Targeted)  # Order equal counts by target name.

# Define target colors.
n_pal <- length(unique(targets$n))
palette_blues <- brewer.pal(min(n_pal, 9), "BuPu")
# Interpolate the palette when there are more than nine targets.
if (n_pal > 9) palette_blues <- colorRampPalette(brewer.pal(9, "BuPu"))(n_pal)

targets$color_group <- as.factor(targets$n)

panel_plot <- ggplot(targets, aes(x = reorder(Gene_Targeted, n), y = n, color = color_group)) + geom_segment(aes(xend = reorder(Gene_Targeted,
    n), y = 0, yend = n), color = "#cccccc", size = 1) + geom_point(size = 4, alpha = 0.85) + scale_color_manual(values = palette_blues) +
    coord_flip() + labs(x = "Target", y = "Number of Drugs", color = "Number of Drugs", title = "Targets with >5 Drugs (Colored by Group)") +
    theme_classic(base_size = 13) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.y = element_text(size = 11,
    color = "grey30"), axis.text.x = element_text(size = 11), legend.position = "right")

print(panel_plot)

# Save PDF and PNG figures.
ggplot2::ggsave(file.path(output_dir, "SupplementaryFigure1B.pdf"), plot = panel_plot,
    width = 11, height = 8, units = "in")
ggplot2::ggsave(file.path(output_dir, "SupplementaryFigure1B.png"), plot = panel_plot,
    width = 11, height = 8, units = "in", dpi = 300, bg = "white")
