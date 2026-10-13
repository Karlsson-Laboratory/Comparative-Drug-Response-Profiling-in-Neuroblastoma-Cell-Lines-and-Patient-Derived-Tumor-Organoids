# Figure 2A: FO5 and Lund/FO4 drug-library overlap using hardcoded counts.

library(eulerr)
library(ggplot2)

data_dir <- "data"
output_dir <- "Figures/Figure2/Figure2A/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fit <- euler(c("FO5" = 97, "Lund/FO4" = 94, "FO5&Lund/FO4" = 431))
panel_plot <- plot(fit, quantities = TRUE)
print(panel_plot)
# Save PDF and PNG figures.
ggsave(file.path(output_dir, "Figure2A.pdf"), plot = panel_plot,
    width = 8, height = 6, units = "in")
ggsave(file.path(output_dir, "Figure2A.png"), plot = panel_plot,
    width = 8, height = 6, units = "in", dpi = 300, bg = "white")
