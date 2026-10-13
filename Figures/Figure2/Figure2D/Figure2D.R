# Figure 2D: Drug-class comparisons of per-drug group mean DSS.

output_dir <- "Figures/Figure2/Figure2D/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(ggpubr)
library(ggplot2)

data_dir <- "data"

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Define model groups.
models_2D <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ", "AMC772T")
models_3D <- c("IMR32", "NB129", "NB067", "LU_NB_1", "LU_NB_2", "LU_NB_3")
models_CTRL <- c("DSS_BM1", "DSS_BM2", "DSS_BM3", "DSS_BM4", "DSS_BM5")
all_models <- c(models_2D, models_3D, models_CTRL)

dss_long <- dss %>%
  select(Drug_Name, Class_Explained, all_of(all_models)) %>%
  mutate(Class_Explained = substr(Class_Explained, 4, nchar(Class_Explained))) %>%
  filter(!is.na(Class_Explained) & Class_Explained != "") %>%
  pivot_longer(
    cols = all_of(all_models),
    names_to = "Model",
    values_to = "DSS"
  ) %>%
  mutate(
    Group = case_when(
      Model %in% models_2D ~ "2D",
      Model %in% models_3D ~ "3D",
      Model %in% models_CTRL ~ "CTRL"
    )
  ) %>%
  filter(!is.na(DSS))

dss_summary <- dss_long %>%
  group_by(Drug_Name, Class_Explained, Group) %>%
  summarise(mean_DSS = mean(DSS, na.rm = TRUE), .groups = "drop")

# Keep drug classes with an unadjusted pairwise Wilcoxon p-value < 0.05.
class_list <- unique(dss_summary$Class_Explained)
comparisons <- list(c("CTRL", "2D"), c("CTRL", "3D"), c("2D", "3D"))

get_p <- function(class_val, group1, group2) {
  test_df <- dss_summary %>% filter(Class_Explained == class_val)
  vals1 <- test_df %>% filter(Group == group1) %>% pull(mean_DSS)
  vals2 <- test_df %>% filter(Group == group2) %>% pull(mean_DSS)
  if(length(vals1) > 0 & length(vals2) > 0) {
    tryCatch(wilcox.test(vals1, vals2)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }
}

signif_classes <- map_lgl(class_list, function(cl) {
  any(sapply(comparisons, function(comp) get_p(cl, comp[1], comp[2]) < 0.05), na.rm = TRUE)
})

classes_to_plot <- class_list[signif_classes]

# Set the control, 2D and 3D display order.
dss_summary <- dss_summary %>%
  mutate(Group = factor(Group, levels = c("CTRL", "2D", "3D")))

# Plot drug means with points behind the boxes.

panel_plot <- ggplot(
  dss_summary %>% dplyr::filter(Class_Explained %in% classes_to_plot),
  aes(x = Group, y = mean_DSS)
) +
  # Draw jittered points behind boxes.
  geom_jitter(
    aes(color = Group),
    width = 0.15, size = 1.8, alpha = 0.6, show.legend = FALSE
  ) +
  # Draw boxes over the points.
  geom_boxplot(
    aes(fill = Group),
    width = 0.6, outlier.shape = NA, color = "black", alpha = 0.95,
    show.legend = FALSE
  ) +
  facet_wrap(~ Class_Explained) +

  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,
    step.increase = 0.28
  ) +

  scale_fill_manual(values = c("CTRL" = "grey40", "2D" = "steelblue", "3D" = "darkorange")) +
  scale_color_manual(values = c("CTRL" = "grey40", "2D" = "steelblue", "3D" = "darkorange")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_classic(base_size = 14) +
  labs(x = "Model Group", y = "Mean DSS")

print(panel_plot)

# Save PDF and PNG figures.
ggplot2::ggsave(file.path(output_dir, "Figure2D.pdf"), plot = panel_plot,
    width = 11, height = 8, units = "in")
ggplot2::ggsave(file.path(output_dir, "Figure2D.png"), plot = panel_plot,
    width = 11, height = 8, units = "in", dpi = 300, bg = "white")
