# Figure 2E: Target enrichment for 3D-versus-2D drug responses.

output_dir <- "Figures/Figure2/Figure2E/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(ggrepel)

data_dir <- "data"

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE, na.strings = 'Na')

# Define model groups.
models_2D <- c("SK.N.SH", "SK.N.BE2", "SK.N.AS", "SH.SY5Y", "SK.N.FI", "SK.N.DZ","AMC772T")
models_3D <- c("NB129", "NB067", "LU_NB_1", "LU_NB_2", "LU_NB_3", "IMR32")
models_CTRL <- c("DSS_BM1", "DSS_BM2", "DSS_BM3", "DSS_BM4", "DSS_BM5")

# Calculate model-group means and their differences.
dss <- dss %>%
  mutate(
    ave_2D = rowMeans(select(., all_of(models_2D)), na.rm = TRUE),
    ave_3D = rowMeans(select(., all_of(models_3D)), na.rm = TRUE),
    Cont_ave = rowMeans(select(., all_of(models_CTRL)), na.rm = TRUE),
    diff_2D_Cont = ave_2D - Cont_ave,
    diff_3D_Cont = ave_3D - Cont_ave,
    diff_3D_2D = ave_3D - ave_2D,
    diff_2D_3D = ave_2D - ave_3D
  )

# Calculate two-sided Fisher target-enrichment tests and an enrichment plot.
fisher_volcano <- function(data, diff_col, hit_threshold, analysis_label) {
  # Split semicolon-separated targets and define hits using the DSS difference threshold.
  data_long <- data %>%
    separate_rows(Gene_Targeted, sep = ";") %>%
    mutate(Gene_Targeted = trimws(Gene_Targeted)) %>%
    filter(!is.na(Gene_Targeted) & Gene_Targeted != "") %>%
    mutate(hit = .data[[diff_col]] > hit_threshold)

  targets <- unique(data_long$Gene_Targeted)
  fisher_results <- map_df(targets, function(target) {
    in_hit <- sum(data_long$Gene_Targeted == target & data_long$hit)
    not_in_hit <- sum(data_long$Gene_Targeted == target & !data_long$hit)
    in_hit_total <- sum(data_long$hit)
    not_in_hit_total <- sum(!data_long$hit)
    not_target_in_hit <- in_hit_total - in_hit
    not_target_not_in_hit <- not_in_hit_total - not_in_hit
    mat <- matrix(c(in_hit, not_in_hit, not_target_in_hit, not_target_not_in_hit), nrow = 2)
    ft <- fisher.test(mat)
    tibble(
      target = target,
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      in_hit = in_hit,
      not_in_hit = not_in_hit
    )
  }) %>%
    mutate(
      log_odds = log10(odds_ratio),
      neglog_p = -log10(p_value),
      size = in_hit + not_in_hit
    )

  # Plot target enrichment.
  log_p_threshold <- 1.3   # Approximately -log10(0.05).
  log_odds_threshold <- 0.3
  fisher_results <- fisher_results %>%
    mutate(
      color = case_when(
        neglog_p > log_p_threshold & log_odds > log_odds_threshold ~ "red",
        neglog_p > log_p_threshold & log_odds < -log_odds_threshold ~ "blue",
        TRUE ~ "gray70"
      )
    )
  label_genes <- fisher_results %>%
    filter(neglog_p > log_p_threshold & abs(log_odds) > log_odds_threshold)

  volcano_plot <- ggplot(fisher_results, aes(x = log_odds, y = neglog_p)) +
    geom_point(aes(size = size, color = color), alpha = 0.7) +
    scale_color_identity() +
    geom_vline(xintercept = c(-log_odds_threshold, log_odds_threshold), linetype = "dotted") +
    geom_hline(yintercept = log_p_threshold, linetype = "dotted") +
    geom_text_repel(
      data = label_genes,
      aes(label = target, color = color),
      size = 4, fontface = "bold",
      segment.color = "black",
      box.padding = 0.5,
      point.padding = 0.4
    ) +
    scale_size_continuous(range = c(2, 8), guide = "none") +
    theme_classic(base_size = 15) +
    labs(
      x = "log10(Odds Ratio)",
      y = "-log10(Fisher Test P-Value)",
      title = paste0(analysis_label)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
    )

  return(list(table = fisher_results, plot = volcano_plot))
}

result_3D_vs_2D <- fisher_volcano(dss, 'diff_3D_2D', 10, '3D vs 2D (3D-2D > 10)')

# Display infinite log odds at finite coordinates; keep the statistical table unchanged.
fisher_table <- result_3D_vs_2D$table
finite_x <- fisher_table$log_odds[is.finite(fisher_table$log_odds)]
pos_inf_x <- max(finite_x) + 0.4
neg_inf_x <- min(finite_x) - 0.4
plot_table <- fisher_table %>% mutate(
    display_x = case_when(log_odds == Inf ~ pos_inf_x, log_odds == -Inf ~ neg_inf_x, TRUE ~ log_odds),
    display_shape = ifelse(is.infinite(log_odds), 17, 16),
    display_color = case_when(neglog_p > 1.3 & log_odds > 0.5 ~ "red",
        neglog_p > 1.3 & log_odds < -0.5 ~ "blue", TRUE ~ "gray70"))
label_genes <- plot_table %>% filter(neglog_p > 1.3, abs(log_odds) > 0.5)
x_breaks <- pretty(range(finite_x), n = 5)
x_breaks <- x_breaks[x_breaks >= min(finite_x) & x_breaks <= max(finite_x)]
x_labels <- as.character(x_breaks)
if (any(fisher_table$log_odds == -Inf, na.rm = TRUE)) { x_breaks <- c(neg_inf_x, x_breaks); x_labels <- c("-Inf", x_labels) }
if (any(fisher_table$log_odds == Inf, na.rm = TRUE)) { x_breaks <- c(x_breaks, pos_inf_x); x_labels <- c(x_labels, "Inf") }

panel_plot <- ggplot(plot_table, aes(x = display_x, y = neglog_p)) +
    geom_point(aes(size = size, color = display_color, shape = display_shape), alpha = 0.7) +
    scale_color_identity() + scale_shape_identity() + scale_size_continuous(range = c(2, 8), guide = "none") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted") + geom_hline(yintercept = 1.3, linetype = "dotted") +
    geom_text_repel(data = label_genes, aes(label = target, color = display_color), fontface = "bold",
        size = 4, seed = 123, box.padding = 0.5, point.padding = 0.4, max.overlaps = Inf) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = expansion(mult = c(0.1, 0.2))) +
    theme_classic(base_size = 15) + labs(x = "log10(Odds Ratio)", y = "-log10(Fisher Test P-Value)")

print(panel_plot)
# Save PDF and PNG figures.
ggsave(file.path(output_dir, "Figure2E.pdf"), plot = panel_plot,
    width = 8, height = 7, units = "in")
ggsave(file.path(output_dir, "Figure2E.png"), plot = panel_plot,
    width = 8, height = 7, units = "in", dpi = 300, bg = "white")
