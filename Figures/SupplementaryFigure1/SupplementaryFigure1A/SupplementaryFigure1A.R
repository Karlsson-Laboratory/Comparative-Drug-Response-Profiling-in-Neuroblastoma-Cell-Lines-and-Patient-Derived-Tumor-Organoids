# Supplementary Figure 1A: Fold growth across seeding densities.

output_dir <- "Figures/SupplementaryFigure1/SupplementaryFigure1A/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)

data_dir <- "data"

# Identify models with both 0-hour and 72-hour luminescence files.
files_0h <- list.files(path = data_dir, pattern = "seeding_density_(.*)_0h_luminescence\\.tsv$")
files_72h <- list.files(path = data_dir, pattern = "seeding_density_(.*)_72h_luminescence\\.tsv$")
cell_lines <- intersect(sub("seeding_density_(.*)_0h_luminescence\\.tsv", "\\1", files_0h), sub("seeding_density_(.*)_72h_luminescence\\.tsv", "\\1",
    files_72h))

# Calculate and plot fold growth for one model.
get_and_plot_fold_growth <- function(cl) {
    file_0h <- paste0("seeding_density_", cl, "_0h_luminescence.tsv")
    file_72h <- paste0("seeding_density_", cl, "_72h_luminescence.tsv")
    df_0h <- read.delim(file.path(data_dir, file_0h), stringsAsFactors = FALSE, header = TRUE)
    df_72h <- read.delim(file.path(data_dir, file_72h), stringsAsFactors = FALSE, header = TRUE)

    # Select seeding-density columns.
    density_cols <- grep("^X[0-9]+", names(df_0h), value = TRUE)
    df_0h <- df_0h %>%
        select(all_of(density_cols))
    df_72h <- df_72h %>%
        select(all_of(density_cols))

    # Calculate fold growth as mean 72-hour luminescence divided by mean 0-hour luminescence.
    df_0h_long <- df_0h %>%
        mutate(timepoint = "0h", replicate = row_number()) %>%
        pivot_longer(-c(timepoint, replicate), names_to = "density", values_to = "value")
    df_72h_long <- df_72h %>%
        mutate(timepoint = "72h", replicate = row_number()) %>%
        pivot_longer(-c(timepoint, replicate), names_to = "density", values_to = "value")
    df_all <- bind_rows(df_0h_long, df_72h_long)

    means <- df_all %>%
        group_by(density, timepoint) %>%
        summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = timepoint, values_from = mean_value) %>%
        mutate(fold_growth = `72h`/`0h`, density_label = as.integer(gsub("^X", "", density))) %>%
        arrange(density_label)

    # Plot fold growth by seeding density.
    p <- ggplot(means, aes(x = factor(density_label, levels = density_label), y = fold_growth)) + geom_col(fill = "white",
        color = "black") + theme_classic(base_size = 14) + labs(x = "Cells per well", y = "Fold Growth",
        title = paste(cl))

    print(p)

    # Save PDF and PNG figures for each model.
    ggsave(file.path(output_dir, paste0("SupplementaryFigure1A_", cl, "_fold_growth.pdf")),
        plot = p, width = 11, height = 8, units = "in")
    ggsave(file.path(output_dir, paste0("SupplementaryFigure1A_", cl, "_fold_growth.png")),
        plot = p, width = 11, height = 8, units = "in", dpi = 300, bg = "white")
}

# Generate figures for all models with paired timepoints.
walk(cell_lines, get_and_plot_fold_growth)
