# Figure 1B: Distribution of drug classes in the FO5 library.

output_dir <- "Figures/Figure1/Figure1B/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(viridis)
library(RColorBrewer)

data_dir <- "data"

# Load drug sensitivity scores.
dss <- read.delim(file.path(data_dir, "fo5_primary_drug_sensitivity_scores.tsv"), stringsAsFactors = FALSE)

# Clean drug-class annotations.
clean_class <- function(x) {
    # Remove the class prefix.
    x <- gsub("^[A-Z]\\._", "", x)
    # Trim whitespace and replace underscores with spaces.
    x <- str_trim(gsub("_", " ", x))
    # Treat empty and 'ssd' class annotations as missing.
    x <- ifelse(x == "" | x == "ssd", NA, x)
    # Standardize class names.
    x <- recode(x, `Kinase inhibitor` = "Kinase Inhibitor", `Conv. Chemo` = "Conventional Chemotherapy",
        Other = "Other", `Kinesin inhibitor` = "Kinesin Inhibitor", `Metabolic modifier` = "Metabolic Modifier",
        `Apoptotic modulator` = "Apoptosis Modulator", `HSP inhibitor` = "HSP Inhibitor", `Differentiating/ epigenetic modifier` = "Epigenetic/Differentiation Modifier",
        `Hormone therapy` = "Hormone Therapy", `Hedgehog inhibitor` = "Hedgehog Inhibitor", `Protease/proteasome inhibitor` = "Proteasome Inhibitor",
        Immunomodulatory = "Immunomodulator", NSAID = "NSAID", Rapalog = "Rapalog")
    x
}

dss$Drug_Class <- clean_class(dss$Class_Explained)

# Remove missing class annotations.
dss <- dss %>%
    filter(!is.na(Drug_Class))

# Calculate class counts and proportions.
df <- dss %>%
    count(Drug_Class) %>%
    mutate(Percentage = n/sum(n), ymax = cumsum(Percentage), ymin = c(0, head(ymax, -1)))

# Plot the class distribution.

panel_plot <- ggplot(df, aes(ymin = ymin, ymax = ymax, xmin = 0, xmax = 1, fill = Drug_Class)) +
    geom_rect(color = NA) + coord_polar(theta = "y") + theme_void(base_size = 14) +
    labs(title = "FO5", fill = "Class") + scale_fill_viridis_d(option = "turbo") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_text(face = "bold"), legend.text = element_text(size = 10))

print(panel_plot)
# Save PDF and PNG figures.
ggsave(file.path(output_dir, "Figure1B.pdf"), plot = panel_plot,
    width = 11, height = 8, units = "in")
ggsave(file.path(output_dir, "Figure1B.png"), plot = panel_plot,
    width = 11, height = 8, units = "in", dpi = 300, bg = "white")
