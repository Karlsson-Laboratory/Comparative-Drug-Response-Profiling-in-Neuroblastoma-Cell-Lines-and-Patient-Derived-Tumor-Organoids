# Supplementary Figure 2F: FO5 and PMC drug-library overlap using hardcoded counts.

library(tidyverse)
library(eulerr)

data_dir <- "data"
output_dir <- "Figures/SupplementaryFigure2/SupplementaryFigure2F/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load FO5 and PMC library inputs.
fo5 <- read.delim(file.path(data_dir, "fo5_primary_drug_sensitivity_scores.tsv"),
    stringsAsFactors = FALSE)
fo5_fo4 <- read.delim(file.path(data_dir, "fo5_fo4_integrated_drug_sensitivity_scores.tsv"),
    stringsAsFactors = FALSE)
pmc_tumoroids <- read.delim(file.path(data_dir, "pmc_tumoroid_drug_response_auc.tsv"),
    stringsAsFactors = FALSE, na.strings = c("NA", "na", "Na", "N/A"))
pmc_cell_lines <- read.delim(file.path(data_dir, "pmc_cell_line_drug_response_auc.tsv"),
    stringsAsFactors = FALSE, na.strings = c("NA", "na", "Na", "N/A"))

# Harmonize drug-name aliases across libraries.
fix_drug_names <- function(x) {
    dplyr::recode(x,
        "Actinomycin D" = "Dactinomycin",
        "6-Mercaptopurine" = "Mercaptopurine",
        "6-Thioguanine" = "Thioguanine",
        "SN-38 (Irinotecan)" = "SN-38",
        "GSK461364" = "GSK-461364",
        "GSK1070916" = "GSK-1070916",
        "GSK269962a" = "GSK269962",
        "S63845" = "S-63845",
        "VER-155008" = "VER_155008",
        "Tubastatin A" = "Tubastatin_A",
        "StemRegenin 1" = "StemRegenin_1",
        "Mitomycin C" = "Mitomycin_C",
        "Valproic Acid" = "Valproic_acid",
        "Sepantronium Bromide" = "Sepantronium_bromide",
        "ENMD2076" = "ENMD-2076",
        "JQ-1 Carboxylic Acid ((+)-JQ-1)" = "JQ1_1",
        "(+)-JQ-1" = "JQ1",
        .default = x
    )
}
fo5_drugs <- unique(fix_drug_names(fo5$Drug_Name))
fo5_fo4_drugs <- unique(fix_drug_names(fo5_fo4$Drug_Name))
pmc_tumoroid_drugs <- unique(fix_drug_names(pmc_tumoroids$Drug_Name))
pmc_cell_line_drugs <- unique(fix_drug_names(pmc_cell_lines$Drug_Name))
pmc_union_drugs <- union(pmc_tumoroid_drugs, pmc_cell_line_drugs)

# Calculate drug memberships from the supplied FO5 and PMC lists.

full_shared_drugs <- intersect(fo5_drugs, pmc_tumoroid_drugs)
full_fo5_only_drugs <- setdiff(fo5_drugs, pmc_tumoroid_drugs)
full_pmc_only_drugs <- setdiff(pmc_tumoroid_drugs, fo5_drugs)
integrated_shared_drugs <- intersect(fo5_fo4_drugs, pmc_tumoroid_drugs)

# Calculate overlap counts for the supplied library combinations.
fo5_lists <- list(fo5_drugs, fo5_fo4_drugs, fo5_drugs, fo5_drugs)
pmc_lists <- list(pmc_tumoroid_drugs, pmc_tumoroid_drugs, pmc_cell_line_drugs, pmc_union_drugs)
computed_counts <- tibble(
    inputs = c("Full FO5 / PMC tumoroids", "FO5-FO4 subset / PMC tumoroids",
        "Full FO5 / PMC cell lines", "Full FO5 / PMC union"),
    FO5_total = lengths(fo5_lists),
    PMC_total = lengths(pmc_lists),
    FO5_only = map2_int(fo5_lists, pmc_lists, ~ length(setdiff(.x, .y))),
    PMC_only = map2_int(fo5_lists, pmc_lists, ~ length(setdiff(.y, .x))),
    shared = map2_int(fo5_lists, pmc_lists, ~ length(intersect(.x, .y)))
)
print(computed_counts)
message("Original counts retained: 294 FO5-only, 72 PMC-only, 137 shared. ",
    "The exact 431-compound FO5 subset/source is unresolved.")

# Plot hardcoded counts: 294 FO5-only, 72 PMC-only and 137 shared.

original_counts <- c("FO5" = 294, "PMC" = 72, "FO5&PMC" = 137)
fit <- euler(original_counts)
panel_plot <- plot(fit,
    fills = list(fill = c("white", "grey90")),
    edges = list(col = "black"),
    quantities = list(cex = 1.4),
    labels = list(cex = 1.4, font = 2))
print(panel_plot)

# Save PDF and PNG figures.
ggplot2::ggsave(file.path(output_dir, "SupplementaryFigure2F.pdf"), plot = panel_plot,
    width = 8, height = 6, units = "in")
ggplot2::ggsave(file.path(output_dir, "SupplementaryFigure2F.png"), plot = panel_plot,
    width = 8, height = 6, units = "in", dpi = 300, bg = "white")
