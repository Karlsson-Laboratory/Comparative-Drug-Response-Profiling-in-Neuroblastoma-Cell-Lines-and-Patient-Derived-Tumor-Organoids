# Supplementary Figure 3F: DSS heatmap with drug-group and model annotations.

library(ComplexHeatmap)
library(circlize)

data_dir <- "data"
output_dir <- "Figures/SupplementaryFigure3/SupplementaryFigure3F/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load drug sensitivity scores.
df <- read.table(file.path(data_dir, "focused_29_drug_topology_validation_dss.csv"),
    sep = ",", header = TRUE)
head(df)
rownames(df) <- df[, 3]
dfm <- as.matrix(df[, 4:18])
head(dfm)

# Generate the exploratory heatmap with default clustering and colors.
ht_exploratory <- Heatmap(dfm, show_row_names = TRUE)
draw(ht_exploratory)
# Save the exploratory heatmap as PDF and PNG.
pdf(file.path(output_dir, "SupplementaryFigure3F_exploratory_clustered.pdf"),
    width = 20, height = 40)
draw(ht_exploratory)
dev.off()
ragg::agg_png(file.path(output_dir, "SupplementaryFigure3F_exploratory_clustered.png"),
    width = 20, height = 40, units = "in", res = 150)
draw(ht_exploratory)
dev.off()

# Map Drug_group categories to the three ordered row sections (10/10/9 drugs).
row_group_labels <- c("2D" = "Effective in 2D", "3D" = "Effective in 3D",
    "Both" = "Effective in both")
row_sections <- factor(unname(row_group_labels[df$Drug_group]),
    levels = c("Effective in 2D", "Effective in 3D", "Effective in both"))
stopifnot(!anyNA(row_sections), identical(as.integer(table(row_sections)), c(10L, 10L, 9L)))

# Derive model prefixes and exact condition suffixes from column names.
# read.table() converts model names to R identifiers such as "SK.N.BE.2.".
column_models <- sub("_[^_]+$", "", colnames(dfm))
column_conditions <- sub("^.*_", "", colnames(dfm))
model_labels <- c("SK.N.BE.2." = "SK-N-BE(2)", "SK.N.DZ" = "SK-N-DZ",
    "SK.N.AS" = "SK-N-AS", "AMC772T" = "AMC772T", "NB067" = "NB067")
model_names <- unname(model_labels[column_models])
column_sections <- factor(model_names, levels = unname(model_labels))

# Pre2D, Pre2DPDO and Pre3D denote first-screen columns.
# Exact suffixes 2D and 3D denote the grown conditions.
condition_labels <- c("Pre2D" = "First screen", "Pre2DPDO" = "First screen",
    "Pre3D" = "First screen", "2D" = "Grown in 2D", "3D" = "Grown in 3D")
screening_condition <- unname(condition_labels[column_conditions])
model_type_labels <- c("SK.N.BE.2." = "Cell line", "SK.N.DZ" = "Cell line",
    "SK.N.AS" = "Cell line", "AMC772T" = "Tumoroid", "NB067" = "Tumoroid")
model_type <- unname(model_type_labels[column_models])
stopifnot(!anyNA(column_sections), !anyNA(screening_condition), !anyNA(model_type),
    identical(as.integer(table(column_sections)), rep(3L, 5)),
    identical(screening_condition, rep(c("First screen", "Grown in 2D", "Grown in 3D"), 5)))

condition_colors <- c("First screen" = "#B2DF8A", "Grown in 2D" = "#4DAF4A",
    "Grown in 3D" = "#235027")
model_type_colors <- c("Cell line" = "#E6EF85", "Tumoroid" = "#66722E")
bottom_annotation <- HeatmapAnnotation(
    `Screening condition` = screening_condition,
    `Model type` = model_type,
    col = list(`Screening condition` = condition_colors, `Model type` = model_type_colors),
    which = "column",
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = grid::unit(3, "mm"),
    gap = grid::unit(0, "mm")
)

# Arrange condition and model-type legend swatches in one row.
condition_legend <- Legend(labels = names(condition_colors),
    legend_gp = grid::gpar(fill = condition_colors), direction = "horizontal", nrow = 1,
    labels_gp = grid::gpar(fontsize = 10))
model_type_legend <- Legend(labels = names(model_type_colors),
    legend_gp = grid::gpar(fill = model_type_colors), direction = "horizontal", nrow = 1,
    labels_gp = grid::gpar(fontsize = 10))
annotation_legends <- packLegend(condition_legend, model_type_legend,
    direction = "horizontal", gap = grid::unit(4, "mm"))

# Map DSS to blue, purple, white and red; color saturation does not change matrix values.

dss_colors <- colorRamp2(c(0, 5, 10, 20), c("#3F4AA5", "#AE8DC3", "white", "#F11920"))
ht <- Heatmap(
    dfm,
    name = "DSS",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    row_split = row_sections,
    column_split = column_sections,
    row_gap = grid::unit(3, "mm"),
    column_gap = grid::unit(3, "mm"),
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = grid::gpar(fontsize = 13),
    column_title = NULL,
    column_labels = model_names,
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 12),
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = 11),
    row_names_max_width = grid::unit(48, "mm"),
    bottom_annotation = bottom_annotation,
    col = dss_colors,
    show_heatmap_legend = FALSE
)
# Place annotation legends above the horizontal DSS legend.
dss_legend <- Legend(title = "DSS", col_fun = dss_colors, at = c(0, 5, 10, 15, 20),
    direction = "horizontal", legend_width = grid::unit(5, "cm"),
    title_position = "topleft", labels_gp = grid::gpar(fontsize = 11))
bottom_legends <- packLegend(annotation_legends, dss_legend,
    direction = "vertical", gap = grid::unit(3, "mm"))
draw(ht, auto_adjust = FALSE, heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom", annotation_legend_list = list(bottom_legends),
    padding = grid::unit(c(6, 6, 6, 6), "mm"))
# Save the nonclustered heatmap as PDF and PNG.
pdf(file.path(output_dir, "SupplementaryFigure3F_nonclustered.pdf"),
    width = 10, height = 10)
draw(ht, auto_adjust = FALSE, heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom", annotation_legend_list = list(bottom_legends),
    padding = grid::unit(c(6, 6, 6, 6), "mm"))
dev.off()
ragg::agg_png(file.path(output_dir, "SupplementaryFigure3F_nonclustered.png"),
    width = 10, height = 10, units = "in", res = 300)
draw(ht, auto_adjust = FALSE, heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom", annotation_legend_list = list(bottom_legends),
    padding = grid::unit(c(6, 6, 6, 6), "mm"))
dev.off()
