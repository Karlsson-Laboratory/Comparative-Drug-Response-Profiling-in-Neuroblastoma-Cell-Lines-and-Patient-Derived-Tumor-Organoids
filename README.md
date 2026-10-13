# Neuroblastoma drug-response figures

Standalone R scripts in `Figures/<figure>/<panel>/` load input tables, calculate panel results, plot them and save PDF/PNG figures. Model groups and panel-specific settings are defined directly in each script.

## Inputs

Keep the following files in `data/`, preserving their headers, encodings and missing-value tokens.

| Filename | Contents |
|---|---|
| `fo5_primary_drug_sensitivity_scores.tsv` | FO5 drug annotations and model DSS values |
| `fo5_fo4_integrated_drug_sensitivity_scores.tsv` | Integrated FO5/FO4 DSS and drug annotations |
| `focused_28_drug_integrated_sensitivity_scores.tsv` | Focused drug-response measurements across conditions |
| `rpmi_nbom_revision_drug_sensitivity_scores.csv` | RPMI/NBOM drug-response measurements |
| `focused_29_drug_topology_validation_dss.csv` | Drug groups and 29-drug, 15-column topology DSS matrix |
| `pmc_tumoroid_drug_response_auc.tsv` | PMC tumoroid drug-response AUC values |
| `pmc_cell_line_drug_response_auc.tsv` | PMC cell-line drug-response AUC values |
| `seeding_density_<model>_<time>_luminescence.tsv` | Replicate seeding-density luminescence |

The 18 density files cover AMC772T, IMR32, NB067, NB129, SHSY5Y, SKNBE2, SKNDZ, SKNFI and SKNSH at `0h` and `72h`. There are 25 input files in total. Scripts use the inputs indicated in their loading statements.

## Run in RStudio

1. Set the working directory to the repository root so `data/` and `Figures/` resolve correctly.
2. Install the packages used by the scripts:

   ```r
   install.packages(c("tidyverse", "ggpubr", "ggrepel", "pheatmap", "eulerr",
       "viridis", "RColorBrewer", "circlize", "ragg", "BiocManager"))
   BiocManager::install("ComplexHeatmap")
   ```

3. Start a fresh R session for each panel, open its script and run the complete file.

Scripts create their output directories and save figures in each panel’s `outputs/` subfolder. Multi-plot scripts use separate filenames for each plot. Generated outputs are excluded from version control by `.gitignore`.
