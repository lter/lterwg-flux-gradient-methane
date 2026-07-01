# lterwg-flux-gradient-methane

Principal Investigators:
- Sparkle Malone (Yale School of the Environment)

## Overview

This repository contains the analysis workflow for NEON CH4 flux-gradient data, including site-level flux processing, ERA5-driven gap-filling, and spatial upscaling to a continental CH4 budget using Random Forest models.

The master script is `workflows/WorkFlow.R`. It sources numbered workflow scripts in order (01 → 17). Run the full pipeline end-to-end by sourcing `WorkFlow.R`, or run individual steps as needed.

## Workflow Scripts

Scripts are named `##_STEP_Description.R` where `##` is the run order and `STEP` is the major phase (`NEON`, `Global`, or `Supp`).

**NEON Data Processing**

- `01_NEON_FilterQC.R` — Applies quality filters to raw NEON gradient flux data and computes QC flags across all sites.
- `02_NEON_SensorHeightPairs.R` — Identifies representative sensor height pairs (RSHP) based on concordance correlation coefficient (CCC ≥ 0.5) for CO2 and H2O.
- `03_NEON_EnsembleGapfill.R` — Ensemble gap-fills gradient fluxes (MBR, AE, WP methods) across all sites.
- `04_NEON_TotalFlux.R` — Combines ensemble gradient flux with storage flux to produce total CH4 flux.

**Site-level Analysis**

- `05_NEON_FluxAnalysis.R` — Computes standardized 30-min flux means, sampling-adjusted daily totals, and annual budgets; classifies sites as sources or sinks.
- `06_NEON_ERA5Gapfill.R` — ERA5-driven half-hourly gap-filling using GAM models fit to temperature and soil moisture covariates.
- `07_NEON_SiteMap.R` — Generates the site map figure used in the methods section.
- `08_NEON_Figures.R` — Produces all publication figures from NEON flux products (run after steps 05–06).
- `09_NEON_FLUXNETComparison.R` — Compares ERA5-gapfilled NEON annual budgets to published FLUXNET, process-model, and chamber reference rates (Delwiche et al. 2021).

**Global Budget / Spatial Upscaling**

- `10_Global_DownloadERA5Land.R` — Downloads ERA5-Land monthly gridded climate inputs from Copernicus CDS.
- `11_Global_DownloadMODIS_WAD2M.R` — Downloads and processes MODIS MCD12C1 land cover and WAD2M inundation data for spatial upscaling.
- `12_Global_SourceProbability.R` — Fits a condition-based upland CH4 source-probability model; validates with leave-one-site-out prediction.
- `13_Global_SpatialUpscalingRF.R` — Monthly Random Forest spatial upscaling using three flux-expression approaches (continuous, threshold, hybrid).
- `14_Global_SpatialUpscalingFiguresRF.R` — Standalone figure script for the RF upscaling pipeline; loads pre-computed outputs and regenerates all publication figures.

**Supplemental**

- `15_Supp_FluxJustification.R` — Addresses reviewer concerns about canopy storage bias, wetland footprint contamination, and within-canopy mixing artifacts.
- `16_Supp_GapfillUncertainty.R` — Sensitivity analysis testing whether source/sink classification is robust to systematic flux-magnitude scaling (0.5×–2.0×).
- `17_Supp_BudgetUncertainty.R` — Quantifies how systematic FG flux bias propagates into the spatially upscaled continental CH4 budget (Tg CH4 yr⁻¹).

**Alternative / Non-RF Upscaling (not run in main workflow)**

- `13alt_Global_SpatialUpscalingMonthly.R` — Older class-balanced source-probability upscaling using ERA5-Land grids (no Random Forest); retained for comparison.
- `14alt_Global_SpatialUpscalingFigures.R` — Figures for the non-RF upscaling approach above.

Deprecated scripts have been moved to `workflows/depreciaded/` and are not part of the active workflow.

## Folder Structure

```
workflows/         # Numbered workflow scripts (01–17) + WorkFlow.R
  check/           # QC audit scripts and diagnostic figures
  depreciaded/     # Archived scripts no longer used in main workflow
  Validation/      # Validation workflow (independent of main pipeline)
functions/         # Shared R functions sourced by workflow scripts
exploratory/       # Exploratory analyses (not part of main workflow)
FIGURES/           # Output figures
```

## Naming Convention

Workflow scripts follow the pattern `##_STEP_Description.R`:
- `##` — two-digit run order (01–17)
- `STEP` — major phase: `NEON`, `Global`, or `Supp`
- `Description` — brief CamelCase description of what the script does

## Contributing Guidelines & Style Guide

When you have a group of people collaborating on a shared project (particularly a code-heavy one), it can be nice to create some guidelines to make sure everyone is contributing in consistent ways. Similarly if your group reaches consensus on a 'style' of file names and/or code it can be good to formalize those rules as well. The standard convention in GitHub is to create a file called "CONTRIBUTING.md" that contains all of this information. If you want some inspiration check out the LTER Scientific Computing team's [CONTRIBUTING.md](https://github.com/lter/scicomp/blob/main/CONTRIBUTING.md) document!

## Related Repositories

- TBD (recommend including user/repository name, link, and brief description)

## Supplementary Resources

LTER Scientific Computing Team [website](https://lter.github.io/scicomp/) & NCEAS' [Resources for Working Groups](https://www.nceas.ucsb.edu/working-group-resources)
