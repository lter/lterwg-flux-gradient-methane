# lterwg-flux-gradient-methane

Principal Investigators:
- Sparkle Malone (Yale School of the Environment)

## Overview

This repository contains the analysis workflow for NEON CH4 flux-gradient data, including site-level flux processing, ERA5-driven gap-filling, and spatial upscaling to a global upland terrestrial CH4 budget using a two-stage Random Forest model (source-probability classification + magnitude regression).

The master script is `workflows/WorkFlow.R`. It sources numbered workflow scripts in order (01 → 17, then 19 — step 18 is reserved and not currently in use). Run the full pipeline end-to-end by sourcing `WorkFlow.R`, or run individual steps as needed.

## Workflow Scripts

Scripts are named `##_STEP_Description.R` where `##` is the run order and `STEP` is the major phase (`NEON`, `Global`, or `Supp`).

**NEON Data Processing**

- `01_NEON_FilterQC.R` — Applies quality filters to raw NEON gradient flux data and computes QC flags across all sites.
- `02_NEON_SensorHeightPairs.R` — Identifies representative sensor height pairs (RSHP) based on concordance correlation coefficient (CCC ≥ 0.5) for CO2 and H2O.
- `03_NEON_EnsembleGapfill.R` — Ensemble gap-fills gradient fluxes (MBR, AE, WP methods) across all sites.
- `04_NEON_TotalFlux.R` — Combines ensemble gradient flux with storage flux to produce total CH4 flux.

**Site-level Analysis**

- `05_NEON_FluxAnalysis.R` — Computes standardized 30-min flux means, sampling-adjusted daily totals, and annual budgets; classifies sites as sources or sinks.
- `06_NEON_ERA5Gapfill.R` — ERA5-driven half-hourly gap-filling using Random Forest models fit to ERA5-Land covariates.
- `07_NEON_SiteMap.R` — Generates the site map figure used in the methods section.
- `08_NEON_Figures.R` — Produces all publication figures from NEON flux products (run after steps 05–06).
- `09_NEON_FLUXNETComparison.R` — Compares ERA5-gapfilled NEON annual budgets to published FLUXNET, process-model, and chamber reference rates (Delwiche et al. 2021). Also classifies NEON sites as climatically "Arid" (aridity_index = MAP / (MAT+10) below a fixed threshold, matching the definition used in steps 12/13) — no NEON site carries a literal "Arid" EcoType label, so this is how the two arid-analog sites (JORN, SRER) are identified for the ecosystem-comparison figures.

**Global Budget / Spatial Upscaling**

- `10_Global_DownloadERA5Land.R` — Downloads ERA5-Land monthly gridded climate inputs from Copernicus CDS.
- `11_Global_DownloadMODIS_WAD2M.R` — Downloads and processes MODIS MCD12C1 land cover (classified into 5 EcoType codes; water, permanent wetlands, urban/built-up, and snow/ice are excluded) and WAD2M inundation data (used to down-weight, not exclude, cell area) for spatial upscaling. The processed MODIS files are a hard requirement for step 13 — the Ecoregions-biome fallback it would otherwise use cannot distinguish cropland/urban/wetland from upland cover.
- `12_SourceProp_MagnitudeModels.R` — Fits the two-stage Random Forest model: Stage 1 is a balanced binary source-probability classifier (P(source)); Stage 2 is a log-magnitude regression, fit separately for weak-sink and weak-source conditions. Climatically "Arid" site-months (JORN, SRER) are excluded from Stage 1 and from Stage 2's source-magnitude model, and are instead routed to a forced-sink prediction — see the script's header docstring for the empirical rationale. Saves a model bundle consumed (not refit) by step 13.
- `12b_Model_FIGURES_RF.R` — Model-performance figures (calibration, Stage 1 classification skill, variable importance, observed-vs-OOB-fitted magnitude) from step 12's output.
- `12c_GeneralizationTest.R` — Leave-one-site-out (LOSO) cross-validation of the Stage 2 magnitude models: refits Stage 1/2 on all NEON sites but one, scores the held-out site, and repeats for every site. Diagnostic only — does not touch the production model bundle that step 13 uses. Slower than the other steps.
- `13_Global_SpatialUpscalingRF.R` — Applies the model bundle fit in step 12 to the global upland grid (does not fit or refit anything itself); produces monthly predictions under three flux-expression approaches (Continuous, Dichotomous, All-Sink). Requires processed MODIS land-cover files from step 11 and will stop with an explicit error if they're missing.
- `14_Global_SpatialUpscalingFiguresRF.R` — Standalone figure script for the RF upscaling pipeline; loads pre-computed outputs from step 13 (and, for one figure, step 12) and regenerates all publication figures without re-running the projection.

**Supplemental**

- `15_Supp_FluxJustification.R` — Addresses reviewer concerns about canopy storage bias, wetland footprint contamination, and within-canopy mixing artifacts.
- `16_Supp_GapfillUncertainty.R` — Sensitivity analysis testing whether source/sink classification is robust to systematic flux-magnitude scaling (0.5×–2.0×).
- `17_Supp_BudgetUncertainty.R` — Quantifies how systematic FG flux bias (±50%, multiplicative) propagates into the spatially upscaled continental CH4 budget (Tg CH4 yr⁻¹) across the three upscaling approaches.
- `19_Supp_NEONRepresentativeness.R` — Dissimilarity-index (DI) / Area-of-Applicability (AOA) analysis comparing NEON's training-data climate space to the global upland upscaling grid, to characterize how representative the NEON network is of the domain it's used to predict over.

**Alternative / Non-RF Upscaling (not run in main workflow)**

- `13alt_Global_SpatialUpscalingMonthly.R` — Older class-balanced source-probability upscaling using ERA5-Land grids (no Random Forest); retained for comparison.
- `14alt_Global_SpatialUpscalingFigures.R` — Figures for the non-RF upscaling approach above.

Deprecated scripts have been moved to `workflows/depreciaded/` and are not part of the active workflow.

## Validation Workflow

The validation pipeline is independent of the main workflow and lives in `workflows/Validation/`. It tests the flux-gradient method against co-located eddy-covariance measurements at three tower sites: SE-Sto (subarctic peat mire, Sweden), SE-Svb (boreal forest with a tall atmospheric profile tower, Sweden), and US-Uaf (boreal black spruce forest over permafrost, Alaska).

Run the validation end-to-end by sourcing `workflows/Validation/Workflow_Validation.R`, which runs steps in order 01 → 02 → 02a → 02b → 05 → 06 → 07. Steps 03/04 (a former diel CH4 analysis and its figures) have been removed from the pipeline; FG-vs-EC agreement is fully covered by steps 05 and 07 without that added modeling cost.

**Validation Scripts**

- `01_VAL_FilterQC.R` — Applies quality filters to raw validation site gradient flux data and computes QC flags.
- `02_VAL_SensorHeightPairs.R` — Identifies representative sensor height pairs (RSHP) for each validation site using a tiered concordance correlation coefficient (CCC) filter on CH4 directly (Tier 1: CCC > 0.1; Tier 2: 0 < CCC ≤ 0.1; Tier 3: best available CCC even if ≤ 0, flagged for caution), applied within candidate pairs allowed by each site's canopy structure.
- `02a` (**not a numbered file in this repo**) — Estimates single-point-approximation storage flux for all three sites. Sourced from the sibling repository `lterwg-flux-gradient` (`workflows/Validation/flow.validation.storage.R` and `VAL_StorageFlux_SinglePoint.R`), not from this repo — a dependency worth knowing about before running the validation pipeline standalone.
- `02b_VAL_TotalFlux.R` — Combines ensemble gradient flux with the single-point storage estimate from 02a to produce total flux (FG_total), falling back to gradient-only flux where no storage estimate is available for a given half-hour.
- `05_VAL_FluxAnalysis.R` — 30-min, daily, and annual flux analysis for both FG and EC fluxes through parallel pipelines for direct comparison.
- `06_VAL_ERA5Gapfill.R` — ERA5-driven half-hourly gap-filling for FG and EC fluxes at validation sites; adapted from `06_NEON_ERA5Gapfill.R`.
- `07_VAL_Supplement.R` — Supplemental figures justifying the FG validation approach (CCC metrics, bias assessment, storage/canopy-height relationship).
- `VAL_CounterGradientFilter.R` — Exploratory counter-gradient QC filter; not part of the main pipeline, retained for reference.

## Folder Structure

```
workflows/
  WorkFlow.R           # Master script — sources steps 01-17, then 19 (18 reserved)
  01-17_*.R, 19_*.R     # Numbered workflow scripts (NEON, Global, Supp)
  12b_*, 12c_*.R        # Model-performance figures / LOSO generalization test (run alongside 12)
  13alt_*, 14alt_*.R    # Alternative non-RF upscaling scripts (not in main workflow)
  Validation/           # Independent validation pipeline (steps 01-02b, 05-07)
  check/                # QC audit scripts and diagnostic figures
  depreciaded/          # Archived scripts no longer used in any workflow
functions/              # Shared R functions sourced by workflow scripts
exploratory/            # Exploratory analyses (not part of main workflow)
FIGURES/                # Output figures
```

## Naming Convention

Workflow scripts follow the pattern `##_STEP_Description.R`:
- `##` — two-digit run order (letter-suffixed steps like `12b`/`12c` run alongside their base-numbered step)
- `STEP` — major phase: `NEON`, `Global`, `Supp`, or `VAL` (validation)
- `Description` — brief CamelCase description of what the script does

## Contributing Guidelines & Style Guide

When you have a group of people collaborating on a shared project (particularly a code-heavy one), it can be nice to create some guidelines to make sure everyone is contributing in consistent ways. Similarly if your group reaches consensus on a 'style' of file names and/or code it can be good to formalize those rules as well. The standard convention in GitHub is to create a file called "CONTRIBUTING.md" that contains all of this information. If you want some inspiration check out the LTER Scientific Computing team's [CONTRIBUTING.md](https://github.com/lter/scicomp/blob/main/CONTRIBUTING.md) document!

## Related Repositories

- `lterwg-flux-gradient` — sibling repository containing shared flux-gradient functions and the validation storage-flux estimation scripts sourced by `Workflow_Validation.R` step 02a.
- `lterwg-flux-gradient-eval` — sibling repository containing shared evaluation/CCC functions (`functions/calc.One2One.CCC_testing.R`) sourced by the validation sensor-height-pair and counter-gradient-filter scripts.

## Supplementary Resources

LTER Scientific Computing Team [website](https://lter.github.io/scicomp/) & NCEAS' [Resources for Working Groups](https://www.nceas.ucsb.edu/working-group-resources)
