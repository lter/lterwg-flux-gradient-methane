# WorkFlow.R
# Master workflow for NEON CH4 flux-gradient analysis and global budget.
# Run scripts in order (01 -> 17). Each step sources a numbered workflow script.

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(sf)

# -------------- Change this stuff -------------
DirRepo.ch4 <- "/Users/sm3466/Library/CloudStorage/Dropbox-YSE/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-methane"
DirRepo.eval <- "/Users/sm3466/Dropbox-YSE/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-eval"

setwd(DirRepo.ch4)
localdir     <- '/Volumes/MaloneLab/Research/FluxGradient'
localdir.ch4 <- '/Volumes/MaloneLab/Research/FluxGradient/METHANE'
DnldFromGoogleDrive <- FALSE # TRUE = grab files from Google Drive; FALSE = use local copies

email <- 'sparklelmalone@gmail.com'
googledrive::drive_auth(email = TRUE)
drive_url    <- googledrive::as_id("https://drive.google.com/drive/folders/1Q99CT77DnqMl2mrUtuikcY47BFpckKw3")
data_folder  <- googledrive::drive_ls(path = drive_url)

metadata  <- read.csv('/Volumes/MaloneLab/Research/FluxGradient/Ameriflux_NEON field-sites.csv')
site.list <- metadata$Site_Id.NEON %>% unique

# ── NEON DATA PROCESSING ──────────────────────────────────────────────────────

# 01 Filter raw NEON gradient flux data and compute QC flags
message('Step 01: Filtering NEON data...')
source(fs::path(DirRepo.ch4, 'workflows/01_NEON_FilterQC.R'))

fileSave <- fs::path(localdir.ch4, "/NEON_GradientFlux_Data_Filter/SITE_DATA_FILTERED_CH4.Rdata")
save(SITE_DATA_FILTERED, file = fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir.ch4, "/NEON_GradientFlux_Data_Filter/SITES_MBR_9min_FILTER_CH4.Rdata")
save(SITES_MBR_9min_FILTER, file = fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir.ch4, "/NEON_GradientFlux_Data_Filter/SITES_AE_9min_FILTER_CH4.Rdata")
save(SITES_AE_9min_FILTER, file = fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir.ch4, "/NEON_GradientFlux_Data_Filter/SITES_WP_9min_FILTER_CH4.Rdata")
save(SITES_WP_9min_FILTER, file = fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

# 02 Identify representative sensor height pairs (RSHP) using CCC for CO2/H2O
message('Step 02: Identifying sensor height pairs...')
source(fs::path(DirRepo.ch4, 'workflows/02_NEON_SensorHeightPairs.R'))

# 03 Ensemble gap-fill the gradient fluxes across methods
message('Step 03: Ensemble gap-filling...')
source(fs::path(DirRepo.ch4, 'workflows/03_NEON_EnsembleGapfill.R'))

# 04 Combine ensemble gradient flux with storage flux to get total flux
message('Step 04: Computing total flux...')
source(fs::path(DirRepo.ch4, 'workflows/04_NEON_TotalFlux.R'))

# ── SITE-LEVEL ANALYSIS ───────────────────────────────────────────────────────

# 05 30-min flux analysis: standardized means, daily totals, annual budgets
message('Step 05: 30-min flux analysis...')
source(fs::path(DirRepo.ch4, 'workflows/05_NEON_FluxAnalysis.R'))

# 06 ERA5-driven half-hourly gap-filling for NEON sites
message('Step 06: ERA5 half-hourly gap-fill...')
source(fs::path(DirRepo.ch4, 'workflows/06_NEON_ERA5Gapfill.R'))

# 07 Site map figure (methods)
message('Step 07: Site map...')
source(fs::path(DirRepo.ch4, 'workflows/07_NEON_SiteMap.R'))

# 08 Publication figures from NEON flux products
message('Step 08: NEON figures...')
source(fs::path(DirRepo.ch4, 'workflows/08_NEON_Figures.R'))

# 09 Compare ERA5-gapfilled NEON fluxes with FLUXNET reference rates
message('Step 09: FLUXNET comparison...')
source(fs::path(DirRepo.ch4, 'workflows/09_NEON_FLUXNETComparison.R'))
## OUTPUT/CH4_flux_medians_by_source_and_behavior.csv

# ── GLOBAL BUDGET / SPATIAL UPSCALING ────────────────────────────────────────

# 10 Download ERA5-Land monthly gridded climate inputs
message('Step 10: Download ERA5-Land...')
source(fs::path(DirRepo.ch4, 'workflows/10_Global_DownloadERA5Land.R'))

# 11 Download and process MODIS land cover and WAD2M inundation inputs
message('Step 11: Download MODIS/WAD2M...')
source(fs::path(DirRepo.ch4, 'workflows/11_Global_DownloadMODIS_WAD2M.R'))

# 12 Condition-based upland CH4 source-probability model
message('Step 12: Source probability model...')
source(fs::path(DirRepo.ch4, 'workflows/12_Global_SourceProbability.R'))

# 13 Random Forest spatial upscaling (monthly, three flux-expression approaches)
message('Step 13: RF spatial upscaling...')
source(fs::path(DirRepo.ch4, 'workflows/13_Global_SpatialUpscalingRF.R'))

# 14 Publication figures for RF spatial upscaling
message('Step 14: RF upscaling figures...')
source(fs::path(DirRepo.ch4, 'workflows/14_Global_SpatialUpscalingFiguresRF.R'))

# ── SUPPLEMENTAL ──────────────────────────────────────────────────────────────

# 15 Supplemental: flux-gradient method justification (canopy storage, footprint)
message('Step 15: Flux justification supplement...')
source(fs::path(DirRepo.ch4, 'workflows/15_Supp_FluxJustification.R'))

# 16 Supplemental: gap-fill flux-magnitude uncertainty analysis
message('Step 16: Gap-fill uncertainty supplement...')
source(fs::path(DirRepo.ch4, 'workflows/16_Supp_GapfillUncertainty.R'))

# 17 Supplemental: spatial budget uncertainty from systematic flux bias
message('Step 17: Budget uncertainty supplement...')
source(fs::path(DirRepo.ch4, 'workflows/17_Supp_BudgetUncertainty.R'))

# ── DEPRECATED (not run in main workflow — kept for reference) ────────────────
# Scripts below have been moved to workflows/depreciaded/
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/NEON.DriveScale.Analysis.R'))
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/NEON.30min.Gapfill.r'))
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/ERA-Upscaling.R'))
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/flow.DIEL.NEON.R'))
# source(fs::path(DirRepo.ch4,'workflows/13alt_Global_SpatialUpscalingMonthly.R'))  # non-RF version
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/NEON.StrongSink.DriverComparison.R'))
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/NEON.DIEL.Analysis2.R'))
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/NEON.MonthlySinkBehavior.Analysis.R'))
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/NEON.Supplementary.ModelDriverPlots.R'))
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/NEON.ConsistencyMagnitude.Analysis.R'))
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/NEON.TotalFlux.AnnualBudget.R'))
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/NEON.DIEL.SynthesisFigure.R'))
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/flow.climate.NEON.R'))
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/flow.soilmoisture.NEON.R'))
# source(fs::path(DirRepo.ch4,'workflows/depreciaded/flow.parms_results.R'))


