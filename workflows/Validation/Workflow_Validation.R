# Workflow_Validation.R
# Master workflow for validation of NEON CH4 flux-gradient methods against
# co-located eddy-covariance towers (SE-Sto, SE-Svb, US-Uaf).
# Run scripts in order (01 -> 02 -> 02a -> 02b -> 05 -> 07). Each step sources
# a numbered workflow script. Steps 03/04 (diel CH4 analysis + its figures)
# have been removed -- see workflows/Validation/03_VAL_DielAnalysis.R for why.
#
# Steps 02a/02b are new: 02a estimates storage flux (single-point
# approximation, all three sites) and 02b adds it to the ensemble gradient
# flux to produce FG_total, which 05/07 now use in place of FG_mean wherever
# they compare against EC_mean (EC inherently includes storage; FG does not
# until 02b has run).

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(sf)

# -------------- Change this stuff -------------
DirRepo      <- "/Users/sm3466/Library/CloudStorage/Dropbox-YSE/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient"
DirRepo.ch4  <- "/Users/sm3466/Library/CloudStorage/Dropbox-YSE/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-methane"
DirRepo.eval <- "/Users/sm3466/Library/CloudStorage/Dropbox-YSE/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-eval"

setwd(DirRepo.ch4)
localdir     <- '/Volumes/MaloneLab/Research/FluxGradient'
localdir.ch4 <- '/Volumes/MaloneLab/Research/FluxGradient/Validation_Sites'
DnldFromGoogleDrive <- FALSE # TRUE = grab files from Google Drive; FALSE = use local copies

email     <- 'sparklelmalone@gmail.com'
googledrive::drive_auth(email = TRUE)
drive_url <- googledrive::as_id("https://drive.google.com/drive/folders/1Q99CT77DnqMl2mrUtuikcY47BFpckKw3")
data_folder <- googledrive::drive_ls(path = drive_url)

site.list <- c('SE-Sto', 'SE-Svb', 'US-Uaf')

# ── VALIDATION DATA PROCESSING ────────────────────────────────────────────────

# 01 Filter raw validation site gradient flux data and compute QC flags
message('Step 01: Filtering validation data...')
source(fs::path(DirRepo.ch4, 'workflows/Validation/01_VAL_FilterQC.R'))

fileSave <- fs::path(localdir.ch4, "SITEval_DATA_FILTERED_CH4.Rdata")
save(SITEval_DATA_FILTERED, file = fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir.ch4, "SITESval_MBR_9min_FILTER_CH4.Rdata")
save(SITESval_MBR_9min_FILTER, file = fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir.ch4, "SITESval_AE_9min_FILTER_CH4.Rdata")
save(SITESval_AE_9min_FILTER, file = fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir.ch4, "SITESval_WP_9min_FILTER_CH4.Rdata")
save(SITESval_WP_9min_FILTER, file = fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

# 02 Identify representative sensor height pairs (RSHP) for validation sites
message('Step 02: Identifying sensor height pairs...')
source(fs::path(DirRepo.ch4, 'workflows/Validation/02_VAL_SensorHeightPairs.R'))

# 02a Estimate storage flux for validation sites (single-point approximation,
# all three sites -- replaces the multi-level column-integration method,
# which never clears its 3-level minimum at SE-Sto/SE-Svb). Also includes
# the SE-Sto raw-unit correction (P_kPa hPa/kPa mislabel; CH4 ppm/ppb
# mislabel) -- see apply_site_unit_corrections() in flow.validation.storage.R.
message('Step 02a: Estimating storage flux (single-point)...')
source(fs::path(DirRepo, 'workflows/Validation/flow.validation.storage.R'))
source(fs::path(DirRepo, 'workflows/Validation/VAL_StorageFlux_SinglePoint.R'))

# 02b Combine ensemble gradient flux with single-point storage flux to
# produce FG_total (fallback to FG_mean alone where no storage estimate is
# available for that half-hour; storage_added flags which case applies).
message('Step 02b: Computing total flux (gradient + storage)...')
source(fs::path(DirRepo.ch4, 'workflows/Validation/02b_VAL_TotalFlux.R'))

# 03/04 REMOVED: diel CH4 analysis (brms Q10/Rref fitting) and its figures.
# FG-vs-EC agreement is fully covered by 05 and 07 without the brms MCMC
# cost. See workflows/Validation/03_VAL_DielAnalysis.R for the rationale and
# how to restore if the diurnal-pattern/Q10 comparison is needed later.

# 05 30-min, daily, and annual flux analysis for FG and EC (parallel pipelines)
message('Step 05: 30-min flux analysis...')
source(fs::path(DirRepo.ch4, 'workflows/Validation/05_VAL_FluxAnalysis.R'))

# 06 ERA5-driven half-hourly gap-filling for FG and EC fluxes at validation sites
message('Step 06: ERA5 half-hourly gap-fill...')
source(fs::path(DirRepo.ch4, 'workflows/Validation/06_VAL_ERA5Gapfill.R'))

# 07 Supplemental figures justifying the FG validation approach
message('Step 07: Validation supplement...')
source(fs::path(DirRepo.ch4, 'workflows/Validation/07_VAL_Supplement.R'))

# ── REFERENCE (not run in main workflow) ──────────────────────────────────────
# VAL_CounterGradientFilter.R — exploratory counter-gradient QC filter; not part
#   of the main validation pipeline but retained in this folder for reference.
# #source(fs::path(DirRepo.ch4,'workflows/07_NEON_SiteMap.R')) # site map if needed




