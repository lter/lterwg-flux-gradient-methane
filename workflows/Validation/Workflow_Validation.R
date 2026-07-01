# Workflow_Validation.R
# Master workflow for validation of NEON CH4 flux-gradient methods against
# co-located eddy-covariance towers (SE-Sto, SE-Svb, US-Uaf).
# Run scripts in order (01 -> 07). Each step sources a numbered workflow script.

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(sf)

# -------------- Change this stuff -------------
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

# 03 Diel CH4 flux analysis by season for validation sites
message('Step 03: Diel analysis...')
source(fs::path(DirRepo.ch4, 'workflows/Validation/03_VAL_DielAnalysis.R'))

# 04 Validation figures comparing FG and EC flux products
message('Step 04: Validation figures...')
source(fs::path(DirRepo.ch4, 'workflows/Validation/04_VAL_Figures.R'))

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




