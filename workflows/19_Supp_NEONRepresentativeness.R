# NEON site representativeness vs. global upland ecosystems.
#
# Question: for the environmental predictors that actually drive the RF
# source-probability and magnitude models (fit in
# 12_SourceProp_MagnitudeModels.R, applied to the global grid in
# 13_Global_SpatialUpscalingRF.R), how well does the NEON upland tower
# network cover the range of conditions
# found across the global upland grid the models are used to predict onto?
# A model can only be trusted where it is asked to extrapolate into
# conditions resembling its training data — this script quantifies that.
#
# Method: a Dissimilarity Index / Area-of-Applicability analysis, adapted
# from Meyer & Pebesma (2021, Methods Ecol Evol) and implemented here from
# first principles (no CAST dependency):
#
#   1. Pull each model's RF variable importance (already written by
#      12_SourceProp_MagnitudeModels.R) and combine the P(source) and
#      magnitude-model importances for the predictors they share (air
#      temperature, soil moisture, MAP, MAT), so continuous predictors are
#      weighted by how much BOTH models actually rely on them.
#   2. Standardize each continuous predictor (scaled by the NEON training
#      data's own SD) and build a variable-importance-weighted Euclidean
#      distance space.
#   3. Within each comparison group — Forest, Grassland, Shrubland, and
#      Arid — compute each global grid cell's Dissimilarity Index (DI): its
#      weighted distance to the nearest NEON reference site-month,
#      normalized by the typical NEON-to-NEON nearest-neighbor distance. A
#      DI of 1 means "as far from training data as a typical NEON site is
#      from its nearest NEON neighbor." Cells with DI beyond the training
#      data's own outlier threshold (Tukey rule: Q3 + 1.5*IQR of training
#      DI) are outside the Area of Applicability (AOA) — climatically novel.
#   4. "Arid" is never a raw EcoType label in NEON's site metadata (only
#      Forest/Grassland/Shrubland/Wetland/Cropland appear); the Arid
#      reference group is instead built from aridity_index < 15 (arid_ai_
#      threshold) — the same threshold 12_SourceProp_MagnitudeModels.R uses
#      to define is_arid for NEON's TRAINING data — regardless of a site's
#      raw EcoType text (this pipeline currently has 2 such sites).
#   5. The global grid's Arid category uses this identical aridity_index
#      threshold, not the biome/MODIS-raster "desert" class that script 13
#      uses on the grid side. That raster-based flag exists there to hard-
#      force P(source)=0 for desert-biome cells in Stage 1 — a flux-forcing
#      rule, not a claim about climate similarity — and for the
#      representative year it happened to tag zero cells, which silently
#      dropped Arid from this analysis entirely until this was corrected.
#      The aridity_index formula (MAP / (MAT+10)) is also numerically
#      degenerate as MAT approaches -10°C from above (denominator -> 0),
#      producing spuriously extreme finite values right at the boundary, not
#      just for MAT <= -10°C. Cells with MAT at or below aridity_mat_floor
#      (set above) are treated as undefined (NA) instead. This guard was
#      validated here first (it collapsed the global aridity_index range
#      from roughly [-5.5M, 698k] to [0.03, 1161]) and has since been
#      applied identically in 12_SourceProp_MagnitudeModels.R (training) and
#      13_Global_SpatialUpscalingRF.R (spatial-projection grid), so training,
#      the grid, and this analysis all use the same well-behaved formula now.
#
# Outputs (written to results_dir / figure_dir — METHANE/OUTPUT and
# METHANE/FIGURES, NOT the Upscaling_Monthly_RF-internal folders):
#   FIGURES/FigS_NEON_Representativeness.png
#   OUTPUT/neon_representativeness_variable_coverage.csv
#   OUTPUT/neon_representativeness_by_ecotype.csv
#   OUTPUT/neon_representativeness_summary.csv
#   OUTPUT/neon_representativeness_summary.txt
#   OUTPUT/neon_representativeness_supplemental_text.md  (manuscript-ready
#     narrative paragraph, numbers regenerated from the tables above on
#     every run — also embedded in the summary.txt)
#
# Prerequisites (written to output_dir, i.e. Upscaling_Monthly_RF/OUTPUT — a
# different folder from results_dir above; run 12 then 13 first):
#   <output_dir>/rf_class_variable_importance.csv       (12_SourceProp_MagnitudeModels.R)
#   <output_dir>/rf_magnitude_variable_importance.csv   (12_SourceProp_MagnitudeModels.R)
#   <output_dir>/era5_template.tif                       (13_Global_SpatialUpscalingRF.R)
#   <output_dir>/mat_map_climatology.tif   (13_Global_SpatialUpscalingRF.R; cached
#                                            MAT/MAP, recomputed here from
#                                            ERA5-Land if absent)
# Raw NEON + spatial inputs (same locations 12_SourceProp_MagnitudeModels.R
# and 13_Global_SpatialUpscalingRF.R read from):
#   30min_site_behavior.csv (also read by 12, for training),
#   NEON_ERA5_30min_site_covariates.csv.gz (covariates-only export used only
#   here — 12 reads the full gapfilled-flux file instead, see the note by
#   neon_covariates_file below),
#   Ecoregions2017.zip (or its already-extracted shapefile),
#   DATA/era5_land_monthly/, DATA/modis_mcd12c1_processed/, DATA/wad2m/
#   (cell membership in the "global upland grid" comes from MODIS land cover
#   alone — Forest/Grassland/Shrubland/Arid, i.e. everything ecotype_lookup
#   keeps after dropping Cropland. WAD2M inundation fraction is used only to
#   scale each surviving cell's area down by its exact upland fraction, not
#   to re-exclude it — see the note by the WAD2M loading block below for why.
#   Proceeds with no area correction at all, with a warning, if DATA/wad2m/
#   is absent)
#
# Requires the RANN package for nearest-neighbor search (fast, kd-tree
# based). Install with install.packages("RANN") if missing.

library(tidyverse)
library(data.table)
library(terra)
library(cowplot)

if (!requireNamespace("RANN", quietly = TRUE))
  stop("Package 'RANN' is required for nearest-neighbor search. Install with install.packages('RANN').")

# ── Paths (mirrors 12_SourceProp_MagnitudeModels.R and 13_Global_SpatialUpscalingRF.R) ──

localdir.ch4 <- Sys.getenv("LOCALDIR_CH4",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/Methane")
spatial_dir <- Sys.getenv("MONTHLY_UPSCALING_DIR",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly")
rf_dir <- Sys.getenv("MONTHLY_RF_DIR",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF")

# output_dir holds the RF-model prerequisites this script READS (variable
# importance CSVs, era5_template.tif) and the MAT/MAP climatology cache it
# shares with 13_Global_SpatialUpscalingRF.R — those live where script 13
# wrote them, under Upscaling_Monthly_RF/OUTPUT, and stay there.
output_dir <- file.path(rf_dir, "OUTPUT")
# This script's OWN results (the representativeness CSVs, summary.txt, and
# supplemental text) are written separately, to the shared METHANE-level
# OUTPUT folder rather than the RF-model-internal one above.
results_dir <- Sys.getenv("METHANE_OUTPUT_DIR",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/OUTPUT")
# Figure lives with the rest of the METHANE-level supplemental figures,
# not under Upscaling_Monthly_RF/FIGURES with the RF-model figures.
figure_dir <- Sys.getenv("METHANE_FIGURE_DIR",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/FIGURES")
data_dir   <- file.path(rf_dir, "DATA")
dir.create(output_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir,  showWarnings = FALSE, recursive = TRUE)

ecoregions_zip     <- file.path(spatial_dir,  "Ecoregions2017.zip")
era5_land_dir      <- file.path(spatial_dir,  "DATA/era5_land_monthly")
modis_ecotype_dir  <- file.path(spatial_dir,  "DATA/modis_mcd12c1_processed")
wad2m_dir          <- file.path(spatial_dir,  "DATA/wad2m")

site_behavior_file <- file.path(localdir.ch4, "OUTPUT/30min_site_behavior.csv")
# NOTE: 12_SourceProp_MagnitudeModels.R reads a "gapfilled_30min" flux file
# for site covariates (it needs the flux values themselves, for training).
# This script only needs the covariates (not flux), so it reads the
# dedicated covariates export directly.
neon_covariates_file <- file.path(localdir.ch4, "OUTPUT/NEON_ERA5_30min_site_covariates.csv.gz")

required_files <- c(
  file.path(output_dir, "rf_class_variable_importance.csv"),
  file.path(output_dir, "rf_magnitude_variable_importance.csv"),
  file.path(output_dir, "era5_template.tif"),
  site_behavior_file, neon_covariates_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0)
  stop("Missing required inputs (run 12_SourceProp_MagnitudeModels.R and ",
       "13_Global_SpatialUpscalingRF.R first if these are RF outputs): ",
       paste(missing_files, collapse = ", "))

# ── Scalar parameters (must match 12_SourceProp_MagnitudeModels.R) ────────────

years_to_process     <- 2000:2025
representative_year  <- max(years_to_process)
neon_upland_ecotypes  <- c("Forest", "Grassland", "Shrubland")  # what NEON actually samples

# arid_ai_threshold matches 12_SourceProp_MagnitudeModels.R's TRAINING-side
# definition: is_arid = aridity_index < 15 (MAP / (MAT+10)). This used to be
# the source of a real bug: 13_Global_SpatialUpscalingRF.R's GRID/prediction
# side originally used a DIFFERENT definition — is_arid = (EcoType ==
# "Arid") from the biome/MODIS raster — and hard-forced P(source) = 0 for
# those cells in Stage 1, which silently zeroed out the newly-real Arid
# category regardless of what the model had learned. That override has
# since been removed: 13_Global_SpatialUpscalingRF.R now classifies Arid the
# identical aridity_index-based way on the grid side and lets the trained
# Stage 1/Stage 2 models' own Arid behavior stand. This script applies the
# same aridity-INDEX threshold on both the NEON and global-grid sides, so
# "Arid" means the same thing everywhere here: apples to apples with how
# NEON's 2 arid-analog sites were identified.
arid_ai_threshold <- 15

# aridity_index = MAP / (MAT + 10) is only numerically well-behaved when
# MAT + 10 is comfortably positive. A plain MAT > -10 guard (used in an
# earlier version of this script) still lets cells just above that boundary
# blow up toward absurd finite values (e.g. an aridity index in the hundreds
# of thousands) as the denominator approaches zero — this floor buffers past
# the boundary, not just past zero. Validated here first (global
# aridity_index range collapsed from roughly [-5.5M, 698k] to [0.03, 1161]),
# then applied identically in 13_Global_SpatialUpscalingRF.R so training,
# the spatial-projection grid, and this analysis all treat the formula the
# same way. Cells with MAT below this floor get aridity_index = NA and are
# excluded via the is.finite() filters downstream.
aridity_mat_floor <- -9  # °C; i.e. require MAT + 10 > 1, not just > 0

# ── Helpers duplicated from 13_Global_SpatialUpscalingRF.R (kept in sync
#    manually; this script is meant to be self-contained) ────────────────────

unzip_if_needed <- function(zip_file, exdir, pattern) {
  if (!dir.exists(exdir)) dir.create(exdir, recursive = TRUE)
  if (length(list.files(exdir, pattern = pattern)) == 0)
    utils::unzip(zip_file, exdir = exdir)
  invisible(list.files(exdir, pattern = pattern, full.names = TRUE))
}

read_era5_land_year <- function(year) {
  file    <- file.path(era5_land_dir, sprintf("era5_land_monthly_%s.nc", year))
  climate <- rast(file); lnames <- names(climate)
  t2m_idx   <- grep("t2m",   lnames, ignore.case = TRUE)
  swvl1_idx <- grep("swvl1", lnames, ignore.case = TRUE)
  swvl2_idx <- grep("swvl2", lnames, ignore.case = TRUE)
  tp_idx    <- grep("(?i)^tp[_=]|^tp$", lnames, perl = TRUE)
  tavg <- climate[[t2m_idx]] - 273.15
  vswc <- (climate[[swvl1_idx]] + climate[[swvl2_idx]]) / 2
  prec <- rast(lapply(seq_along(tp_idx), function(i)
    climate[[tp_idx[i]]] * 1000 * lubridate::days_in_month(lubridate::make_date(year, i, 1))))
  names(tavg) <- sprintf("tavg_%02d", 1:12)
  names(vswc) <- sprintf("vswc_%02d", 1:12)
  names(prec) <- sprintf("prec_%02d", 1:12)
  list(tavg = tavg, vswc = vswc, prec = prec)
}

# ── Static rasters: template, MAT/MAP climatology, ecotype, WAD2M inundation ──

template <- rast(file.path(output_dir, "era5_template.tif"))

# WAD2M inundation fraction — identical loading logic to
# 13_Global_SpatialUpscalingRF.R's get_inundation_fraction(). Used only to
# weight each surviving cell's area down by however much of it is wetland,
# matching how the upscaling model itself only ever predicts onto the UPLAND
# fraction of a cell, not the whole cell. It is NOT used to decide whether a
# cell belongs in the grid at all — that membership call is MODIS land
# cover's (Forest/Grassland/Shrubland/Arid-i.e.-barren, via eco_r/
# ecotype_lookup below). WAD2M is a coarse (0.5deg) wetland-extent model and
# can mark a cell majority-inundated — e.g. floodplain grassland/shrubland,
# scattered wetlands within a forest matrix — even where MODIS's finer-
# resolution land cover correctly shows upland vegetation; hard-excluding on
# WAD2M alone dropped real upland area from earlier versions of this script.
# So: trust land cover for membership, use WAD2M only for the area
# correction (a cell that's biome-classified "Grassland" but 40% seasonally
# inundated contributes only 60% of its area as upland-representative).
wad2m_files <- list.files(wad2m_dir, pattern = "\\.nc$", full.names = TRUE)
use_wad2m   <- length(wad2m_files) > 0
wad2m       <- if (use_wad2m) rast(wad2m_files[1]) else NULL
wad2m_yrs   <- if (use_wad2m) 2000:2020 else integer()
if (!use_wad2m) message("No WAD2M inundation data found at ", wad2m_dir,
  " — proceeding without an inundation correction (every cell treated as fully upland).")

get_inundation_fraction <- function(year, month) {
  if (!use_wad2m) { r <- template; values(r) <- 0; names(r) <- "inundation_fraction"; return(r) }
  if (year %in% wad2m_yrs) {
    r <- wad2m[[(year - min(wad2m_yrs)) * 12 + month]]
  } else {
    r <- mean(wad2m[[seq(month, nlyr(wad2m), by = 12)]], na.rm = TRUE)
  }
  if (!compareGeom(r, template, stopOnError = FALSE)) r <- resample(r, template, method = "bilinear")
  names(r) <- "inundation_fraction"; r
}

mat_map_cache <- file.path(output_dir, "mat_map_climatology.tif")
if (file.exists(mat_map_cache)) {
  message("Loading cached MAT/MAP climatology...")
  mat_map <- rast(mat_map_cache)
  mat <- mat_map[["MAT"]]; map <- mat_map[["MAP"]]
} else {
  message("No cached MAT/MAP climatology found — computing from ", min(years_to_process),
          "-", max(years_to_process), " ERA5-Land (one-time cost; re-run ",
          "13_Global_SpatialUpscalingRF.R to cache this for next time)...")
  mat_sum <- map_sum <- template
  values(mat_sum) <- 0; values(map_sum) <- 0; n_layers <- 0
  for (yr in years_to_process) {
    cy <- read_era5_land_year(yr)
    mat_sum <- mat_sum + sum(cy$tavg, na.rm = TRUE)
    map_sum <- map_sum + sum(cy$prec, na.rm = TRUE)
    n_layers <- n_layers + 12
  }
  mat <- mat_sum / n_layers; names(mat) <- "MAT"
  map <- map_sum / length(years_to_process); names(map) <- "MAP"
  writeRaster(c(mat, map), mat_map_cache, overwrite = TRUE)
}

ecoregion_files <- unzip_if_needed(ecoregions_zip,
  file.path(data_dir, "ecoregions2017"), "Ecoregions2017\\.shp$")
ecoregions   <- vect(ecoregion_files[1])
biome_field  <- intersect(c("BIOME_NUM","BIOME"), names(ecoregions))[1]
biome_raster <- rasterize(ecoregions, template, field = biome_field, touches = TRUE)
ecotype_raster <- subst(biome_raster, from = 1:14,
  to = c(2,2,2,2,2,2, 3,3,NA,3, 4,4,4,NA))
names(ecotype_raster) <- "ecotype_code"

ecotype_lookup <- tibble(
  ecotype_code = c(1,2,3,4,5),
  EcoType      = c("Cropland","Forest","Grassland","Shrubland","Arid")
) %>% filter(EcoType != "Cropland")

modis_files <- list.files(modis_ecotype_dir,
  pattern = "^MODIS_MCD12C1_ecotype_[0-9]{4}\\.tif$", full.names = TRUE)
if (length(modis_files) > 0) {
  modis_years <- as.integer(str_match(basename(modis_files), "([0-9]{4})\\.tif$")[,2])
  sel   <- modis_years[which.min(abs(modis_years - representative_year))]
  eco_r <- rast(file.path(modis_ecotype_dir, sprintf("MODIS_MCD12C1_ecotype_%s.tif", sel)))
  if (!same.crs(eco_r, template)) crs(eco_r) <- crs(template)
  if (!compareGeom(eco_r, template, stopOnError = FALSE)) eco_r <- resample(eco_r, template, method = "near")
  names(eco_r) <- "ecotype_code"
} else {
  eco_r <- ecotype_raster
}

cell_area_mha <- cellSize(template, unit = "m") / 1e10
names(cell_area_mha) <- "area_mha"

# ── Global upland grid: representative-year annual-mean covariates ────────────
# Same predictor set used to train the models (EcoType, mean air temperature,
# mean soil moisture, MAP, MAT, aridity index), built for one representative
# year rather than the full 2000-2025 x 12-month stack that
# 13_Global_SpatialUpscalingRF.R predicts onto — sufficient to characterize
# the shape of the global predictor space without re-running the full loop.

climate_year <- read_era5_land_year(representative_year)
tavg_annual  <- mean(climate_year$tavg, na.rm = TRUE); names(tavg_annual) <- "mean_ERA5_Tair_C"
vswc_annual  <- mean(climate_year$vswc, na.rm = TRUE); names(vswc_annual) <- "mean_ERA5_VSWC"

# Annual-mean inundation fraction for the representative year (mirrors the
# annual-averaging already used for Tair/VSWC above).
inund_annual <- mean(rast(lapply(1:12, function(m) get_inundation_fraction(representative_year, m))),
                     na.rm = TRUE)
names(inund_annual) <- "inundation_fraction"
# WAD2M doesn't cover every terrestrial cell (e.g. no data at that pixel in
# any month) — those come out of the mean() above as NA/NaN. Left as NA,
# they'd silently drop the whole cell below via as.data.frame(..., na.rm =
# TRUE), even when MODIS land cover (eco_r) confirms it's real upland
# vegetation — the exact bug this is fixing. Missing WAD2M coverage is not
# evidence of inundation, so treat it as 0 (no correction), matching how
# 13_Global_SpatialUpscalingRF.R's grid loop already handles this
# (`inund[is.na(inund)] <- 0`).
inund_annual[is.na(inund_annual)] <- 0

stk <- c(eco_r, tavg_annual, vswc_annual, mat, map, cell_area_mha, inund_annual)
global_grid <- as.data.frame(stk, xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
  mutate(ecotype_code = as.integer(ecotype_code)) %>%
  inner_join(ecotype_lookup, by = "ecotype_code") %>%
  mutate(
    # See aridity_mat_floor above: MAT below this floor is treated as
    # undefined (NA) rather than left to blow up toward a spurious extreme
    # finite value as MAT + 10 -> 0.
    aridity_index         = if_else(MAT > aridity_mat_floor, MAP / (MAT + 10), NA_real_),
    inundation_fraction   = pmin(pmax(inundation_fraction, 0), 1),
    upland_area_fraction  = 1 - inundation_fraction,
    # Same weighting 13_Global_SpatialUpscalingRF.R applies before predicting:
    # a cell that's 40% seasonally inundated only contributes 60% of its
    # area as "upland." This is a down-WEIGHT only — membership in the grid
    # already came from MODIS land cover via the inner_join above (Forest/
    # Grassland/Shrubland/Arid), so WAD2M inundation never re-excludes a
    # cell here (see the note by the WAD2M loading block for why).
    area_mha              = area_mha * upland_area_fraction,
    # Arid defined via the SAME aridity_index < 15 threshold used for NEON
    # training (see arid_ai_threshold note above) and, since the Arid-
    # classification fix, by 13_Global_SpatialUpscalingRF.R's own grid-side
    # classification too — not from the biome/MODIS raster's own "desert"
    # class (eco_r/ecotype_lookup), which is used only to decide grid CELL
    # MEMBERSHIP (is this cell upland at all), not whether it counts as
    # climatically Arid. A biome-classified desert cell that doesn't cross
    # this aridity_index threshold folds into Shrubland here, mirroring
    # 13's own EcoType remap for its spatial-projection grid — neither
    # script hard-forces P(source)=0 for these cells anymore; the models'
    # own trained Arid behavior is used instead.
    is_arid                = as.integer(!is.na(aridity_index) & aridity_index < arid_ai_threshold),
    EcoType                = if_else(is_arid == 1, "Arid",
                                     if_else(EcoType == "Arid", "Shrubland", EcoType))
  ) %>%
  filter(area_mha > 0, is.finite(aridity_index),
         is.finite(mean_ERA5_Tair_C), is.finite(mean_ERA5_VSWC))

message(sprintf(
  "Global upland grid: %d cells (representative year %d), %.1f Mha total (MODIS-land-cover upland cells, WAD2M-area-weighted, %s); %d cells (%.1f Mha) classified Arid via aridity_index < %d",
  nrow(global_grid), representative_year, sum(global_grid$area_mha),
  if (use_wad2m) "WAD2M available" else "WAD2M NOT found, no inundation correction applied",
  sum(global_grid$is_arid == 1), sum(global_grid$area_mha[global_grid$is_arid == 1]), arid_ai_threshold))

# ── NEON training predictor space (site-months, upland only) ──────────────────
# Built directly from the covariates export rather than replaying
# 12_SourceProp_MagnitudeModels.R's flux-labeling logic, since only the
# predictor values (not the flux/weak_source label) matter here.

site_climate_normals <- read.csv(site_behavior_file) %>%
  mutate(SITE_ID = as.character(SITE_ID)) %>%
  distinct(SITE_ID, MAP, MAT)

neon_covariates <- data.table::fread(neon_covariates_file) %>% as_tibble() %>%
  mutate(SITE_ID = as.character(SITE_ID), Year = as.integer(Year), month = as.integer(month),
         ERA5_Tair_C = as.numeric(ERA5_Tair_C), ERA5_VSWC = as.numeric(ERA5_VSWC)) %>%
  filter(is.finite(Year), is.finite(month), is.finite(ERA5_Tair_C), is.finite(ERA5_VSWC),
         # Matches the EcoType levels the RF models were actually trained
         # with (12_SourceProp_MagnitudeModels.R filters to these three raw
         # EcoType levels before the Arid recode; NEON's 2 Cropland and 2
         # Wetland sites are excluded here for the same reason they drop
         # out of training).
         EcoType %in% neon_upland_ecotypes)

neon_predictors <- neon_covariates %>%
  reframe(.by = c(SITE_ID, EcoType, Year, month),
    mean_ERA5_Tair_C = mean(ERA5_Tair_C, na.rm = TRUE),
    mean_ERA5_VSWC   = mean(ERA5_VSWC,   na.rm = TRUE)) %>%
  inner_join(site_climate_normals, by = "SITE_ID") %>%
  # Same aridity_mat_floor guard applied to the global grid above — a no-op
  # for NEON in practice (no site is anywhere near that cold), kept for
  # consistency since both sides should treat the formula identically.
  mutate(aridity_index = if_else(MAT > aridity_mat_floor, MAP / (MAT + 10), NA_real_)) %>%
  filter(is.finite(MAP), is.finite(MAT), is.finite(aridity_index))

# "Arid" is never a raw EcoType label in NEON's site metadata (only
# Forest/Grassland/Shrubland/Wetland/Cropland) — but 13_Global_SpatialUpscalingRF.R
# doesn't define aridity from that label either. It flags a cell/site as
# arid from the SAME aridity_index threshold used here: aridity_index < 15
# (arid_ai_threshold, set above), independent of the EcoType text label. A
# NEON site can be labeled "Shrubland" and still cross that threshold.
# Comparing against global Arid cells therefore has to use the aridity-
# index-derived flag, not the (always-absent) text label — otherwise every
# Arid grid cell looks like it has zero NEON analogs, which is wrong.
neon_predictors <- neon_predictors %>%
  mutate(
    is_arid          = as.integer(!is.na(aridity_index) & aridity_index < arid_ai_threshold),
    comparison_group = if_else(is_arid == 1, "Arid", EcoType)
  )

n_neon_sites      <- n_distinct(neon_predictors$SITE_ID)
n_arid_neon_sites <- neon_predictors %>% filter(is_arid == 1) %>% pull(SITE_ID) %>% n_distinct()
message(sprintf(
  "NEON training predictor space: %d site-months across %d upland sites (%s); %d of those sites cross the aridity_index < %d threshold and are treated as the Arid reference group",
  nrow(neon_predictors), n_neon_sites, paste(neon_upland_ecotypes, collapse = ", "),
  n_arid_neon_sites, arid_ai_threshold))

# ── Combined variable-importance weights ───────────────────────────────────────
# The P(source) model's importance and the magnitude models' importance
# (Weak-sink + Weak-source, whose predictors are the same variables
# z-standardized) are normalized within each model to sum to 1, then
# averaged per shared raw variable, then renormalized — so a variable that
# matters a lot to any of the three RF fits pulls more weight in the
# distance metric below.

rf_class_importance <- read.csv(file.path(output_dir, "rf_class_variable_importance.csv"))
rf_mag_importance   <- read.csv(file.path(output_dir, "rf_magnitude_variable_importance.csv"))

predictor_to_raw <- c(
  "mean_ERA5_Tair_C" = "Tair", "z_Tair" = "Tair",
  "mean_ERA5_VSWC"   = "VSWC", "z_VSWC" = "VSWC",
  "MAP" = "MAP", "z_MAP" = "MAP",
  "MAT" = "MAT", "z_MAT" = "MAT",
  "aridity_index" = "aridity_index",
  "EcoType" = "EcoType"
  # source_probability and is_arid are model-derived quantities, not raw
  # site conditions, and are excluded from this comparison.
)

var_lookup <- tibble(
  raw_variable = c("Tair", "VSWC", "MAP", "MAT", "aridity_index", "EcoType"),
  column       = c("mean_ERA5_Tair_C", "mean_ERA5_VSWC", "MAP", "MAT", "aridity_index", "EcoType"),
  pretty_label = c("Air temperature (°C)", "Soil moisture (m³ m⁻³)",
                   "Mean annual precip. (mm)", "Mean annual temp. (°C)",
                   "Aridity index (MAP / (MAT+10))", "Ecosystem type")
)

normalize_importance <- function(df) {
  df %>% group_by(model) %>% mutate(importance_norm = importance / sum(importance)) %>% ungroup()
}

combined_importance <- bind_rows(
    normalize_importance(rf_class_importance),
    normalize_importance(rf_mag_importance)
  ) %>%
  mutate(raw_variable = predictor_to_raw[predictor]) %>%
  filter(!is.na(raw_variable)) %>%
  group_by(raw_variable) %>%
  summarise(weight = mean(importance_norm), .groups = "drop") %>%
  mutate(weight = weight / sum(weight)) %>%
  left_join(var_lookup, by = "raw_variable") %>%
  arrange(desc(weight))

message("Combined, renormalized variable weights (P(source) + both magnitude models):")
print(combined_importance %>% select(raw_variable, pretty_label, weight))

continuous_vars <- combined_importance$raw_variable[combined_importance$raw_variable != "EcoType"]
z_cols <- paste0("z_", continuous_vars)
w <- combined_importance$weight[match(continuous_vars, combined_importance$raw_variable)]
w <- w / sum(w)  # renormalized among continuous predictors only (EcoType handled by stratification below)

# ── Standardize continuous predictors on the NEON training data's own SD ──────

var_cols <- var_lookup$column[match(continuous_vars, var_lookup$raw_variable)]
std_stats <- neon_predictors %>%
  summarise(across(all_of(var_cols), list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))))

standardize <- function(df) {
  for (i in seq_along(continuous_vars)) {
    col <- var_cols[i]; zc <- z_cols[i]
    m <- std_stats[[paste0(col, "_mean")]]; s <- std_stats[[paste0(col, "_sd")]]
    df[[zc]] <- (df[[col]] - m) / s
  }
  df
}

neon_z   <- standardize(neon_predictors)
global_z <- standardize(global_grid)

# ── Dissimilarity Index / Area of Applicability, stratified by comparison
#    group (Forest / Grassland / Shrubland / Arid) ─────────────────────────
# (see header for the method summary and the Meyer & Pebesma 2021 reference)
#
# The reference (NEON) pool for each group is neon_z$comparison_group —
# Arid pulls from NEON site-months that cross the aridity_index threshold,
# regardless of their raw EcoType label (see note above). The query (global
# grid) pool is global_z$EcoType, which already applies the same
# aridity_index-based Arid remap used in the global_grid construction above
# (a biome-classified desert cell that doesn't cross the threshold folds
# into Shrubland; one that does becomes Arid) — this mirrors
# 13_Global_SpatialUpscalingRF.R's own grid-side classification, so
# training, the spatial-projection grid, and this analysis all agree on
# what counts as Arid.

comparison_groups <- c("Forest", "Grassland", "Shrubland", "Arid")

compute_di_for_group <- function(grp) {
  ref <- neon_z  %>% filter(comparison_group == grp)
  qry <- global_z %>% filter(EcoType == grp)
  if (nrow(qry) == 0) return(NULL)
  if (nrow(ref) < 2) {
    message(sprintf("  %s: only %d NEON reference site-month(s) — DI left undefined for this group.",
      grp, nrow(ref)))
    return(qry %>% mutate(DI = NA_real_, aoa_threshold = NA_real_, inside_AOA = FALSE,
                           ecotype_has_neon_analog = FALSE))
  }

  ref_mat <- as.matrix(ref[, z_cols])  %*% diag(sqrt(w), nrow = length(w))
  qry_mat <- as.matrix(qry[, z_cols])  %*% diag(sqrt(w), nrow = length(w))

  nn_train   <- RANN::nn2(ref_mat, ref_mat, k = 2)  # k=2: nearest OTHER training point
  d_train_nn <- nn_train$nn.dists[, 2]
  norm_const <- mean(d_train_nn)
  di_train   <- d_train_nn / norm_const
  threshold_val <- unname(quantile(di_train, 0.75) + 1.5 * IQR(di_train))

  nn_query <- RANN::nn2(ref_mat, qry_mat, k = 1)
  qry %>% mutate(
    DI                     = as.numeric(nn_query$nn.dists) / norm_const,
    aoa_threshold          = threshold_val,
    inside_AOA             = DI <= threshold_val,
    ecotype_has_neon_analog = TRUE
  )
}

global_grid_di <- map(comparison_groups, compute_di_for_group) %>% list_rbind() %>%
  mutate(representativeness = case_when(
    !ecotype_has_neon_analog ~ "No NEON EcoType analog",
    inside_AOA               ~ "Within AOA",
    TRUE                     ~ "Outside AOA / novel"
  ))

# ── Summary tables ─────────────────────────────────────────────────────────────

representativeness_by_ecotype <- global_grid_di %>%
  group_by(EcoType) %>%
  summarise(
    total_area_mha             = sum(area_mha),
    area_well_represented_mha  = sum(area_mha[ecotype_has_neon_analog & inside_AOA]),
    area_outside_aoa_mha       = sum(area_mha[ecotype_has_neon_analog & !inside_AOA]),
    area_no_ecotype_analog_mha = sum(area_mha[!ecotype_has_neon_analog]),
    pct_well_represented       = 100 * area_well_represented_mha / total_area_mha,
    .groups = "drop"
  ) %>%
  left_join(
    neon_predictors %>% group_by(comparison_group) %>%
      summarise(n_neon_sites = n_distinct(SITE_ID), .groups = "drop") %>%
      rename(EcoType = comparison_group),
    by = "EcoType"
  ) %>%
  mutate(n_neon_sites = replace_na(n_neon_sites, 0L)) %>%
  arrange(desc(total_area_mha))

overall_summary <- tibble(
  representative_year   = representative_year,
  n_neon_sites           = n_neon_sites,
  n_neon_site_months     = nrow(neon_predictors),
  n_global_cells         = nrow(global_grid_di),
  total_upland_area_mha  = sum(global_grid_di$area_mha),
  pct_area_well_represented = 100 * sum(global_grid_di$area_mha[global_grid_di$ecotype_has_neon_analog & global_grid_di$inside_AOA]) /
                               sum(global_grid_di$area_mha),
  pct_area_outside_aoa_same_ecotype = 100 * sum(global_grid_di$area_mha[global_grid_di$ecotype_has_neon_analog & !global_grid_di$inside_AOA]) /
                               sum(global_grid_di$area_mha),
  pct_area_no_ecotype_analog = 100 * sum(global_grid_di$area_mha[!global_grid_di$ecotype_has_neon_analog]) /
                               sum(global_grid_di$area_mha)
)

variable_coverage <- map(continuous_vars, function(v) {
  col <- var_lookup$column[var_lookup$raw_variable == v]
  neon_vals   <- neon_predictors[[col]]
  neon_range  <- range(neon_vals, na.rm = TRUE)
  neon_p05_95 <- quantile(neon_vals, c(0.05, 0.95), na.rm = TRUE)
  gv <- global_grid_di[[col]]; ga <- global_grid_di$area_mha
  tibble(
    raw_variable  = v,
    pretty_label  = var_lookup$pretty_label[var_lookup$raw_variable == v],
    weight        = combined_importance$weight[combined_importance$raw_variable == v],
    neon_min = neon_range[1], neon_max = neon_range[2],
    neon_p05 = unname(neon_p05_95[1]), neon_p95 = unname(neon_p05_95[2]),
    global_min = min(gv, na.rm = TRUE), global_max = max(gv, na.rm = TRUE),
    pct_area_within_neon_range  = 100 * sum(ga[gv >= neon_range[1]  & gv <= neon_range[2]],  na.rm = TRUE) / sum(ga),
    pct_area_within_neon_p05_95 = 100 * sum(ga[gv >= neon_p05_95[1] & gv <= neon_p05_95[2]], na.rm = TRUE) / sum(ga)
  )
}) %>% list_rbind() %>% arrange(desc(weight))

write.csv(variable_coverage,               file.path(results_dir, "neon_representativeness_variable_coverage.csv"), row.names = FALSE)
write.csv(representativeness_by_ecotype,   file.path(results_dir, "neon_representativeness_by_ecotype.csv"),        row.names = FALSE)
write.csv(overall_summary,                 file.path(results_dir, "neon_representativeness_summary.csv"),           row.names = FALSE)

# ── Supplemental text (narrative, regenerated from the tables above every run) ─
# A manuscript-ready paragraph with the key numbers filled in from
# overall_summary / representativeness_by_ecotype / variable_coverage, so it
# never drifts out of sync with the CSVs/figure — re-running this script
# regenerates the numbers in place. Follows the same "Supplemental Text —"
# convention as 16_Supp_GapfillUncertainty.R / 17_Supp_BudgetUncertainty.R,
# but (unlike those) is also written to disk, not just printed to console.

best_ecotype  <- representativeness_by_ecotype %>% arrange(desc(pct_well_represented)) %>% slice(1)
worst_ecotype <- representativeness_by_ecotype %>% arrange(pct_well_represented)       %>% slice(1)
top_dim       <- combined_importance %>% arrange(desc(weight)) %>% slice(1)
top_continuous_var <- variable_coverage %>% slice(1)  # variable_coverage already arrange(desc(weight))

supp_text <- paste(strwrap(sprintf(
"Supplemental Text — Representativeness of the NEON upland tower network relative to global upland ecosystems.

To evaluate how well the NEON CH4 flux-gradient network represents the range of
environmental conditions the spatially upscaled Random Forest models (source-
probability and magnitude, see Methods) are asked to extrapolate onto, we computed
a Dissimilarity Index (DI) and Area of Applicability (AOA) for the global upland grid
following Meyer & Pebesma (2021), using a Euclidean distance in predictor space
standardized by the NEON training data and weighted by each predictor's combined
(source-probability + magnitude model) variable importance. %s was the single most
influential predictor (weight = %.2f), followed by %s (weight = %.2f).

Of %.0f Mha of global upland area (%d grid cells, representative year %d), %.1f%%
fell within the AOA of the %d NEON sites (%d site-months across Forest, Grassland,
Shrubland, and aridity-index-defined Arid reference groups) and is therefore
considered climatically well represented by the training data; the remaining %.1f%%
fell outside the AOA within an ecosystem type NEON does sample (climatically novel
conditions within a familiar biome), and %.1f%% had no NEON analog for that
ecosystem type at all.

Representativeness varied by ecosystem type: %s was best represented (%.1f%% of its
%.0f Mha within the AOA, n = %d NEON sites), while %s was least represented (%.1f%%
of its %.0f Mha within the AOA, n = %d NEON sites). Despite this joint (multivariate)
representativeness gap, individual predictors showed much higher marginal coverage
in isolation: %s, the highest-weighted continuous predictor (weight = %.2f), had
%.0f%% of global upland area falling within NEON's univariate range. This gap between
high marginal (single-variable) coverage and low joint (multivariate) representativeness
is the reason a per-variable range check is not a substitute for the DI/AOA analysis,
and indicates that global upscaled predictions should be interpreted with the greatest
caution in regions and combinations of conditions falling outside the AOA, even where
any one driver individually resembles a NEON site.",
  top_dim$pretty_label[1], top_dim$weight[1],
  combined_importance$pretty_label[2], combined_importance$weight[2],
  overall_summary$total_upland_area_mha, overall_summary$n_global_cells, overall_summary$representative_year,
  overall_summary$pct_area_well_represented,
  overall_summary$n_neon_sites, overall_summary$n_neon_site_months,
  overall_summary$pct_area_outside_aoa_same_ecotype,
  overall_summary$pct_area_no_ecotype_analog,
  best_ecotype$EcoType[1],  best_ecotype$pct_well_represented[1],  best_ecotype$total_area_mha[1],  best_ecotype$n_neon_sites[1],
  worst_ecotype$EcoType[1], worst_ecotype$pct_well_represented[1], worst_ecotype$total_area_mha[1], worst_ecotype$n_neon_sites[1],
  top_continuous_var$pretty_label[1], top_continuous_var$weight[1], top_continuous_var$pct_area_within_neon_range[1]
), width = 90), collapse = "\n")

writeLines(supp_text, file.path(results_dir, "neon_representativeness_supplemental_text.md"))
cat("\n", supp_text, "\n\n", sep = "")

capture.output({
  cat("NEON site representativeness vs. global upland ecosystems\n\n")
  cat(sprintf("Representative year: %d | NEON sites: %d (%d site-months) | Global cells: %d\n\n",
    representative_year, n_neon_sites, nrow(neon_predictors), nrow(global_grid_di)))
  cat("--- Combined variable weights ---\n"); print(combined_importance %>% select(raw_variable, pretty_label, weight))
  cat("\n--- Univariate coverage ---\n"); print(variable_coverage)
  cat("\n--- Representativeness by EcoType ---\n"); print(representativeness_by_ecotype)
  cat("\n--- Overall summary ---\n"); print(overall_summary)
  cat("\n--- Supplemental Text ---\n\n"); cat(supp_text, "\n")
}, file = file.path(results_dir, "neon_representativeness_summary.txt"))

# ── Figures ─────────────────────────────────────────────────────────────────

# plot.title.position = "plot" anchors each title to the full plot width
# (including the axis-label gutter), not just the panel area — without it,
# titles on narrower panels (e.g. A's coord_flip bar chart, whose y-axis
# text eats into the panel width) get right-clipped by the plot device.
fig_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 12),
        plot.title.position = "plot", plot.margin = margin(t = 5.5, r = 10, b = 5.5, l = 5.5),
        strip.background = element_rect(fill = "grey92"), strip.text = element_text(size = 10),
        axis.title = element_text(size = 11), axis.text = element_text(size = 9))

# aoa_colors (panels B/D) vs. dist_colors (panels A/C) are deliberately
# built from two non-overlapping hue families — blue/orange/dark-red for
# AOA membership, purple/teal for NEON-vs-global predictor comparisons — so
# a color can never mean "AOA status" in one panel and "data source" in
# another within the same figure.
aoa_colors  <- c("Within AOA" = "#2166AC", "Outside AOA / novel" = "#D95F02", "No NEON EcoType analog" = "#7A1216")
dist_colors <- c("NEON site-months" = "#762A83", "Global upland grid" = "#1B9E77")

# Note on object names vs. printed panel letters: the R objects below are
# named pA/pB/pC/pD for what each one computes (importance / predictor
# coverage / AOA map / EcoType bars), but the LETTER PRINTED ON THE FIGURE
# doesn't follow that same order — it follows physical position in
# fig_top/fig_repres below (top-left, top-right, then down), so the object
# called pD is the one titled "B." and the object called pB is titled "C."
# Keep the panel letter in each labs(title=...)/draw_label(...) string in
# sync with where plot_grid actually places it, not with the object name.

# A: combined variable importance weights driving the distance metric. Uses
# dist_colors' NEON site-months color for the continuous predictors, since
# this panel and the one titled "C" are two views of the same
# NEON-training-data story (what drives the distance metric vs. how NEON's
# own values are distributed); EcoType is called out separately in grey
# because it enters the model as a stratifying factor, not a continuous
# distance dimension.
pA <- combined_importance %>%
  mutate(pretty_label = fct_reorder(pretty_label, weight)) %>%
  ggplot(aes(x = pretty_label, y = weight, fill = raw_variable == "EcoType")) +
  geom_col(width = 0.65, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "grey45", `FALSE` = unname(dist_colors["NEON site-months"]))) +
  coord_flip() +
  labs(title = "A.",
       x = NULL, y = "Normalized weight (P(source) + magnitude models)") +
  fig_theme

# pB object -> printed as panel "C": NEON vs. global (area-weighted)
# distributions, ordered by importance. Makes the point the whole figure is
# built around: NEON's univariate (single-variable) coverage of the global
# range can be high — printed directly on each panel as "% of global area
# within NEON range", from variable_coverage's pct_area_within_neon_range —
# even where panel D's joint, multivariate AOA membership says most area is
# still climatically novel. High marginal overlap is fully compatible with
# low joint representativeness; that gap is the reason a per-variable range
# check isn't a substitute for the DI/AOA analysis in D.
zoom_vars <- c("aridity_index", "MAP")  # long right tails; zoom to see the bulk of the distribution

# Simple weighted-quantile via linear interpolation on the cumulative weight
# — good enough for setting a plot axis limit, not used for any reported stat.
weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w)
  x <- x[keep]; w <- w[keep]
  o <- order(x); x <- x[o]; w <- w[o]
  cum_w <- cumsum(w) / sum(w)
  approx(cum_w, x, xout = probs, rule = 2)$y
}

dist_compare <- map(continuous_vars, function(v) {
  col <- var_lookup$column[var_lookup$raw_variable == v]
  bind_rows(
    tibble(value = neon_predictors[[col]], weight = 1,
           source = "NEON site-months", raw_variable = v),
    tibble(value = global_grid_di[[col]], weight = global_grid_di$area_mha,
           source = "Global upland grid", raw_variable = v)
  )
}) %>% list_rbind() %>%
  left_join(var_lookup %>% select(raw_variable, pretty_label), by = "raw_variable") %>%
  mutate(pretty_label = factor(pretty_label,
    levels = var_lookup$pretty_label[match(continuous_vars, var_lookup$raw_variable)]))

# Built as separate small-multiple panels (rather than one facet_wrap) so
# aridity_index and MAP can each get their own zoomed x-axis via
# coord_cartesian — which clips the viewport only, after the density is
# estimated over the FULL data (unlike scale_x_continuous(limits=...), which
# would drop the tail before the KDE runs and bias the visible shape) — and
# so each panel can carry its own %-of-area-covered text.
make_dist_panel <- function(v) {
  col         <- var_lookup$column[var_lookup$raw_variable == v]
  label       <- var_lookup$pretty_label[var_lookup$raw_variable == v]
  d           <- dist_compare %>% filter(raw_variable == v)
  pct_covered <- variable_coverage$pct_area_within_neon_range[variable_coverage$raw_variable == v]

  p <- ggplot(d, aes(x = value, weight = weight, fill = source, color = source)) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    scale_fill_manual(values = dist_colors, breaks = names(dist_colors), name = NULL, drop = TRUE) +
    scale_color_manual(values = dist_colors, breaks = names(dist_colors), name = NULL, drop = TRUE) +
    labs(title = label, x = NULL, y = "Density") +
    fig_theme + theme(plot.title = element_text(face = "plain", size = 10))

  if (v %in% zoom_vars) {
    neon_q <- quantile(neon_predictors[[col]], c(0.025, 0.975), na.rm = TRUE)
    glob_q <- weighted_quantile(global_grid_di[[col]], global_grid_di$area_mha, c(0.025, 0.975))
    p <- p + coord_cartesian(xlim = range(c(neon_q, glob_q)))
  }

  p + annotate("text", x = Inf, y = Inf, hjust = 1.2, vjust = 1.4, size = 3.5, fontface = "bold",
               label = sprintf("%.0f%%", pct_covered))
}

dist_panels  <- map(continuous_vars, make_dist_panel)
legend_b     <- get_legend(dist_panels[[1]] + theme(legend.position = "top",
                                                     legend.background = element_rect(fill = "transparent", color = NA),
                                                     legend.key        = element_rect(fill = "transparent", color = NA)))
dist_panels  <- map(dist_panels, ~ . + theme(legend.position = "none"))

pB <- plot_grid(
  ggdraw() + draw_label("C. ",
                         fontface = "bold", size = 12, x = 0.01, hjust = 0),
  legend_b,
  plot_grid(plotlist = dist_panels, ncol = 3),
  ncol = 1, rel_heights = c(0.06, 0.08, 1)
)

# pC object -> printed as panel "D": DI / Area-of-Applicability map
map_df <- global_grid_di  # representativeness classification already attached above

pC <- map_df %>%
  ggplot(aes(x = x, y = y, fill = representativeness)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  coord_equal(expand = FALSE) +
  scale_fill_manual(values = aoa_colors, name = NULL) +
  labs(title = "D.",
       x = "Longitude", y = "Latitude") +
  fig_theme + theme(panel.grid = element_blank(), legend.position = "bottom") +
  geom_label(data = overall_summary,
             aes(x = -Inf, y = Inf, label = sprintf("%.1f%%", pct_area_well_represented)),
             inherit.aes = FALSE, hjust = -0.15, vjust = 1.4, size = 4.5, fontface = "bold",
             color = "black", fill = "white", label.size = 0, alpha = 0.85)

# pD object -> printed as panel "B": area-weighted representativeness by Ecosystem Type
ecotype_pct <- representativeness_by_ecotype %>%
  mutate(
    pct_no_analog   = 100 * area_no_ecotype_analog_mha / total_area_mha,
    pct_outside_aoa = 100 * area_outside_aoa_mha       / total_area_mha,
    pct_represented = 100 * area_well_represented_mha  / total_area_mha
  )

pD <- ecotype_pct %>%
  select(EcoType, `No NEON EcoType analog` = pct_no_analog,
         `Outside AOA / novel` = pct_outside_aoa, `Within AOA` = pct_represented) %>%
  pivot_longer(-EcoType, names_to = "representativeness", values_to = "pct_area") %>%
  mutate(representativeness = factor(representativeness, levels = names(aoa_colors))) %>%
  ggplot(aes(x = EcoType, y = pct_area, fill = representativeness)) +
  geom_col(width = 0.6) +
  # % "Within AOA" (the NEON-coverage number) printed just above each bar,
  # since that segment is too thin to hold a legible in-bar label.
  geom_text(data = ecotype_pct, aes(x = EcoType, y = 100, label = sprintf("%.1f%%", pct_represented)),
            inherit.aes = FALSE, vjust = -0.4, size = 3.4, fontface = "bold") +
  scale_y_continuous(limits = c(0, 108), breaks = c(0, 25, 50, 75, 100)) +
  scale_fill_manual(values = aoa_colors, name = NULL) +
  labs(title = "B. ", x = "Ecosystem Type", y = "% of Ecosystem Type area") +
  fig_theme + theme(legend.position = "none")

# Physical layout: A next to pD (printed "B"), then pB (printed "C"), then
# pC (printed "D") — reading order top-left -> top-right -> down -> down
# already matches A/B/C/D once the labs()/draw_label() titles above are read,
# even though the object names below don't.
fig_top    <- plot_grid(pA, pD, ncol = 2, rel_widths = c(1.1, 1))
fig_repres <- plot_grid(fig_top, pB, pC, ncol = 1, rel_heights = c(0.8, 1.1, 1.2))

ggsave(file.path(figure_dir, "FigS_NEON_Representativeness.png"),
  fig_repres, width = 9, height = 12.5, units = "in", dpi = 300, bg = "white")

message("NEON representativeness supplement complete. Outputs in: ", rf_dir)
