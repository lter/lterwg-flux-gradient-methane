# Monthly condition-based spatial CH4 upscaling — Random Forest, three approaches.
#
# Spatial projection ONLY: loads the source-probability and magnitude models
# already fit and tested by 12_SourceProp_MagnitudeModels.R (via
# source_magnitude_model_bundle.rds) and applies them to the global upland
# grid, month by month, 2000-2025. Does not fit or evaluate any model
# itself — see 12_SourceProp_MagnitudeModels.R for training, OOB skill, and
# the FLUXNET-CH4 external validation.
#
# Three flux expression approaches, sharing the same Stage 1/Stage 2 models:
#
#   Approach 1 — Continuous:
#     flux = P(source) × source_mag + (1 − P(source)) × sink_mag
#     No hard threshold; probability-weighted blend.
#
#   Approach 2 — Dichotomous:
#     Hard threshold at P(source) ≥ 0.5; assigns cell entirely to source or sink flux.
#
#   Approach 3 — All-Sink (theoretical maximum sink):
#     Every upland cell is assigned the sink magnitude regardless of
#     P(source). This represents the most extreme possible terrestrial sink
#     estimate derivable from ecosystem-scale observations.
#     A GMB threshold sensitivity sweep (0.30–0.99) is retained as a diagnostic
#     showing that even forcing all cells to sinks never recovers the GMB target.
#
# "Arid" is classified here the SAME way 12_SourceProp_MagnitudeModels.R
# classified it for training — aridity_index < arid_ai_threshold (loaded
# from the model bundle, not redefined here) — not the biome/MODIS raster's
# own "desert" class. P(source) IS hard-forced to 0 for Arid cells: the
# only 2 NEON sites crossing this threshold (JORN, SRER) are empirically
# 100% weak-source, an unrepresentative sample for learning a genuine
# P(source) pattern, so 12_SourceProp_MagnitudeModels.R excludes Arid from
# Stage 1 and Stage 2's source model entirely and this script asserts that
# design (model_bundle$arid_forced_sink) rather than assuming it. Arid
# cells are instead routed deterministically through the sink magnitude
# model, which extrapolates from Forest/Grassland/Shrubland sink behavior
# via the shared continuous covariates.
#
# Outputs: /Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF/OUTPUT/

library(tidyverse)
library(terra)
library(ranger)

# ── Paths ─────────────────────────────────────────────────────────────────────

spatial_dir <- Sys.getenv("MONTHLY_UPSCALING_DIR",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly")
rf_dir <- Sys.getenv("MONTHLY_RF_DIR",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF")

if (!dir.exists(spatial_dir)) stop("GLM directory not found (needed for comparison): ", spatial_dir)

output_dir <- file.path(rf_dir, "OUTPUT")
figure_dir <- file.path(rf_dir, "FIGURES")
data_dir   <- file.path(rf_dir, "DATA")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir,   showWarnings = FALSE, recursive = TRUE)

ecoregions_zip     <- file.path(spatial_dir,  "Ecoregions2017.zip")
era5_land_dir      <- file.path(spatial_dir,  "DATA/era5_land_monthly")
modis_ecotype_dir  <- file.path(spatial_dir,  "DATA/modis_mcd12c1_processed")
wad2m_dir          <- file.path(spatial_dir,  "DATA/wad2m")

ecoregions_extracted_dir <- file.path(data_dir, "ecoregions2017")
ecoregions_already_extracted <- length(list.files(ecoregions_extracted_dir,
  pattern = "Ecoregions2017\\.shp$")) > 0

model_bundle_file <- file.path(output_dir, "source_magnitude_model_bundle.rds")

required_files <- c(model_bundle_file)
if (!ecoregions_already_extracted) required_files <- c(required_files, ecoregions_zip)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0)
  stop("Missing: ", paste(missing_files, collapse = ", "),
       "\n(", model_bundle_file, " is written by 12_SourceProp_MagnitudeModels.R — run that first.)")

# MODIS land-cover files are required, not optional. The Ecoregions2017
# biome-based fallback (built below regardless, for defensiveness) has no
# concept of cropland, urban, or wetland at all -- it classifies purely by
# ecoregion/biome, so a city or farm field sitting inside a temperate
# forest ecoregion would be silently counted as Forest upland. The MODIS
# path is what actually excludes those land-cover classes (see script
# 11_Global_DownloadMODIS_WAD2M.R's classify_modis_to_ecotype(): IGBP Water,
# Permanent Wetlands, Urban and Built-up, and Snow/Ice all map to NA;
# Croplands map to a code ecotype_lookup below filters out). Silently
# falling back to the biome-only method would silently reintroduce
# cropland/urban/wetland cells into the upland grid, so this stops instead.
modis_files_check <- list.files(modis_ecotype_dir,
  pattern = "^MODIS_MCD12C1_ecotype_[0-9]{4}\\.tif$", full.names = TRUE)
if (length(modis_files_check) == 0)
  stop("No MODIS land-cover files found in ", modis_ecotype_dir,
       ". These are required to exclude cropland/urban/wetland cells from ",
       "the upland grid (the Ecoregions-biome fallback cannot do this). ",
       "Run 11_Global_DownloadMODIS_WAD2M.R first.")

# ── Load fitted models (12_SourceProp_MagnitudeModels.R) ──────────────────────
# Loading scalar thresholds from the bundle too (not redefining them here) is
# deliberate: arid_ai_threshold/aridity_mat_floor/binary_threshold drifting
# apart between the fitting script and this one is exactly how the
# Arid-EcoType bug happened before.

message("Loading fitted models from ", model_bundle_file, " ...")
model_bundle <- readRDS(model_bundle_file)

rf_model_A              <- model_bundle$rf_model_A
iso_cal_A                <- model_bundle$iso_cal_A
sink_mag_model            <- model_bundle$sink_mag_model
source_mag_model          <- model_bundle$source_mag_model
magnitude_standardizers    <- model_bundle$magnitude_standardizers
training_vswc_range        <- model_bundle$training_vswc_range
training_prec_range        <- model_bundle$training_prec_range
ecotype_levels              <- model_bundle$ecotype_levels
arid_ai_threshold           <- model_bundle$arid_ai_threshold
aridity_mat_floor           <- model_bundle$aridity_mat_floor
binary_threshold            <- model_bundle$binary_threshold
arid_forced_sink            <- model_bundle$arid_forced_sink

# Asserted rather than silently assumed: if 12_SourceProp_MagnitudeModels.R
# ever changes how it handles Arid, this script's hard-coded P(source) = 0
# override below needs to change with it, not silently drift out of sync
# (the same class of bug arid_ai_threshold/aridity_mat_floor had before).
if (!isTRUE(arid_forced_sink))
  stop("Model bundle does not have arid_forced_sink = TRUE -- this script's ",
       "Arid handling (P(source) hard-forced to 0, sink model only) assumes ",
       "it does. Re-run 12_SourceProp_MagnitudeModels.R or update this script ",
       "to match its current Arid design.")

# ── Scalar parameters (spatial-projection only) ────────────────────────────────

gmb_soil_sink_tg_ch4_yr    <- -35           # GMB upland soil sink target
arid_shrubland_fill_temp_threshold <- 15
years_to_process            <- 2000:2025
# Threshold grid for GMB-constrained search (fine steps to allow interpolation)
gmb_threshold_grid          <- seq(0.30, 0.99, by = 0.01)

# ── Helper functions ──────────────────────────────────────────────────────────

gC_m2_yr_to_tg_ch4 <- function(flux, area_mha) flux * area_mha * 0.0133333333333333

apply_standardizers <- function(dat, std) {
  dat %>% mutate(
    z_Tair = (mean_ERA5_Tair_C - std$Tair$center) / std$Tair$scale,
    z_VSWC = (mean_ERA5_VSWC   - std$VSWC$center) / std$VSWC$scale,
    z_MAP  = (MAP               - std$MAP$center)  / std$MAP$scale,
    z_MAT  = (MAT               - std$MAT$center)  / std$MAT$scale
  )
}

unzip_if_needed <- function(zip_file, exdir, pattern) {
  if (!dir.exists(exdir)) dir.create(exdir, recursive = TRUE)
  if (length(list.files(exdir, pattern = pattern)) == 0)
    utils::unzip(zip_file, exdir = exdir)
  invisible(list.files(exdir, pattern = pattern, full.names = TRUE))
}

# Duplicated from 12_SourceProp_MagnitudeModels.R (kept in sync manually,
# same pattern used between 13/19) — this script only APPLIES the already-
# fitted models loaded above, it never fits anything itself.
predict_rf_prob <- function(rf_model, newdata) {
  as.numeric(predict(rf_model, data = newdata)$predictions[, "1"])
}

calibrate_isotonic <- function(iso_cal, new_probs) {
  pmin(pmax(approx(iso_cal$x, iso_cal$y, xout = new_probs,
                   method = "linear", rule = 2, ties = "ordered")$y, 0), 1)
}

predict_rf_magnitude <- function(fit, newdata) {
  mag <- exp(predict(fit$model, data = newdata)$predictions)
  if (identical(fit$state, "Weak-sink")) -mag else mag
}

# ── Spatial data infrastructure ────────────────────────────────────────────────

ecoregion_files <- unzip_if_needed(ecoregions_zip,
  file.path(data_dir, "ecoregions2017"), "Ecoregions2017\\.shp$")

era5_land_files <- file.path(era5_land_dir,
  sprintf("era5_land_monthly_%s.nc", years_to_process))
use_era5_land <- all(file.exists(era5_land_files))

read_era5_land_year <- function(year) {
  file    <- file.path(era5_land_dir, sprintf("era5_land_monthly_%s.nc", year))
  climate <- rast(file); lnames <- names(climate)
  t2m_idx   <- grep("t2m",   lnames, ignore.case = TRUE)
  swvl1_idx <- grep("swvl1", lnames, ignore.case = TRUE)
  swvl2_idx <- grep("swvl2", lnames, ignore.case = TRUE)
  tp_idx    <- grep("(?i)^tp[_=]|^tp$", lnames, perl = TRUE)
  counts    <- lengths(list(t2m_idx, swvl1_idx, swvl2_idx, tp_idx))
  if (!all(counts == 12L)) stop(sprintf(
    "Expected 12 layers each for t2m/swvl1/swvl2/tp in %s; found %s.",
    basename(file), paste(c("t2m","swvl1","swvl2","tp"), counts, sep="=", collapse=", ")))
  tavg <- climate[[t2m_idx]] - 273.15
  vswc <- (climate[[swvl1_idx]] + climate[[swvl2_idx]]) / 2
  prec <- rast(lapply(seq_along(tp_idx), function(i)
    climate[[tp_idx[i]]] * 1000 * lubridate::days_in_month(lubridate::make_date(year, i, 1))))
  names(tavg) <- sprintf("tavg_%02d", 1:12)
  names(vswc) <- sprintf("vswc_%02d", 1:12)
  names(prec) <- sprintf("prec_%02d", 1:12)
  list(tavg = tavg, vswc = vswc, prec = prec)
}

if (use_era5_land) {
  message("Using ERA5-Land monthly grids from ", era5_land_dir)
  first_era5 <- read_era5_land_year(years_to_process[1])
  tavg <- first_era5$tavg; prec <- first_era5$prec; template <- tavg[[1]]
  mat_sum <- map_sum <- template
  values(mat_sum) <- 0; values(map_sum) <- 0; n_mat_layers <- 0
  for (yr in years_to_process) {
    cy <- read_era5_land_year(yr)
    mat_sum <- mat_sum + sum(cy$tavg, na.rm = TRUE)
    map_sum <- map_sum + sum(cy$prec, na.rm = TRUE)
    n_mat_layers <- n_mat_layers + 12
  }
  mat <- mat_sum / n_mat_layers; names(mat) <- "MAT"
  map <- map_sum / length(years_to_process); names(map) <- "MAP"
} else {
  wc_t <- file.path(spatial_dir, "wc2.1_10m_tavg.zip")
  wc_p <- file.path(spatial_dir, "wc2.1_10m_prec.zip")
  if (!file.exists(wc_t) || !file.exists(wc_p))
    stop("ERA5-Land absent and WorldClim fallback zips not found.")
  message("ERA5-Land not found; using WorldClim fallback.")
  tavg <- rast(sort(unzip_if_needed(wc_t, file.path(data_dir,"worldclim_tavg"), "tavg_.*\\.tif$")))
  prec <- rast(sort(unzip_if_needed(wc_p, file.path(data_dir,"worldclim_prec"), "prec_.*\\.tif$")))
  names(tavg) <- sprintf("tavg_%02d", 1:12); names(prec) <- sprintf("prec_%02d", 1:12)
  template <- tavg[[1]]
  mat <- mean(tavg, na.rm = TRUE); names(mat) <- "MAT"
  map <- sum(prec, na.rm = TRUE);  names(map) <- "MAP"
}

# Cache the MAT/MAP climatology so downstream scripts (e.g. the NEON
# representativeness supplement) don't have to re-read all 26 years of
# ERA5-Land/WorldClim files just to reproduce these two static rasters.
writeRaster(c(mat, map), file.path(output_dir, "mat_map_climatology.tif"), overwrite = TRUE)

cell_area_mha <- cellSize(template, unit = "m") / 1e10
names(cell_area_mha) <- "area_mha"

ecoregions    <- vect(ecoregion_files[1])
biome_field   <- intersect(c("BIOME_NUM","BIOME"), names(ecoregions))[1]
if (is.na(biome_field)) stop("Cannot find BIOME_NUM/BIOME in ecoregions shapefile.")
biome_raster  <- rasterize(ecoregions, template, field = biome_field, touches = TRUE)
names(biome_raster) <- "biome_num"
ecotype_raster <- subst(biome_raster, from = 1:14,
  to = c(2,2,2,2,2,2, 3,3,NA,3, 4,4,4,NA))
names(ecotype_raster) <- "ecotype_code"

ecotype_lookup <- tibble(
  ecotype_code = c(1,2,3,4,5),
  EcoType      = c("Cropland","Forest","Grassland","Shrubland","Arid")
) %>% filter(EcoType != "Cropland")

modis_files <- list.files(modis_ecotype_dir,
  pattern = "^MODIS_MCD12C1_ecotype_[0-9]{4}\\.tif$", full.names = TRUE)
use_modis   <- length(modis_files) > 0
modis_years <- if (use_modis)
  as.integer(str_match(basename(modis_files), "([0-9]{4})\\.tif$")[,2]) else integer()

get_ecotype_raster <- function(year) {
  if (!use_modis) return(ecotype_raster)
  sel <- modis_years[which.min(abs(modis_years - year))]
  r   <- rast(file.path(modis_ecotype_dir, sprintf("MODIS_MCD12C1_ecotype_%s.tif", sel)))
  if (!same.crs(r, template)) crs(r) <- crs(template)
  if (!compareGeom(r, template, stopOnError = FALSE)) r <- resample(r, template, method = "near")
  names(r) <- "ecotype_code"; r
}

wad2m_files <- list.files(wad2m_dir, pattern = "\\.nc$", full.names = TRUE)
use_wad2m   <- length(wad2m_files) > 0
wad2m       <- if (use_wad2m) rast(wad2m_files[1]) else NULL
wad2m_yrs   <- if (use_wad2m) 2000:2020 else integer()

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

# ── Spatial projection loop ───────────────────────────────────────────────────
# Both Stage 1/Stage 2 models predict in the same pass.
# Stores per cell: source_prob_A, sink_flux, source_flux, area.

message("Starting spatial loop (", length(years_to_process), " years × 12 months)...")
monthly_cell_predictions <- vector("list", length(years_to_process) * 12)
idx <- 1L

for (year in years_to_process) {
  climate_year <- if (use_era5_land) read_era5_land_year(year) else
    list(tavg = tavg, vswc = NULL, prec = prec)

  for (m in 1:12) {
    message("  Year ", year, " month ", m)
    eco_r  <- get_ecotype_raster(year)
    inund  <- get_inundation_fraction(year, m)
    inund[is.na(inund)] <- 0

    stk <- c(eco_r, inund, climate_year$tavg[[m]], climate_year$prec[[m]], mat, map, cell_area_mha)
    names(stk) <- c("ecotype_code","inundation_fraction","tavg","prec","MAT","MAP","area_mha")

    if (use_era5_land) {
      vl <- climate_year$vswc[[m]]
      vl[is.na(vl) & !is.na(eco_r) & eco_r == 5 &
         !is.na(climate_year$tavg[[m]]) & climate_year$tavg[[m]] > arid_shrubland_fill_temp_threshold] <- 0
      stk <- c(stk, vl); names(stk)[nlyr(stk)] <- "vswc"
    }

    dat <- as.data.frame(stk, xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
      filter(!is.na(ecotype_code), !is.na(tavg), !is.na(prec), !is.na(MAT), !is.na(MAP), !is.na(area_mha)) %>%
      mutate(ecotype_code = as.integer(ecotype_code)) %>%
      inner_join(ecotype_lookup, by = "ecotype_code") %>%
      mutate(
        Year                 = year,
        month                = m,
        inundation_fraction  = pmin(pmax(inundation_fraction, 0), 1),
        upland_area_fraction = 1 - inundation_fraction,
        area_mha             = area_mha * upland_area_fraction,
        mean_ERA5_Tair_C = tavg,
        mean_ERA5_VSWC   = if (use_era5_land) vswc else
          scales::rescale(pmin(pmax(prec, training_prec_range[1]), training_prec_range[2]),
                          to = training_vswc_range, from = training_prec_range),
        # See aridity_mat_floor (loaded from model_bundle above).
        aridity_index = if_else(MAT > aridity_mat_floor, MAP / (MAT + 10), NA_real_),
        # Classified the SAME way 12_SourceProp_MagnitudeModels.R classified
        # training data: aridity_index < arid_ai_threshold, not the
        # biome/MODIS raster's own "Arid" class from ecotype_lookup/eco_r
        # (checked here via the pre-remap EcoType, right before it's
        # overwritten below). A biome-classified desert cell that does NOT
        # cross the aridity_index threshold folds into Shrubland (the
        # raster's next-best guess); any cell that DOES cross the threshold
        # becomes "Arid" and is handled by the models' own trained Arid
        # behavior (see below), not a hard override.
        is_arid = as.integer(!is.na(aridity_index) & aridity_index < arid_ai_threshold),
        EcoType = factor(if_else(is_arid == 1, "Arid",
                                 if_else(EcoType == "Arid", "Shrubland", as.character(EcoType))),
                         levels = ecotype_levels)
      ) %>%
      filter(area_mha > 0, is.finite(aridity_index))

    dat <- apply_standardizers(dat, magnitude_standardizers)

    # ── Stage 1: predict Model A ──────────────────────────────────────────────
    # Arid cells ARE hard-forced to P(source) = 0. rf_model_A never trained
    # on Arid data at all (12_SourceProp_MagnitudeModels.R excludes it from
    # Stage 1 entirely), because the only 2 NEON sites crossing the
    # aridity_index threshold (JORN, SRER) are empirically 100% weak-source
    # -- the opposite of the standard "arid soils are net CH4 sinks" prior --
    # so there's no trustworthy basis for a learned Arid P(source). Rather
    # than let rf_model_A extrapolate onto a level it never saw, Arid is
    # routed deterministically to the sink magnitude model instead, which
    # itself excludes Arid's own site-months and predicts Arid cells by
    # extrapolating from Forest/Grassland/Shrubland sink behavior via the
    # shared continuous covariates. source_mag_model is still evaluated
    # below (dat$predicted_source_flux) for every cell for code simplicity,
    # but it's multiplied by source_prob_A = 0 for Arid cells in the
    # Continuous approach below, and Approach 2's P >= binary_threshold
    # dichotomous rule never fires for Arid either, so Arid never actually
    # receives a source-model flux value in any approach.
    dat$source_prob_A_raw <- predict_rf_prob(rf_model_A, dat)
    dat$source_prob_A     <- if_else(dat$is_arid == 1, 0,
                                     calibrate_isotonic(iso_cal_A, dat$source_prob_A_raw))

    # Stage 2 uses Model A probability as predictor
    dat <- dat %>% mutate(source_probability = source_prob_A)
    dat$predicted_sink_flux   <- pmin(predict_rf_magnitude(sink_mag_model,   dat), 0)
    dat$predicted_source_flux <- pmax(predict_rf_magnitude(source_mag_model, dat), 0)

    monthly_cell_predictions[[idx]] <- dat %>%
      select(Year, cell, x, y, month, EcoType, is_arid, area_mha,
             inundation_fraction, upland_area_fraction,
             source_prob_A,
             predicted_sink_flux, predicted_source_flux)
    idx <- idx + 1L
  }
}

cell_preds <- bind_rows(monthly_cell_predictions)

# ── Approach 1: Continuous ────────────────────────────────────────────────────
# flux = P_A × source_mag + (1 − P_A) × sink_mag (no hard threshold)

cell_preds <- cell_preds %>%
  mutate(
    flux_continuous  = source_prob_A * predicted_source_flux +
                       (1 - source_prob_A) * predicted_sink_flux,
    tg_continuous    = gC_m2_yr_to_tg_ch4(flux_continuous, area_mha),

    # Approach 2: Dichotomous — same Model A probabilities, hard threshold at 0.5
    flux_balanced    = if_else(source_prob_A >= binary_threshold,
                               predicted_source_flux, predicted_sink_flux),
    tg_balanced      = gC_m2_yr_to_tg_ch4(flux_balanced, area_mha)
  )

# ── Approach 3: All-Sink (theoretical maximum sink) ───────────────────────────
# Every cell is assigned the sink magnitude flux, regardless of P(source).
# This is the most extreme possible sink estimate from the model.

cell_preds <- cell_preds %>%
  mutate(
    flux_constrained = predicted_sink_flux,
    tg_constrained   = gC_m2_yr_to_tg_ch4(flux_constrained, area_mha)
  )

# ── GMB threshold sensitivity (diagnostic, decoupled from Approach 3) ─────────
# Sweep P(source) thresholds 0.30–0.99 to show the budget never reaches −35.

message("Running GMB threshold sensitivity sweep (diagnostic)...")
gmb_sensitivity <- map_dfr(gmb_threshold_grid, function(p) {
  cell_preds %>%
    mutate(tg = gC_m2_yr_to_tg_ch4(
      if_else(source_prob_A >= p, predicted_source_flux, predicted_sink_flux),
      area_mha)) %>%
    group_by(Year) %>%
    summarise(annual_tg = sum(tg, na.rm = TRUE), .groups = "drop") %>%
    summarise(mean_tg = mean(annual_tg), .groups = "drop") %>%
    mutate(threshold = p)
})

# Record the minimum achievable budget from the sweep (floor at P* = 0.99)
P_star_nearest <- gmb_sensitivity$threshold[which.min(abs(gmb_sensitivity$mean_tg - gmb_soil_sink_tg_ch4_yr))]
budget_floor   <- min(gmb_sensitivity$mean_tg, na.rm = TRUE)
message(sprintf("GMB sensitivity sweep: budget floor = %.1f Tg/yr at P* = %.2f (target −35 never reached)",
                budget_floor, P_star_nearest))

# ── Annual budget summaries ───────────────────────────────────────────────────

annual_budget <- cell_preds %>%
  group_by(Year) %>%
  summarise(
    continuous_tg_ch4_yr   = sum(tg_continuous,  na.rm = TRUE),
    balanced_tg_ch4_yr     = sum(tg_balanced,     na.rm = TRUE),
    constrained_tg_ch4_yr  = sum(tg_constrained,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(gmb_reference_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr)

annual_budget_long <- annual_budget %>%
  pivot_longer(c(continuous_tg_ch4_yr, balanced_tg_ch4_yr, constrained_tg_ch4_yr),
               names_to = "approach", values_to = "annual_tg_ch4_yr") %>%
  mutate(approach = recode(approach,
    continuous_tg_ch4_yr  = "Continuous",
    balanced_tg_ch4_yr    = "Dichotomous",
    constrained_tg_ch4_yr = "All-Sink"
  ))

budget_summary <- annual_budget_long %>%
  group_by(approach) %>%
  summarise(
    mean_tg_ch4_yr  = mean(annual_tg_ch4_yr, na.rm = TRUE),
    sd_tg_ch4_yr    = sd(annual_tg_ch4_yr,   na.rm = TRUE),
    min_tg_ch4_yr   = min(annual_tg_ch4_yr,  na.rm = TRUE),
    max_tg_ch4_yr   = max(annual_tg_ch4_yr,  na.rm = TRUE),
    pct_of_gmb      = 100 * mean_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr),
    .groups = "drop"
  )

comparison_budget <- budget_summary

# ── Write outputs ─────────────────────────────────────────────────────────────

write.csv(annual_budget,          file.path(output_dir, "annual_budget_three_approaches.csv"),    row.names = FALSE)
write.csv(annual_budget_long,     file.path(output_dir, "annual_budget_long.csv"),                row.names = FALSE)
write.csv(budget_summary,         file.path(output_dir, "budget_summary_three_approaches.csv"),   row.names = FALSE)
write.csv(comparison_budget,      file.path(output_dir, "comparison_budget_all_approaches.csv"),  row.names = FALSE)
write.csv(gmb_sensitivity,        file.path(output_dir, "gmb_threshold_sensitivity.csv"),         row.names = FALSE)
write.csv(data.frame(
  P_star                   = P_star_nearest,   # nearest threshold to GMB target; budget floor never reached
  budget_floor_tg_ch4_yr  = budget_floor,
  binary_threshold         = binary_threshold,
  gmb_target_tg_ch4_yr    = gmb_soil_sink_tg_ch4_yr,
  model_bundle_source       = model_bundle_file
), file.path(output_dir, "spatial_projection_parameters.csv"), row.names = FALSE)

saveRDS(cell_preds, file.path(output_dir, "monthly_cell_predictions_2000_2025.rds"))
terra::writeRaster(template, file.path(output_dir, "era5_template.tif"), overwrite = TRUE)

# ── Summary text ──────────────────────────────────────────────────────────────

capture.output({
  cat("RF spatial upscaling — three-approach comparison\n\n")
  cat("Models loaded from: ", model_bundle_file, "\n\n", sep = "")
  cat("Approach 1 — Continuous: P(source)-weighted flux blend (no hard threshold)\n")
  cat(sprintf("Approach 2 — Dichotomous: hard threshold at P(source) >= %.2f\n", binary_threshold))
  cat(sprintf("Approach 3 — All-Sink: all cells assigned sink flux (theoretical maximum sink)\n"))
  cat(sprintf("  GMB sensitivity floor: %.1f Tg/yr at P* = %.2f; GMB target (%.0f) never reached\n\n",
    budget_floor, P_star_nearest, gmb_soil_sink_tg_ch4_yr))

  cat("─── Annual Budget Summary — Three Approaches ───\n")
  print(budget_summary)

  cat("\n─── Full Budget Comparison ───\n")
  print(comparison_budget)
}, file = file.path(output_dir, "spatial_trial_summary_RF.txt"))

# ── Figures ───────────────────────────────────────────────────────────────────

fig_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = "bottom", strip.background = element_rect(fill = "grey92"),
        axis.title = element_text(size = 11), axis.text = element_text(size = 10))

approach_colors <- c(
  "Continuous"     = "#009688",  # material teal
  "Dichotomous"    = "#9C27B0",  # material purple
  "GMB-Dichotomous" = "#E91E63"   # material pink
)
gmb_line_color <- "#2166AC"

# Fig 1: Annual budget time series — three RF approaches
p1 <- annual_budget_long %>%
  ggplot(aes(x = Year, y = annual_tg_ch4_yr, color = approach)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -34,
           fill = gmb_line_color, alpha = 0.10) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = gmb_line_color,
             linetype = "dashed", linewidth = 0.7) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_color_manual(values = approach_colors) +
  labs(title = "Annual Net Exchange — Three Approaches",
       x = "Year", y = "Net exchange (Tg CH₄ yr⁻¹)", color = NULL) +
  fig_theme

ggsave(file.path(figure_dir, "Fig1_budget_three_approaches.png"),
  p1, width = 9, height = 4.5, units = "in", dpi = 300, bg = "white")

# Fig 2: Budget bar chart — mean ± range
p2 <- comparison_budget %>%
  mutate(approach = str_wrap(approach, 30)) %>%
  ggplot(aes(x = approach, y = mean_tg_ch4_yr, fill = approach)) +
  geom_col(width = 0.65, alpha = 0.85, show.legend = FALSE) +
  geom_errorbar(aes(ymin = min_tg_ch4_yr, ymax = max_tg_ch4_yr),
                width = 0.2, linewidth = 0.6) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = gmb_line_color,
             linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  annotate("text", x = Inf, y = gmb_soil_sink_tg_ch4_yr + 1.5,
           label = "GMB reference (−35)", hjust = 1.05, size = 3.2, color = gmb_line_color) +
  scale_fill_manual(values = approach_colors) +
  labs(title = "Mean Annual Budget (2000–2025) — All Approaches",
       x = NULL, y = "Net exchange (Tg CH₄ yr⁻¹)") +
  coord_flip() +
  fig_theme

ggsave(file.path(figure_dir, "Fig2_budget_bar.png"),
  p2, width = 8, height = 4, units = "in", dpi = 300, bg = "white")

# Fig 4: GMB threshold sensitivity curve with P* marked
p4 <- gmb_sensitivity %>%
  ggplot(aes(x = threshold, y = mean_tg)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -34,
           fill = gmb_line_color, alpha = 0.10) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = gmb_line_color,
             linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_line(color = "#4A148C", linewidth = 1.1) +
  geom_point(color = "#4A148C", size = 1.5) +
  list(
    geom_vline(xintercept = P_star_nearest, color = "#4A148C", linetype = "dotted", linewidth = 0.9),
    annotate("label", x = P_star_nearest, y = Inf,
             label = sprintf("Budget floor: %.1f Tg\n(P = %.2f)", budget_floor, P_star_nearest),
             hjust = -0.1, vjust = 1.4, size = 3.2, color = "#4A148C",
             fill = "white", label.size = 0)
  ) +
  labs(title = "Approach 3: GMB-Constrained Threshold Sensitivity (P(source))",
       x = "P(source) threshold", y = "Mean annual budget (Tg CH₄ yr⁻¹)") +
  fig_theme

ggsave(file.path(figure_dir, "Fig4_gmb_threshold_sensitivity.png"),
  p4, width = 7, height = 4, units = "in", dpi = 300, bg = "white")

message("Three-approach RF spatial upscaling complete. Outputs in: ", rf_dir)
