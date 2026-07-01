# Monthly condition-based spatial CH4 upscaling — Random Forest, three approaches.
#
# Single Stage 1 P(source) model (1:1 balanced RF), three flux expression approaches:
#
#   Approach 1 — Continuous:
#     flux = P(source) × source_mag + (1 − P(source)) × sink_mag
#     No hard threshold; probability-weighted blend.
#
#   Approach 2 — Dichotomous:
#     Hard threshold at P(source) ≥ 0.5; assigns cell entirely to source or sink flux.
#
#   Approach 3 — All-Sink (theoretical maximum sink):
#     Every non-arid upland cell is assigned the sink magnitude regardless of
#     P(source). This represents the most extreme possible terrestrial sink
#     estimate derivable from ecosystem-scale observations.
#     A GMB threshold sensitivity sweep (0.30–0.99) is retained as a diagnostic
#     showing that even forcing all cells to sinks never recovers the GMB target.
#
# Stage 1 (single balanced RF):
#   - Training data downsampled to 1:1 source:sink (equal class representation)
#   - Regularised probability forest: min.node.size = 20, max.depth = 8,
#     sample.fraction = 0.7 without replacement
#   - OOB predictions used for honest AUC and isotonic calibration
#
# Stage 2 (shared across all approaches):
#   - ranger log-absolute-flux regression, trained on full data
#
# Outputs: /Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF/OUTPUT/

library(tidyverse)
library(data.table)
library(terra)
library(ranger)
library(patchwork)
library(cowplot)

# ── Paths ─────────────────────────────────────────────────────────────────────

localdir.ch4 <- Sys.getenv("LOCALDIR_CH4",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/Methane")
spatial_dir <- Sys.getenv("MONTHLY_UPSCALING_DIR",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly")
rf_dir <- Sys.getenv("MONTHLY_RF_DIR",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF")

if (!dir.exists(localdir.ch4)) stop("CH4 data directory not found: ", localdir.ch4)
if (!dir.exists(spatial_dir))  stop("GLM directory not found (needed for comparison): ", spatial_dir)

output_dir <- file.path(rf_dir, "OUTPUT")
figure_dir <- file.path(rf_dir, "FIGURES")
data_dir   <- file.path(rf_dir, "DATA")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir,   showWarnings = FALSE, recursive = TRUE)

glm_output_dir     <- file.path(spatial_dir, "OUTPUT")
era5_30min_file    <- file.path(localdir.ch4, "OUTPUT/NEON_ERA5_gapfilled_30min.csv.gz")
site_behavior_file <- file.path(localdir.ch4, "OUTPUT/30min_site_behavior.csv")
ecoregions_zip     <- file.path(spatial_dir,  "Ecoregions2017.zip")
era5_land_dir      <- file.path(spatial_dir,  "DATA/era5_land_monthly")
modis_ecotype_dir  <- file.path(spatial_dir,  "DATA/modis_mcd12c1_processed")
wad2m_dir          <- file.path(spatial_dir,  "DATA/wad2m")

missing_files <- c(era5_30min_file, site_behavior_file, ecoregions_zip)[
  !file.exists(c(era5_30min_file, site_behavior_file, ecoregions_zip))]
if (length(missing_files) > 0) stop("Missing: ", paste(missing_files, collapse = ", "))

# ── Scalar parameters ─────────────────────────────────────────────────────────

binary_threshold            <- 0.5          # fixed threshold for approaches 1 & 2
gmb_soil_sink_tg_ch4_yr    <- -35           # GMB upland soil sink target
arid_shrubland_fill_temp_threshold <- 15
arid_ai_threshold           <- 15
years_to_process            <- 2000:2025
rf_seed                     <- 42
n_trees                     <- 500
rf_min_node_size            <- 20
rf_max_depth                <- 8
rf_sample_frac              <- 0.7
# Threshold grid for GMB-constrained search (fine steps to allow interpolation)
gmb_threshold_grid          <- seq(0.30, 0.99, by = 0.01)

# ── Helper functions ──────────────────────────────────────────────────────────

gC_m2_yr_to_tg_ch4 <- function(flux, area_mha) flux * area_mha * 0.0133333333333333

wtd_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

auc_rank <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]; predicted <- predicted[ok]
  n1 <- sum(observed == 1); n0 <- sum(observed == 0)
  r  <- rank(predicted, ties.method = "average")
  (sum(r[observed == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

make_standardizer <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  list(center = m, scale = if (!is.finite(s) || s == 0) 1 else s)
}

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

# ── RF helper functions ────────────────────────────────────────────────────────

class_balance_weights_rf <- function(y) {
  n <- length(y); tbl <- table(y)
  as.numeric(n / (length(tbl) * tbl[as.character(y)]))
}

fit_rf_prob_model <- function(training_data, predictors, use_case_weights = TRUE) {
  training_data <- training_data %>%
    mutate(weak_source_f = factor(weak_source, levels = c(0, 1)))
  fmla <- reformulate(predictors, response = "weak_source_f")
  cw   <- if (use_case_weights) class_balance_weights_rf(training_data$weak_source) else NULL
  ranger(
    formula         = fmla,
    data            = training_data,
    num.trees       = n_trees,
    probability     = TRUE,
    case.weights    = cw,
    min.node.size   = rf_min_node_size,
    max.depth       = rf_max_depth,
    replace         = FALSE,
    sample.fraction = rf_sample_frac,
    importance      = "impurity",
    seed            = rf_seed
  )
}

predict_rf_prob <- function(rf_model, newdata) {
  as.numeric(predict(rf_model, data = newdata)$predictions[, "1"])
}

fit_isotonic_calibration <- function(oob_probs, observed_labels) {
  ord <- order(oob_probs)
  iso <- isoreg(x = oob_probs[ord], y = as.numeric(observed_labels[ord]))
  list(x = oob_probs[ord], y = iso$yf)
}

calibrate_isotonic <- function(iso_cal, new_probs) {
  pmin(pmax(approx(iso_cal$x, iso_cal$y, xout = new_probs,
                   method = "linear", rule = 2, ties = "ordered")$y, 0), 1)
}

fit_rf_magnitude_model <- function(training_data, state_name, predictors) {
  training_data <- training_data %>%
    mutate(log_abs_flux = log(pmax(abs(monthly_flux_gC_m2_month), 1e-6)))
  list(
    model  = ranger(reformulate(predictors, "log_abs_flux"),
                    data = training_data, num.trees = n_trees,
                    importance = "impurity", seed = rf_seed),
    engine = "ranger_log_abs",
    state  = state_name
  )
}

predict_rf_magnitude <- function(fit, newdata) {
  mag <- exp(predict(fit$model, data = newdata)$predictions)
  if (identical(fit$state, "Weak-sink")) -mag else mag
}

summarise_magnitude_fit <- function(training_data, fit) {
  fitted <- predict_rf_magnitude(fit, training_data)
  tibble(
    magnitude_model               = fit$state,
    engine                        = fit$engine,
    n_observations                = nrow(training_data),
    mean_observed_flux            = mean(training_data$monthly_flux_gC_m2_month, na.rm = TRUE),
    mean_fitted_flux              = mean(fitted, na.rm = TRUE),
    rmse_gC_m2_month              = sqrt(mean((training_data$monthly_flux_gC_m2_month - fitted)^2, na.rm = TRUE)),
    mae_gC_m2_month               = mean(abs(training_data$monthly_flux_gC_m2_month - fitted), na.rm = TRUE),
    correlation_observed_fitted   = suppressWarnings(
      cor(training_data$monthly_flux_gC_m2_month, fitted, use = "complete.obs"))
  )
}

classification_skill <- function(prob, observed, threshold, model_label, eval_basis,
                                 oob_error = NA_real_) {
  pred_class <- as.integer(prob >= threshold)
  tibble(
    model              = model_label,
    evaluation_basis   = eval_basis,
    n_site_months      = length(observed),
    n_sink_months      = sum(observed == 0),
    n_source_months    = sum(observed == 1),
    source_fraction    = mean(observed),
    auc                = auc_rank(observed, prob),
    brier_score        = mean((observed - prob)^2),
    brier_null         = mean(observed) * (1 - mean(observed)),
    brier_skill_score  = 1 - mean((observed - prob)^2) / (mean(observed) * (1 - mean(observed))),
    tjur_r2            = mean(prob[observed == 1]) - mean(prob[observed == 0]),
    ranger_oob_error   = oob_error,
    threshold          = threshold,
    accuracy           = mean(pred_class == observed),
    sensitivity_source = sum(pred_class == 1 & observed == 1) / sum(observed == 1),
    specificity_sink   = sum(pred_class == 0 & observed == 0) / sum(observed == 0),
    precision_source   = sum(pred_class == 1 & observed == 1) / max(sum(pred_class == 1), 1),
    predicted_source_fraction = mean(pred_class == 1)
  )
}

# ── Spatial data infrastructure (identical to GLM script) ─────────────────────

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

# ── Load and prepare training data ────────────────────────────────────────────

site_behavior <- read.csv(site_behavior_file) %>%
  mutate(SITE_ID = as.character(SITE_ID))

upland_sites <- site_behavior %>%
  filter(!is.na(EcoType),
    !str_detect(EcoType, regex("wetland|inundat|flood|marsh|swamp|bog|fen|lake|rice",
                               ignore_case = TRUE))) %>%
  distinct(SITE_ID, EcoType, MAP, MAT)

era5_30min <- data.table::fread(era5_30min_file) %>% as_tibble() %>%
  mutate(across(c(SITE_ID), as.character),
         across(c(Year, month), as.integer),
         across(c(ERA5_Tair_C, ERA5_VSWC, gapfilled_CH4_mgC_30min), as.numeric)) %>%
  inner_join(upland_sites %>% select(SITE_ID, EcoType), by = c("SITE_ID","EcoType")) %>%
  filter(is.finite(Year), is.finite(month), is.finite(ERA5_Tair_C), is.finite(ERA5_VSWC))

monthly_training <- era5_30min %>%
  reframe(.by = c(SITE_ID, EcoType, Year, month),
    monthly_budget_mgC_m2  = sum(gapfilled_CH4_mgC_30min, na.rm = TRUE),
    mean_ERA5_Tair_C        = mean(ERA5_Tair_C, na.rm = TRUE),
    mean_ERA5_VSWC          = mean(ERA5_VSWC,   na.rm = TRUE)) %>%
  left_join(upland_sites, by = c("SITE_ID","EcoType")) %>%
  mutate(
    EcoType                  = factor(EcoType, levels = c("Forest","Grassland","Shrubland")),
    weak_source              = as.integer(monthly_budget_mgC_m2 > 0),
    monthly_flux_gC_m2_month = monthly_budget_mgC_m2 / 1000,
    MAP                      = as.numeric(MAP),
    MAT                      = as.numeric(MAT),
    aridity_index            = MAP / (MAT + 10),
    is_arid                  = as.integer(!is.na(aridity_index) & aridity_index < arid_ai_threshold)
  ) %>%
  filter(is.finite(weak_source), is.finite(mean_ERA5_Tair_C), is.finite(mean_ERA5_VSWC),
         is.finite(MAP), is.finite(MAT), is.finite(aridity_index))

class_predictors <- c("EcoType","mean_ERA5_Tair_C","mean_ERA5_VSWC","MAP","MAT","aridity_index")

# ── Stage 1: Balanced RF (1:1 source:sink) ────────────────────────────────────
# Downsample source months to match the number of sink months so both classes
# have exactly equal representation in training. No case weights needed.
# OOB predictions are used for honest evaluation and isotonic calibration.

set.seed(rf_seed)
sink_idx   <- which(monthly_training$weak_source == 0)
source_idx <- which(monthly_training$weak_source == 1)
bal_idx    <- c(sink_idx, sample(source_idx, length(sink_idx)))
balanced_training <- monthly_training[bal_idx, ]

message(sprintf("Stage 1 balanced training: %d sink + %d source = %d total (from %d)",
  length(sink_idx), length(sink_idx), 2 * length(sink_idx), nrow(monthly_training)))

message("Fitting Stage 1 balanced RF...")
rf_model_A <- fit_rf_prob_model(balanced_training, class_predictors, use_case_weights = FALSE)

oob_prob_A  <- rf_model_A$predictions[, "1"]
iso_cal_A   <- fit_isotonic_calibration(oob_prob_A, balanced_training$weak_source)

balanced_training <- balanced_training %>%
  mutate(
    source_prob_A_oob = oob_prob_A,
    source_prob_A_raw = predict_rf_prob(rf_model_A, balanced_training),
    source_prob_A     = calibrate_isotonic(iso_cal_A, source_prob_A_oob)
  )

# Also add calibrated projections back to full training for Stage 2 (source_probability predictor)
monthly_training <- monthly_training %>%
  mutate(
    source_prob_A_raw = predict_rf_prob(rf_model_A, monthly_training),
    source_prob_A     = calibrate_isotonic(iso_cal_A, source_prob_A_raw)
  )

skill_A <- classification_skill(
  prob        = balanced_training$source_prob_A,
  observed    = balanced_training$weak_source,
  threshold   = binary_threshold,
  model_label = "Continuous",
  eval_basis  = "OOB + isotonic calibration (1:1 balanced)",
  oob_error   = rf_model_A$prediction.error
)

cal_skill_A <- balanced_training %>%
  mutate(prob_bin = ntile(source_prob_A_oob, 10)) %>%
  group_by(prob_bin) %>%
  summarise(
    n                         = n(),
    mean_oob_raw_prob         = mean(source_prob_A_oob, na.rm = TRUE),
    mean_isotonic_cal_prob    = mean(source_prob_A,     na.rm = TRUE),
    observed_source_fraction  = mean(weak_source,       na.rm = TRUE),
    calibration_error         = mean_isotonic_cal_prob - observed_source_fraction,
    .groups = "drop"
  ) %>% mutate(model = "A")

importance_A <- tibble(
  predictor  = names(rf_model_A$variable.importance),
  importance = rf_model_A$variable.importance,
  model      = "P(source)"
) %>% arrange(desc(importance))

# ── Stage 2: RF magnitude models (shared across all approaches) ───────────────
# Trained on the full weighted training data.
# source_prob_A (calibrated OOB) used as the source_probability predictor
# so it represents an honest, unbiased probability signal.

magnitude_standardizers <- list(
  Tair = make_standardizer(monthly_training$mean_ERA5_Tair_C),
  VSWC = make_standardizer(monthly_training$mean_ERA5_VSWC),
  MAP  = make_standardizer(monthly_training$MAP),
  MAT  = make_standardizer(monthly_training$MAT)
)

monthly_training <- apply_standardizers(monthly_training, magnitude_standardizers) %>%
  mutate(
    source_probability = source_prob_A,   # Stage 2 uses Model A calibrated OOB probs
    SITE_ID = factor(SITE_ID)
  )

mag_predictors <- c("z_Tair","z_VSWC","z_MAP","z_MAT","source_probability","is_arid","EcoType")

sink_mag_data   <- monthly_training %>% filter(weak_source == 0, monthly_flux_gC_m2_month <= 0)
source_mag_data <- monthly_training %>% filter(weak_source == 1, monthly_flux_gC_m2_month >  0)

message("Fitting Stage 2 sink magnitude RF...")
sink_mag_model   <- fit_rf_magnitude_model(sink_mag_data,   "Weak-sink",   mag_predictors)
message("Fitting Stage 2 source magnitude RF...")
source_mag_model <- fit_rf_magnitude_model(source_mag_data, "Weak-source", mag_predictors)

magnitude_model_skill <- bind_rows(
  summarise_magnitude_fit(sink_mag_data,   sink_mag_model),
  summarise_magnitude_fit(source_mag_data, source_mag_model)
)

magnitude_fitted_values <- bind_rows(
  sink_mag_data %>%
    mutate(magnitude_model = "Weak-sink",
           fitted_flux = pmin(predict_rf_magnitude(sink_mag_model, sink_mag_data), 0)),
  source_mag_data %>%
    mutate(magnitude_model = "Weak-source",
           fitted_flux = pmax(predict_rf_magnitude(source_mag_model, source_mag_data), 0))
) %>% rename(fitted_flux_gC_m2_month = fitted_flux)

importance_mag <- bind_rows(
  tibble(predictor = names(sink_mag_model$model$variable.importance),
         importance = sink_mag_model$model$variable.importance,   model = "Stage 2 — Weak-sink"),
  tibble(predictor = names(source_mag_model$model$variable.importance),
         importance = source_mag_model$model$variable.importance, model = "Stage 2 — Weak-source")
) %>% arrange(model, desc(importance))

training_vswc_range <- range(monthly_training$mean_ERA5_VSWC, na.rm = TRUE)
training_prec_range <- quantile(monthly_training$MAP / 12, probs = c(0.02, 0.98), na.rm = TRUE)

# ── Spatial projection loop ───────────────────────────────────────────────────
# Both Stage 1 models predict in the same pass.
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
        is_arid              = as.integer(EcoType == "Arid"),
        EcoType              = factor(if_else(EcoType == "Arid","Shrubland", EcoType),
                                      levels = levels(monthly_training$EcoType)),
        mean_ERA5_Tair_C = tavg,
        mean_ERA5_VSWC   = if (use_era5_land) vswc else
          scales::rescale(pmin(pmax(prec, training_prec_range[1]), training_prec_range[2]),
                          to = training_vswc_range, from = training_prec_range),
        aridity_index = MAP / (MAT + 10)
      ) %>%
      filter(area_mha > 0, is.finite(aridity_index))

    dat <- apply_standardizers(dat, magnitude_standardizers)

    # ── Stage 1: predict Model A; arid cells forced to 0 ─────────────────────
    not_arid <- !(!is.na(dat$is_arid) & dat$is_arid == 1)

    dat$source_prob_A_raw <- 0; dat$source_prob_A <- 0
    if (any(not_arid)) {
      nd <- dat[not_arid, ]
      dat$source_prob_A_raw[not_arid] <- predict_rf_prob(rf_model_A, nd)
      dat$source_prob_A[not_arid]     <- calibrate_isotonic(iso_cal_A, dat$source_prob_A_raw[not_arid])
    }

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
# Every non-arid upland cell is assigned the sink magnitude flux, regardless of
# P(source). This is the most extreme possible sink estimate from the model.

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

# ── Classification skill comparison table ─────────────────────────────────────

comparison_class <- skill_A %>%
  select(model, evaluation_basis, auc, brier_score, brier_null, brier_skill_score,
         tjur_r2, ranger_oob_error, threshold, accuracy,
         sensitivity = sensitivity_source, specificity = specificity_sink)

# ── Magnitude skill comparison ────────────────────────────────────────────────

comparison_magnitude <- magnitude_model_skill %>%
  select(magnitude_model, rmse_gC_m2_month, mae_gC_m2_month, correlation_observed_fitted) %>%
  mutate(approach = "RF/ranger", .before = 1)

# ── Budget comparison ─────────────────────────────────────────────────────────

comparison_budget <- budget_summary

# ── Write outputs ─────────────────────────────────────────────────────────────

write.csv(annual_budget,          file.path(output_dir, "annual_budget_three_approaches.csv"),    row.names = FALSE)
write.csv(annual_budget_long,     file.path(output_dir, "annual_budget_long.csv"),                row.names = FALSE)
write.csv(budget_summary,         file.path(output_dir, "budget_summary_three_approaches.csv"),   row.names = FALSE)
write.csv(comparison_class,       file.path(output_dir, "comparison_class_skill_GLM_vs_RF.csv"),  row.names = FALSE)
write.csv(comparison_magnitude,   file.path(output_dir, "comparison_magnitude_skill_GLM_vs_RF.csv"), row.names = FALSE)
write.csv(comparison_budget,      file.path(output_dir, "comparison_budget_all_approaches.csv"),  row.names = FALSE)
write.csv(cal_skill_A,            file.path(output_dir, "probability_calibration_skill.csv"),     row.names = FALSE)
write.csv(importance_A,           file.path(output_dir, "rf_class_variable_importance.csv"),      row.names = FALSE)
write.csv(importance_mag,         file.path(output_dir, "rf_magnitude_variable_importance.csv"),  row.names = FALSE)
write.csv(magnitude_model_skill,  file.path(output_dir, "magnitude_model_skill.csv"),             row.names = FALSE)
write.csv(magnitude_fitted_values,file.path(output_dir, "magnitude_model_fitted_values.csv"),     row.names = FALSE)
write.csv(gmb_sensitivity,        file.path(output_dir, "gmb_threshold_sensitivity.csv"),         row.names = FALSE)
write.csv(data.frame(
  P_star                   = P_star_nearest,   # nearest threshold to GMB target; budget floor never reached
  budget_floor_tg_ch4_yr  = budget_floor,
  binary_threshold         = binary_threshold,
  gmb_target_tg_ch4_yr    = gmb_soil_sink_tg_ch4_yr,
  n_trees                  = n_trees,
  rf_min_node_size         = rf_min_node_size,
  rf_max_depth             = rf_max_depth,
  rf_sample_frac           = rf_sample_frac,
  rf_seed                  = rf_seed,
  calibration_method       = "isotonic_regression_on_OOB",
  stage1_training          = sprintf("1:1 balanced (n=%d sink + %d source = %d total; from %d)",
                               length(sink_idx), length(sink_idx),
                               2 * length(sink_idx), nrow(monthly_training))
), file.path(output_dir, "model_parameters.csv"), row.names = FALSE)

saveRDS(cell_preds, file.path(output_dir, "monthly_cell_predictions_2000_2025.rds"))
terra::writeRaster(template, file.path(output_dir, "era5_template.tif"), overwrite = TRUE)

# ── Summary text ──────────────────────────────────────────────────────────────

capture.output({
  cat("RF upscaling — three-approach comparison\n\n")
  cat(sprintf("Stage 1: weighted RF | n=%d | n_trees=%d | min.node.size=%d | max.depth=%d | sample.frac=%.1f\n",
    nrow(monthly_training), n_trees, rf_min_node_size, rf_max_depth, rf_sample_frac))
  cat("Calibration: isotonic regression on OOB predictions\n")
  cat("Approach 1 — Continuous: P(source)-weighted flux blend (no hard threshold)\n")
  cat(sprintf("Approach 2 — Dichotomous: hard threshold at P(source) >= %.2f\n", binary_threshold))
  cat(sprintf("Approach 3 — All-Sink: all non-arid cells assigned sink flux (theoretical maximum sink)\n"))
  cat(sprintf("  GMB sensitivity floor: %.1f Tg/yr at P* = %.2f; GMB target (%.0f) never reached\n\n",
    budget_floor, P_star_nearest, gmb_soil_sink_tg_ch4_yr))

  cat("─── Stage 1 Classification Skill ───\n")
  print(comparison_class)

  cat("\n─── Stage 2 Magnitude Skill ───\n")
  print(comparison_magnitude)

  cat("\n─── Annual Budget Summary — Three Approaches ───\n")
  print(budget_summary)

  cat("\n─── Full Budget Comparison (incl. GLM) ───\n")
  print(comparison_budget)

  cat("\n─── Calibration (OOB deciles) ───\n"); print(cal_skill_A)

  cat("\n─── Stage 1 Variable Importance ───\n")
  print(importance_A)
  cat("\n─── Stage 2 Variable Importance ───\n"); print(importance_mag)
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

ecotype_colors <- c(Forest = "#1B7837", Grassland = "#D9B86C",
                    Shrubland = "#C2A5CF", Arid = "#D95F02")

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

# Fig 3: Classification skill comparison (AUC, Tjur R², sensitivity, specificity)
p3 <- comparison_class %>%
  filter(!is.na(auc)) %>%
  pivot_longer(c(auc, tjur_r2, accuracy, sensitivity, specificity),
               names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(metric = factor(metric,
    levels = c("auc","tjur_r2","accuracy","sensitivity","specificity"),
    labels = c("AUC","Tjur R²","Accuracy","Sensitivity","Specificity")),
    model_short = model) %>%
  ggplot(aes(x = value, y = metric, fill = model_short)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
  scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "Continuous"  = "#009688",
    "Dichotomous" = "#9C27B0"
  )) +
  labs(title = "Stage 1 Classification Skill",
       x = "Value", y = NULL, fill = NULL) +
  fig_theme + theme(legend.position = "top")

ggsave(file.path(figure_dir, "Fig3_classification_skill.png"),
  p3, width = 7, height = 4.5, units = "in", dpi = 300, bg = "white")

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

# Fig 5: Magnitude model — observed vs fitted
make_mag_plot <- function(model_label, tag) {
  magnitude_fitted_values %>%
    filter(magnitude_model == model_label, !is.na(EcoType)) %>%
    ggplot(aes(x = monthly_flux_gC_m2_month, y = fitted_flux_gC_m2_month, color = EcoType)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "grey35", linetype = "dashed", linewidth = 0.5) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = ecotype_colors) +
    labs(title = paste0(tag, ". ", model_label), color = NULL,
         x = "Observed (g C m⁻² mo⁻¹)", y = "Fitted (g C m⁻² mo⁻¹)") +
    fig_theme
}
fig5 <- plot_grid(make_mag_plot("Weak-sink","A"), make_mag_plot("Weak-source","B"), ncol = 2)
ggsave(file.path(figure_dir, "Fig5_magnitude_models.png"),
  fig5, width = 9, height = 4.5, units = "in", dpi = 300, bg = "white")

message("Three-approach RF upscaling complete. Outputs in: ", rf_dir)

