# Monthly condition-based spatial CH4 upscaling trial.
#
# This is a spatial implementation of the monthly class-balanced source
# probability model. When available, it uses downloaded ERA5-Land monthly grids
# plus the 2017 terrestrial ecoregions map stored in:
#   /Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly
#
# If ERA5-Land NetCDF files are absent, the script falls back to WorldClim
# climatology as a proof-of-concept mode.

library(tidyverse)
library(data.table)
library(terra)
library(lme4)
library(patchwork)
library(cowplot)

localdir.ch4 <- Sys.getenv(
  "LOCALDIR_CH4",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/Methane"
)

spatial_dir <- Sys.getenv(
  "MONTHLY_UPSCALING_DIR",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly"
)

if (!dir.exists(localdir.ch4)) {
  stop("CH4 data directory does not exist: ", localdir.ch4)
}
if (!dir.exists(spatial_dir)) {
  dir.create(spatial_dir, recursive = TRUE)
}

output_dir <- file.path(spatial_dir, "OUTPUT")
figure_dir <- file.path(spatial_dir, "FIGURES")
data_dir <- file.path(spatial_dir, "DATA")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

era5_30min_file <- file.path(localdir.ch4, "OUTPUT/NEON_ERA5_gapfilled_30min.csv.gz")
# Static site attributes (SITE_ID, EcoType, MAP, MAT) from 05_NEON_FluxAnalysis.R's
# actively-regenerated summary, not the orphaned OUTPUT/30min_site_behavior.csv.
site_attributes_file <- file.path(localdir.ch4, "OUTPUT/NEON_scale_annual_budget_summary.csv")

worldclim_tavg_zip <- file.path(spatial_dir, "wc2.1_10m_tavg.zip")
worldclim_prec_zip <- file.path(spatial_dir, "wc2.1_10m_prec.zip")
ecoregions_zip <- file.path(spatial_dir, "Ecoregions2017.zip")
era5_land_dir <- file.path(data_dir, "era5_land_monthly")
modis_ecotype_dir <- file.path(data_dir, "modis_mcd12c1_processed")
wad2m_dir <- file.path(data_dir, "wad2m")

required_files <- c(era5_30min_file, site_attributes_file, ecoregions_zip)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

# Threshold is derived from the training distribution after model fitting (see below).
# This percentile is taken over observed source months so the threshold is always
# within the range of values the GLM actually predicts, avoiding extrapolation artefacts.
source_probability_threshold_percentile <- 0.85
source_probability_thresholds <- seq(0.50, 0.95, by = 0.05)  # retained for sensitivity analysis
# No hard inundation threshold; upland area is scaled continuously by (1 - inundation_fraction).
# Cells with 100% inundation are excluded downstream via filter(area_mha > 0).
arid_shrubland_fill_temp_threshold <- 15  # °C monthly mean; NA VSWC filled with 0 for Arid cells above this
arid_ai_threshold <- 15  # de Martonne AI below which NEON training sites are flagged is_arid = 1
                          # (captures JORN ~8.9, SRER ~13.0, MOAB ~10.5; excludes wetter shrublands)
gmb_soil_sink_tg_ch4_yr <- -35
magnitude_constraint_weight <- 0.65
years_to_repeat <- 2000:2025
years_to_process <- years_to_repeat

gC_m2_yr_to_tg_ch4 <- function(flux_gC_m2_yr, area_mha) {
  flux_gC_m2_yr * area_mha * 0.0133333333333333
}

wtd_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

add_class_balance_weights <- function(dat) {
  n_total <- nrow(dat)
  n_class <- dat %>% count(weak_source, name = "n_class")
  dat %>%
    left_join(n_class, by = "weak_source") %>%
    mutate(class_balance_weight = n_total / (length(unique(weak_source)) * n_class)) %>%
    dplyr::select(-n_class)
}

fit_monthly_balanced_model <- function(training_data, formula) {
  training_data <- add_class_balance_weights(training_data)
  glm(
    formula,
    data = training_data,
    family = quasibinomial(),
    weights = class_balance_weight
  )
}

safe_logit <- function(p) {
  qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
}

fit_probability_calibration <- function(training_data) {
  calibration_data <- add_class_balance_weights(training_data) %>%
    mutate(source_probability_logit = safe_logit(source_probability_raw))

  glm(
    weak_source ~ source_probability_logit,
    data = calibration_data,
    family = quasibinomial(),
    weights = class_balance_weight
  )
}

calibrate_source_probability <- function(probability_model, p) {
  as.numeric(predict(
    probability_model,
    newdata = tibble(source_probability_logit = safe_logit(p)),
    type = "response"
  ))
}

auc_rank <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  n_source <- sum(observed == 1)
  n_sink <- sum(observed == 0)
  prediction_rank <- rank(predicted, ties.method = "average")
  (sum(prediction_rank[observed == 1]) - n_source * (n_source + 1) / 2) /
    (n_source * n_sink)
}

make_standardizer <- function(x) {
  center <- mean(x, na.rm = TRUE)
  scale_value <- sd(x, na.rm = TRUE)
  if (!is.finite(scale_value) || scale_value == 0) {
    scale_value <- 1
  }
  list(center = center, scale = scale_value)
}

apply_standardizers <- function(dat, standardizers) {
  dat %>%
    mutate(
      z_Tair = (mean_ERA5_Tair_C - standardizers$Tair$center) / standardizers$Tair$scale,
      z_VSWC = (mean_ERA5_VSWC - standardizers$VSWC$center) / standardizers$VSWC$scale,
      z_MAP = (MAP - standardizers$MAP$center) / standardizers$MAP$scale,
      z_MAT = (MAT - standardizers$MAT$center) / standardizers$MAT$scale
    )
}

fit_flux_magnitude_model <- function(training_data, state_name) {
  training_data <- training_data %>%
    mutate(
      SITE_ID = factor(SITE_ID),
      log_abs_flux = log(pmax(abs(monthly_flux_gC_m2_month), 1e-6))
    )

  # is_arid: fixed effect offset for arid NEON towers (AI < threshold); allows
  # desert observations to calibrate a separate baseline sink rate without
  # requiring Arid as a random-effect level (spatial Arid cells use allow.new.levels=TRUE)
  lmer_formula <- log_abs_flux ~
    z_Tair + z_VSWC + z_MAP + z_MAT + source_probability + is_arid +
    (1 + z_Tair + z_VSWC | EcoType) + (1 | SITE_ID)
  lmer_intercept_formula <- log_abs_flux ~
    z_Tair + z_VSWC + z_MAP + z_MAT + source_probability + is_arid +
    (1 | EcoType) + (1 | SITE_ID)
  lm_formula <- log_abs_flux ~
    z_Tair + z_VSWC + z_MAP + z_MAT + source_probability + is_arid + EcoType

  fit <- tryCatch(
    lme4::lmer(
      lmer_formula,
      data = training_data,
      REML = TRUE,
      control = lme4::lmerControl(
        check.conv.singular = "ignore",
        check.conv.grad = "ignore"
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || isTRUE(lme4::isSingular(fit, tol = 1e-4))) {
    fit <- tryCatch(
      lme4::lmer(
        lmer_intercept_formula,
        data = training_data,
        REML = TRUE,
        control = lme4::lmerControl(
          check.conv.singular = "ignore",
          check.conv.grad = "ignore"
        )
      ),
      error = function(e) NULL
    )
  }

  if (!is.null(fit)) {
    return(list(model = fit, engine = "lmer_log_abs", state = state_name))
  }

  list(
    model = lm(lm_formula, data = training_data),
    engine = "lm",
    state = state_name
  )
}

predict_flux_magnitude <- function(fit, newdata) {
  log_abs_prediction <- if (identical(fit$engine, "lmer_log_abs")) {
    if (!"SITE_ID" %in% names(newdata)) {
      newdata$SITE_ID <- factor("SPATIAL")
    }
    as.numeric(predict(fit$model, newdata = newdata, allow.new.levels = TRUE))
  } else {
    as.numeric(predict(fit$model, newdata = newdata))
  }

  signed_magnitude <- exp(log_abs_prediction)
  if (identical(fit$state, "Weak sink")) {
    -signed_magnitude
  } else {
    signed_magnitude
  }
}

summarise_magnitude_fit <- function(training_data, fit, state_name) {
  fitted_flux <- predict_flux_magnitude(fit, training_data)
  tibble(
    magnitude_model = state_name,
    engine = fit$engine,
    n_observations = nrow(training_data),
    mean_observed_flux_gC_m2_month = mean(training_data$monthly_flux_gC_m2_month, na.rm = TRUE),
    mean_fitted_flux_gC_m2_month = mean(fitted_flux, na.rm = TRUE),
    rmse_gC_m2_month = sqrt(mean((training_data$monthly_flux_gC_m2_month - fitted_flux)^2, na.rm = TRUE)),
    mae_gC_m2_month = mean(abs(training_data$monthly_flux_gC_m2_month - fitted_flux), na.rm = TRUE),
    correlation_observed_fitted = suppressWarnings(cor(training_data$monthly_flux_gC_m2_month, fitted_flux, use = "complete.obs"))
  )
}

unzip_if_needed <- function(zip_file, exdir, expected_pattern) {
  if (!dir.exists(exdir)) {
    dir.create(exdir, recursive = TRUE)
  }
  existing <- list.files(exdir, pattern = expected_pattern, full.names = TRUE)
  if (length(existing) == 0) {
    utils::unzip(zip_file, exdir = exdir)
  }
  invisible(list.files(exdir, pattern = expected_pattern, full.names = TRUE))
}

ecoregion_files <- unzip_if_needed(ecoregions_zip, file.path(data_dir, "ecoregions2017"), "Ecoregions2017\\.shp$")

era5_land_files <- file.path(era5_land_dir, sprintf("era5_land_monthly_%s.nc", years_to_process))
use_era5_land <- all(file.exists(era5_land_files))

read_era5_land_year <- function(year) {
  file <- file.path(era5_land_dir, sprintf("era5_land_monthly_%s.nc", year))
  climate <- rast(file)
  lnames <- names(climate)

  t2m_idx   <- grep("t2m",   lnames, ignore.case = TRUE)
  swvl1_idx <- grep("swvl1", lnames, ignore.case = TRUE)
  swvl2_idx <- grep("swvl2", lnames, ignore.case = TRUE)
  tp_idx    <- grep("(?i)^tp[_=]|^tp$", lnames, perl = TRUE)

  counts <- lengths(list(t2m_idx, swvl1_idx, swvl2_idx, tp_idx))
  if (!all(counts == 12L)) {
    stop(sprintf(
      "Expected 12 layers each for t2m/swvl1/swvl2/tp in %s; found %s.\nAll layer names: %s",
      basename(file),
      paste(c("t2m", "swvl1", "swvl2", "tp"), counts, sep = "=", collapse = ", "),
      paste(lnames, collapse = ", ")
    ))
  }

  tavg_year <- climate[[t2m_idx]] - 273.15
  vswc_year <- (climate[[swvl1_idx]] + climate[[swvl2_idx]]) / 2
  prec_year <- rast(lapply(seq_along(tp_idx), function(i) {
    climate[[tp_idx[i]]] * 1000 *
      lubridate::days_in_month(lubridate::make_date(year, i, 1))
  }))
  names(tavg_year) <- sprintf("tavg_%02d", 1:12)
  names(vswc_year) <- sprintf("vswc_%02d", 1:12)
  names(prec_year) <- sprintf("prec_%02d", 1:12)
  list(tavg = tavg_year, vswc = vswc_year, prec = prec_year)
}

if (use_era5_land) {
  message("Using ERA5-Land monthly grids from ", era5_land_dir)
  first_era5 <- read_era5_land_year(years_to_process[1])
  tavg <- first_era5$tavg
  prec <- first_era5$prec
  template <- tavg[[1]]

  mat_sum <- template
  map_sum <- template
  values(mat_sum) <- 0
  values(map_sum) <- 0
  n_mat_layers <- 0

  for (year in years_to_process) {
    climate_year <- read_era5_land_year(year)
    mat_sum <- mat_sum + sum(climate_year$tavg, na.rm = TRUE)
    map_sum <- map_sum + sum(climate_year$prec, na.rm = TRUE)
    n_mat_layers <- n_mat_layers + 12
  }
  mat <- mat_sum / n_mat_layers
  map <- map_sum / length(years_to_process)
  names(mat) <- "MAT"
  names(map) <- "MAP"
} else {
  if (!file.exists(worldclim_tavg_zip) || !file.exists(worldclim_prec_zip)) {
    stop("ERA5-Land files are missing and WorldClim fallback zip files are not available.")
  }
  message("ERA5-Land files not found; using WorldClim climatology fallback.")
  tavg_files <- unzip_if_needed(worldclim_tavg_zip, file.path(data_dir, "worldclim_tavg"), "tavg_.*\\.tif$")
  prec_files <- unzip_if_needed(worldclim_prec_zip, file.path(data_dir, "worldclim_prec"), "prec_.*\\.tif$")
  tavg <- rast(sort(tavg_files))
  prec <- rast(sort(prec_files))
  names(tavg) <- sprintf("tavg_%02d", 1:12)
  names(prec) <- sprintf("prec_%02d", 1:12)
  template <- tavg[[1]]
  mat <- mean(tavg, na.rm = TRUE)
  names(mat) <- "MAT"
  map <- sum(prec, na.rm = TRUE)
  names(map) <- "MAP"
}

cell_area_mha <- cellSize(template, unit = "m") / 1e10
names(cell_area_mha) <- "area_mha"

ecoregions <- vect(ecoregion_files[1])
ecoregion_names <- names(ecoregions)
biome_field <- case_when(
  "BIOME_NUM" %in% ecoregion_names ~ "BIOME_NUM",
  "BIOME" %in% ecoregion_names ~ "BIOME",
  TRUE ~ NA_character_
)
if (is.na(biome_field)) {
  stop("Could not find BIOME_NUM/BIOME field in Ecoregions2017 shapefile.")
}

biome_raster <- rasterize(ecoregions, template, field = biome_field, touches = TRUE)
names(biome_raster) <- "biome_num"

ecotype_raster <- subst(
  biome_raster,
  from = 1:14,
  to = c(
    2, 2, 2, 2, 2, 2, # forest biomes
    3, 3, NA, 3,      # grassland/savanna; flooded grasslands excluded
    4, 4, 4, NA       # tundra/scrub/desert as shrubland; mangroves excluded
  )
)
names(ecotype_raster) <- "ecotype_code"

ecotype_lookup <- tibble(
  ecotype_code = c(1, 2, 3, 4, 5),
  EcoType = c("Cropland", "Forest", "Grassland", "Shrubland", "Arid")
) %>%
  filter(EcoType != "Cropland")

modis_ecotype_files <- list.files(
  modis_ecotype_dir,
  pattern = "^MODIS_MCD12C1_ecotype_[0-9]{4}\\.tif$",
  full.names = TRUE
)
use_modis_landcover <- length(modis_ecotype_files) > 0
modis_years <- if (use_modis_landcover) {
  as.integer(str_match(basename(modis_ecotype_files), "([0-9]{4})\\.tif$")[, 2])
} else {
  integer()
}

get_ecotype_raster <- function(year) {
  if (use_modis_landcover) {
    selected_year <- modis_years[which.min(abs(modis_years - year))]
    modis_file <- file.path(modis_ecotype_dir, sprintf("MODIS_MCD12C1_ecotype_%s.tif", selected_year))
    r <- rast(modis_file)
    if (!same.crs(r, template)) {
      crs(r) <- crs(template)
    }
    if (!compareGeom(r, template, stopOnError = FALSE)) {
      r <- resample(r, template, method = "near")
    }
    names(r) <- "ecotype_code"
    return(r)
  }
  ecotype_raster
}

wad2m_files <- list.files(wad2m_dir, pattern = "\\.nc$", full.names = TRUE)
use_wad2m_inundation <- length(wad2m_files) > 0
wad2m <- if (use_wad2m_inundation) rast(wad2m_files[1]) else NULL
wad2m_available_years <- if (use_wad2m_inundation) 2000:2020 else integer()

get_inundation_fraction <- function(year, month) {
  if (!use_wad2m_inundation) {
    r <- template
    values(r) <- 0
    names(r) <- "inundation_fraction"
    return(r)
  }
  selected_year <- if (year %in% wad2m_available_years) year else NA_integer_
  if (is.na(selected_year)) {
    selected_layers <- seq(month, nlyr(wad2m), by = 12)
    r <- mean(wad2m[[selected_layers]], na.rm = TRUE)
  } else {
    layer_index <- (selected_year - min(wad2m_available_years)) * 12 + month
    r <- wad2m[[layer_index]]
  }
  if (!compareGeom(r, template, stopOnError = FALSE)) {
    r <- resample(r, template, method = "bilinear")
  }
  names(r) <- "inundation_fraction"
  r
}

site_attributes <- read.csv(site_attributes_file) %>%
  mutate(SITE_ID = as.character(SITE_ID))

upland_sites <- site_attributes %>%
  filter(
    !is.na(EcoType),
    !str_detect(EcoType, regex("wetland|inundat|flood|marsh|swamp|bog|fen|lake|rice|crop|agri", ignore_case = TRUE))
  ) %>%
  distinct(SITE_ID, EcoType, MAP, MAT)

era5_30min <- data.table::fread(era5_30min_file) %>%
  as_tibble() %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    Year = as.integer(Year),
    month = as.integer(month),
    ERA5_Tair_C = as.numeric(ERA5_Tair_C),
    ERA5_VSWC = as.numeric(ERA5_VSWC),
    gapfilled_CH4_mgC_30min = as.numeric(gapfilled_CH4_mgC_30min)
  ) %>%
  inner_join(upland_sites %>% dplyr::select(SITE_ID, EcoType), by = c("SITE_ID", "EcoType")) %>%
  filter(is.finite(Year), is.finite(month), is.finite(ERA5_Tair_C), is.finite(ERA5_VSWC))

monthly_training <- era5_30min %>%
  reframe(
    .by = c(SITE_ID, EcoType, Year, month),
    monthly_budget_mgC_m2 = sum(gapfilled_CH4_mgC_30min, na.rm = TRUE),
    mean_ERA5_Tair_C = mean(ERA5_Tair_C, na.rm = TRUE),
    mean_ERA5_VSWC = mean(ERA5_VSWC, na.rm = TRUE)
  ) %>%
  left_join(upland_sites, by = c("SITE_ID", "EcoType")) %>%
  mutate(
    EcoType = factor(EcoType, levels = c("Forest", "Grassland", "Shrubland")),
    weak_source = as.integer(monthly_budget_mgC_m2 > 0),
    monthly_flux_gC_m2_month = monthly_budget_mgC_m2 / 1000,
    MAP = as.numeric(MAP),
    MAT = as.numeric(MAT),
    # de Martonne aridity index: low values = hyper-arid; captures moisture-temperature
    # balance that the additive Tair + VSWC terms cannot express independently
    aridity_index = MAP / (MAT + 10),
    # Binary indicator for arid NEON towers (AI < threshold); used as fixed effect
    # in the sink magnitude lmer so desert observations calibrate arid sink rates
    is_arid = as.integer(!is.na(aridity_index) & aridity_index < arid_ai_threshold)
  ) %>%
  filter(
    is.finite(weak_source),
    is.finite(mean_ERA5_Tair_C),
    is.finite(mean_ERA5_VSWC),
    is.finite(MAP),
    is.finite(MAT),
    is.finite(aridity_index)
  )

# Tair * VSWC interaction: temperature drives methanogenesis only when moisture
# is present; the interaction suppresses P(source) in hot-dry cells.
# Aridity index (MAP / (MAT + 10)): explicitly penalises hyper-arid regimes
# where neither soil moisture nor long-term climate supports source behaviour.
monthly_formula <- weak_source ~ EcoType +
  scale(mean_ERA5_Tair_C) * scale(mean_ERA5_VSWC) +
  scale(MAP) + scale(MAT) + scale(aridity_index)

monthly_model <- fit_monthly_balanced_model(monthly_training, monthly_formula)

monthly_training <- monthly_training %>%
  mutate(source_probability_raw = as.numeric(predict(monthly_model, newdata = monthly_training, type = "response")))

probability_calibration_model <- fit_probability_calibration(monthly_training)

monthly_training <- monthly_training %>%
  mutate(
    source_probability = calibrate_source_probability(probability_calibration_model, source_probability_raw)
  )

# Derive threshold from the training distribution: Xth percentile of predicted
# P(source) among observed source months. This guarantees the threshold lies within
# the model's output range and that some training source months exceed it.
source_probability_threshold <- quantile(
  monthly_training$source_probability[monthly_training$weak_source == 1],
  probs = source_probability_threshold_percentile,
  na.rm = TRUE
)
message(sprintf(
  "Source probability threshold (%.0fth pct of source months): %.4f",
  source_probability_threshold_percentile * 100,
  source_probability_threshold
))

magnitude_standardizers <- list(
  Tair = make_standardizer(monthly_training$mean_ERA5_Tair_C),
  VSWC = make_standardizer(monthly_training$mean_ERA5_VSWC),
  MAP = make_standardizer(monthly_training$MAP),
  MAT = make_standardizer(monthly_training$MAT)
)

monthly_training <- apply_standardizers(monthly_training, magnitude_standardizers) %>%
  mutate(SITE_ID = factor(SITE_ID))

probability_calibration_skill <- monthly_training %>%
  mutate(probability_bin = ntile(source_probability, 10)) %>%
  group_by(probability_bin) %>%
  summarise(
    n = dplyr::n(),
    mean_raw_source_probability = mean(source_probability_raw, na.rm = TRUE),
    mean_calibrated_source_probability = mean(source_probability, na.rm = TRUE),
    observed_source_fraction = mean(weak_source, na.rm = TRUE),
    .groups = "drop"
  )

probability_threshold_training_class <- as.integer(monthly_training$source_probability >= source_probability_threshold)
class_probability_model_skill <- tibble(
  model = "No-season monthly source-probability GLM",
  n_site_months = nrow(monthly_training),
  n_sink_months = sum(monthly_training$weak_source == 0),
  n_source_months = sum(monthly_training$weak_source == 1),
  source_fraction = mean(monthly_training$weak_source),
  auc = auc_rank(monthly_training$weak_source, monthly_training$source_probability),
  brier_score = mean((monthly_training$weak_source - monthly_training$source_probability)^2),
  tjur_r2 = mean(monthly_training$source_probability[monthly_training$weak_source == 1]) -
    mean(monthly_training$source_probability[monthly_training$weak_source == 0]),
  null_deviance = monthly_model$null.deviance,
  residual_deviance = monthly_model$deviance,
  deviance_explained = 1 - monthly_model$deviance / monthly_model$null.deviance,
  threshold = source_probability_threshold,
  accuracy = mean(probability_threshold_training_class == monthly_training$weak_source),
  sensitivity_source = sum(probability_threshold_training_class == 1 & monthly_training$weak_source == 1) /
    sum(monthly_training$weak_source == 1),
  specificity_sink = sum(probability_threshold_training_class == 0 & monthly_training$weak_source == 0) /
    sum(monthly_training$weak_source == 0),
  precision_source = sum(probability_threshold_training_class == 1 & monthly_training$weak_source == 1) /
    sum(probability_threshold_training_class == 1),
  predicted_source_fraction = mean(probability_threshold_training_class == 1)
)

sink_magnitude_training <- monthly_training %>%
  filter(weak_source == 0, monthly_flux_gC_m2_month <= 0)

source_magnitude_training <- monthly_training %>%
  filter(weak_source == 1, monthly_flux_gC_m2_month > 0)

sink_magnitude_model <- fit_flux_magnitude_model(sink_magnitude_training, "Weak sink")
source_magnitude_model <- fit_flux_magnitude_model(source_magnitude_training, "Weak source")

magnitude_model_skill <- bind_rows(
  summarise_magnitude_fit(sink_magnitude_training, sink_magnitude_model, "Weak sink"),
  summarise_magnitude_fit(source_magnitude_training, source_magnitude_model, "Weak source")
)

magnitude_fitted_values <- bind_rows(
  sink_magnitude_training %>%
    mutate(
      magnitude_model = "Weak sink",
      fitted_flux_gC_m2_month = predict_flux_magnitude(sink_magnitude_model, sink_magnitude_training)
    ),
  source_magnitude_training %>%
    mutate(
      magnitude_model = "Weak source",
      fitted_flux_gC_m2_month = predict_flux_magnitude(source_magnitude_model, source_magnitude_training)
    )
) %>%
  mutate(
    fitted_flux_gC_m2_month = if_else(
      magnitude_model == "Weak sink",
      pmin(fitted_flux_gC_m2_month, 0),
      pmax(fitted_flux_gC_m2_month, 0)
    )
  )

# ── State–magnitude consistency diagnostic ─────────────────────────────────────
# Apply BOTH magnitude models to every training observation. Compare GLM-predicted
# state against the sign of each model's fitted flux to flag site-months where
# the state model and magnitude model disagree.
state_magnitude_diagnostic <- monthly_training %>%
  mutate(
    glm_predicted_state    = if_else(source_probability >= source_probability_threshold,
                                     "Weak source", "Weak sink"),
    observed_state         = if_else(weak_source == 1, "Weak source", "Weak sink"),
    state_agreement        = glm_predicted_state == observed_state,
    fitted_sink_flux       = pmin(predict_flux_magnitude(sink_magnitude_model,   .),  0),
    fitted_source_flux     = pmax(predict_flux_magnitude(source_magnitude_model, .),  0),
    # For each GLM-predicted state, does the matching magnitude model give the right sign?
    magnitude_sign_correct = case_when(
      glm_predicted_state == "Weak sink"   & fitted_sink_flux   < 0 ~ TRUE,
      glm_predicted_state == "Weak source" & fitted_source_flux > 0 ~ TRUE,
      TRUE                                                           ~ FALSE
    ),
    # Full agreement: GLM state matches observed AND magnitude sign matches GLM state
    full_agreement         = state_agreement & magnitude_sign_correct
  )

state_magnitude_agreement_summary <- state_magnitude_diagnostic %>%
  group_by(EcoType) %>%
  summarise(
    n                          = dplyr::n(),
    pct_state_agreement        = round(100 * mean(state_agreement,        na.rm = TRUE), 1),
    pct_magnitude_sign_correct = round(100 * mean(magnitude_sign_correct, na.rm = TRUE), 1),
    pct_full_agreement         = round(100 * mean(full_agreement,         na.rm = TRUE), 1),
    n_glm_source               = sum(glm_predicted_state == "Weak source"),
    n_glm_sink                 = sum(glm_predicted_state == "Weak sink"),
    n_obs_source               = sum(observed_state      == "Weak source"),
    n_obs_sink                 = sum(observed_state      == "Weak sink"),
    .groups = "drop"
  )

state_magnitude_flagged <- state_magnitude_diagnostic %>%
  filter(!full_agreement) %>%
  dplyr::select(
    SITE_ID, EcoType, Year, month,
    monthly_flux_gC_m2_month, observed_state, glm_predicted_state,
    source_probability, state_agreement, magnitude_sign_correct,
    fitted_sink_flux, fitted_source_flux,
    mean_ERA5_Tair_C, mean_ERA5_VSWC, MAP, MAT, aridity_index
  ) %>%
  arrange(EcoType, SITE_ID, Year, month)

write.csv(state_magnitude_agreement_summary,
          file.path(output_dir, "state_magnitude_agreement_summary.csv"), row.names = FALSE)
write.csv(state_magnitude_flagged,
          file.path(output_dir, "state_magnitude_flagged_disagreements.csv"), row.names = FALSE)
# ───────────────────────────────────────────────────────────────────────────────

training_vswc_range <- range(monthly_training$mean_ERA5_VSWC, na.rm = TRUE)
training_prec_range <- quantile(monthly_training$MAP / 12, probs = c(0.02, 0.98), na.rm = TRUE)


monthly_components <- list()
monthly_cell_predictions <- list()
monthly_threshold_components <- list()
monthly_expected_components <- list()
monthly_continuous_components <- list()
component_index <- 1

for (year in years_to_process) {
  if (use_era5_land) {
    climate_year <- read_era5_land_year(year)
  } else {
    climate_year <- list(tavg = tavg, vswc = NULL, prec = prec)
  }

  for (m in 1:12) {
    message("Predicting year ", year, " month ", m)
    year_ecotype_raster <- get_ecotype_raster(year)
    inundation_fraction <- get_inundation_fraction(year, m)
    # WAD2M uses -9999 as fill for non-wetland land cells; terra reads those as NA.
    # For any cell with NA inundation, treat as 0 (no inundation).
    # Ocean cells are still excluded downstream by the !is.na(ecotype_code) filter.
    inundation_fraction[is.na(inundation_fraction)] <- 0
    month_stack <- c(
      year_ecotype_raster,
      inundation_fraction,
      climate_year$tavg[[m]],
      climate_year$prec[[m]],
      mat,
      map,
      cell_area_mha
    )
    names(month_stack) <- c("ecotype_code", "inundation_fraction", "tavg", "prec", "MAT", "MAP", "area_mha")

    if (use_era5_land) {
      vswc_layer <- climate_year$vswc[[m]]
      tavg_layer <- climate_year$tavg[[m]]
      # Fill NA VSWC with 0 for Arid cells (ecotype_code 5) where monthly
      # temperature confirms an arid/hot regime. Hyper-arid areas (Sahara,
      # Arabian Peninsula, etc.) often have missing VSWC in ERA5-Land because
      # the land surface model does not converge for bare rock; temperature
      # check prevents applying the fill to cold barren/tundra cells.
      shrubland_arid_mask <- !is.na(year_ecotype_raster) &
        year_ecotype_raster == 5 &
        !is.na(tavg_layer) &
        tavg_layer > arid_shrubland_fill_temp_threshold
      vswc_layer[is.na(vswc_layer) & shrubland_arid_mask] <- 0
      month_stack <- c(month_stack, vswc_layer)
      names(month_stack)[nlyr(month_stack)] <- "vswc"
    }

    dat_raw <- as.data.frame(month_stack, xy = TRUE, cells = TRUE, na.rm = TRUE)
    # Diagnostic: report ecotype code distribution to track where Arid cells survive.
    if (year == min(years) && m == 1) {
      ecotype_counts <- table(dat_raw$ecotype_code, useNA = "always")
      message(sprintf(
        "[DIAG] Year %d Month %d — raw raster cells by ecotype code: %s",
        year, m,
        paste(names(ecotype_counts), ecotype_counts, sep = "=", collapse = ", ")
      ))
    }

    dat <- dat_raw %>%
      filter(!is.na(ecotype_code), !is.na(tavg), !is.na(prec), !is.na(MAT), !is.na(MAP), !is.na(area_mha)) %>%
      mutate(ecotype_code = as.integer(ecotype_code)) %>%
      inner_join(ecotype_lookup, by = "ecotype_code") %>%
      mutate(
        Year = year,
        month = m,
        inundation_fraction = pmin(pmax(inundation_fraction, 0), 1),
        upland_area_fraction = 1 - inundation_fraction,
        area_mha = area_mha * upland_area_fraction,
        # Compute is_arid from character EcoType before factor conversion;
        # then remap Arid → Shrubland so lmer random effects have a valid level.
        # The is_arid fixed effect in the sink model carries the arid offset.
        is_arid = as.integer(EcoType == "Arid"),
        EcoType = factor(
          if_else(EcoType == "Arid", "Shrubland", EcoType),
          levels = levels(monthly_training$EcoType)
        ),
        mean_ERA5_Tair_C = tavg,
        mean_ERA5_VSWC = if (use_era5_land) {
          vswc
        } else {
          scales::rescale(
            pmin(pmax(prec, training_prec_range[1]), training_prec_range[2]),
            to = training_vswc_range,
            from = training_prec_range
          )
        },
        aridity_index = MAP / (MAT + 10)
      ) %>%
      filter(area_mha > 0, is.finite(aridity_index))

    # Diagnostic: confirm Arid cells survive the full assembly pipeline.
    if (year == min(years) && m == 1) {
      message(sprintf(
        "[DIAG] Year %d Month %d — dat rows after join+filter: %d total, %d Arid (is_arid==1)",
        year, m, nrow(dat), sum(dat$is_arid == 1, na.rm = TRUE)
      ))
    }

    dat <- apply_standardizers(dat, magnitude_standardizers) %>%
      mutate(
        SITE_ID = factor("SPATIAL")
      )

    # Source probability: run GLM only for Forest/Grassland/Shrubland cells.
    # Arid cells are forced to source_probability = 0 (always Weak sink).
    arid_mask <- !is.na(dat$is_arid) & dat$is_arid == 1
    dat$source_probability_raw <- 0
    if (any(!arid_mask)) {
      dat$source_probability_raw[!arid_mask] <- as.numeric(
        predict(monthly_model, newdata = dat[!arid_mask, ], type = "response")
      )
    }
    dat$source_probability <- 0
    if (any(!arid_mask)) {
      dat$source_probability[!arid_mask] <- calibrate_source_probability(
        probability_calibration_model,
        dat$source_probability_raw[!arid_mask]
      )
    }
    dat$predicted_sink_flux_gC_m2_month <- pmin(
      predict_flux_magnitude(sink_magnitude_model, dat),
      0
    )
    dat$predicted_source_flux_gC_m2_month <- pmax(
      predict_flux_magnitude(source_magnitude_model, dat),
      0
    )

    dat <- dat %>%
      mutate(
        selected_exchange_class = if_else(
          source_probability >= source_probability_threshold,
          "Weak source",
          "Weak sink"
        ),
        expected_flux_gC_m2_month = if_else(
          selected_exchange_class == "Weak source",
          predicted_source_flux_gC_m2_month,
          predicted_sink_flux_gC_m2_month
        ),
        expected_flux_continuous_gC_m2_month = source_probability * predicted_source_flux_gC_m2_month +
          (1 - source_probability) * predicted_sink_flux_gC_m2_month
      )

    monthly_cell_predictions[[component_index]] <- dat %>%
      dplyr::select(
        Year,
        cell,
        x,
        y,
        month,
        EcoType,
        is_arid,
        area_mha,
        inundation_fraction,
        upland_area_fraction,
        source_probability_raw,
        source_probability,
        predicted_sink_flux_gC_m2_month,
        predicted_source_flux_gC_m2_month,
        expected_flux_gC_m2_month,
        selected_exchange_class
      )

    monthly_expected_components[[component_index]] <- dat %>%
      mutate(monthly_tg_ch4 = gC_m2_yr_to_tg_ch4(expected_flux_gC_m2_month, area_mha)) %>%
      group_by(Year, month, EcoType, is_arid) %>%
      summarise(
        area_mha = sum(area_mha, na.rm = TRUE),
        monthly_tg_ch4 = sum(monthly_tg_ch4, na.rm = TRUE),
        mean_monthly_rate_gC_m2_month = wtd_mean(expected_flux_gC_m2_month, .data$area_mha),
        mean_source_probability = wtd_mean(source_probability, .data$area_mha),
        mean_sink_flux_gC_m2_month = wtd_mean(predicted_sink_flux_gC_m2_month, .data$area_mha),
        mean_source_flux_gC_m2_month = wtd_mean(predicted_source_flux_gC_m2_month, .data$area_mha),
        mean_tavg_C = wtd_mean(tavg, .data$area_mha),
        mean_prec_mm = wtd_mean(prec, .data$area_mha),
        mean_vswc = wtd_mean(mean_ERA5_VSWC, .data$area_mha),
        mean_inundation_fraction = wtd_mean(inundation_fraction, .data$area_mha),
        .groups = "drop"
      )

    monthly_continuous_components[[component_index]] <- dat %>%
      mutate(monthly_tg_ch4 = gC_m2_yr_to_tg_ch4(expected_flux_continuous_gC_m2_month, area_mha)) %>%
      group_by(Year, month, EcoType) %>%
      summarise(
        area_mha = sum(area_mha, na.rm = TRUE),
        monthly_tg_ch4 = sum(monthly_tg_ch4, na.rm = TRUE),
        .groups = "drop"
      )

    monthly_components[[component_index]] <- dat %>%
      group_by(Year, month, EcoType, is_arid, selected_exchange_class) %>%
      summarise(
        total_area_mha = sum(area_mha, na.rm = TRUE),
        mean_source_probability = wtd_mean(source_probability, .data$area_mha),
        mean_tavg_C = wtd_mean(tavg, .data$area_mha),
        mean_prec_mm = wtd_mean(prec, .data$area_mha),
        mean_vswc = wtd_mean(mean_ERA5_VSWC, .data$area_mha),
        mean_inundation_fraction = wtd_mean(inundation_fraction, .data$area_mha),
        .groups = "drop"
      ) %>%
      rename(area_mha = total_area_mha)

    component_index <- component_index + 1
  }
}

monthly_area_by_class <- bind_rows(monthly_components)
monthly_cell_predictions <- bind_rows(monthly_cell_predictions)
monthly_expected_flux_components <- bind_rows(monthly_expected_components)
monthly_continuous_flux_components <- bind_rows(monthly_continuous_components)

annual_source_sink_area <- monthly_area_by_class %>%
  group_by(Year, month, selected_exchange_class) %>%
  summarise(monthly_area_mha = sum(area_mha, na.rm = TRUE), .groups = "drop") %>%
  group_by(Year, selected_exchange_class) %>%
  summarise(mean_monthly_area_mha = mean(monthly_area_mha, na.rm = TRUE), .groups = "drop") %>%
  mutate(input_climate_note = if_else(
    use_era5_land,
    "ERA5-Land monthly temperature and soil moisture vary by year.",
    "WorldClim monthly climatology repeated for each year; area series is flat by construction."
  ))

annual_cell_class_2025 <- monthly_cell_predictions %>%
  filter(Year == 2025) %>%
  group_by(cell, x, y, EcoType) %>%
  summarise(
    area_mha = first(area_mha),
    source_months = sum(selected_exchange_class == "Weak source", na.rm = TRUE),
    mean_source_probability = mean(source_probability, na.rm = TRUE),
    annual_exchange_class = if_else(source_months >= 6, "Weak source", "Weak sink"),
    .groups = "drop"
  )

monthly_cell_class_change_summary <- monthly_cell_predictions %>%
  group_by(cell, x, y, EcoType) %>%
  summarise(
    mean_monthly_area_mha = mean(area_mha, na.rm = TRUE),
    n_months = dplyr::n(),
    source_months = sum(selected_exchange_class == "Weak source", na.rm = TRUE),
    sink_months = sum(selected_exchange_class == "Weak sink", na.rm = TRUE),
    mean_source_probability = mean(source_probability, na.rm = TRUE),
    min_source_probability = min(source_probability, na.rm = TRUE),
    max_source_probability = max(source_probability, na.rm = TRUE),
    changed_monthly_class = n_distinct(selected_exchange_class) > 1,
    .groups = "drop"
  )

annual_cell_class_all_years <- monthly_cell_predictions %>%
  group_by(Year, cell, x, y, EcoType) %>%
  summarise(
    mean_monthly_area_mha = mean(area_mha, na.rm = TRUE),
    source_months = sum(selected_exchange_class == "Weak source", na.rm = TRUE),
    mean_source_probability = mean(source_probability, na.rm = TRUE),
    annual_exchange_class = if_else(source_months >= 6, "Weak source", "Weak sink"),
    .groups = "drop"
  )

annual_cell_class_change_summary <- annual_cell_class_all_years %>%
  group_by(cell, x, y, EcoType) %>%
  summarise(
    mean_monthly_area_mha = mean(mean_monthly_area_mha, na.rm = TRUE),
    n_years = dplyr::n(),
    source_years = sum(annual_exchange_class == "Weak source", na.rm = TRUE),
    sink_years = sum(annual_exchange_class == "Weak sink", na.rm = TRUE),
    mean_source_probability = mean(mean_source_probability, na.rm = TRUE),
    changed_annual_class = n_distinct(annual_exchange_class) > 1,
    .groups = "drop"
  )

cell_class_change_totals <- bind_rows(
  monthly_cell_class_change_summary %>%
    summarise(
      class_change_scale = "Monthly over 2000-2025",
      n_cells = dplyr::n(),
      n_cells_changed_class = sum(changed_monthly_class, na.rm = TRUE),
      percent_cells_changed_class = 100 * n_cells_changed_class / n_cells,
      total_area_mha = sum(mean_monthly_area_mha, na.rm = TRUE),
      changed_area_mha = sum(mean_monthly_area_mha[changed_monthly_class], na.rm = TRUE),
      percent_area_changed_class = 100 * changed_area_mha / total_area_mha
    ),
  annual_cell_class_change_summary %>%
    summarise(
      class_change_scale = "Annual dominant class over 2000-2025",
      n_cells = dplyr::n(),
      n_cells_changed_class = sum(changed_annual_class, na.rm = TRUE),
      percent_cells_changed_class = 100 * n_cells_changed_class / n_cells,
      total_area_mha = sum(mean_monthly_area_mha, na.rm = TRUE),
      changed_area_mha = sum(mean_monthly_area_mha[changed_annual_class], na.rm = TRUE),
      percent_area_changed_class = 100 * changed_area_mha / total_area_mha
    )
  )

annual_hybrid_flux_2000_2025 <- monthly_expected_flux_components %>%
  group_by(Year) %>%
  summarise(
    n_monthly_ecotype_components = dplyr::n(),
    annual_net_exchange_tg_ch4_yr = sum(monthly_tg_ch4, na.rm = TRUE),
    global_methane_budget_soil_sink_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_global_soil_sink_magnitude = 100 * annual_net_exchange_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr),
    .groups = "drop"
  ) %>%
  mutate(magnitude_model = "Binary: hard-threshold + lmer conditional magnitude")

annual_continuous_flux_2000_2025 <- monthly_continuous_flux_components %>%
  group_by(Year) %>%
  summarise(
    n_monthly_ecotype_components = dplyr::n(),
    annual_net_exchange_tg_ch4_yr = sum(monthly_tg_ch4, na.rm = TRUE),
    global_methane_budget_soil_sink_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_global_soil_sink_magnitude = 100 * annual_net_exchange_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr),
    .groups = "drop"
  ) %>%
  mutate(magnitude_model = "Continuous: P(source) × source mag + (1-P) × sink mag")

annual_expected_flux_2000_2025 <- bind_rows(
  annual_hybrid_flux_2000_2025,
  annual_continuous_flux_2000_2025
) %>%
  mutate(
    input_climate_note = if_else(
      use_era5_land,
      "ERA5-Land monthly temperature and soil moisture vary by year.",
      "WorldClim monthly climatology repeated for each year; no interannual ERA5 variability included."
    )
  )

input_notes <- tribble(
  ~item, ~note,
  "Spatial data directory", spatial_dir,
  "Temperature input", if_else(use_era5_land, "ERA5-Land monthly 2 m temperature, converted from K to degrees C.", "WorldClim v2.1 10-minute monthly average temperature climatology."),
  "Precipitation input", if_else(use_era5_land, "ERA5-Land total precipitation converted to monthly mm and used for long-term MAP.", "WorldClim v2.1 10-minute monthly precipitation climatology."),
  "Land-cover/ecosystem input", if_else(use_modis_landcover, "Annual MODIS MCD12C1 processed ecosystem rasters; nearest available year used for missing endpoints.", "Ecoregions2017 terrestrial biome polygons rasterized to the climate grid."),
  "Wetland/inundation exclusion", if_else(use_wad2m_inundation, "Monthly WAD2M inundation fraction used to continuously scale upland area by (1 - inundation_fraction); no hard threshold applied. Cells with 100% inundation are excluded (area_mha = 0). Hyper-arid Shrubland cells with NA ERA5-Land VSWC and monthly mean temperature > 15 C are assigned VSWC = 0 to recover bare-soil cells missing from the ERA5 land surface model.", "Flooded grasslands/savannas and mangroves excluded from the ecoregion biome map; no dynamic inundation fraction file was available."),
  "Cropland exclusion", "Cropland excluded from upscaling. NEON cropland sites are likely not representative of global dryland cropland methane dynamics, and the source probability GLM classified ~91% of global cropland area as Weak source. Analysis covers natural upland ecosystems: Forest, Grassland, Shrubland only.",
  "Soil moisture input", if_else(use_era5_land, "ERA5-Land volumetric soil water layers 1 and 2 averaged as the monthly VSWC predictor.", "ERA5/ERA5-Land monthly soil moisture was not available locally; monthly precipitation was rescaled to the NEON training VSWC range as a temporary proxy."),
  "Magnitude model", "Class probability is estimated with the balanced monthly probability GLM (EcoType, ERA5-Land Tair, VSWC, MAP, MAT), then calibrated. Separate lmer log-absolute-magnitude models for weak-sink and weak-source months predict signed flux from ERA5-Land Tair, VSWC, MAP, MAT, calibrated source probability, and EcoType/site random effects. Expected flux = P(source) * source_magnitude + (1-P(source)) * sink_magnitude, applied cell-by-cell. No static rates used.",
  "Interannual treatment", if_else(use_era5_land, "Monthly climate predictors vary by year from 2000-2025.", "The same monthly climatology is repeated for 2000-2025; outputs are not year-specific ERA5 estimates."),
  "Interpretation", if_else(use_era5_land, "ERA5-Land spatial upscaling trial; land-cover and inundation are dynamic only when processed MODIS/WAD2M files are present.", "Proof-of-concept spatial mechanics only; not a final global methane budget estimate.")
)

write.csv(monthly_area_by_class, file.path(output_dir, "monthly_area_by_ecotype_class.csv"), row.names = FALSE)
write.csv(monthly_expected_flux_components, file.path(output_dir, "monthly_expected_flux_components.csv"), row.names = FALSE)
write.csv(annual_expected_flux_2000_2025, file.path(output_dir, "annual_expected_flux_2000_2025.csv"), row.names = FALSE)
write.csv(magnitude_model_skill, file.path(output_dir, "magnitude_model_skill.csv"), row.names = FALSE)
write.csv(magnitude_fitted_values, file.path(output_dir, "magnitude_model_fitted_values.csv"), row.names = FALSE)
write.csv(class_probability_model_skill, file.path(output_dir, "class_probability_model_skill.csv"), row.names = FALSE)
write.csv(probability_calibration_skill, file.path(output_dir, "probability_calibration_skill.csv"), row.names = FALSE)
write.csv(monthly_cell_class_change_summary, file.path(output_dir, "monthly_cell_class_change_summary_2000_2025.csv"), row.names = FALSE)
write.csv(annual_cell_class_all_years, file.path(output_dir, "annual_cell_class_all_years_2000_2025.csv"), row.names = FALSE)
write.csv(annual_cell_class_change_summary, file.path(output_dir, "annual_cell_class_change_summary_2000_2025.csv"), row.names = FALSE)
write.csv(cell_class_change_totals, file.path(output_dir, "cell_class_change_totals_2000_2025.csv"), row.names = FALSE)
write.csv(annual_source_sink_area, file.path(output_dir, "annual_source_sink_area_2000_2025.csv"), row.names = FALSE)
write.csv(annual_cell_class_2025, file.path(output_dir, "annual_source_sink_map_cells_2025.csv"), row.names = FALSE)
write.csv(input_notes, file.path(output_dir, "spatial_trial_input_notes.csv"), row.names = FALSE)
# Save monthly_cell_predictions as RDS (large; needed by ERA-SpatialUpscaling-Figures.R)
saveRDS(monthly_cell_predictions,
        file.path(output_dir, "monthly_cell_predictions_2000_2025.rds"))
# Save ERA5 template raster so the figure script can reconstruct maps without reloading ERA5
terra::writeRaster(template,
                   file.path(output_dir, "era5_template.tif"),
                   overwrite = TRUE)
# Save scalar parameters needed by the figure script
write.csv(
  data.frame(
    source_probability_threshold            = source_probability_threshold,
    source_probability_threshold_percentile = source_probability_threshold_percentile,
    gmb_soil_sink_tg_ch4_yr                = gmb_soil_sink_tg_ch4_yr
  ),
  file.path(output_dir, "model_parameters.csv"),
  row.names = FALSE
)

capture.output(
  {
    cat("Monthly spatial upscaling — lmer magnitude models\n\n")
    cat("Class-probability model formula:\n")
    print(monthly_formula)
    cat("\nBalanced monthly model summary:\n")
    print(summary(monthly_model))
    cat("\nMagnitude model skill:\n")
    print(magnitude_model_skill)
    cat("\nClass-probability model skill:\n")
    print(class_probability_model_skill)
    cat("\nProbability calibration skill:\n")
    print(probability_calibration_skill)
    cat("\nState-magnitude consistency diagnostic (by EcoType):\n")
    print(state_magnitude_agreement_summary)
    cat("\nFlagged disagreements (n):", nrow(state_magnitude_flagged), "\n")
    cat("\nAnnual expected-flux estimate (lmer magnitude models):\n")
    print(annual_expected_flux_2000_2025)
    cat("\nInput notes:\n")
    print(input_notes)
  },
  file = file.path(output_dir, "spatial_trial_summary.txt")
)

# ─────────────────────────────────────────────────────────────────────────────
# Figures
# All figures are produced by ERA-SpatialUpscaling-Figures.R, which loads the
# saved outputs written above (monthly_cell_predictions_2000_2025.rds,
# era5_template.tif, model_parameters.csv, and the CSV files) and re-generates
# Figures 1–4 without re-fitting any models or re-projecting any rasters.
# Run that script independently to iterate on figures.
#
# For convenience, the figure script is sourced here so a full pipeline run
# still produces figures.  Set SKIP_FIGURES=true to skip this step.
# ─────────────────────────────────────────────────────────────────────────────
if (!identical(tolower(Sys.getenv("SKIP_FIGURES", unset = "false")), "true")) {
  figures_script <- file.path(dirname(sys.frame(1)$ofile), "ERA-SpatialUpscaling-Figures.R")
  if (!file.exists(figures_script)) {
    # Fall back: look for the script relative to the main script's location
    figures_script <- file.path(
      dirname(normalizePath(sys.call()[[2]], mustWork = FALSE)),
      "ERA-SpatialUpscaling-Figures.R"
    )
  }
  if (file.exists(figures_script)) {
    message("Sourcing figure script: ", figures_script)
    source(figures_script, local = FALSE)
  } else {
    message("Figure script not found; skipping. Run ERA-SpatialUpscaling-Figures.R manually.")
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Legacy palette / theme block kept here only for in-session use
# (e.g. quick diagnostic plots during a live run).
# ─────────────────────────────────────────────────────────────────────────────
fig_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor    = element_blank(),
    plot.title          = element_text(face = "bold", size = 12),
    plot.tag            = element_text(face = "bold", size = 12),
    legend.position     = "bottom",
    legend.title        = element_text(size = 11),
    legend.text         = element_text(size = 10),
    strip.background    = element_rect(fill = "grey92"),
    strip.text          = element_text(size = 11),
    axis.title          = element_text(size = 11),
    axis.text           = element_text(size = 10)
  )

map_theme <- fig_theme +
  theme(panel.grid = element_blank(), legend.key.width = unit(1.2, "cm"))

# Colour helpers
sink_col   <- "#2166AC"
source_col <- "#B2182B"
flux_scale <- function(lim, name = "g C m⁻² yr⁻¹") {
  scale_fill_gradient2(
    low = sink_col, mid = "#D4C08A", high = source_col,
    midpoint = 0, limits = c(-lim, lim),
    breaks = pretty(c(-lim, lim), n = 5),
    name = name, na.value = "transparent"
  )
}
prob_scale <- scale_fill_gradient2(
  low = sink_col, mid = "#D4C08A", high = source_col,
  midpoint = 0.5, limits = c(0, 1),
  breaks = c(0, 0.25, 0.5, 0.75, 1),
  name = "P(source)", na.value = "transparent"
)

make_flux_map <- function(data, flux_col, title_label, tag_label,
                          flux_abs_max = flux_map_abs_max) {
  data %>%
    ggplot(aes(x = x, y = y, fill = .data[[flux_col]])) +
    geom_tile(width = res(template)[1], height = res(template)[2]) +
    coord_equal(expand = FALSE) +
    flux_scale(flux_abs_max) +
    labs(title = title_label, x = "Longitude", y = "Latitude", tag = tag_label) +
    map_theme
}

# ── Shared data prep ──────────────────────────────────────────────────────────
mcp <- monthly_cell_predictions %>%
  mutate(EcoType_plot = if_else(is_arid == 1L, "Arid", as.character(EcoType)))

# Annual source/sink area (binary classification)
plot_annual_source_sink_area <- annual_source_sink_area %>%
  ggplot(aes(x = Year, y = mean_monthly_area_mha, color = selected_exchange_class)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.6) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_color_manual(values = c("Weak sink" = sink_col, "Weak source" = source_col)) +
  labs(x = "Year", y = "Area (Mha)", color = NULL) +
  fig_theme

# ─────────────────────────────────────────────────────────────────────────────
# Shared data prep for figures
# ─────────────────────────────────────────────────────────────────────────────

# Annual classification raster (binary, 2025)
annual_class_raster <- template
names(annual_class_raster) <- "annual_exchange_class"
values(annual_class_raster) <- NA_real_
annual_class_raster[annual_cell_class_2025$cell] <- if_else(
  annual_cell_class_2025$annual_exchange_class == "Weak source", 2, 1
)
annual_class_map_df <- as.data.frame(annual_class_raster, xy = TRUE, na.rm = TRUE) %>%
  mutate(annual_exchange_class = factor(annual_exchange_class,
    levels = c(1, 2), labels = c("Weak sink", "Weak source")))

# Cell-level fluxes 2025 (binary + continuous)
cell_flux_2025 <- monthly_cell_predictions %>%
  filter(Year == 2025) %>%
  mutate(flux_cont = source_probability * predicted_source_flux_gC_m2_month +
           (1 - source_probability) * predicted_sink_flux_gC_m2_month) %>%
  group_by(cell, x, y) %>%
  summarise(
    annual_binary_gC_m2_yr     = sum(expected_flux_gC_m2_month, na.rm = TRUE),
    annual_continuous_gC_m2_yr = sum(flux_cont,                 na.rm = TRUE),
    .groups = "drop"
  )
flux_map_abs_max <- max(abs(range(
  c(cell_flux_2025$annual_binary_gC_m2_yr, cell_flux_2025$annual_continuous_gC_m2_yr),
  na.rm = TRUE)))

# Cell-level flux at high threshold (0.95)
high_threshold <- 0.95
cell_flux_2025_ht <- monthly_cell_predictions %>%
  filter(Year == 2025) %>%
  mutate(flux_ht = if_else(source_probability >= high_threshold,
                            predicted_source_flux_gC_m2_month,
                            predicted_sink_flux_gC_m2_month)) %>%
  group_by(cell, x, y) %>%
  summarise(annual_ht_gC_m2_yr = sum(flux_ht, na.rm = TRUE), .groups = "drop")
flux_ht_abs_max <- max(abs(range(cell_flux_2025_ht$annual_ht_gC_m2_yr, na.rm = TRUE)))

# Seasonal P(source) 2025
prob_seasonal_2025 <- monthly_cell_predictions %>%
  filter(Year == 2025) %>%
  mutate(season = case_when(
    month %in% c(12, 1, 2) ~ "DJF",
    month %in% c(6, 7, 8)  ~ "JJA",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(season)) %>%
  group_by(cell, x, y, season) %>%
  summarise(mean_source_probability = mean(source_probability, na.rm = TRUE), .groups = "drop")

# Threshold sensitivity
threshold_sensitivity <- map_dfr(source_probability_thresholds, function(thresh) {
  monthly_cell_predictions %>%
    mutate(
      flux_t  = if_else(source_probability >= thresh,
                         predicted_source_flux_gC_m2_month, predicted_sink_flux_gC_m2_month),
      tg_t    = gC_m2_yr_to_tg_ch4(flux_t, area_mha)
    ) %>%
    group_by(Year) %>%
    summarise(annual_tg = sum(tg_t, na.rm = TRUE), .groups = "drop") %>%
    mutate(threshold = thresh)
})
thresh_summary <- threshold_sensitivity %>%
  group_by(threshold) %>%
  summarise(mean_tg = mean(annual_tg), min_tg = min(annual_tg), max_tg = max(annual_tg),
            .groups = "drop")

# Budget by EcoType
budget_by_ecotype_year <- mcp %>%
  mutate(tg_cell = gC_m2_yr_to_tg_ch4(expected_flux_gC_m2_month, area_mha)) %>%
  group_by(Year, EcoType_plot) %>%
  summarise(annual_tg = sum(tg_cell, na.rm = TRUE), .groups = "drop")

# Seasonal flux by EcoType
flux_seasonal <- mcp %>%
  mutate(
    flux_bin  = expected_flux_gC_m2_month,
    flux_cont = source_probability * predicted_source_flux_gC_m2_month +
      (1 - source_probability) * predicted_sink_flux_gC_m2_month
  ) %>%
  group_by(month, EcoType_plot) %>%
  summarise(
    Binary     = weighted.mean(flux_bin,  area_mha, na.rm = TRUE),
    Continuous = weighted.mean(flux_cont, area_mha, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(Binary, Continuous), names_to = "approach", values_to = "flux")

# P(source) time series by EcoType
prob_ecotype_ts <- mcp %>%
  group_by(Year, EcoType_plot) %>%
  summarise(mean_prob = weighted.mean(source_probability, area_mha, na.rm = TRUE),
            .groups = "drop")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1  Source-probability model  (6.5 × 8.0 in)
# A: P(source) distribution by EcoType
# B: Annual mean P(source) map, 2025
# C: Annual source/sink area time series, 2000-2025
# ─────────────────────────────────────────────────────────────────────────────

p1a <- mcp %>%
  group_by(cell, EcoType_plot) %>%
  summarise(mean_prob = mean(source_probability, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = mean_prob, fill = EcoType_plot, color = EcoType_plot)) +
  geom_density(alpha = 0.28, linewidth = 0.7) +
  geom_vline(xintercept = source_probability_threshold, color = "grey30",
             linetype = "dotted", linewidth = 0.8) +
  scale_fill_manual(values  = ecotype_colors) +
  scale_color_manual(values = ecotype_colors) +
  labs(title = "A. P(source) by Ecosystem Type",
       x = "Mean P(source)", y = "Density",
       fill = NULL, color = NULL) +
  fig_theme +
  theme(legend.position      = c(0.18, 0.75),
        legend.background    = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size      = unit(0.35, "cm"),
        legend.text          = element_text(size = 6.5))

p1b <- annual_cell_class_2025 %>%
  ggplot(aes(x = x, y = y, fill = mean_source_probability)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  coord_equal(expand = FALSE) +
  prob_scale +
  labs(title = "B. Continuous P(source), 2025", x = "Longitude", y = "Latitude") +
  map_theme +
  theme(legend.position   = c(0.08, 0.28),
        legend.direction  = "vertical",
        legend.key.width  = unit(0.4, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.text       = element_text(size = 6.5),
        legend.title      = element_text(size = 7),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA))

p1c <- annual_class_map_df %>%
  ggplot(aes(x = x, y = y, fill = annual_exchange_class)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  coord_equal(expand = FALSE) +
  scale_fill_manual(values = c("Weak sink" = sink_col, "Weak source" = source_col),
                    labels = c("Weak sink" = "Weak-sink", "Weak source" = "Weak-source"),
                    na.value = "transparent") +
  labs(title = "C. Binary P(source), 2025", x = "Longitude", y = "Latitude",
       fill = NULL) +
  map_theme +
  theme(plot.margin      = margin(0, 0, 0, 0),
        legend.position  = "top",
        legend.direction = "horizontal",
        legend.key.size  = unit(0.4, "cm"),
        legend.text      = element_text(size = 7))

p1d <- annual_source_sink_area %>%
  ggplot(aes(x = Year, y = mean_monthly_area_mha / 100, color = selected_exchange_class)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_color_manual(values = c("Weak sink" = sink_col, "Weak source" = source_col),
                     labels = c("Weak sink" = "Weak-sink", "Weak source" = "Weak-source")) +
  labs(title = "D. Annual Source/Sink Area, 2000-2025",
       x = "Year", y = "Area (100 Mha)", color = NULL) +
  fig_theme

# Embed D (time series) as inset in bottom-left of C (exchange class map)
p1c_inset <- ggdraw() +
  draw_plot(p1c) +
  draw_plot(p1d + theme(plot.title = element_blank(),
                        axis.title = element_text(size = 7),
                        axis.text  = element_text(size = 6),
                        legend.position = "none",
                        plot.background = element_rect(fill = "transparent", color = NA),
                        panel.background = element_rect(fill = "transparent")),
            x = 0.13, y = 0.08, width = 0.234, height = 0.263) +
  draw_label("D.", x = 0.13, y = 0.345, hjust = 0, vjust = 0, size = 8, fontface = "bold")

fig1 <- plot_grid(
  plot_grid(p1a, p1b, ncol = 2, rel_widths = c(0.7, 1.3)),
  p1c_inset,
  ncol = 1, rel_heights = c(1, 1.1)
)

ggsave(file.path(figure_dir, "Fig1_source_probability.png"),
  fig1, width = 8.0, height = 7.0, units = "in", dpi = 300, bg = "white")


# ─────────────────────────────────────────────────────────────────────────────
# Supplement Fig S1  Probability model performance  (6.5 × 7.0 in, 2×2)
# A: Calibration curve
# B: Classification performance metrics bar chart
# C: Annual mean P(source) map, 2025
# D: P(source) distribution by EcoType
# ─────────────────────────────────────────────────────────────────────────────

ps <- class_probability_model_skill

pS1a <- probability_calibration_skill %>%
  ggplot(aes(x = mean_calibrated_source_probability, y = observed_source_fraction)) +
  geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = "dashed",
              linewidth = 0.5) +
  geom_point(aes(size = n), color = source_col, alpha = 0.85) +
  geom_line(color = source_col, linewidth = 0.8) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_size_continuous(range = c(1.5, 5)) +
  labs(title = "A. Probability Calibration",
       x = "Calibrated P(source)", y = "Observed fraction", size = "n months") +
  fig_theme

pS1b <- tibble(
  metric = c("AUC", "Tjur R²", "Dev. explained", "Accuracy",
             "Sensitivity", "Specificity"),
  value  = c(ps$auc, ps$tjur_r2, ps$deviance_explained,
             ps$accuracy, ps$sensitivity_source, ps$specificity_sink),
  group  = c("Discrimination", "Discrimination", "Discrimination",
             "Classification", "Classification", "Classification")
) %>%
  mutate(metric = factor(metric,
    levels = c("AUC", "Tjur R²", "Dev. explained",
               "Accuracy", "Sensitivity", "Specificity"))) %>%
  ggplot(aes(x = value, y = metric, fill = group)) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f", value)), hjust = -0.1, size = 3.2) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_x_continuous(limits = c(0, 1.12), expand = c(0, 0)) +
  scale_fill_manual(values = c("Discrimination" = "#7570B3",
                                "Classification" = "#1B9E77")) +
  labs(title = sprintf("B. Model Performance  (threshold = %.2f)", ps$threshold),
       x = "Value", y = NULL, fill = NULL) +
  fig_theme + theme(legend.position = "bottom")

pS1c <- annual_cell_class_2025 %>%
  ggplot(aes(x = x, y = y, fill = mean_source_probability)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  coord_equal(expand = FALSE) +
  prob_scale +
  labs(title = "C. Annual Mean P(source), 2025", x = "Longitude", y = "Latitude") +
  map_theme

pS1d <- mcp %>%
  group_by(cell, EcoType_plot) %>%
  summarise(mean_prob = mean(source_probability, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = mean_prob, fill = EcoType_plot, color = EcoType_plot)) +
  geom_density(alpha = 0.28, linewidth = 0.7) +
  geom_vline(xintercept = source_probability_threshold, color = "grey30",
             linetype = "dotted", linewidth = 0.8) +
  annotate("text", x = min(source_probability_threshold + 0.06, 0.90), y = Inf,
           label = sprintf("threshold = %.2f", source_probability_threshold),
           hjust = 0, vjust = 1.4, size = 3.0, color = "grey30") +
  scale_fill_manual(values  = ecotype_colors) +
  scale_color_manual(values = ecotype_colors) +
  labs(title = "D. P(source) Distribution by Ecosystem Type",
       x = "Mean P(source)", y = "Density", fill = NULL, color = NULL) +
  fig_theme

figS1 <- plot_grid(
  plot_grid(pS1a, pS1b, ncol = 2),
  plot_grid(pS1c, pS1d, ncol = 2),
  ncol = 1, rel_heights = c(1, 1)
)

ggsave(file.path(figure_dir, "FigS1_probability_model_performance.png"),
  figS1, width = 6.5, height = 7.0, units = "in", dpi = 300, bg = "white")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2  Magnitude models  (6.5 × 8.0 in)
# A: Sink model | B: Source model
# C: Seasonal flux cycle by EcoType
# ─────────────────────────────────────────────────────────────────────────────

make_mag_plot <- function(model_label, panel_tag) {
  n_obs <- sum(magnitude_fitted_values$magnitude_model == model_label)
  magnitude_fitted_values %>%
    filter(magnitude_model == model_label, !is.na(EcoType)) %>%
    ggplot(aes(x = monthly_flux_gC_m2_month, y = fitted_flux_gC_m2_month, color = EcoType)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "grey35", linetype = "dashed", linewidth = 0.5) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = ecotype_colors) +
    labs(title = paste0(panel_tag, ". ", model_label, " model  (n = ", n_obs, ")"),
         x = "Observed (g C m⁻² mo⁻¹)",
         y = "Fitted (g C m⁻² mo⁻¹)",
         color = NULL) +
    fig_theme
}

p2a <- make_mag_plot("Weak sink",   "A")
p2b <- make_mag_plot("Weak source", "B")

p2c <- flux_seasonal %>%
  ggplot(aes(x = month, y = flux, color = EcoType_plot, linetype = approach)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 1.3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = ecotype_colors) +
  scale_linetype_manual(values = c("Binary" = "solid", "Continuous" = "dashed")) +
  facet_wrap(~EcoType_plot, scales = "free_y", ncol = 2) +
  labs(title = "C. Seasonal Flux Cycle by Ecosystem Type",
       x = "Month", y = "Flux (g C m⁻² mo⁻¹)",
       color = NULL, linetype = "Approach") +
  fig_theme +
  guides(color = "none")

# Shared EcoType legend extracted from p2a
shared_ecotype_legend <- get_legend(
  p2a + theme(legend.position = "top",
              legend.key.size = unit(0.4, "cm"),
              legend.text = element_text(size = 8))
)

fig2 <- plot_grid(
  shared_ecotype_legend,
  plot_grid(p2a + theme(legend.position = "none"),
            p2b + theme(legend.position = "none"),
            ncol = 2),
  p2c,
  ncol = 1, rel_heights = c(0.07, 1, 1.4)
)

ggsave(file.path(figure_dir, "Fig2_magnitude_models.png"),
  fig2, width = 8.0, height = 8.0, units = "in", dpi = 300, bg = "white")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 3  Budget  (6.5 × 8.5 in)
# A: Annual net exchange time series, both approaches
# B: Budget by EcoType (stacked)
# C: Threshold sensitivity with data-driven and 0.95 marked
# ─────────────────────────────────────────────────────────────────────────────

p3a <- annual_expected_flux_2000_2025 %>%
  ggplot(aes(x = Year, y = annual_net_exchange_tg_ch4_yr, color = magnitude_model)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -35,
           fill = sink_col, alpha = 0.10) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = sink_col,
             linetype = "dashed", linewidth = 0.7) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_color_manual(values = model_colors,
                     labels = c("Binary: hard-threshold + lmer conditional magnitude" = "Binary",
                                "Continuous: P(source) × source mag + (1-P) × sink mag" = "Continuous")) +
  labs(title = "A. Annual Net Exchange, 2000–2025",
       x = "Year", y = "Net exchange (Tg CH₄ yr⁻¹)", color = NULL) +
  fig_theme

p3a_bar <- annual_expected_flux_2000_2025 %>%
  filter(Year == 2025) %>%
  mutate(approach = c("Binary: hard-threshold + lmer conditional magnitude" = "Binary",
                      "Continuous: P(source) × source mag + (1-P) × sink mag" = "Continuous")[magnitude_model]) %>%
  ggplot(aes(x = approach, y = annual_net_exchange_tg_ch4_yr, fill = magnitude_model)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = sink_col,
             linetype = "dashed", linewidth = 0.7) +
  scale_fill_manual(values = model_colors) +
  labs(title = "2025 Budget", x = NULL, y = NULL) +
  fig_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7))

p3b <- budget_by_ecotype_year %>%
  ggplot(aes(x = Year, y = annual_tg, fill = EcoType_plot)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_area(alpha = 0.82, position = "stack") +
  geom_line(aes(color = EcoType_plot), position = "stack",
            linewidth = 0.35, show.legend = FALSE) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_fill_manual(values  = ecotype_colors) +
  scale_color_manual(values = ecotype_colors) +
  labs(title = "B. Binary Budget by Ecosystem Type",
       x = "Year", y = "Net exchange (Tg CH₄ yr⁻¹)", fill = NULL) +
  fig_theme

p3c <- threshold_sensitivity %>%
  ggplot(aes(x = threshold, y = annual_tg, group = Year)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -35,
           fill = sink_col, alpha = 0.10) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = sink_col,
             linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = source_probability_threshold, color = "grey40",
             linetype = "dotted", linewidth = 0.9) +
  geom_vline(xintercept = high_threshold, color = "grey40",
             linetype = "longdash", linewidth = 0.7) +
  geom_line(alpha = 0.15, linewidth = 0.4, color = as.character(binary_color)) +
  geom_ribbon(data = thresh_summary,
              aes(y = mean_tg, ymin = min_tg, ymax = max_tg, group = 1),
              fill = as.character(binary_color), alpha = 0.17, color = NA) +
  geom_line(data = thresh_summary, aes(y = mean_tg, group = 1),
            color = as.character(binary_color), linewidth = 1.1) +
  geom_point(data = thresh_summary, aes(y = mean_tg, group = 1),
             color = as.character(binary_color), size = 2) +
  annotate("text", x = source_probability_threshold + 0.005, y = Inf,
           label = sprintf("data-driven\n%.2f", source_probability_threshold),
           hjust = 0, vjust = 1.2, size = 3.1, color = "grey30") +
  annotate("text", x = high_threshold + 0.005, y = Inf,
           label = "high\n0.95", hjust = 0, vjust = 1.2, size = 3.1, color = "grey30") +
  scale_x_continuous(breaks = source_probability_thresholds) +
  labs(title = "C. Binary Approach: Threshold Sensitivity",
       x = "P(source) threshold", y = "Net exchange (Tg CH₄ yr⁻¹)") +
  fig_theme

fig3 <- plot_grid(
  plot_grid(p3a, p3a_bar, ncol = 2, rel_widths = c(1.8, 1)),
  p3b,
  p3c,
  ncol = 1, rel_heights = c(1, 1, 1)
)

ggsave(file.path(figure_dir, "Fig3_budget.png"),
  fig3, width = 7.5, height = 8.5, units = "in", dpi = 300, bg = "white")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4  Spatial flux maps  (6.5 × 9.0 in)
# A: Annual exchange class 2025 (binary)
# B: Annual flux, binary approach
# C: Annual flux, continuous approach
# ─────────────────────────────────────────────────────────────────────────────

make_flux_map <- function(data, flux_col, title_label, flux_abs_max = flux_map_abs_max) {
  data %>%
    ggplot(aes(x = x, y = y, fill = .data[[flux_col]])) +
    geom_tile(width = res(template)[1], height = res(template)[2]) +
    coord_equal(expand = FALSE) +
    flux_scale(flux_abs_max) +
    labs(title = title_label, x = "Longitude", y = "Latitude") +
    map_theme
}

p4a <- annual_class_map_df %>%
  ggplot(aes(x = x, y = y, fill = annual_exchange_class)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  coord_equal(expand = FALSE) +
  scale_fill_manual(values = c("Weak sink" = sink_col, "Weak source" = source_col),
                    na.value = "transparent") +
  labs(title = "A. Annual Exchange Class, 2025  (binary)",
       x = "Longitude", y = "Latitude", fill = NULL) +
  map_theme

p4b <- make_flux_map(cell_flux_2025, "annual_binary_gC_m2_yr",
                     "A. Annual Net Flux, 2025  (binary)")

p4c <- make_flux_map(cell_flux_2025, "annual_continuous_gC_m2_yr",
                     "B. Annual Net Flux, 2025  (continuous)")

fig4 <- plot_grid(p4b, p4c, ncol = 1, rel_heights = c(1, 1))

ggsave(file.path(figure_dir, "Fig4_spatial_maps.png"),
  fig4, width = 6.5, height = 6.5, units = "in", dpi = 300, bg = "white")

# Also write individual maps for separate use
ggsave(file.path(figure_dir, "annual_expected_flux_time_series_2000_2025.png"),
  p3a + labs(title = NULL), width = 6.5, height = 3.5, units = "in", dpi = 300)
ggsave(file.path(figure_dir, "ERA5_source_probability_multipanel.png"),
  fig1, width = 6.5, height = 8.5, units = "in", dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "ERA5_magnitude_models_multipanel.png"),
  fig2, width = 6.5, height = 8.0, units = "in", dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "ERA5_budget_comparison_multipanel.png"),
  fig3, width = 6.5, height = 8.5, units = "in", dpi = 300, bg = "white")

message("Wrote monthly spatial upscaling trial outputs to ", spatial_dir)
