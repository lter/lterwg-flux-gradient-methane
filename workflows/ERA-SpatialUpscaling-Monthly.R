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
site_behavior_file <- file.path(localdir.ch4, "OUTPUT/30min_site_behavior.csv")
rate_file <- file.path(localdir.ch4, "OUTPUT/ERA_Upscaling_stage2_flux_rates_rate_scenarios.csv")
soil_chamber_file <- file.path(localdir.ch4, "OUTPUT/soil_chamber_CH4_flux_reference_values.csv")

worldclim_tavg_zip <- file.path(spatial_dir, "wc2.1_10m_tavg.zip")
worldclim_prec_zip <- file.path(spatial_dir, "wc2.1_10m_prec.zip")
ecoregions_zip <- file.path(spatial_dir, "Ecoregions2017.zip")
era5_land_dir <- file.path(data_dir, "era5_land_monthly")
modis_ecotype_dir <- file.path(data_dir, "modis_mcd12c1_processed")
wad2m_dir <- file.path(data_dir, "wad2m")

required_files <- c(era5_30min_file, site_behavior_file, rate_file, ecoregions_zip)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

source_probability_threshold <- 0.80
source_probability_thresholds <- seq(0.50, 0.95, by = 0.05)
inundation_exclusion_threshold <- 0.05
gmb_soil_sink_tg_ch4_yr <- -35
magnitude_constraint_weight <- 0.65
years_to_repeat <- 2000:2025
years_to_process <- years_to_repeat
rate_scenario_colors <- c(
  "Chamber/process sink + NEON source" = "#1B9E77",
  "NEON ERA5 rates" = "#D95F02"
)

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

  lmer_formula <- log_abs_flux ~
    z_Tair + z_VSWC + z_MAP + z_MAT + source_probability +
    (1 + z_Tair + z_VSWC | EcoType) + (1 | SITE_ID)
  lmer_intercept_formula <- log_abs_flux ~
    z_Tair + z_VSWC + z_MAP + z_MAT + source_probability +
    (1 | EcoType) + (1 | SITE_ID)
  lm_formula <- log_abs_flux ~
    z_Tair + z_VSWC + z_MAP + z_MAT + source_probability + EcoType

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
  names(climate) <- c(
    sprintf("t2m_%02d", 1:12),
    sprintf("swvl1_%02d", 1:12),
    sprintf("swvl2_%02d", 1:12),
    sprintf("tp_%02d", 1:12)
  )
  tavg_year <- climate[[sprintf("t2m_%02d", 1:12)]] - 273.15
  vswc_year <- (climate[[sprintf("swvl1_%02d", 1:12)]] + climate[[sprintf("swvl2_%02d", 1:12)]]) / 2
  prec_year <- rast(lapply(1:12, function(month) {
    climate[[sprintf("tp_%02d", month)]] * 1000 * lubridate::days_in_month(lubridate::make_date(year, month, 1))
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
  ecotype_code = c(1, 2, 3, 4),
  EcoType = c("Cropland", "Forest", "Grassland", "Shrubland")
)

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

site_behavior <- read.csv(site_behavior_file) %>%
  mutate(SITE_ID = as.character(SITE_ID))

upland_sites <- site_behavior %>%
  filter(
    !is.na(EcoType),
    !str_detect(EcoType, regex("wetland|inundat|flood|marsh|swamp|bog|fen|lake|rice", ignore_case = TRUE))
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
    EcoType = factor(EcoType, levels = c("Cropland", "Forest", "Grassland", "Shrubland")),
    weak_source = as.integer(monthly_budget_mgC_m2 > 0),
    monthly_flux_gC_m2_month = monthly_budget_mgC_m2 / 1000,
    MAP = as.numeric(MAP),
    MAT = as.numeric(MAT)
  ) %>%
  filter(
    is.finite(weak_source),
    is.finite(mean_ERA5_Tair_C),
    is.finite(mean_ERA5_VSWC),
    is.finite(MAP),
    is.finite(MAT)
  )

monthly_formula <- weak_source ~ EcoType +
  scale(mean_ERA5_Tair_C) + scale(mean_ERA5_VSWC) + scale(MAP) + scale(MAT)

monthly_model <- fit_monthly_balanced_model(monthly_training, monthly_formula)

monthly_training <- monthly_training %>%
  mutate(source_probability_raw = as.numeric(predict(monthly_model, newdata = monthly_training, type = "response")))

probability_calibration_model <- fit_probability_calibration(monthly_training)

monthly_training <- monthly_training %>%
  mutate(
    source_probability = calibrate_source_probability(probability_calibration_model, source_probability_raw)
  )

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

training_vswc_range <- range(monthly_training$mean_ERA5_VSWC, na.rm = TRUE)
training_prec_range <- quantile(monthly_training$MAP / 12, probs = c(0.02, 0.98), na.rm = TRUE)

rate_scenarios <- read.csv(rate_file) %>%
  mutate(
    EcoType = as.character(EcoType),
    exchange_class = as.character(exchange_class),
    calibrated_rate_gC_m2_yr = as.numeric(calibrated_rate_gC_m2_yr),
    calibrated_rate_se_gC_m2_yr = as.numeric(calibrated_rate_se_gC_m2_yr)
  ) %>%
  filter(rate_scenario %in% c("NEON ERA5 rates", "Chamber/process sink + NEON source"))

upland_chamber_source_cap_gC_m2_yr <- if (file.exists(soil_chamber_file)) {
  chamber_refs <- read.csv(soil_chamber_file) %>%
    mutate(
      ecosystem_class = as.character(ecosystem_class),
      daily_high_mgC_m2_day = as.numeric(daily_high_mgC_m2_day)
    ) %>%
    filter(
      !str_detect(ecosystem_class, regex("wetland|inundat|flood|marsh|swamp|bog|fen|lake|rice", ignore_case = TRUE)),
      is.finite(daily_high_mgC_m2_day),
      daily_high_mgC_m2_day > 0
    )
  if (nrow(chamber_refs) > 0) {
    max(chamber_refs$daily_high_mgC_m2_day, na.rm = TRUE) * 365 / 1000
  } else {
    NA_real_
  }
} else {
  NA_real_
}

magnitude_constraints <- rate_scenarios %>%
  filter(rate_scenario == "Chamber/process sink + NEON source") %>%
  dplyr::select(EcoType, exchange_class, calibrated_rate_gC_m2_yr) %>%
  pivot_wider(names_from = exchange_class, values_from = calibrated_rate_gC_m2_yr) %>%
  transmute(
    EcoType = factor(EcoType, levels = levels(monthly_training$EcoType)),
    chamber_process_sink_floor_gC_m2_month = `Weak sink` / 12,
    upland_neon_source_cap_gC_m2_month = `Weak source` / 12,
    upland_chamber_source_cap_gC_m2_month = upland_chamber_source_cap_gC_m2_yr / 12
  )

monthly_components <- list()
monthly_cell_predictions <- list()
monthly_threshold_components <- list()
monthly_expected_components <- list()
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
      month_stack <- c(month_stack, climate_year$vswc[[m]])
      names(month_stack)[nlyr(month_stack)] <- "vswc"
    }

    dat <- as.data.frame(month_stack, xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
      filter(!is.na(ecotype_code), !is.na(inundation_fraction), !is.na(tavg), !is.na(prec), !is.na(MAT), !is.na(MAP), !is.na(area_mha)) %>%
      mutate(ecotype_code = as.integer(ecotype_code)) %>%
      inner_join(ecotype_lookup, by = "ecotype_code") %>%
      mutate(
        Year = year,
        month = m,
        inundation_fraction = pmin(pmax(inundation_fraction, 0), 1),
        upland_area_fraction = if_else(inundation_fraction <= inundation_exclusion_threshold, 1 - inundation_fraction, 0),
        area_mha = area_mha * upland_area_fraction,
        EcoType = factor(EcoType, levels = levels(monthly_training$EcoType)),
        mean_ERA5_Tair_C = tavg,
        mean_ERA5_VSWC = if (use_era5_land) {
          vswc
        } else {
          scales::rescale(
            pmin(pmax(prec, training_prec_range[1]), training_prec_range[2]),
            to = training_vswc_range,
            from = training_prec_range
          )
        }
      ) %>%
      filter(area_mha > 0)

    dat <- apply_standardizers(dat, magnitude_standardizers) %>%
      mutate(SITE_ID = factor("SPATIAL"))

    dat$source_probability_raw <- as.numeric(predict(monthly_model, newdata = dat, type = "response"))
    dat$source_probability <- calibrate_source_probability(
      probability_calibration_model,
      dat$source_probability_raw
    )
    dat$predicted_sink_flux_gC_m2_month <- pmin(
      predict_flux_magnitude(sink_magnitude_model, dat),
      0
    )
    dat$predicted_source_flux_gC_m2_month <- pmax(
      predict_flux_magnitude(source_magnitude_model, dat),
      0
    )

    dat <- dat %>%
      left_join(magnitude_constraints, by = "EcoType") %>%
      mutate(
        constrained_sink_flux_gC_m2_month = pmin(
          0,
          (1 - magnitude_constraint_weight) * predicted_sink_flux_gC_m2_month +
            magnitude_constraint_weight * chamber_process_sink_floor_gC_m2_month
        ),
        constrained_source_flux_gC_m2_month = pmax(
          0,
          if_else(
            predicted_source_flux_gC_m2_month > coalesce(upland_chamber_source_cap_gC_m2_month, upland_neon_source_cap_gC_m2_month),
            (1 - magnitude_constraint_weight) * predicted_source_flux_gC_m2_month +
              magnitude_constraint_weight * coalesce(upland_chamber_source_cap_gC_m2_month, upland_neon_source_cap_gC_m2_month),
            predicted_source_flux_gC_m2_month
          )
        ),
        expected_flux_gC_m2_month = source_probability * predicted_source_flux_gC_m2_month +
          (1 - source_probability) * predicted_sink_flux_gC_m2_month,
        constrained_expected_flux_gC_m2_month = source_probability * constrained_source_flux_gC_m2_month +
          (1 - source_probability) * constrained_sink_flux_gC_m2_month
      )

    monthly_threshold_components[[component_index]] <- purrr::map_dfr(
      source_probability_thresholds,
      function(threshold) {
        dat %>%
          mutate(
            source_probability_threshold = threshold,
            selected_exchange_class = if_else(
              source_probability >= threshold,
              "Weak source",
              "Weak sink"
            )
          ) %>%
          group_by(Year, source_probability_threshold, month, EcoType, selected_exchange_class) %>%
          summarise(area_mha = sum(area_mha, na.rm = TRUE), .groups = "drop")
      }
    )

    dat <- dat %>%
      mutate(
        selected_exchange_class = if_else(
          source_probability >= source_probability_threshold,
          "Weak source",
          "Weak sink"
        )
      )

    monthly_cell_predictions[[component_index]] <- dat %>%
      dplyr::select(
        Year,
        cell,
        x,
        y,
        month,
        EcoType,
        area_mha,
        inundation_fraction,
        upland_area_fraction,
        source_probability_raw,
        source_probability,
        predicted_sink_flux_gC_m2_month,
        predicted_source_flux_gC_m2_month,
        expected_flux_gC_m2_month,
        constrained_sink_flux_gC_m2_month,
        constrained_source_flux_gC_m2_month,
        constrained_expected_flux_gC_m2_month,
        selected_exchange_class
      )

    monthly_expected_components[[component_index]] <- dat %>%
      dplyr::select(
        Year,
        month,
        EcoType,
        area_mha,
        source_probability_raw,
        source_probability,
        predicted_sink_flux_gC_m2_month,
        predicted_source_flux_gC_m2_month,
        expected_flux_gC_m2_month,
        constrained_sink_flux_gC_m2_month,
        constrained_source_flux_gC_m2_month,
        constrained_expected_flux_gC_m2_month,
        tavg,
        prec,
        mean_ERA5_VSWC,
        inundation_fraction
      ) %>%
      pivot_longer(
        cols = c(expected_flux_gC_m2_month, constrained_expected_flux_gC_m2_month),
        names_to = "magnitude_scenario",
        values_to = "monthly_rate_gC_m2_month"
      ) %>%
      mutate(
        magnitude_scenario = recode(
          magnitude_scenario,
          expected_flux_gC_m2_month = "Continuous expected flux, condition-only magnitude",
          constrained_expected_flux_gC_m2_month = "Continuous expected flux, chamber-constrained magnitude"
        ),
        monthly_tg_ch4 = gC_m2_yr_to_tg_ch4(monthly_rate_gC_m2_month, area_mha)
      ) %>%
      group_by(Year, month, EcoType, magnitude_scenario) %>%
      summarise(
        total_area_mha = sum(area_mha, na.rm = TRUE),
        monthly_tg_ch4 = sum(monthly_tg_ch4, na.rm = TRUE),
        mean_monthly_rate_gC_m2_month = wtd_mean(monthly_rate_gC_m2_month, .data$area_mha),
        mean_source_probability = wtd_mean(source_probability, .data$area_mha),
        mean_sink_flux_gC_m2_month = wtd_mean(predicted_sink_flux_gC_m2_month, .data$area_mha),
        mean_source_flux_gC_m2_month = wtd_mean(predicted_source_flux_gC_m2_month, .data$area_mha),
        mean_constrained_sink_flux_gC_m2_month = wtd_mean(constrained_sink_flux_gC_m2_month, .data$area_mha),
        mean_constrained_source_flux_gC_m2_month = wtd_mean(constrained_source_flux_gC_m2_month, .data$area_mha),
        mean_tavg_C = wtd_mean(tavg, .data$area_mha),
        mean_prec_mm = wtd_mean(prec, .data$area_mha),
        mean_vswc = wtd_mean(mean_ERA5_VSWC, .data$area_mha),
        mean_inundation_fraction = wtd_mean(inundation_fraction, .data$area_mha),
        .groups = "drop"
      ) %>%
      rename(area_mha = total_area_mha)

    monthly_components[[component_index]] <- dat %>%
      group_by(Year, month, EcoType, selected_exchange_class) %>%
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
monthly_threshold_area_by_class <- bind_rows(monthly_threshold_components)
monthly_expected_flux_components <- bind_rows(monthly_expected_components)

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

monthly_upscaled_components <- monthly_area_by_class %>%
  left_join(
    rate_scenarios %>%
      dplyr::select(rate_scenario, EcoType, exchange_class, calibrated_rate_gC_m2_yr, calibrated_rate_se_gC_m2_yr),
    by = c("EcoType", "selected_exchange_class" = "exchange_class"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    monthly_rate_gC_m2_month = calibrated_rate_gC_m2_yr / 12,
    monthly_rate_se_gC_m2_month = calibrated_rate_se_gC_m2_yr / 12,
    monthly_tg_ch4 = gC_m2_yr_to_tg_ch4(monthly_rate_gC_m2_month, area_mha),
    monthly_tg_ch4_se_component = abs(gC_m2_yr_to_tg_ch4(monthly_rate_se_gC_m2_month, area_mha))
  )

threshold_sensitivity_components <- monthly_threshold_area_by_class %>%
  left_join(
    rate_scenarios %>%
      dplyr::select(rate_scenario, EcoType, exchange_class, calibrated_rate_gC_m2_yr, calibrated_rate_se_gC_m2_yr),
    by = c("EcoType", "selected_exchange_class" = "exchange_class"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    monthly_rate_gC_m2_month = calibrated_rate_gC_m2_yr / 12,
    monthly_rate_se_gC_m2_month = calibrated_rate_se_gC_m2_yr / 12,
    monthly_tg_ch4 = gC_m2_yr_to_tg_ch4(monthly_rate_gC_m2_month, area_mha),
    monthly_tg_ch4_se_component = abs(gC_m2_yr_to_tg_ch4(monthly_rate_se_gC_m2_month, area_mha))
  )

threshold_sensitivity <- threshold_sensitivity_components %>%
  group_by(Year, source_probability_threshold, rate_scenario) %>%
  summarise(
    annual_net_exchange_tg_ch4_yr = sum(monthly_tg_ch4, na.rm = TRUE),
    approximate_se_tg_ch4_yr = sqrt(sum(monthly_tg_ch4_se_component^2, na.rm = TRUE)),
    mean_monthly_source_area_mha = sum(area_mha[selected_exchange_class == "Weak source"], na.rm = TRUE) / 12,
    mean_monthly_sink_area_mha = sum(area_mha[selected_exchange_class == "Weak sink"], na.rm = TRUE) / 12,
    .groups = "drop"
  )

annual_expected_flux_2000_2025 <- monthly_expected_flux_components %>%
  group_by(Year, magnitude_scenario) %>%
  summarise(
    n_monthly_ecotype_components = dplyr::n(),
    annual_net_exchange_tg_ch4_yr = sum(monthly_tg_ch4, na.rm = TRUE),
    mean_source_probability = wtd_mean(mean_source_probability, area_mha),
    global_methane_budget_soil_sink_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_global_soil_sink_magnitude = 100 * annual_net_exchange_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr),
    .groups = "drop"
  ) %>%
  mutate(input_climate_note = if_else(
    use_era5_land,
    "ERA5-Land monthly temperature and soil moisture vary by year.",
    "WorldClim monthly climatology repeated for each year; no interannual ERA5 variability included."
  ))

threshold_sensitivity_summary <- threshold_sensitivity %>%
  group_by(source_probability_threshold, rate_scenario) %>%
  summarise(
    mean_annual_net_exchange_tg_ch4_yr = mean(annual_net_exchange_tg_ch4_yr, na.rm = TRUE),
    min_annual_net_exchange_tg_ch4_yr = min(annual_net_exchange_tg_ch4_yr, na.rm = TRUE),
    max_annual_net_exchange_tg_ch4_yr = max(annual_net_exchange_tg_ch4_yr, na.rm = TRUE),
    mean_monthly_source_area_mha = mean(mean_monthly_source_area_mha, na.rm = TRUE),
    mean_monthly_sink_area_mha = mean(mean_monthly_sink_area_mha, na.rm = TRUE),
    .groups = "drop"
  )

climatological_annual_estimate <- monthly_upscaled_components %>%
  group_by(Year, rate_scenario) %>%
  summarise(
    source_probability_threshold = source_probability_threshold,
    n_monthly_ecotype_components = dplyr::n(),
    annual_net_exchange_tg_ch4_yr = sum(monthly_tg_ch4, na.rm = TRUE),
    approximate_se_tg_ch4_yr = sqrt(sum(monthly_tg_ch4_se_component^2, na.rm = TRUE)),
    global_methane_budget_soil_sink_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_global_soil_sink_magnitude = 100 * annual_net_exchange_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr),
    .groups = "drop"
  )

annual_estimates_2000_2025 <- climatological_annual_estimate %>%
  mutate(input_climate_note = if_else(
    use_era5_land,
    "ERA5-Land monthly temperature and soil moisture vary by year.",
    "WorldClim monthly climatology repeated for each year; no interannual ERA5 variability included."
  ))

input_notes <- tribble(
  ~item, ~note,
  "Spatial data directory", spatial_dir,
  "Temperature input", if_else(use_era5_land, "ERA5-Land monthly 2 m temperature, converted from K to degrees C.", "WorldClim v2.1 10-minute monthly average temperature climatology."),
  "Precipitation input", if_else(use_era5_land, "ERA5-Land total precipitation converted to monthly mm and used for long-term MAP.", "WorldClim v2.1 10-minute monthly precipitation climatology."),
  "Land-cover/ecosystem input", if_else(use_modis_landcover, "Annual MODIS MCD12C1 processed ecosystem rasters; nearest available year used for missing endpoints.", "Ecoregions2017 terrestrial biome polygons rasterized to the climate grid."),
  "Wetland/inundation exclusion", if_else(use_wad2m_inundation, paste0("Monthly WAD2M inundation fraction used to remove cells above ", inundation_exclusion_threshold, " inundated fraction and down-weight remaining upland area by 1 - inundation_fraction."), "Flooded grasslands/savannas and mangroves excluded from the ecoregion biome map; no dynamic inundation fraction file was available."),
  "Cropland limitation", if_else(use_modis_landcover, "MODIS cropland and cropland/natural vegetation mosaic classes are mapped to Cropland.", "Ecoregions2017 does not map croplands explicitly, so this run includes forest, grassland, and shrubland/desert/tundra-like uplands only."),
  "Soil moisture input", if_else(use_era5_land, "ERA5-Land volumetric soil water layers 1 and 2 averaged as the monthly VSWC predictor.", "ERA5/ERA5-Land monthly soil moisture was not available locally; monthly precipitation was rescaled to the NEON training VSWC range as a temporary proxy."),
  "Continuous magnitude model", "Class probability is estimated with the balanced monthly probability model using ecosystem type, ERA5-Land temperature, ERA5-Land soil moisture, MAP, and MAT, then calibrated before use in the expected-flux equation. Separate weak-sink and weak-source log-absolute-magnitude models predict signed flux from ERA5-Land temperature, ERA5-Land soil moisture, MAP, MAT, calibrated source probability, EcoType partial pooling, and site random effects during training. Season is intentionally excluded from both class-probability and magnitude models.",
  "Magnitude constraints", paste0("The constrained continuous scenario applies soft shrinkage, not hard clipping, toward chamber/process sink rates and the maximum positive non-wetland upland chamber source bound available in the local chamber reference table. Constraint weight = ", magnitude_constraint_weight, ". FLUXNET-CH4 source classes are not used for source bounds."),
  "Interannual treatment", if_else(use_era5_land, "Monthly climate predictors vary by year from 2000-2025.", "The same monthly climatology is repeated for 2000-2025; outputs are not year-specific ERA5 estimates."),
  "Interpretation", if_else(use_era5_land, "ERA5-Land spatial upscaling trial; land-cover and inundation are dynamic only when processed MODIS/WAD2M files are present.", "Proof-of-concept spatial mechanics only; not a final global methane budget estimate.")
)

write.csv(monthly_area_by_class, file.path(output_dir, "monthly_area_by_ecotype_class.csv"), row.names = FALSE)
write.csv(monthly_upscaled_components, file.path(output_dir, "monthly_upscaled_components.csv"), row.names = FALSE)
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
write.csv(threshold_sensitivity, file.path(output_dir, "spatial_threshold_sensitivity_by_year.csv"), row.names = FALSE)
write.csv(threshold_sensitivity_summary, file.path(output_dir, "spatial_threshold_sensitivity.csv"), row.names = FALSE)
write.csv(climatological_annual_estimate, file.path(output_dir, "annual_estimate_by_year.csv"), row.names = FALSE)
write.csv(annual_estimates_2000_2025, file.path(output_dir, "annual_estimates_2000_2025.csv"), row.names = FALSE)
write.csv(input_notes, file.path(output_dir, "spatial_trial_input_notes.csv"), row.names = FALSE)

capture.output(
  {
    cat("Monthly spatial upscaling trial\n\n")
    cat("Model formula:\n")
    print(monthly_formula)
    cat("\nBalanced monthly model summary:\n")
    print(summary(monthly_model))
    cat("\nMagnitude model skill:\n")
    print(magnitude_model_skill)
    cat("\nClass-probability model skill:\n")
    print(class_probability_model_skill)
    cat("\nProbability calibration skill:\n")
    print(probability_calibration_skill)
    cat("\nContinuous expected annual flux estimate:\n")
    print(annual_expected_flux_2000_2025)
    cat("\nClimatological annual estimate:\n")
    print(climatological_annual_estimate)
    cat("\nInput notes:\n")
    print(input_notes)
  },
  file = file.path(output_dir, "spatial_trial_summary.txt")
)

plot_monthly_components <- monthly_upscaled_components %>%
  group_by(rate_scenario, month, EcoType, selected_exchange_class) %>%
  summarise(monthly_tg_ch4 = mean(tapply(monthly_tg_ch4, Year, sum), na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = month, y = monthly_tg_ch4, fill = selected_exchange_class)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_col(width = 0.82) +
  facet_grid(rate_scenario ~ EcoType, scales = "free_y") +
  scale_fill_manual(values = c("Weak sink" = "#2166AC", "Weak source" = "#B2182B")) +
  scale_x_continuous(breaks = 1:12) +
  labs(
    title = "Monthly Spatial Upscaling Trial Components",
    subtitle = if_else(use_era5_land, "Mean monthly component across 2000-2025 using ERA5-Land climate.", "WorldClim/ecoregion climatological proof-of-concept; not final ERA5 2000-2025 spatial estimate."),
    x = "Month",
    y = "Monthly exchange (Tg CH4)",
    fill = "Assigned class"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  file.path(figure_dir, "monthly_spatial_components.png"),
  plot_monthly_components,
  width = 12,
  height = 7,
  units = "in",
  dpi = 300
)

plot_annual_estimate <- climatological_annual_estimate %>%
  group_by(rate_scenario) %>%
  summarise(
    annual_net_exchange_tg_ch4_yr = mean(annual_net_exchange_tg_ch4_yr, na.rm = TRUE),
    approximate_se_tg_ch4_yr = mean(approximate_se_tg_ch4_yr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = rate_scenario, y = annual_net_exchange_tg_ch4_yr, fill = rate_scenario)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -35, fill = "#2166AC", alpha = 0.12) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_errorbar(
    aes(
      ymin = annual_net_exchange_tg_ch4_yr - approximate_se_tg_ch4_yr,
      ymax = annual_net_exchange_tg_ch4_yr + approximate_se_tg_ch4_yr
    ),
    width = 0.12
  ) +
  coord_flip() +
  scale_fill_manual(values = rate_scenario_colors) +
  labs(
    title = "Mean Annual Spatial Upscaling Estimate",
    subtitle = if_else(use_era5_land, "Mean 2000-2025; dashed blue line = GMB soil sink.", "Climatological spatial proxy; dashed blue line = GMB soil sink."),
    x = NULL,
    y = "Net exchange (Tg CH4 yr-1)",
    fill = "Rate scenario"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(4, 8, 4, 4)
  )

ggsave(
  file.path(figure_dir, "annual_spatial_trial_estimate.png"),
  plot_annual_estimate,
  width = 9,
  height = 4.8,
  units = "in",
  dpi = 300
)

plot_annual_time_series <- annual_estimates_2000_2025 %>%
  mutate(
    ymin = annual_net_exchange_tg_ch4_yr - approximate_se_tg_ch4_yr,
    ymax = annual_net_exchange_tg_ch4_yr + approximate_se_tg_ch4_yr
  ) %>%
  ggplot(aes(x = Year, y = annual_net_exchange_tg_ch4_yr, color = rate_scenario, fill = rate_scenario)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -35, fill = "#2166AC", alpha = 0.12) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.16, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.7) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5), minor_breaks = 2000:2025) +
  scale_color_manual(values = rate_scenario_colors) +
  scale_fill_manual(values = rate_scenario_colors) +
  labs(
    title = "Annual Spatial Upscaling Trial, 2000-2025",
    subtitle = if_else(use_era5_land, "ERA5-Land temperature and soil moisture drive annual variation.", "Repeated monthly climatology; annual inputs not available."),
    x = "Year",
    y = "Net exchange (Tg CH4 yr-1)",
    color = "Rate scenario",
    fill = "Rate scenario"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5),
    plot.margin = margin(4, 4, 4, 8),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "annual_spatial_trial_time_series_2000_2025.png"),
  plot_annual_time_series,
  width = 9.5,
  height = 5.2,
  units = "in",
  dpi = 300
)

plot_continuous_expected_flux <- annual_expected_flux_2000_2025 %>%
  mutate(
    magnitude_scenario_label = recode(
      magnitude_scenario,
      "Continuous expected flux, condition-only magnitude" = "Condition-only magnitude",
      "Continuous expected flux, chamber-constrained magnitude" = "Chamber-constrained magnitude"
    )
  ) %>%
  ggplot(aes(x = Year, y = annual_net_exchange_tg_ch4_yr, color = magnitude_scenario_label)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -35, fill = "#2166AC", alpha = 0.12) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.7) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5), minor_breaks = 2000:2025) +
  scale_color_manual(values = c(
    "Condition-only magnitude" = "#7570B3",
    "Chamber-constrained magnitude" = "#1B9E77"
  )) +
  labs(
    title = "Continuous Expected-Flux Spatial Upscaling, 2000-2025",
    subtitle = "P(source) weights source and sink magnitudes; no hard threshold.",
    x = "Year",
    y = "Net exchange (Tg CH4 yr-1)",
    color = "Magnitude scenario"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5),
    plot.margin = margin(4, 8, 4, 4),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "annual_expected_flux_time_series_2000_2025.png"),
  plot_continuous_expected_flux,
  width = 9.5,
  height = 5.2,
  units = "in",
  dpi = 300
)

plot_magnitude_model_skill <- magnitude_fitted_values %>%
  ggplot(aes(x = monthly_flux_gC_m2_month, y = fitted_flux_gC_m2_month, color = EcoType)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, color = "grey35", linetype = "dashed", linewidth = 0.5) +
  geom_point(alpha = 0.45, size = 1.7) +
  facet_wrap(~ magnitude_model, scales = "free") +
  labs(
    title = "Conditional Magnitude Model Fit",
    subtitle = "Separate monthly source and sink models with EcoType partial pooling.",
    x = "Observed monthly flux (g C m-2 month-1)",
    y = "Fitted monthly flux (g C m-2 month-1)",
    color = "Ecosystem type"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "magnitude_model_observed_vs_fitted.png"),
  plot_magnitude_model_skill,
  width = 9.5,
  height = 5.2,
  units = "in",
  dpi = 300
)

plot_probability_calibration <- probability_calibration_skill %>%
  ggplot(aes(x = mean_calibrated_source_probability, y = observed_source_fraction)) +
  geom_abline(slope = 1, intercept = 0, color = "grey35", linetype = "dashed", linewidth = 0.6) +
  geom_point(aes(size = n), color = "#B2182B", alpha = 0.85) +
  geom_line(color = "#B2182B", linewidth = 0.8) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "Monthly Source-Probability Calibration",
    subtitle = "Binned check of calibrated balanced probabilities.",
    x = "Mean calibrated source probability",
    y = "Observed source fraction",
    size = "Months"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5),
    plot.margin = margin(4, 4, 4, 8),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "monthly_source_probability_calibration.png"),
  plot_probability_calibration,
  width = 6.5,
  height = 5.4,
  units = "in",
  dpi = 300
)

plot_continuous_annual_budget_comparison <- annual_expected_flux_2000_2025 %>%
  mutate(
    magnitude_scenario_label = recode(
      magnitude_scenario,
      "Continuous expected flux, condition-only magnitude" = "Condition-only magnitude",
      "Continuous expected flux, chamber-constrained magnitude" = "Chamber-constrained magnitude"
    )
  ) %>%
  group_by(magnitude_scenario_label) %>%
  summarise(
    mean_annual_net_exchange_tg_ch4_yr = mean(annual_net_exchange_tg_ch4_yr, na.rm = TRUE),
    min_annual_net_exchange_tg_ch4_yr = min(annual_net_exchange_tg_ch4_yr, na.rm = TRUE),
    max_annual_net_exchange_tg_ch4_yr = max(annual_net_exchange_tg_ch4_yr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(
    x = magnitude_scenario_label,
    y = mean_annual_net_exchange_tg_ch4_yr,
    fill = magnitude_scenario_label
  )) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -35, fill = "#2166AC", alpha = 0.12) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_col(width = 0.62) +
  geom_errorbar(
    aes(
      ymin = min_annual_net_exchange_tg_ch4_yr,
      ymax = max_annual_net_exchange_tg_ch4_yr
    ),
    width = 0.12
  ) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Condition-only magnitude" = "#7570B3",
    "Chamber-constrained magnitude" = "#1B9E77"
  )) +
  labs(
    title = "Annual Budget Comparison",
    subtitle = "Means and annual min-max, 2000-2025.",
    x = NULL,
    y = "Net exchange (Tg CH4 yr-1)",
    fill = "Magnitude scenario"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

plot_annual_source_sink_area <- annual_source_sink_area %>%
  ggplot(aes(x = Year, y = mean_monthly_area_mha, color = selected_exchange_class)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.7) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5), minor_breaks = 2000:2025) +
  scale_color_manual(values = c("Weak sink" = "#2166AC", "Weak source" = "#B2182B")) +
  labs(
    title = "Annual Source/Sink Area, 2000-2025",
    subtitle = if_else(use_era5_land, "Area is mean monthly classified area from year-specific ERA5-Land climate.", "Area is mean monthly classified area; flat lines reflect repeated monthly climatology."),
    x = "Year",
    y = "Mean monthly area (Mha)",
    color = "Assigned class"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "annual_source_sink_area_time_series_2000_2025.png"),
  plot_annual_source_sink_area,
  width = 9.5,
  height = 5.2,
  units = "in",
  dpi = 300
)

annual_class_raster <- template
names(annual_class_raster) <- "annual_exchange_class"
values(annual_class_raster) <- NA_real_
annual_class_raster[annual_cell_class_2025$cell] <- if_else(
  annual_cell_class_2025$annual_exchange_class == "Weak source",
  2,
  1
)

annual_class_map <- as.data.frame(annual_class_raster, xy = TRUE, na.rm = TRUE) %>%
  mutate(
    annual_exchange_class = factor(
      annual_exchange_class,
      levels = c(1, 2),
      labels = c("Weak sink", "Weak source")
    )
  )

plot_annual_source_sink_map <- annual_class_map %>%
  ggplot(aes(x = x, y = y, fill = annual_exchange_class)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  coord_equal(expand = FALSE) +
  scale_fill_manual(values = c("Weak sink" = "#2166AC", "Weak source" = "#B2182B"), na.value = "transparent") +
  labs(
    title = "Mapped Source/Sink Classification, 2025",
    subtitle = if_else(use_era5_land, "Dominant 2025 class; weak source if >=6 months exceed 0.80.", "Dominant repeated-climatology class; weak source if >=6 months exceed 0.80."),
    x = "Longitude",
    y = "Latitude",
    fill = "Annual class"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "source_sink_map_2025.png"),
  plot_annual_source_sink_map,
  width = 11,
  height = 5.8,
  units = "in",
  dpi = 300
)

plot_threshold_sensitivity <- threshold_sensitivity %>%
  ggplot(aes(x = source_probability_threshold, y = annual_net_exchange_tg_ch4_yr, color = rate_scenario, fill = rate_scenario, group = interaction(rate_scenario, Year))) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -35, fill = "#2166AC", alpha = 0.12) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = source_probability_threshold, color = "grey40", linetype = "dotted", linewidth = 0.6) +
  geom_line(alpha = 0.18, linewidth = 0.45) +
  geom_ribbon(
    data = threshold_sensitivity_summary,
    aes(
      y = mean_annual_net_exchange_tg_ch4_yr,
      ymin = min_annual_net_exchange_tg_ch4_yr,
      ymax = max_annual_net_exchange_tg_ch4_yr,
      group = rate_scenario
    ),
    alpha = 0.14,
    color = NA
  ) +
  geom_line(
    data = threshold_sensitivity_summary,
    aes(y = mean_annual_net_exchange_tg_ch4_yr, group = rate_scenario),
    linewidth = 1
  ) +
  geom_point(
    data = threshold_sensitivity_summary,
    aes(y = mean_annual_net_exchange_tg_ch4_yr, group = rate_scenario),
    size = 1.8
  ) +
  scale_x_continuous(breaks = source_probability_thresholds) +
  scale_color_manual(values = rate_scenario_colors) +
  scale_fill_manual(values = rate_scenario_colors) +
  labs(
    title = "Spatial Source-Threshold Sensitivity",
    subtitle = "Higher thresholds require stronger evidence for weak-source classification.",
    x = "Weak-source probability threshold",
    y = "Net exchange (Tg CH4 yr-1)",
    color = "Rate scenario",
    fill = "Rate scenario"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "spatial_source_threshold_sensitivity.png"),
  plot_threshold_sensitivity,
  width = 9.5,
  height = 5.4,
  units = "in",
  dpi = 300
)

plot_annual_source_sink_area_inset <- plot_annual_source_sink_area +
  labs(
    title = "D. Annual Source/Sink Area",
    subtitle = NULL,
    x = NULL,
    y = NULL,
    color = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA),
    plot.margin = margin(3, 3, 3, 3)
  )

plot_source_sink_map_with_area_inset <- (
  plot_annual_source_sink_map +
    labs(title = "C. Mapped Source/Sink Classification, 2025") +
    theme(plot.margin = margin(4, 4, 4, 4)) +
    inset_element(
      plot_annual_source_sink_area_inset,
      left = 0.03,
      bottom = 0.08,
      right = 0.26,
      top = 0.40,
      align_to = "panel",
      on_top = TRUE
    )
)

plot_hard_threshold_rate_legend <- get_legend(
  plot_annual_time_series +
    guides(color = guide_legend(nrow = 1), fill = "none") +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
)

plot_hard_threshold_top_panels <- plot_grid(
  plot_annual_estimate +
    labs(tag = "A") +
    guides(fill = "none") +
    theme(
      legend.position = "none",
      plot.tag = element_text(face = "bold", size = 13)
    ),
  plot_annual_time_series +
    labs(tag = "B") +
    guides(color = "none", fill = "none") +
    theme(
      legend.position = "none",
      plot.tag = element_text(face = "bold", size = 13)
    ),
  ncol = 2,
  rel_widths = c(1.06, 1)
)

plot_hard_threshold_top_row <- plot_grid(
  plot_hard_threshold_top_panels,
  plot_hard_threshold_rate_legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

plot_hard_threshold_multipanel <- (
  plot_hard_threshold_top_row /
    plot_source_sink_map_with_area_inset /
    (plot_threshold_sensitivity + labs(title = "E. Spatial Source-Threshold Sensitivity"))
) +
  plot_layout(heights = c(1.02, 1.55, 1.05))

ggsave(
  file.path(figure_dir, "ERA5_hard_threshold_multipanel.png"),
  plot_hard_threshold_multipanel,
  width = 13,
  height = 15.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

plot_continuous_scenario_legend <- get_legend(
  plot_continuous_expected_flux +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_continuous_top_panels <- plot_grid(
  plot_continuous_expected_flux +
    labs(tag = "A") +
    theme(
      legend.position = "none",
      plot.tag = element_text(face = "bold", size = 13)
    ),
  plot_continuous_annual_budget_comparison +
    labs(tag = "B") +
    guides(fill = "none") +
    theme(
      legend.position = "none",
      plot.tag = element_text(face = "bold", size = 13)
    ),
  ncol = 2,
  rel_widths = c(1.08, 1)
)

plot_continuous_top_row <- plot_grid(
  plot_continuous_top_panels,
  plot_continuous_scenario_legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

plot_continuous_bottom_row <- plot_grid(
  plot_magnitude_model_skill +
    labs(tag = "C") +
    theme(plot.tag = element_text(face = "bold", size = 13)),
  plot_probability_calibration +
    labs(tag = "D") +
    theme(plot.tag = element_text(face = "bold", size = 13)),
  ncol = 2,
  rel_widths = c(1.18, 1)
)

plot_continuous_multipanel <- plot_grid(
  plot_continuous_top_row,
  plot_continuous_bottom_row,
  ncol = 1,
  rel_heights = c(0.95, 1.1)
)

ggsave(
  file.path(figure_dir, "ERA5_continuous_expected_flux_multipanel.png"),
  plot_continuous_multipanel,
  width = 13,
  height = 10.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Wrote monthly spatial upscaling trial outputs to ", spatial_dir)
