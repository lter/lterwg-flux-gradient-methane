# ERA5-driven half-hourly CH4 gap-fill for validation towers (SE-Sto, SE-Svb, US-Uaf)
#
# Adapted from NEON.ERA5.HalfHourlyGapfill.R.
# Trains two random-forest gap-fill models for each of the three validation sites:
#   (1) FG flux (flux-gradient method, from flow.30min.analysis.VAL.R)
#   (2) EC flux (eddy covariance reference, parallel pipeline)
# Both gap-fill products are then used to compute annual budgets and a direct
# FG vs. EC comparison at the 30-min and annual scale.
#
# Requires:
#   OUTPUT/VAL_30min_ch4_model_data.csv      (from flow.30min.analysis.VAL.R)
#   OUTPUT/VAL_30min_gapfill_annual_budgets.csv  (from Val.30min.Gapfill.R)
#
# Outputs (VAL_ prefix throughout):
#   VAL_ERA5_hourly_site_covariates.csv.gz
#   VAL_ERA5_30min_site_covariates.csv.gz
#   VAL_ERA5_gapfilled_30min.csv.gz            — both FG and EC filled columns
#   VAL_ERA5_gapfilled_annual_budget_by_year.csv
#   VAL_ERA5_mean_annual_budget.csv
#   VAL_ERA5_FG_vs_EC_30min_stats.csv
#   VAL_ERA5_FG_vs_EC_annual.csv
#   VAL_ERA5_rf_oob_fit_metrics.csv
#   VAL_ERA5_rf_sign_sensitivity.csv
#   VAL_ERA5_halfhour_gapfill_model_summary.txt

library(tidyverse)
library(ranger)
library(slider)
library(jsonlite)
library(data.table)
library(lubridate)

`%||%` <- function(x, y) if (is.null(x)) y else x

# ── Directories ───────────────────────────────────────────────────────────────
localdir.val <- Sys.getenv(
  "LOCALDIR_VAL",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/Validation_Sites"
)

if (!dir.exists(localdir.val)) stop("Validation data directory not found: ", localdir.val)
setwd(localdir.val)
dir.create("OUTPUT",  showWarnings = FALSE, recursive = TRUE)
dir.create("FIGURES", showWarnings = FALSE, recursive = TRUE)

# ── Input files ───────────────────────────────────────────────────────────────
model_data_file   <- "OUTPUT/VAL_30min_ch4_model_data.csv"
gam_budgets_file  <- "OUTPUT/VAL_30min_gapfill_annual_budgets.csv"   # optional, for comparison
era5_hourly_file  <- "OUTPUT/VAL_ERA5_hourly_site_covariates.csv.gz"
era5_halfhour_file <- "OUTPUT/VAL_ERA5_30min_site_covariates.csv.gz"

if (!file.exists(model_data_file))
  stop("Missing: ", model_data_file, "\nRun flow.30min.analysis.VAL.R first.")

# ── Constants ─────────────────────────────────────────────────────────────────
flux_to_mgC_30min       <- 2 * 0.0000288872 * 1000
seconds_per_30min       <- 30 * 60
ug_c_per_umol_c         <- 12.011
min_training_obs        <- 50L    # minimum obs to include a site in RF training
B_sign                  <- 500L   # bootstrap draws for sign sensitivity
behavior_levels         <- c("Consistent sink", "Fluctuating", "Consistent source")
season_levels           <- c("Winter", "Spring", "Summer", "Autumn")

mg_c_30min_to_umol_c_s  <- function(x) x * 1000 / ug_c_per_umol_c / seconds_per_30min

# ── Site metadata ─────────────────────────────────────────────────────────────
# Coordinates extracted from site attr files.
site_metadata <- tibble(
  SITE_ID   = c("SE-Sto", "SE-Svb", "US-Uaf"),
  latitude  = c(64.18203, 64.25611, 64.8663),
  longitude = c(19.55654, 19.7745, -147.8555),
  EcoType   = c("Wetland", "Forest", "Wetland")
)

# ── Load 30-min flux data ─────────────────────────────────────────────────────
ch4_30min <- read.csv(model_data_file) %>%
  mutate(
    SITE_ID      = as.character(SITE_ID),
    time.rounded = lubridate::ymd_hms(time.rounded, tz = "UTC", truncated = 3),
    Date         = as.Date(time.rounded),
    Year         = as.integer(format(time.rounded, "%Y")),
    month        = as.integer(format(time.rounded, "%m")),
    doy          = as.integer(format(time.rounded, "%j")),
    hour_num     = as.numeric(format(time.rounded, "%H")) +
                   as.numeric(format(time.rounded, "%M")) / 60,
    season       = factor(season, levels = season_levels),
    CH4_FG_mgC_30min = as.numeric(CH4_FG_mgC_30min),
    CH4_EC_mgC_30min = as.numeric(CH4_EC_mgC_30min)
  ) %>%
  filter(!is.na(SITE_ID), !is.na(time.rounded))

site_years <- ch4_30min %>% distinct(SITE_ID, Year)
year_min   <- min(site_years$Year, na.rm = TRUE)
year_max   <- max(site_years$Year, na.rm = TRUE)
start_date <- as.Date(paste0(year_min, "-01-01"))
end_date   <- as.Date(paste0(year_max, "-12-31"))

message("Data spans ", start_date, " to ", end_date,
        " for ", n_distinct(ch4_30min$SITE_ID), " sites.")

# ── ERA5 fetch (Open-Meteo Archive API) ───────────────────────────────────────
fetch_open_meteo_era5_site <- function(site, latitude, longitude,
                                        start_date, end_date) {
  hourly_vars <- paste(
    "temperature_2m", "soil_moisture_0_to_7cm", "soil_temperature_0_to_7cm",
    "soil_moisture_7_to_28cm", "precipitation", "shortwave_radiation",
    "relative_humidity_2m", sep = ","
  )
  query <- list(
    latitude   = latitude,  longitude  = longitude,
    start_date = as.character(start_date),
    end_date   = as.character(end_date),
    hourly     = hourly_vars, timezone = "UTC", models = "era5"
  )
  url <- paste0(
    "https://archive-api.open-meteo.com/v1/archive?",
    paste(paste0(names(query), "=",
                 vapply(query, function(x) utils::URLencode(as.character(x), reserved = TRUE),
                        character(1))),
          collapse = "&")
  )
  parsed <- jsonlite::fromJSON(url)
  tibble(
    SITE_ID        = site,
    time_hour      = as.POSIXct(parsed$hourly$time, format = "%Y-%m-%dT%H:%M", tz = "UTC"),
    ERA5_Tair_C    = as.numeric(parsed$hourly$temperature_2m),
    ERA5_VSWC      = as.numeric(parsed$hourly$soil_moisture_0_to_7cm),
    ERA5_Tsoil_C   = as.numeric(parsed$hourly$soil_temperature_0_to_7cm   %||% NA_real_),
    ERA5_VSWC_deep = as.numeric(parsed$hourly$soil_moisture_7_to_28cm     %||% NA_real_),
    ERA5_precip_mm = as.numeric(parsed$hourly$precipitation                %||% NA_real_),
    ERA5_SW_Wm2    = as.numeric(parsed$hourly$shortwave_radiation          %||% NA_real_),
    ERA5_RH_pct    = as.numeric(parsed$hourly$relative_humidity_2m        %||% NA_real_)
  )
}

# ── Cache check / fetch ───────────────────────────────────────────────────────
required_era5_cols <- c("ERA5_Tair_C", "ERA5_VSWC", "ERA5_Tsoil_C",
                         "ERA5_VSWC_deep", "ERA5_precip_mm", "ERA5_SW_Wm2", "ERA5_RH_pct")

era5_cache_valid <- FALSE
if (file.exists(era5_hourly_file)) {
  era5_hourly <- data.table::fread(era5_hourly_file) %>%
    as_tibble() %>%
    mutate(time_hour = as.POSIXct(time_hour, tz = "UTC"))
  cached_sites <- unique(era5_hourly$SITE_ID)
  cache_start  <- as.Date(min(era5_hourly$time_hour, na.rm = TRUE))
  cache_end    <- as.Date(max(era5_hourly$time_hour, na.rm = TRUE))
  missing_sites <- setdiff(site_metadata$SITE_ID, cached_sites)
  missing_cols  <- setdiff(required_era5_cols, names(era5_hourly))

  if (cache_start <= start_date && cache_end >= end_date &&
      length(missing_sites) == 0 && length(missing_cols) == 0) {
    message("ERA5 cache valid — skipping download.")
    era5_cache_valid <- TRUE
  } else {
    message("ERA5 cache incomplete — re-fetching.")
    era5_hourly <- NULL
  }
}

if (!era5_cache_valid) {
  message("Fetching ERA5 covariates for ", nrow(site_metadata), " validation sites...")
  fetch_results <- vector("list", nrow(site_metadata))
  for (i in seq_len(nrow(site_metadata))) {
    row <- site_metadata[i, ]
    message("  [", i, "/", nrow(site_metadata), "] ", row$SITE_ID)
    fetch_results[[i]] <- tryCatch(
      fetch_open_meteo_era5_site(row$SITE_ID, row$latitude, row$longitude,
                                  start_date, end_date),
      error = function(e) {
        stop("ERA5 fetch failed for ", row$SITE_ID, ": ", conditionMessage(e))
      }
    )
    Sys.sleep(0.5)
  }
  era5_hourly <- bind_rows(fetch_results)
  data.table::fwrite(era5_hourly, era5_hourly_file)
  message("ERA5 hourly data cached to ", era5_hourly_file)
}

# ── Interpolate hourly ERA5 → 30-min site-year grid ──────────────────────────
era5_halfhour_cache_valid <- FALSE
if (file.exists(era5_halfhour_file)) {
  era5_halfhour <- data.table::fread(era5_halfhour_file) %>%
    as_tibble() %>%
    mutate(
      time.rounded = as.POSIXct(time.rounded, tz = "UTC"),
      season  = factor(season, levels = season_levels),
      EcoType = factor(EcoType)
    )
  cached_sites  <- unique(era5_halfhour$SITE_ID)
  missing_sites <- setdiff(site_metadata$SITE_ID, cached_sites)
  missing_cols  <- setdiff(required_era5_cols, names(era5_halfhour))
  hh_start <- as.Date(min(era5_halfhour$time.rounded, na.rm = TRUE))
  hh_end   <- as.Date(max(era5_halfhour$time.rounded, na.rm = TRUE))

  if (hh_start <= start_date && hh_end >= end_date &&
      length(missing_sites) == 0 && length(missing_cols) == 0) {
    message("ERA5 30-min cache valid — skipping interpolation.")
    era5_halfhour_cache_valid <- TRUE
  } else {
    message("ERA5 30-min cache incomplete — re-interpolating.")
    era5_halfhour <- NULL
  }
}

if (!era5_halfhour_cache_valid) {
  # Full 30-min grid covering all site-years
  site_year_grid <- site_years %>%
    mutate(
      start_t = as.POSIXct(paste0(Year, "-01-01 00:00:00"), tz = "UTC"),
      end_t   = as.POSIXct(paste0(Year, "-12-31 23:30:00"), tz = "UTC")
    ) %>%
    mutate(time.rounded = purrr::map2(start_t, end_t, ~ seq(.x, .y, by = "30 min"))) %>%
    dplyr::select(SITE_ID, Year, time.rounded) %>%
    unnest(time.rounded) %>%
    mutate(
      Date     = as.Date(time.rounded),
      month    = as.integer(format(time.rounded, "%m")),
      doy      = as.integer(format(time.rounded, "%j")),
      hour_num = as.numeric(format(time.rounded, "%H")) +
                 as.numeric(format(time.rounded, "%M")) / 60,
      season   = factor(
        case_when(
          month %in% c(12, 1, 2) ~ "Winter",
          month %in% 3:5          ~ "Spring",
          month %in% 6:8          ~ "Summer",
          TRUE                    ~ "Autumn"
        ), levels = season_levels
      )
    ) %>%
    left_join(site_metadata %>% dplyr::select(SITE_ID, EcoType), by = "SITE_ID") %>%
    mutate(EcoType = factor(EcoType))

  # Linear interpolation helper
  era5_interp <- function(floor_v, ceil_v, w)
    if_else(is.finite(floor_v) & is.finite(ceil_v),
            floor_v + w * (ceil_v - floor_v),
            coalesce(floor_v, ceil_v))

  era5_halfhour <- site_year_grid %>%
    mutate(
      time_floor   = as.POSIXct(floor(as.numeric(time.rounded) / 3600) * 3600,
                                 origin = "1970-01-01", tz = "UTC"),
      time_ceiling = time_floor + hours(1),
      hw           = as.numeric(difftime(time.rounded, time_floor, units = "hours"))
    ) %>%
    left_join(era5_hourly %>% rename(
      time_floor = time_hour,
      Tf = ERA5_Tair_C, Vf = ERA5_VSWC, Tsf = ERA5_Tsoil_C, Vdf = ERA5_VSWC_deep,
      Pf = ERA5_precip_mm, Sf = ERA5_SW_Wm2, Rf = ERA5_RH_pct),
      by = c("SITE_ID", "time_floor")) %>%
    left_join(era5_hourly %>% rename(
      time_ceiling = time_hour,
      Tc = ERA5_Tair_C, Vc = ERA5_VSWC, Tsc = ERA5_Tsoil_C, Vdc = ERA5_VSWC_deep,
      Pc = ERA5_precip_mm, Sc = ERA5_SW_Wm2, Rc = ERA5_RH_pct),
      by = c("SITE_ID", "time_ceiling")) %>%
    mutate(
      ERA5_Tair_C    = era5_interp(Tf,  Tc,  hw),
      ERA5_VSWC      = era5_interp(Vf,  Vc,  hw),
      ERA5_Tsoil_C   = era5_interp(Tsf, Tsc, hw),
      ERA5_VSWC_deep = era5_interp(Vdf, Vdc, hw),
      ERA5_SW_Wm2    = era5_interp(Sf,  Sc,  hw),
      ERA5_RH_pct    = era5_interp(Rf,  Rc,  hw),
      ERA5_precip_mm = coalesce(Pf, Pc)   # precipitation: floor-hour rate
    ) %>%
    dplyr::select(SITE_ID, Year, time.rounded, Date, month, doy, hour_num, season, EcoType,
                  ERA5_Tair_C, ERA5_VSWC, ERA5_Tsoil_C, ERA5_VSWC_deep,
                  ERA5_precip_mm, ERA5_SW_Wm2, ERA5_RH_pct)

  data.table::fwrite(era5_halfhour, era5_halfhour_file)
  message("ERA5 30-min covariates cached to ", era5_halfhour_file)
}

# ── Climatological prior (per site × month × hour) ───────────────────────────
# Separately for FG and EC — each gets its own prior for training.

make_clim_prior <- function(ch4_30min, flux_col) {
  obs <- ch4_30min %>%
    filter(is.finite(.data[[flux_col]])) %>%
    mutate(season = as.character(season))

  smh  <- obs %>% reframe(.by = c(SITE_ID, month, hour_num),
                            clim_smh = mean(.data[[flux_col]], na.rm = TRUE))
  ssh  <- obs %>% reframe(.by = c(SITE_ID, season, hour_num),
                            clim_ssh = mean(.data[[flux_col]], na.rm = TRUE))
  sh   <- obs %>% reframe(.by = c(SITE_ID, hour_num),
                            clim_sh  = mean(.data[[flux_col]], na.rm = TRUE))
  gmh  <- obs %>% reframe(.by = c(month, hour_num),
                            clim_gmh = mean(.data[[flux_col]], na.rm = TRUE))
  list(smh = smh, ssh = ssh, sh = sh, gmh = gmh)
}

clim_FG <- make_clim_prior(ch4_30min, "CH4_FG_mgC_30min")
clim_EC <- make_clim_prior(ch4_30min, "CH4_EC_mgC_30min")

add_clim_prior <- function(df, clim) {
  df %>%
    left_join(clim$smh, by = c("SITE_ID", "month", "hour_num")) %>%
    left_join(clim$ssh %>% mutate(season = factor(season, levels = season_levels)),
              by = c("SITE_ID", "season", "hour_num")) %>%
    left_join(clim$sh,  by = c("SITE_ID", "hour_num")) %>%
    left_join(clim$gmh, by = c("month",   "hour_num")) %>%
    mutate(climatological_CH4_mgC_30min = coalesce(clim_smh, clim_ssh, clim_sh, clim_gmh)) %>%
    dplyr::select(-clim_smh, -clim_ssh, -clim_sh, -clim_gmh)
}

# ── Join ERA5 covariates with observed flux ───────────────────────────────────
base_data <- era5_halfhour %>%
  left_join(
    ch4_30min %>% dplyr::select(SITE_ID, time.rounded,
                                 CH4_FG_mgC_30min, CH4_EC_mgC_30min),
    by = c("SITE_ID", "time.rounded")
  ) %>%
  arrange(SITE_ID, time.rounded) %>%
  group_by(SITE_ID) %>%
  mutate(
    tair_roll7   = slider::slide_dbl(ERA5_Tair_C,    mean, .before = 335L,
                                     .complete = FALSE, na_rm = TRUE),
    precip_roll7 = slider::slide_dbl(ERA5_precip_mm, sum,  .before = 335L,
                                     .complete = FALSE, na_rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    ERA5_VPD_kPa = 0.6108 * exp(17.27 * ERA5_Tair_C / (ERA5_Tair_C + 237.3)) *
                   pmax(0, 1 - ERA5_RH_pct / 100),
    Tair_VSWC    = ERA5_Tair_C * ERA5_VSWC,
    sin_hour     = sin(2 * pi * hour_num / 24),
    cos_hour     = cos(2 * pi * hour_num / 24),
    season       = factor(season, levels = season_levels),
    EcoType      = factor(EcoType),
    SITE_ID      = factor(SITE_ID)   # factor for RF
  )

# Add climatological priors
gapfill_data_FG <- add_clim_prior(base_data, clim_FG) %>%
  mutate(observed_FG = is.finite(CH4_FG_mgC_30min))

gapfill_data_EC <- add_clim_prior(base_data, clim_EC) %>%
  mutate(observed_EC = is.finite(CH4_EC_mgC_30min))

# ── RF predictors ─────────────────────────────────────────────────────────────
# With only 3 sites, SITE_ID is included as a predictor in a single "site" RF.
# A population RF (without SITE_ID) is also retained for sites / periods with
# no training data.

rf_predictors_pop <- c(
  "ERA5_Tair_C", "ERA5_Tsoil_C", "ERA5_VSWC", "ERA5_VSWC_deep",
  "ERA5_precip_mm", "ERA5_SW_Wm2", "ERA5_VPD_kPa",
  "tair_roll7", "precip_roll7", "Tair_VSWC",
  "climatological_CH4_mgC_30min",
  "sin_hour", "cos_hour", "doy", "month", "season", "EcoType"
)
rf_predictors_site <- c(rf_predictors_pop, "SITE_ID")

# ── Fit RF for a given flux column ─────────────────────────────────────────────
fit_val_rf <- function(gapfill_data, flux_col, observed_col, flux_label) {
  training <- gapfill_data %>%
    filter(.data[[observed_col]],
           is.finite(.data[[flux_col]]),
           is.finite(ERA5_Tair_C), is.finite(ERA5_VSWC),
           is.finite(hour_num), is.finite(doy),
           is.finite(climatological_CH4_mgC_30min),
           !is.na(season), !is.na(EcoType)) %>%
    rename(response = all_of(flux_col)) %>%
    droplevels()

  message(flux_label, " training: ", nrow(training), " observations across ",
          n_distinct(training$SITE_ID), " sites.")

  if (nrow(training) < min_training_obs)
    stop("Too few training obs for ", flux_label, ": ", nrow(training))

  training_rf <- training %>%
    mutate(
      SITE_ID = factor(as.character(SITE_ID)),
      season  = factor(season, levels = season_levels),
      EcoType = factor(EcoType)
    )

  # mtry tuning on population predictors
  mtry_candidates <- unique(c(4L, 6L, 9L, floor(length(rf_predictors_pop) / 2L)))
  mtry_oob <- vapply(mtry_candidates, function(m) {
    rf_tmp <- ranger(
      formula   = response ~ .,
      data      = training_rf[, c("response", rf_predictors_pop)],
      num.trees = 100L, mtry = m, min.node.size = 5, sample.fraction = 0.7, seed = 42
    )
    rf_tmp$r.squared
  }, numeric(1))
  best_mtry_pop  <- mtry_candidates[which.max(mtry_oob)]
  best_mtry_site <- min(best_mtry_pop + 1L, length(rf_predictors_site))

  message(sprintf("  %s mtry sweep: best pop=%d (OOB R²=%.4f)",
                  flux_label, best_mtry_pop, max(mtry_oob)))

  set.seed(42)
  rf_site <- ranger(
    formula       = response ~ .,
    data          = training_rf[, c("response", rf_predictors_site)],
    num.trees     = 500, mtry = best_mtry_site, min.node.size = 5,
    sample.fraction = 0.7, importance = "permutation", seed = 42
  )

  set.seed(42)
  rf_pop <- ranger(
    formula       = response ~ .,
    data          = training_rf[, c("response", rf_predictors_pop)],
    num.trees     = 500, mtry = best_mtry_pop, min.node.size = 5,
    sample.fraction = 0.7, importance = "permutation", seed = 42
  )

  list(
    rf_site    = rf_site,
    rf_pop     = rf_pop,
    training   = training_rf,
    best_mtry_site = best_mtry_site,
    best_mtry_pop  = best_mtry_pop,
    flux_label = flux_label
  )
}

message("Fitting RF for FG flux...")
rf_FG <- fit_val_rf(gapfill_data_FG, "CH4_FG_mgC_30min", "observed_FG", "FG")

message("Fitting RF for EC flux...")
rf_EC <- fit_val_rf(gapfill_data_EC, "CH4_EC_mgC_30min", "observed_EC", "EC")

save(rf_FG, rf_EC, file = "OUTPUT/VAL_ERA5_rf_models.Rdata")

# ── OOB fit metrics ───────────────────────────────────────────────────────────
oob_metrics <- function(rf_obj) {
  tr <- rf_obj$training %>%
    mutate(
      oob_site = as.numeric(rf_obj$rf_site$predictions),
      oob_pop  = as.numeric(rf_obj$rf_pop$predictions),
      res_site = response - oob_site,
      res_pop  = response - oob_pop
    )

  site_metrics <- tr %>%
    group_by(SITE_ID) %>%
    summarise(
      n_obs                = n(),
      oob_rmse_site        = sqrt(mean(res_site^2, na.rm = TRUE)),
      oob_r_site           = suppressWarnings(cor(response, oob_site, use = "complete.obs")),
      sign_accuracy_site   = mean(sign(response) == sign(oob_site), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(flux = rf_obj$flux_label)

  global <- tr %>%
    summarise(
      SITE_ID            = "ALL",
      n_obs              = n(),
      oob_rmse_site      = sqrt(mean(res_site^2, na.rm = TRUE)),
      oob_r_site         = suppressWarnings(cor(response, oob_site, use = "complete.obs")),
      sign_accuracy_site = mean(sign(response) == sign(oob_site), na.rm = TRUE),
      flux               = rf_obj$flux_label
    )

  bind_rows(site_metrics, global)
}

oob_all <- bind_rows(oob_metrics(rf_FG), oob_metrics(rf_EC))
write.csv(oob_all, "OUTPUT/VAL_ERA5_rf_oob_fit_metrics.csv", row.names = FALSE)

# ── Predict on full grid ──────────────────────────────────────────────────────
apply_rf <- function(gapfill_data, rf_obj, flux_col, observed_col) {
  site_levels <- levels(rf_obj$training$SITE_ID)

  pred_data <- gapfill_data %>%
    mutate(
      SITE_ID_rf = factor(
        if_else(as.character(SITE_ID) %in% site_levels, as.character(SITE_ID), site_levels[1]),
        levels = site_levels
      ),
      season  = factor(season,  levels = season_levels),
      EcoType = factor(EcoType, levels = levels(rf_obj$training$EcoType))
    )

  pred_pop <- predict(rf_obj$rf_pop,
                       data = pred_data %>% dplyr::select(all_of(rf_predictors_pop)))$predictions
  pred_site <- predict(rf_obj$rf_site,
                        data = pred_data %>%
                          mutate(SITE_ID = SITE_ID_rf) %>%
                          dplyr::select(all_of(rf_predictors_site)))$predictions

  gapfill_data %>%
    mutate(
      pred_pop_flux   = as.numeric(pred_pop),
      pred_site_flux  = as.numeric(pred_site),
      # Use site RF (includes SITE_ID) as primary prediction
      pred_flux       = pred_site_flux,
      observed_flux   = .data[[observed_col]],
      CH4_obs         = .data[[flux_col]],
      gapfilled_flux  = if_else(observed_flux, CH4_obs, pred_flux),
      prediction_type = if_else(observed_flux, "observed", "site_rf_gapfill")
    )
}

pred_FG <- apply_rf(gapfill_data_FG, rf_FG, "CH4_FG_mgC_30min", "observed_FG")
pred_EC <- apply_rf(gapfill_data_EC, rf_EC, "CH4_EC_mgC_30min", "observed_EC")

# Merge FG and EC predictions into one wide table
era5_gapfilled_30min <- pred_FG %>%
  dplyr::select(SITE_ID, Year, time.rounded, Date, month, doy, hour_num, season, EcoType,
                ERA5_Tair_C, ERA5_Tsoil_C, ERA5_VSWC, ERA5_VSWC_deep,
                ERA5_precip_mm, ERA5_SW_Wm2, ERA5_VPD_kPa,
                tair_roll7, precip_roll7, climatological_CH4_mgC_30min,
                CH4_FG_obs     = CH4_obs,
                FG_pred        = pred_flux,
                FG_gapfilled   = gapfilled_flux,
                FG_obs_flag    = observed_flux,
                FG_pred_type   = prediction_type) %>%
  left_join(
    pred_EC %>% dplyr::select(SITE_ID, time.rounded,
                               CH4_EC_obs   = CH4_obs,
                               EC_pred      = pred_flux,
                               EC_gapfilled = gapfilled_flux,
                               EC_obs_flag  = observed_flux,
                               EC_pred_type = prediction_type),
    by = c("SITE_ID", "time.rounded")
  ) %>%
  mutate(SITE_ID = as.character(SITE_ID))

data.table::fwrite(era5_gapfilled_30min, "OUTPUT/VAL_ERA5_gapfilled_30min.csv.gz")
message("30-min gapfilled output written.")

# ── Annual budgets ────────────────────────────────────────────────────────────
annual_by_year <- era5_gapfilled_30min %>%
  group_by(SITE_ID, Year) %>%
  summarise(
    n_30min            = n(),
    n_FG_obs           = sum(FG_obs_flag, na.rm = TRUE),
    n_EC_obs           = sum(EC_obs_flag, na.rm = TRUE),
    FG_coverage        = n_FG_obs / n_30min,
    EC_coverage        = n_EC_obs / n_30min,
    FG_annual_gC_m2_yr = sum(FG_gapfilled, na.rm = TRUE) / 1000,
    EC_annual_gC_m2_yr = sum(EC_gapfilled, na.rm = TRUE) / 1000,
    FG_obs_sum_gC      = sum(CH4_FG_obs[FG_obs_flag],  na.rm = TRUE) / 1000,
    EC_obs_sum_gC      = sum(CH4_EC_obs[EC_obs_flag],  na.rm = TRUE) / 1000,
    .groups = "drop"
  ) %>%
  mutate(
    FG_minus_EC_annual = FG_annual_gC_m2_yr - EC_annual_gC_m2_yr,
    sign_agree_annual  = sign(FG_annual_gC_m2_yr) == sign(EC_annual_gC_m2_yr)
  )

mean_annual <- annual_by_year %>%
  group_by(SITE_ID) %>%
  summarise(
    n_years                    = n(),
    FG_mean_annual_gC_m2_yr    = mean(FG_annual_gC_m2_yr, na.rm = TRUE),
    FG_sd_annual_gC_m2_yr      = sd(FG_annual_gC_m2_yr,   na.rm = TRUE),
    EC_mean_annual_gC_m2_yr    = mean(EC_annual_gC_m2_yr, na.rm = TRUE),
    EC_sd_annual_gC_m2_yr      = sd(EC_annual_gC_m2_yr,   na.rm = TRUE),
    FG_behavior = case_when(
      mean(FG_annual_gC_m2_yr > 0) >= 0.75 ~ "Consistent source",
      mean(FG_annual_gC_m2_yr > 0) <= 0.25 ~ "Consistent sink",
      TRUE ~ "Fluctuating"
    ),
    EC_behavior = case_when(
      mean(EC_annual_gC_m2_yr > 0) >= 0.75 ~ "Consistent source",
      mean(EC_annual_gC_m2_yr > 0) <= 0.25 ~ "Consistent sink",
      TRUE ~ "Fluctuating"
    ),
    mean_FG_minus_EC = FG_mean_annual_gC_m2_yr - EC_mean_annual_gC_m2_yr,
    sign_agree       = sign(FG_mean_annual_gC_m2_yr) == sign(EC_mean_annual_gC_m2_yr),
    .groups = "drop"
  )

write.csv(annual_by_year, "OUTPUT/VAL_ERA5_gapfilled_annual_budget_by_year.csv", row.names = FALSE)
write.csv(mean_annual,    "OUTPUT/VAL_ERA5_mean_annual_budget.csv",              row.names = FALSE)

# ── 30-min FG vs EC stats (ERA5-filled) ──────────────────────────────────────
fg_ec_30min_stats <- era5_gapfilled_30min %>%
  filter(FG_obs_flag, EC_obs_flag,   # only paired observed values
         is.finite(CH4_FG_obs), is.finite(CH4_EC_obs)) %>%
  group_by(SITE_ID) %>%
  summarise(
    n_pairs      = n(),
    r_pearson    = suppressWarnings(cor(CH4_FG_obs, CH4_EC_obs, use = "complete.obs")),
    rmse         = sqrt(mean((CH4_FG_obs - CH4_EC_obs)^2, na.rm = TRUE)),
    mae          = mean(abs(CH4_FG_obs - CH4_EC_obs),     na.rm = TRUE),
    mean_bias    = mean(CH4_FG_obs - CH4_EC_obs,          na.rm = TRUE),
    sign_accuracy = mean(sign(CH4_FG_obs) == sign(CH4_EC_obs), na.rm = TRUE),
    .groups = "drop"
  )

write.csv(fg_ec_30min_stats, "OUTPUT/VAL_ERA5_FG_vs_EC_30min_stats.csv", row.names = FALSE)

# ── Sign sensitivity (bootstrap) ─────────────────────────────────────────────
# Same approach as NEON.ERA5.HalfHourlyGapfill.R Layer 3.

per_site_rmse_FG <- oob_all %>%
  filter(flux == "FG", SITE_ID != "ALL") %>%
  dplyr::select(SITE_ID, oob_rmse_site)

per_site_rmse_EC <- oob_all %>%
  filter(flux == "EC", SITE_ID != "ALL") %>%
  dplyr::select(SITE_ID, oob_rmse_site)

sign_sensitivity <- annual_by_year %>%
  left_join(per_site_rmse_FG %>% rename(rmse_FG = oob_rmse_site), by = "SITE_ID") %>%
  left_join(per_site_rmse_EC %>% rename(rmse_EC = oob_rmse_site), by = "SITE_ID") %>%
  mutate(
    rmse_FG = coalesce(rmse_FG, median(per_site_rmse_FG$oob_rmse_site, na.rm = TRUE)),
    rmse_EC = coalesce(rmse_EC, median(per_site_rmse_EC$oob_rmse_site, na.rm = TRUE)),
    FG_gap_sd = rmse_FG * sqrt(pmax(n_30min - n_FG_obs, 0L)) / 1000,
    EC_gap_sd = rmse_EC * sqrt(pmax(n_30min - n_EC_obs, 0L)) / 1000
  ) %>%
  rowwise() %>%
  mutate(
    FG_boot = list(FG_obs_sum_gC + rnorm(B_sign, FG_annual_gC_m2_yr - FG_obs_sum_gC, FG_gap_sd)),
    EC_boot = list(EC_obs_sum_gC + rnorm(B_sign, EC_annual_gC_m2_yr - EC_obs_sum_gC, EC_gap_sd)),
    FG_p_source   = mean(unlist(FG_boot) > 0),
    EC_p_source   = mean(unlist(EC_boot) > 0),
    FG_ci_lwr     = quantile(unlist(FG_boot), 0.025),
    FG_ci_upr     = quantile(unlist(FG_boot), 0.975),
    EC_ci_lwr     = quantile(unlist(EC_boot), 0.025),
    EC_ci_upr     = quantile(unlist(EC_boot), 0.975),
    FG_sign_stable = (FG_ci_lwr > 0) | (FG_ci_upr < 0),
    EC_sign_stable = (EC_ci_lwr > 0) | (EC_ci_upr < 0)
  ) %>%
  ungroup() %>%
  dplyr::select(-FG_boot, -EC_boot)

write.csv(sign_sensitivity, "OUTPUT/VAL_ERA5_rf_sign_sensitivity.csv", row.names = FALSE)

# ── FG vs EC annual comparison table ─────────────────────────────────────────
fg_ec_annual <- mean_annual %>%
  dplyr::select(SITE_ID, FG_mean_annual_gC_m2_yr, FG_sd_annual_gC_m2_yr, FG_behavior,
                EC_mean_annual_gC_m2_yr, EC_sd_annual_gC_m2_yr, EC_behavior,
                mean_FG_minus_EC, sign_agree)

write.csv(fg_ec_annual, "OUTPUT/VAL_ERA5_FG_vs_EC_annual.csv", row.names = FALSE)

# ── Summary report ────────────────────────────────────────────────────────────
capture.output(
  {
    cat("# Validation ERA5 Half-Hourly Gap-fill — FG and EC comparison\n")
    cat("Generated:", format(Sys.time()), "\n\n")

    cat("## RF OOB performance (site-level)\n")
    print(oob_all)

    cat("\n## ERA5 annual budgets (g C m-2 yr-1)\n")
    print(mean_annual)

    cat("\n## 30-min FG vs EC agreement (observed half-hours only)\n")
    print(fg_ec_30min_stats)

    cat("\n## Sign stability (bootstrap)\n")
    print(sign_sensitivity %>%
            dplyr::select(SITE_ID, Year, FG_annual_gC_m2_yr, EC_annual_gC_m2_yr,
                          FG_p_source, EC_p_source, FG_sign_stable, EC_sign_stable))
  },
  file = "OUTPUT/VAL_ERA5_halfhour_gapfill_model_summary.txt"
)

message("\nVal.ERA5.HalfHourlyGapfill.R complete. Outputs in ", file.path(localdir.val, "OUTPUT"))
