# ERA5-driven half-hourly total-flux CH4 gap filling for NEON sites.
#lo
# This workflow obtains hourly ERA5 point covariates for NEON tower sites,
# interpolates them to 30-minute timestamps, fits a total-flux GAM using ERA5
# temperature and soil volumetric water content, and compares the resulting
# annual budgets to the model-standardized 30-minute annual budget from
# NEON.30min.Gapfill.r and to the annual categories from flow.30min.analysis.R.

library(tidyverse)
library(ggplot2)
library(ranger)
library(slider)
library(patchwork)
library(jsonlite)
library(data.table)

localdir <- "/Volumes/MaloneLab/Research/FluxGradient"
localdir.ch4 <- "/Volumes/MaloneLab/Research/FluxGradient/Methane"
setwd(localdir.ch4)

dir.create("OUTPUT", showWarnings = FALSE)
dir.create("FIGURES", showWarnings = FALSE)

model_data_file <- "OUTPUT/30min_ch4_model_data.csv"
site_standardized_flux_file <- "OUTPUT/30min_site_standardized_flux.csv"
daily_flux_summary_file <- "OUTPUT/NEON_scale_daily_flux_summary.csv"
scaled_annual_budget_file <- "OUTPUT/NEON_scale_annual_budget_summary.csv"
model_standardized_budget_file <- "OUTPUT/NEON_30min_gapfill_annual_budgets.csv"
reference_annual_category_file <- "OUTPUT/NEON_scale_annual_budget_summary.csv"
metadata_file <- file.path(localdir, "Ameriflux_NEON field-sites.csv")

required_files <- c(
  model_data_file, site_standardized_flux_file, daily_flux_summary_file,
  scaled_annual_budget_file, model_standardized_budget_file,
  reference_annual_category_file, metadata_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

behavior_levels <- c("Weak-sink", "Fluctuating", "Weak-source")
# Color convention: blue = sink (uptake), grey = fluctuating, red = source (emission)
behavior_colors <- c(
  "Weak-sink"   = "#2166AC",
  "Fluctuating"       = "#4D4D4D",
  "Weak-source" = "#B2182B"
)
ecotype_colors <- c(
  "Cropland" = "#E69F00",
  "Forest" = "#009E73",
  "Grassland" = "#F0E442",
  "Shrubland" = "#CC79A7",
  "Wetland" = "#56B4E9"
)

era5_hourly_file <- "OUTPUT/NEON_ERA5_hourly_site_covariates.csv.gz"
era5_halfhour_file <- "OUTPUT/NEON_ERA5_30min_site_covariates.csv.gz"
era5_fetch_log_file <- "OUTPUT/NEON_ERA5_fetch_log.csv"

flux_to_mgC_30min <- 2 * 0.0000288872 * 1000
seconds_per_30min <- 30 * 60
ug_c_per_umol_c <- 12.011
min_training_obs_per_site <- 50

mg_c_30min_to_umol_c_s <- function(x) {
  x * 1000 / ug_c_per_umol_c / seconds_per_30min
}

ch4_30min <- read.csv(model_data_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    time.rounded = lubridate::ymd_hms(time.rounded, tz = "UTC", truncated = 3),
    Date = as.Date(time.rounded),
    Year = as.integer(format(time.rounded, "%Y")),
    month = as.integer(format(time.rounded, "%m")),
    doy = as.integer(format(time.rounded, "%j")),
    hour_num = as.numeric(format(time.rounded, "%H")) + as.numeric(format(time.rounded, "%M")) / 60,
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    CH4_mgC_30min = as.numeric(CH4_mgC_30min)
  ) %>%
  filter(!is.na(SITE_ID), !is.na(time.rounded), !is.na(Year))

reference_annual_behavior <- read.csv(reference_annual_category_file) %>%
  transmute(
    SITE_ID = as.character(SITE_ID),
    reference_annual_behavior = factor(behavior_annual_scaled, levels = behavior_levels),
    reference_prop_source_years = prop_positive_scaled_daily_annual,
    reference_annual_budget_gC_m2_yr = flux_scaled_daily_annual_gC_m2_yr
  )

site_years <- ch4_30min %>%
  distinct(SITE_ID, Year)

year_min <- min(site_years$Year, na.rm = TRUE)
year_max <- max(site_years$Year, na.rm = TRUE)
start_date <- as.Date(paste0(year_min, "-01-01"))
end_date <- as.Date(paste0(year_max, "-12-31"))

site_metadata <- read.csv(metadata_file) %>%
  transmute(
    SITE_ID = Site_Id.NEON,
    latitude = Latitude..degrees.,
    longitude = Longitude..degrees.,
    EcoType = case_when(
      Vegetation.Abbreviation..IGBP. %in% c("ENF", "DBF", "MF", "EBF", "SAV") ~ "Forest",
      Vegetation.Abbreviation..IGBP. == "WET" ~ "Wetland",
      Vegetation.Abbreviation..IGBP. == "GRA" ~ "Grassland",
      Vegetation.Abbreviation..IGBP. %in% c("CVM", "CRO") ~ "Cropland",
      Vegetation.Abbreviation..IGBP. == "OSH" ~ "Shrubland",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(SITE_ID %in% unique(ch4_30min$SITE_ID), is.finite(latitude), is.finite(longitude))

if (nrow(site_metadata) == 0) {
  stop("No site coordinates available for sites in 30-minute model data.")
}

fetch_open_meteo_era5_site <- function(site, latitude, longitude, start_date, end_date) {
  hourly_vars <- paste(
    "temperature_2m",
    "soil_moisture_0_to_7cm",
    "soil_temperature_0_to_7cm",
    "soil_moisture_7_to_28cm",
    "precipitation",
    "shortwave_radiation",
    "relative_humidity_2m",
    sep = ","
  )
  query <- list(
    latitude = latitude,
    longitude = longitude,
    start_date = as.character(start_date),
    end_date = as.character(end_date),
    hourly = hourly_vars,
    timezone = "UTC",
    models = "era5"
  )
  url <- paste0(
    "https://archive-api.open-meteo.com/v1/archive?",
    paste(paste0(names(query), "=", utils::URLencode(as.character(query), reserved = TRUE)), collapse = "&")
  )

  parsed <- jsonlite::fromJSON(url)

  required_vars <- c("time", "temperature_2m", "soil_moisture_0_to_7cm")
  missing_vars  <- required_vars[!required_vars %in% names(parsed$hourly)]
  if (length(missing_vars) > 0) {
    stop("Open-Meteo response missing variables for ", site, ": ", paste(missing_vars, collapse = ", "))
  }

  tibble(
    SITE_ID          = site,
    time_hour        = as.POSIXct(parsed$hourly$time, format = "%Y-%m-%dT%H:%M", tz = "UTC"),
    ERA5_Tair_C      = as.numeric(parsed$hourly$temperature_2m),
    ERA5_VSWC        = as.numeric(parsed$hourly$soil_moisture_0_to_7cm),
    ERA5_Tsoil_C     = as.numeric(parsed$hourly$soil_temperature_0_to_7cm   %||% NA_real_),
    ERA5_VSWC_deep   = as.numeric(parsed$hourly$soil_moisture_7_to_28cm     %||% NA_real_),
    ERA5_precip_mm   = as.numeric(parsed$hourly$precipitation                %||% NA_real_),
    ERA5_SW_Wm2      = as.numeric(parsed$hourly$shortwave_radiation          %||% NA_real_),
    ERA5_RH_pct      = as.numeric(parsed$hourly$relative_humidity_2m        %||% NA_real_),
    ERA5_latitude    = parsed$latitude,
    ERA5_longitude   = parsed$longitude,
    ERA5_elevation   = parsed$elevation %||% NA_real_
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

era5_cache_valid <- FALSE
if (file.exists(era5_hourly_file)) {
  message("Checking cached ERA5 hourly covariates: ", era5_hourly_file)
  era5_hourly <- data.table::fread(era5_hourly_file) %>%
    as_tibble() %>%
    mutate(time_hour = as.POSIXct(time_hour, tz = "UTC"))

  cached_start <- as.Date(min(era5_hourly$time_hour, na.rm = TRUE))
  cached_end   <- as.Date(max(era5_hourly$time_hour, na.rm = TRUE))
  cached_sites <- unique(era5_hourly$SITE_ID)
  required_sites <- unique(site_metadata$SITE_ID)
  missing_sites  <- setdiff(required_sites, cached_sites)

  required_hourly_cols <- c("ERA5_Tair_C", "ERA5_VSWC", "ERA5_Tsoil_C",
                             "ERA5_VSWC_deep", "ERA5_precip_mm", "ERA5_SW_Wm2", "ERA5_RH_pct")
  missing_cols <- setdiff(required_hourly_cols, names(era5_hourly))

  if (cached_start <= start_date && cached_end >= end_date &&
      length(missing_sites) == 0 && length(missing_cols) == 0) {
    message(
      "Cache covers ", cached_start, " to ", cached_end,
      " for all ", length(cached_sites), " sites — skipping download."
    )
    era5_cache_valid <- TRUE
  } else {
    if (cached_start > start_date || cached_end < end_date) {
      message(
        "Cache date range (", cached_start, " to ", cached_end,
        ") does not cover required range (", start_date, " to ", end_date,
        "). Re-fetching."
      )
    }
    if (length(missing_sites) > 0) {
      message(
        "Cache is missing ", length(missing_sites), " site(s): ",
        paste(missing_sites, collapse = ", "), ". Re-fetching."
      )
    }
    if (length(missing_cols) > 0) {
      message(
        "Cache missing ERA5 columns: ", paste(missing_cols, collapse = ", "), ". Re-fetching."
      )
    }
    era5_hourly <- NULL
  }
}

if (!era5_cache_valid) {
  message("Fetching ERA5 covariates from Open-Meteo Archive API for ", nrow(site_metadata), " sites.")

  fetch_results <- vector("list", nrow(site_metadata))
  fetch_log <- vector("list", nrow(site_metadata))

  for (i in seq_len(nrow(site_metadata))) {
    site_row <- site_metadata[i, ]
    message("[", i, "/", nrow(site_metadata), "] Fetching ", site_row$SITE_ID)
    result <- tryCatch(
      {
        dat <- fetch_open_meteo_era5_site(
          site = site_row$SITE_ID,
          latitude = site_row$latitude,
          longitude = site_row$longitude,
          start_date = start_date,
          end_date = end_date
        )
        list(data = dat, status = "ok", message = NA_character_)
      },
      error = function(e) list(data = NULL, status = "error", message = conditionMessage(e))
    )

    fetch_results[[i]] <- result$data
    fetch_log[[i]] <- tibble(
      SITE_ID = site_row$SITE_ID,
      latitude = site_row$latitude,
      longitude = site_row$longitude,
      status = result$status,
      message = result$message
    )
    Sys.sleep(0.4)
  }

  fetch_log <- bind_rows(fetch_log)
  write.csv(fetch_log, era5_fetch_log_file, row.names = FALSE)

  if (any(fetch_log$status != "ok")) {
    stop(
      "ERA5 fetch failed for ",
      sum(fetch_log$status != "ok"),
      " sites. See ", era5_fetch_log_file,
      ". First error: ", fetch_log$message[fetch_log$status != "ok"][1]
    )
  }

  era5_hourly <- bind_rows(fetch_results)
  data.table::fwrite(era5_hourly, era5_hourly_file)
}

site_year_grid <- site_years %>%
  mutate(
    start_time = as.POSIXct(paste0(Year, "-01-01 00:00:00"), tz = "UTC"),
    end_time = as.POSIXct(paste0(Year, "-12-31 23:30:00"), tz = "UTC")
  ) %>%
  mutate(time.rounded = purrr::map2(start_time, end_time, ~ seq(.x, .y, by = "30 min"))) %>%
  dplyr::select(SITE_ID, Year, time.rounded) %>%
  unnest(time.rounded) %>%
  mutate(
    Date = as.Date(time.rounded),
    month = as.integer(format(time.rounded, "%m")),
    doy = as.integer(format(time.rounded, "%j")),
    hour_num = as.numeric(format(time.rounded, "%H")) + as.numeric(format(time.rounded, "%M")) / 60,
    season = factor(
      case_when(
        month %in% c(12, 1, 2) ~ "Winter",
        month %in% c(3, 4, 5) ~ "Spring",
        month %in% c(6, 7, 8) ~ "Summer",
        TRUE ~ "Autumn"
      ),
      levels = c("Winter", "Spring", "Summer", "Autumn")
    )
  ) %>%
  left_join(site_metadata %>% dplyr::select(SITE_ID, EcoType), by = "SITE_ID") %>%
  mutate(EcoType = factor(EcoType))

era5_halfhour_cache_valid <- FALSE
if (file.exists(era5_halfhour_file)) {
  message("Checking cached ERA5 half-hour covariates: ", era5_halfhour_file)
  era5_halfhour <- data.table::fread(era5_halfhour_file) %>%
    as_tibble() %>%
    mutate(
      time.rounded = as.POSIXct(time.rounded, tz = "UTC"),
      Date = as.Date(Date),
      season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
      EcoType = factor(EcoType)
    )

  cached_hh_start  <- as.Date(min(era5_halfhour$time.rounded, na.rm = TRUE))
  cached_hh_end    <- as.Date(max(era5_halfhour$time.rounded, na.rm = TRUE))
  cached_hh_sites  <- unique(era5_halfhour$SITE_ID)
  required_sites   <- unique(site_metadata$SITE_ID)
  missing_hh_sites <- setdiff(required_sites, cached_hh_sites)

  required_hh_cols <- c("ERA5_Tair_C", "ERA5_VSWC", "ERA5_Tsoil_C",
                         "ERA5_VSWC_deep", "ERA5_precip_mm", "ERA5_SW_Wm2", "ERA5_RH_pct")
  missing_hh_cols  <- setdiff(required_hh_cols, names(era5_halfhour))

  if (cached_hh_start <= start_date && cached_hh_end >= end_date &&
      length(missing_hh_sites) == 0 && length(missing_hh_cols) == 0) {
    message(
      "Half-hour cache covers ", cached_hh_start, " to ", cached_hh_end,
      " for all ", length(cached_hh_sites), " sites — skipping interpolation."
    )
    era5_halfhour_cache_valid <- TRUE
  } else {
    if (cached_hh_start > start_date || cached_hh_end < end_date) {
      message(
        "Half-hour cache date range (", cached_hh_start, " to ", cached_hh_end,
        ") does not cover required range (", start_date, " to ", end_date,
        "). Re-interpolating."
      )
    }
    if (length(missing_hh_sites) > 0) {
      message(
        "Half-hour cache is missing ", length(missing_hh_sites), " site(s): ",
        paste(missing_hh_sites, collapse = ", "), ". Re-interpolating."
      )
    }
    if (length(missing_hh_cols) > 0) {
      message(
        "Half-hour cache missing columns: ", paste(missing_hh_cols, collapse = ", "), ". Re-interpolating."
      )
    }
    era5_halfhour <- NULL
  }
}

if (!era5_halfhour_cache_valid) {
  message("Interpolating hourly ERA5 covariates to 30-minute site-year grid.")
  # Linear interpolation with fallback to available endpoint
  era5_interp <- function(floor, ceiling, w)
    if_else(is.finite(floor) & is.finite(ceiling),
            floor + w * (ceiling - floor),
            coalesce(floor, ceiling))

  era5_halfhour <- site_year_grid %>%
    mutate(
      time_floor_hour = as.POSIXct(
        floor(as.numeric(time.rounded) / 3600) * 3600,
        origin = "1970-01-01",
        tz = "UTC"
      ),
      time_ceiling_hour = time_floor_hour + lubridate::hours(1),
      hour_weight = as.numeric(difftime(time.rounded, time_floor_hour, units = "hours"))
    ) %>%
    left_join(
      era5_hourly %>%
        dplyr::select(SITE_ID, time_floor_hour = time_hour,
                      ERA5_Tair_C_floor    = ERA5_Tair_C,
                      ERA5_VSWC_floor      = ERA5_VSWC,
                      ERA5_Tsoil_C_floor   = ERA5_Tsoil_C,
                      ERA5_VSWC_deep_floor = ERA5_VSWC_deep,
                      ERA5_precip_mm_floor = ERA5_precip_mm,
                      ERA5_SW_Wm2_floor    = ERA5_SW_Wm2,
                      ERA5_RH_pct_floor    = ERA5_RH_pct),
      by = c("SITE_ID", "time_floor_hour")
    ) %>%
    left_join(
      era5_hourly %>%
        dplyr::select(SITE_ID, time_ceiling_hour = time_hour,
                      ERA5_Tair_C_ceiling    = ERA5_Tair_C,
                      ERA5_VSWC_ceiling      = ERA5_VSWC,
                      ERA5_Tsoil_C_ceiling   = ERA5_Tsoil_C,
                      ERA5_VSWC_deep_ceiling = ERA5_VSWC_deep,
                      ERA5_precip_mm_ceiling = ERA5_precip_mm,
                      ERA5_SW_Wm2_ceiling    = ERA5_SW_Wm2,
                      ERA5_RH_pct_ceiling    = ERA5_RH_pct),
      by = c("SITE_ID", "time_ceiling_hour")
    ) %>%
    mutate(
      ERA5_Tair_C    = era5_interp(ERA5_Tair_C_floor,    ERA5_Tair_C_ceiling,    hour_weight),
      ERA5_VSWC      = era5_interp(ERA5_VSWC_floor,      ERA5_VSWC_ceiling,      hour_weight),
      ERA5_Tsoil_C   = era5_interp(ERA5_Tsoil_C_floor,   ERA5_Tsoil_C_ceiling,   hour_weight),
      ERA5_VSWC_deep = era5_interp(ERA5_VSWC_deep_floor, ERA5_VSWC_deep_ceiling, hour_weight),
      ERA5_SW_Wm2    = era5_interp(ERA5_SW_Wm2_floor,    ERA5_SW_Wm2_ceiling,    hour_weight),
      ERA5_RH_pct    = era5_interp(ERA5_RH_pct_floor,    ERA5_RH_pct_ceiling,    hour_weight),
      # Precipitation is a rate (mm/hr) — take floor-hour value (no interpolation)
      ERA5_precip_mm = coalesce(ERA5_precip_mm_floor, ERA5_precip_mm_ceiling)
    ) %>%
    dplyr::select(SITE_ID, Year, time.rounded, Date, month, doy, hour_num, season, EcoType,
                  ERA5_Tair_C, ERA5_VSWC, ERA5_Tsoil_C, ERA5_VSWC_deep,
                  ERA5_precip_mm, ERA5_SW_Wm2, ERA5_RH_pct)

  data.table::fwrite(era5_halfhour, era5_halfhour_file)
}

# ── Climatological flux lookup tables ────────────────────────────────────────
# Computed from all observed CH4 half-hours. Three levels of specificity;
# joined in cascade so the finest available value is used as a predictor.
# For gap periods the lookup is out-of-sample (no data exists for that slot
# by definition). For observed training steps the mean includes the target
# observation (1/n dilution) — acceptable for gap-filling applications.

obs_ch4 <- ch4_30min %>%
  filter(is.finite(CH4_mgC_30min)) %>%
  mutate(season = as.character(season))

clim_site_month_hour <- obs_ch4 %>%
  reframe(
    .by = c(SITE_ID, month, hour_num),
    clim_site_month_hour  = mean(CH4_mgC_30min, na.rm = TRUE),
    n_clim_site_month_hour = n()
  )

clim_site_season_hour <- obs_ch4 %>%
  reframe(
    .by = c(SITE_ID, season, hour_num),
    clim_site_season_hour = mean(CH4_mgC_30min, na.rm = TRUE)
  )

clim_site_hour <- obs_ch4 %>%
  reframe(
    .by = c(SITE_ID, hour_num),
    clim_site_hour = mean(CH4_mgC_30min, na.rm = TRUE)
  )

# Global fallback (month × hour) for sites with no lookup at all
clim_global_month_hour <- obs_ch4 %>%
  reframe(
    .by = c(month, hour_num),
    clim_global_month_hour = mean(CH4_mgC_30min, na.rm = TRUE)
  )

message("Climatological lookup tables computed (",
        nrow(clim_site_month_hour), " site-month-hour cells, ",
        n_distinct(clim_site_month_hour$SITE_ID), " sites).")

gapfill_data <- era5_halfhour %>%
  left_join(
    ch4_30min %>%
      dplyr::select(SITE_ID, time.rounded, CH4_mgC_30min),
    by = c("SITE_ID", "time.rounded")
  ) %>%
  # ── Antecedent / derived features (computed per site in temporal order) ──────
  arrange(SITE_ID, time.rounded) %>%
  group_by(SITE_ID) %>%
  mutate(
    # 7-day rolling mean air temperature (336 half-hour steps = 7 × 48)
    tair_roll7   = slider::slide_dbl(ERA5_Tair_C,    mean, .before = 335L, .complete = FALSE, na_rm = TRUE),
    # 7-day rolling sum precipitation (mm accumulated over previous 7 days)
    precip_roll7 = slider::slide_dbl(ERA5_precip_mm, sum,  .before = 335L, .complete = FALSE, na_rm = TRUE)
  ) %>%
  ungroup() %>%
  # ── Climatological lookup joins (cascade: finest available resolution) ────────
  left_join(
    clim_site_month_hour %>% dplyr::select(-n_clim_site_month_hour),
    by = c("SITE_ID", "month", "hour_num")
  ) %>%
  left_join(
    clim_site_season_hour %>%
      mutate(season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))),
    by = c("SITE_ID", "season", "hour_num")
  ) %>%
  left_join(clim_site_hour,         by = c("SITE_ID", "hour_num")) %>%
  left_join(clim_global_month_hour, by = c("month",   "hour_num")) %>%
  mutate(
    # Vapour pressure deficit (kPa) from T and RH
    ERA5_VPD_kPa  = 0.6108 * exp(17.27 * ERA5_Tair_C / (ERA5_Tair_C + 237.3)) *
                    pmax(0, 1 - ERA5_RH_pct / 100),
    # Temperature × soil moisture interaction (explicit for RF)
    Tair_VSWC     = ERA5_Tair_C * ERA5_VSWC,
    # Climatological prior: site × time-of-day mean — finest available resolution wins
    climatological_CH4_mgC_30min = coalesce(
      clim_site_month_hour, clim_site_season_hour,
      clim_site_hour, clim_global_month_hour
    ),
    observed_flux = is.finite(CH4_mgC_30min),
    sin_hour  = sin(2 * pi * hour_num / 24),
    cos_hour  = cos(2 * pi * hour_num / 24),
    SITE_ID   = factor(SITE_ID),
    season    = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    EcoType   = factor(EcoType),
    month     = as.integer(month)
  )

training_data <- gapfill_data %>%
  filter(
    observed_flux,
    is.finite(CH4_mgC_30min),
    is.finite(ERA5_Tair_C),
    is.finite(ERA5_VSWC),
    is.finite(hour_num),
    is.finite(doy),
    is.finite(climatological_CH4_mgC_30min),
    !is.na(season),
    !is.na(EcoType)
  ) %>%
  droplevels()

if (nrow(training_data) < 500) {
  stop("Too few observations with ERA5 covariates to fit model: ", nrow(training_data))
}

site_training_counts <- training_data %>%
  count(SITE_ID, name = "n_training_obs")

# ── Gap-fill models ───────────────────────────────────────────────────────────
#
# Two models are trained:
#   rf_site  : RF with SITE_ID as a predictor (site-level offsets).
#              Used for sites with >= min_training_obs_per_site observations.
#   rf_pop   : RF without SITE_ID (population-level, used for all gap periods).
# OOB predictions used for ALL evaluation (no data leakage).

rf_predictors_pop <- c(
  "ERA5_Tair_C", "ERA5_Tsoil_C", "ERA5_VSWC", "ERA5_VSWC_deep",
  "ERA5_precip_mm", "ERA5_SW_Wm2", "ERA5_VPD_kPa",
  "tair_roll7", "precip_roll7", "Tair_VSWC",
  "climatological_CH4_mgC_30min",
  "sin_hour", "cos_hour", "doy", "month", "season", "EcoType"
)
rf_predictors_site <- c(rf_predictors_pop, "SITE_ID")

training_rf <- training_data %>%
  mutate(
    SITE_ID = factor(as.character(SITE_ID)),
    season  = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    EcoType = factor(EcoType)
  )

# ── mtry tuning (population model) ───────────────────────────────────────────
message("Tuning mtry for population RF (100-tree sweep)...")
mtry_candidates <- unique(c(4L, 6L, 9L, floor(length(rf_predictors_pop) / 2L)))
mtry_oob <- vapply(mtry_candidates, function(m) {
  rf_tmp <- ranger(
    formula = CH4_mgC_30min ~ .,
    data    = training_rf[, c("CH4_mgC_30min", rf_predictors_pop)],
    num.trees = 100L, mtry = m, min.node.size = 5,
    sample.fraction = 0.7, seed = 42
  )
  rf_tmp$r.squared
}, numeric(1))
best_mtry_pop  <- mtry_candidates[which.max(mtry_oob)]
best_mtry_site <- min(best_mtry_pop + 1L, length(rf_predictors_site))
message(sprintf("mtry sweep: %s | best pop mtry = %d (OOB R² = %.4f)",
                paste(mtry_candidates, round(mtry_oob, 4), sep = "→", collapse = "; "),
                best_mtry_pop, max(mtry_oob)))

# ── RF site model ─────────────────────────────────────────────────────────────
message("Fitting site-level RF (", nrow(training_rf), " obs, mtry = ", best_mtry_site, ")...")
set.seed(42)
era5_gapfill_rf_site <- ranger(
  formula         = CH4_mgC_30min ~ .,
  data            = training_rf[, c("CH4_mgC_30min", rf_predictors_site)],
  num.trees       = 500,
  mtry            = best_mtry_site,
  min.node.size   = 5,
  sample.fraction = 0.7,
  importance      = "permutation",
  seed            = 42
)

# ── RF population model ───────────────────────────────────────────────────────
message("Fitting population RF (mtry = ", best_mtry_pop, ")...")
set.seed(42)
era5_gapfill_rf_pop <- ranger(
  formula         = CH4_mgC_30min ~ .,
  data            = training_rf[, c("CH4_mgC_30min", rf_predictors_pop)],
  num.trees       = 500,
  mtry            = best_mtry_pop,
  min.node.size   = 5,
  sample.fraction = 0.7,
  importance      = "permutation",
  seed            = 42
)

capture.output(
  {
    cat("NEON ERA5-driven half-hourly total-flux gapfill — Random Forest\n")
    cat("Training observations:", nrow(training_rf), "\n")
    cat("ERA5 hourly file:", era5_hourly_file, "\n")
    cat("ERA5 half-hour file:", era5_halfhour_file, "\n\n")

    cat("── Site-level RF (with SITE_ID, mtry =", best_mtry_site, ") ──\n")
    cat("OOB R²:  ", round(era5_gapfill_rf_site$r.squared, 4), "\n")
    cat("OOB MSE: ", round(era5_gapfill_rf_site$prediction.error, 6), "\n\n")

    cat("── Population RF (no SITE_ID, mtry =", best_mtry_pop, ") ──\n")
    cat("OOB R²:  ", round(era5_gapfill_rf_pop$r.squared, 4), "\n")
    cat("OOB MSE: ", round(era5_gapfill_rf_pop$prediction.error, 6), "\n\n")

    cat("── Variable importance (site-level RF) ──\n")
    imp <- sort(era5_gapfill_rf_site$variable.importance, decreasing = TRUE)
    print(round(imp, 6))
  },
  file = "OUTPUT/NEON_ERA5_halfhour_gapfill_model_summary.txt"
)

# ── OOB fit metrics (honest — no data leakage) ───────────────────────────────

training_fit <- training_rf %>%
  mutate(
    oob_pred_site          = as.numeric(era5_gapfill_rf_site$predictions),
    oob_pred_pop           = as.numeric(era5_gapfill_rf_pop$predictions),
    residual_site          = CH4_mgC_30min - oob_pred_site,
    residual_pop           = CH4_mgC_30min - oob_pred_pop,
    # convenience aliases for downstream code
    fitted_CH4_mgC_30min   = oob_pred_site,
    residual_CH4_mgC_30min = residual_site
  )

fit_metrics <- training_fit %>%
  summarise(
    n_training               = dplyr::n(),
    # Site RF (OOB)
    rmse_site_mgC_m2_30min  = sqrt(mean(residual_site^2,  na.rm = TRUE)),
    mae_site_mgC_m2_30min   = mean(abs(residual_site),    na.rm = TRUE),
    bias_site_mgC_m2_30min  = mean(residual_site,         na.rm = TRUE),
    r_site                  = suppressWarnings(cor(CH4_mgC_30min, oob_pred_site, use = "complete.obs")),
    sign_accuracy_site_30min = mean(sign(CH4_mgC_30min) == sign(oob_pred_site), na.rm = TRUE),
    oob_r2_site             = era5_gapfill_rf_site$r.squared,
    # Population RF (OOB)
    rmse_pop_mgC_m2_30min   = sqrt(mean(residual_pop^2,   na.rm = TRUE)),
    mae_pop_mgC_m2_30min    = mean(abs(residual_pop),     na.rm = TRUE),
    bias_pop_mgC_m2_30min   = mean(residual_pop,          na.rm = TRUE),
    r_pop                   = suppressWarnings(cor(CH4_mgC_30min, oob_pred_pop, use = "complete.obs")),
    sign_accuracy_pop_30min = mean(sign(CH4_mgC_30min) == sign(oob_pred_pop),  na.rm = TRUE),
    oob_r2_pop              = era5_gapfill_rf_pop$r.squared,
    observed_sd_mgC_m2_30min = sd(CH4_mgC_30min, na.rm = TRUE)
  )

# Legacy aliases for downstream summary text and flow.plots.R
fit_metrics$rmse_mgC_m2_30min          <- fit_metrics$rmse_site_mgC_m2_30min
fit_metrics$mae_mgC_m2_30min           <- fit_metrics$mae_site_mgC_m2_30min
fit_metrics$bias_mgC_m2_30min          <- fit_metrics$bias_site_mgC_m2_30min
fit_metrics$correlation_observed_fitted <- fit_metrics$r_site

write.csv(fit_metrics, "OUTPUT/NEON_ERA5_halfhour_gapfill_fit_metrics.csv", row.names = FALSE)

# Per-site OOB metrics (used in sign sensitivity)
fit_metrics_site <- training_fit %>%
  group_by(SITE_ID) %>%
  summarise(
    n_obs                   = dplyr::n(),
    oob_rmse_site           = sqrt(mean(residual_site^2,  na.rm = TRUE)),
    oob_r_site              = suppressWarnings(cor(CH4_mgC_30min, oob_pred_site, use = "complete.obs")),
    sign_accuracy_30min     = mean(sign(CH4_mgC_30min) == sign(oob_pred_site), na.rm = TRUE),
    .groups = "drop"
  )

set.seed(20260519)
fit_plot_data <- training_fit %>%
  slice_sample(n = min(10000, nrow(training_fit)))

write.csv(fit_plot_data, "OUTPUT/NEON_ERA5_halfhour_gapfill_fit_plot_data.csv", row.names = FALSE)


# ── Partial dependence (marginal effects) for the population RF ───────────────
# Each driver is varied across its observed range while all other predictors are
# held at reference values. The population RF (no SITE_ID) is used so predictions
# reflect the marginal effect without site-level confounding.

reference_values <- training_rf %>%
  summarise(
    hour_num       = 12,
    doy            = 182,
    month          = 7L,
    ERA5_Tair_C    = median(ERA5_Tair_C,    na.rm = TRUE),
    ERA5_VSWC      = median(ERA5_VSWC,      na.rm = TRUE),
    ERA5_Tsoil_C   = median(ERA5_Tsoil_C,   na.rm = TRUE),
    ERA5_VSWC_deep = median(ERA5_VSWC_deep, na.rm = TRUE),
    ERA5_precip_mm = median(ERA5_precip_mm, na.rm = TRUE),
    ERA5_SW_Wm2    = median(ERA5_SW_Wm2,    na.rm = TRUE),
    ERA5_VPD_kPa   = median(ERA5_VPD_kPa,   na.rm = TRUE),
    tair_roll7     = median(tair_roll7,      na.rm = TRUE),
    precip_roll7   = median(precip_roll7,    na.rm = TRUE),
    Tair_VSWC      = median(Tair_VSWC,       na.rm = TRUE),
    climatological_CH4_mgC_30min = median(climatological_CH4_mgC_30min, na.rm = TRUE),
    season         = names(sort(table(season),  decreasing = TRUE))[1],
    EcoType        = names(sort(table(EcoType), decreasing = TRUE))[1]
  )

make_effect_grid <- function(driver, values) {
  tibble(
    driver         = driver,
    driver_value   = values,
    hour_num       = reference_values$hour_num,
    doy            = reference_values$doy,
    month          = reference_values$month,
    ERA5_Tair_C    = reference_values$ERA5_Tair_C,
    ERA5_VSWC      = reference_values$ERA5_VSWC,
    ERA5_Tsoil_C   = reference_values$ERA5_Tsoil_C,
    ERA5_VSWC_deep = reference_values$ERA5_VSWC_deep,
    ERA5_precip_mm = reference_values$ERA5_precip_mm,
    ERA5_SW_Wm2    = reference_values$ERA5_SW_Wm2,
    ERA5_VPD_kPa   = reference_values$ERA5_VPD_kPa,
    tair_roll7     = reference_values$tair_roll7,
    precip_roll7   = reference_values$precip_roll7,
    Tair_VSWC      = reference_values$Tair_VSWC,
    climatological_CH4_mgC_30min = reference_values$climatological_CH4_mgC_30min,
    season         = reference_values$season,
    EcoType        = reference_values$EcoType
  ) %>%
    mutate(
      sin_hour       = sin(2 * pi * hour_num / 24),
      cos_hour       = cos(2 * pi * hour_num / 24),
      hour_num       = if_else(driver == "hour_num",       driver_value, hour_num),
      doy            = if_else(driver == "doy",            driver_value, doy),
      ERA5_Tair_C    = if_else(driver == "ERA5_Tair_C",    driver_value, ERA5_Tair_C),
      ERA5_VSWC      = if_else(driver == "ERA5_VSWC",      driver_value, ERA5_VSWC),
      ERA5_Tsoil_C   = if_else(driver == "ERA5_Tsoil_C",   driver_value, ERA5_Tsoil_C),
      ERA5_VSWC_deep = if_else(driver == "ERA5_VSWC_deep", driver_value, ERA5_VSWC_deep),
      ERA5_precip_mm = if_else(driver == "ERA5_precip_mm", driver_value, ERA5_precip_mm),
      ERA5_SW_Wm2    = if_else(driver == "ERA5_SW_Wm2",    driver_value, ERA5_SW_Wm2),
      ERA5_VPD_kPa   = if_else(driver == "ERA5_VPD_kPa",   driver_value, ERA5_VPD_kPa),
      tair_roll7     = if_else(driver == "tair_roll7",     driver_value, tair_roll7),
      precip_roll7   = if_else(driver == "precip_roll7",   driver_value, precip_roll7),
      climatological_CH4_mgC_30min = if_else(driver == "climatological_CH4_mgC_30min",
                                             driver_value, climatological_CH4_mgC_30min),
      Tair_VSWC      = ERA5_Tair_C * ERA5_VSWC,   # keep consistent with training
      sin_hour       = sin(2 * pi * hour_num / 24),
      cos_hour       = cos(2 * pi * hour_num / 24),
      season         = factor(season,  levels = levels(training_rf$season)),
      EcoType        = factor(EcoType, levels = levels(training_rf$EcoType))
    )
}

make_range <- function(col, lo = 0.02, hi = 0.98, n = 100)
  seq(quantile(training_rf[[col]], lo, na.rm = TRUE),
      quantile(training_rf[[col]], hi, na.rm = TRUE),
      length.out = n)

effect_grid <- bind_rows(
  make_effect_grid("hour_num",       seq(0, 23.5, by = 0.5)),
  make_effect_grid("doy",            seq(1, 366, length.out = 100)),
  make_effect_grid("ERA5_Tair_C",    make_range("ERA5_Tair_C")),
  make_effect_grid("ERA5_VSWC",      make_range("ERA5_VSWC")),
  make_effect_grid("ERA5_Tsoil_C",   make_range("ERA5_Tsoil_C")),
  make_effect_grid("ERA5_VPD_kPa",   make_range("ERA5_VPD_kPa")),
  make_effect_grid("ERA5_SW_Wm2",    make_range("ERA5_SW_Wm2")),
  make_effect_grid("ERA5_precip_mm", make_range("ERA5_precip_mm")),
  make_effect_grid("tair_roll7",     make_range("tair_roll7")),
  make_effect_grid("precip_roll7",   make_range("precip_roll7")),
  make_effect_grid("climatological_CH4_mgC_30min", make_range("climatological_CH4_mgC_30min"))
)

effect_preds <- predict(era5_gapfill_rf_pop, data = effect_grid)$predictions

driver_labels <- c(
  hour_num       = "Hour of day",
  doy            = "Day of year",
  ERA5_Tair_C    = "ERA5 air temperature (°C)",
  ERA5_VSWC      = "ERA5 soil moisture 0–7 cm (m³/m³)",
  ERA5_Tsoil_C   = "ERA5 soil temperature 0–7 cm (°C)",
  ERA5_VPD_kPa   = "ERA5 vapour pressure deficit (kPa)",
  ERA5_SW_Wm2    = "ERA5 shortwave radiation (W m⁻²)",
  ERA5_precip_mm = "ERA5 precipitation (mm hr⁻¹)",
  tair_roll7     = "7-day rolling mean air temp (°C)",
  precip_roll7   = "7-day cumulative precip (mm)"
)

effect_grid <- effect_grid %>%
  mutate(
    pred_CH4_mgC_30min  = as.numeric(effect_preds),
    se_CH4_mgC_30min    = NA_real_,
    lower_CH4_mgC_30min = NA_real_,
    upper_CH4_mgC_30min = NA_real_,
    driver_label = driver_labels[driver]
  )

write.csv(effect_grid, "OUTPUT/NEON_ERA5_halfhour_gapfill_model_effects.csv", row.names = FALSE)



# ── Prediction: apply RF models to full site-year grid ───────────────────────
# Sites with >= min_training_obs_per_site use the site-level RF (SITE_ID known).
# Sites with too few training obs use the population RF (no SITE_ID).

model_site_levels <- levels(training_rf$SITE_ID)

prediction_data <- gapfill_data %>%
  mutate(SITE_ID_original = as.character(SITE_ID)) %>%
  left_join(site_training_counts %>% mutate(SITE_ID = as.character(SITE_ID)), by = "SITE_ID") %>%
  mutate(
    n_training_obs  = replace_na(n_training_obs, 0L),
    use_site_model  = SITE_ID_original %in% model_site_levels &
                      n_training_obs >= min_training_obs_per_site,
    # Ensure SITE_ID factor has the same levels as training for the site RF
    SITE_ID_rf = factor(
      if_else(SITE_ID_original %in% model_site_levels, SITE_ID_original, model_site_levels[1]),
      levels = model_site_levels
    ),
    season  = factor(season,  levels = c("Winter", "Spring", "Summer", "Autumn")),
    EcoType = factor(EcoType, levels = levels(training_rf$EcoType))
  )

# Population predictions (no SITE_ID) — population RF
pred_pop_all <- predict(
  era5_gapfill_rf_pop,
  data = prediction_data %>% dplyr::select(all_of(rf_predictors_pop))
)$predictions

# Site predictions only for rows where the site was adequately trained
pred_site_all <- rep(NA_real_, nrow(prediction_data))
site_rows <- which(prediction_data$use_site_model)

if (length(site_rows) > 0) {
  pred_site_all[site_rows] <- predict(
    era5_gapfill_rf_site,
    data = prediction_data[site_rows, ] %>%
      mutate(SITE_ID = SITE_ID_rf) %>%
      dplyr::select(all_of(rf_predictors_site))
  )$predictions
}

era5_gapfilled_30min <- prediction_data %>%
  mutate(
    pred_CH4_mgC_30min_population = as.numeric(pred_pop_all),
    pred_CH4_mgC_30min_site       = as.numeric(pred_site_all),
    pred_CH4_mgC_30min = if_else(
      use_site_model,
      pred_CH4_mgC_30min_site,
      pred_CH4_mgC_30min_population
    ),
    gapfilled_CH4_mgC_30min = if_else(observed_flux, CH4_mgC_30min, pred_CH4_mgC_30min),
    prediction_type = case_when(
      observed_flux   ~ "observed",
      use_site_model  ~ "site_rf_gapfill",
      TRUE            ~ "population_rf_gapfill"
    ),
    SITE_ID = SITE_ID_original
  ) %>%
  dplyr::select(
    SITE_ID, Year, time.rounded, Date, month, doy, hour_num, season, EcoType,
    ERA5_Tair_C, ERA5_Tsoil_C, ERA5_VSWC, ERA5_VSWC_deep,
    ERA5_precip_mm, ERA5_SW_Wm2, ERA5_VPD_kPa,
    tair_roll7, precip_roll7, Tair_VSWC, climatological_CH4_mgC_30min,
    CH4_mgC_30min, observed_flux,
    pred_CH4_mgC_30min, gapfilled_CH4_mgC_30min, prediction_type, n_training_obs
  )

data.table::fwrite(era5_gapfilled_30min, "OUTPUT/NEON_ERA5_gapfilled_30min.csv.gz")

# ── Sign sensitivity analysis ─────────────────────────────────────────────────
#
# Three layers of sign uncertainty quantification:
#
#  Layer 1 — 30-min OOB sign accuracy
#    For each observed half-hour, does the OOB prediction have the correct sign?
#    Computed per site so we know which sites the model is most uncertain about.
#
#  Layer 2 — Annual budget sign from OOB predictions
#    Treating all observations as if they were gap-filled (using OOB predictions),
#    does the annual budget have the same sign as the observed annual budget?
#    Sign accuracy at the annual scale is the relevant metric for the paper.
#
#  Layer 3 — Bootstrap sign probability for actual gap-filled budgets
#    For each site-year, the actual annual budget mixes observed (fixed) and
#    gap-filled (uncertain) observations.  The gap-fill uncertainty is
#    approximated as N(0, oob_rmse × sqrt(n_gapfilled)) per site (conservative:
#    assumes independent errors across time steps).  Drawing B samples from this
#    distribution gives a distribution of plausible annual budgets, and
#    P(source) = fraction of draws that are positive.  Sites with P(source)
#    near 0.5 have sign-uncertain budgets.

# Layer 1: OOB sign accuracy (already in fit_metrics_site)
write.csv(fit_metrics_site, "OUTPUT/NEON_ERA5_rf_oob_sign_accuracy_by_site.csv", row.names = FALSE)

# Layer 2: Annual budget sign from OOB predictions vs. observed
oob_annual_budget <- training_rf %>%
  mutate(
    oob_pred_site = as.numeric(era5_gapfill_rf_site$predictions),
    Year = as.integer(format(as.Date(substr(as.character(SITE_ID), 1, 4), "%Y"), "%Y"))
  )

# Re-extract Year from training_data (which has it)
oob_annual_budget <- training_data %>%
  mutate(
    oob_pred_site = as.numeric(era5_gapfill_rf_site$predictions),
    oob_pred_pop  = as.numeric(era5_gapfill_rf_pop$predictions)
  ) %>%
  group_by(SITE_ID, Year) %>%
  summarise(
    n_obs                      = dplyr::n(),
    annual_obs_gC_m2_yr        = sum(CH4_mgC_30min,  na.rm = TRUE) / 1000,
    annual_oob_site_gC_m2_yr   = sum(oob_pred_site,  na.rm = TRUE) / 1000,
    annual_oob_pop_gC_m2_yr    = sum(oob_pred_pop,   na.rm = TRUE) / 1000,
    .groups = "drop"
  ) %>%
  mutate(
    sign_agree_site = sign(annual_obs_gC_m2_yr) == sign(annual_oob_site_gC_m2_yr),
    sign_agree_pop  = sign(annual_obs_gC_m2_yr) == sign(annual_oob_pop_gC_m2_yr)
  )

oob_sign_summary <- oob_annual_budget %>%
  summarise(
    n_site_years         = dplyr::n(),
    sign_accuracy_site   = mean(sign_agree_site, na.rm = TRUE),
    sign_accuracy_pop    = mean(sign_agree_pop,  na.rm = TRUE)
  )

write.csv(oob_annual_budget,  "OUTPUT/NEON_ERA5_rf_oob_annual_budget_sign.csv",    row.names = FALSE)
write.csv(oob_sign_summary,   "OUTPUT/NEON_ERA5_rf_oob_sign_summary.csv",          row.names = FALSE)

# Layer 3: Bootstrap P(source) for gap-filled annual budgets
# Uses the per-site OOB RMSE as the per-observation uncertainty.
# Gap-fill contribution uncertainty: SD = oob_rmse × sqrt(n_gapfilled)
# (conservative: independent errors; actual uncertainty is smaller if errors are
# correlated in time, as is typical for gap-filling.)

B_sign <- 500L   # bootstrap draws per site-year

set.seed(42)
era5_annual_budget_pre <- era5_gapfilled_30min %>%
  group_by(SITE_ID, Year) %>%
  summarise(
    n_30min                         = dplyr::n(),
    n_observed                      = sum(observed_flux,    na.rm = TRUE),
    n_gapfilled                     = sum(!observed_flux,   na.rm = TRUE),
    observed_coverage               = mean(observed_flux,   na.rm = TRUE),
    n_training_obs                  = max(n_training_obs,   na.rm = TRUE),
    annual_budget_gC_m2_yr          = sum(gapfilled_CH4_mgC_30min, na.rm = TRUE) / 1000,
    model_only_annual_budget_gC_m2_yr = sum(pred_CH4_mgC_30min,    na.rm = TRUE) / 1000,
    observed_partial_gC_m2          = sum(CH4_mgC_30min[observed_flux],  na.rm = TRUE) / 1000,
    gapfilled_partial_gC_m2         = sum(gapfilled_CH4_mgC_30min[!observed_flux], na.rm = TRUE) / 1000,
    .groups = "drop"
  )

sign_sensitivity <- era5_annual_budget_pre %>%
  left_join(
    fit_metrics_site %>% dplyr::select(SITE_ID = SITE_ID, oob_rmse_site),
    by = "SITE_ID"
  ) %>%
  mutate(
    # Fall back to global median RMSE for sites without OOB metrics
    oob_rmse_site = coalesce(
      oob_rmse_site,
      median(fit_metrics_site$oob_rmse_site, na.rm = TRUE)
    ),
    # SD of gap-fill contribution (assuming independent errors across time steps)
    gap_fill_sd_gC_m2_yr = oob_rmse_site * sqrt(pmax(n_gapfilled, 0L)) / 1000,
    # Sign margin: how far the budget is from zero relative to gap-fill uncertainty
    sign_margin_sd = if_else(
      gap_fill_sd_gC_m2_yr > 0,
      annual_budget_gC_m2_yr / gap_fill_sd_gC_m2_yr,
      NA_real_
    )
  ) %>%
  rowwise() %>%
  mutate(
    # Bootstrap: draw B samples of the gap-fill contribution from N(mean, SD)
    boot_budgets  = list(
      observed_partial_gC_m2 +
        rnorm(B_sign, mean = gapfilled_partial_gC_m2, sd = gap_fill_sd_gC_m2_yr)
    ),
    p_source      = mean(unlist(boot_budgets) > 0),
    budget_ci_lwr = quantile(unlist(boot_budgets), 0.025),
    budget_ci_upr = quantile(unlist(boot_budgets), 0.975),
    # Stable = 95% CI does not straddle zero
    sign_stable   = (budget_ci_lwr > 0) | (budget_ci_upr < 0),
    sign_category = case_when(
      p_source >= 0.95 ~ "Stable source",
      p_source <= 0.05 ~ "Stable sink",
      p_source >  0.50 ~ "Probable source (uncertain)",
      TRUE             ~ "Probable sink (uncertain)"
    )
  ) %>%
  ungroup() %>%
  dplyr::select(-boot_budgets)

write.csv(sign_sensitivity, "OUTPUT/NEON_ERA5_rf_sign_sensitivity.csv", row.names = FALSE)

# Site-level sign sensitivity summary (aggregate across years)
sign_sensitivity_site <- sign_sensitivity %>%
  group_by(SITE_ID) %>%
  summarise(
    n_years              = dplyr::n(),
    mean_annual_budget   = mean(annual_budget_gC_m2_yr,  na.rm = TRUE),
    mean_p_source        = mean(p_source,                 na.rm = TRUE),
    frac_sign_stable     = mean(sign_stable,              na.rm = TRUE),
    frac_stable_source   = mean(sign_category == "Stable source",               na.rm = TRUE),
    frac_stable_sink     = mean(sign_category == "Stable sink",                 na.rm = TRUE),
    frac_uncertain       = mean(grepl("uncertain", sign_category),               na.rm = TRUE),
    sign_sensitivity_class = case_when(
      mean_p_source >= 0.95 ~ "Stable source",
      mean_p_source <= 0.05 ~ "Stable sink",
      mean_p_source >  0.50 ~ "Probable source (sign-uncertain)",
      TRUE                  ~ "Probable sink (sign-uncertain)"
    ),
    .groups = "drop"
  )

write.csv(sign_sensitivity_site, "OUTPUT/NEON_ERA5_rf_sign_sensitivity_by_site.csv", row.names = FALSE)

message(sprintf(
  "Sign sensitivity: %d/%d site-years have stable sign; %d sites sign-uncertain at site level.",
  sum(sign_sensitivity$sign_stable, na.rm = TRUE),
  nrow(sign_sensitivity),
  sum(grepl("uncertain", sign_sensitivity_site$sign_sensitivity_class))
))

era5_annual_budget <- era5_annual_budget_pre %>% arrange(SITE_ID, Year)

era5_mean_annual_budget <- era5_annual_budget %>%
  reframe(
    .by = SITE_ID,
    n_years = dplyr::n(),
    era5_source_years = sum(annual_budget_gC_m2_yr > 0, na.rm = TRUE),
    era5_prop_source_years = mean(annual_budget_gC_m2_yr > 0, na.rm = TRUE),
    mean_observed_coverage = mean(observed_coverage, na.rm = TRUE),
    n_training_obs = max(n_training_obs, na.rm = TRUE),
    mean_era5_gapfilled_annual_budget_gC_m2_yr = mean(annual_budget_gC_m2_yr, na.rm = TRUE),
    median_era5_gapfilled_annual_budget_gC_m2_yr = median(annual_budget_gC_m2_yr, na.rm = TRUE),
    sd_era5_gapfilled_annual_budget_gC_m2_yr = sd(annual_budget_gC_m2_yr, na.rm = TRUE),
    mean_era5_model_only_annual_budget_gC_m2_yr = mean(model_only_annual_budget_gC_m2_yr, na.rm = TRUE)
  ) %>%
  mutate(
    era5_annual_behavior = case_when(
      era5_prop_source_years >= 1 ~ "Weak-source",
      era5_prop_source_years <= 0 ~ "Weak-sink",
      TRUE ~ "Fluctuating"
    ),
    era5_annual_behavior = factor(era5_annual_behavior, levels = behavior_levels)
  )

model_standardized_budget <- read.csv(model_standardized_budget_file) %>%
  transmute(
    SITE_ID = as.character(SITE_ID),
    model_standardized_annual_budget_gC_m2_yr = annual_budget_mean_gC_m2_yr,
    model_standardized_lwr_gC_m2_yr = annual_budget_lwr_gC_m2_yr,
    model_standardized_upr_gC_m2_yr = annual_budget_upr_gC_m2_yr,
    model_standardized_annual_behavior = factor(annual_behavior, levels = behavior_levels),
    model_standardized_prob_annual_source = prob_annual_source
  )

budget_comparison <- era5_mean_annual_budget %>%
  inner_join(model_standardized_budget, by = "SITE_ID") %>%
  left_join(reference_annual_behavior, by = "SITE_ID") %>%
  mutate(
    budget_difference_era5_minus_model_standardized_gC_m2_yr =
      mean_era5_gapfilled_annual_budget_gC_m2_yr - model_standardized_annual_budget_gC_m2_yr,
    absolute_difference_gC_m2_yr = abs(budget_difference_era5_minus_model_standardized_gC_m2_yr),
    sign_agreement = sign(mean_era5_gapfilled_annual_budget_gC_m2_yr) == sign(model_standardized_annual_budget_gC_m2_yr),
    era5_behavior = case_when(
      mean_era5_gapfilled_annual_budget_gC_m2_yr < 0 ~ "Annual sink",
      TRUE ~ "Annual source"
    )
  ) %>%
  arrange(desc(absolute_difference_gC_m2_yr))

era5_annual_class_change_summary <- budget_comparison %>%
  count(reference_annual_behavior, era5_annual_behavior, name = "n_sites") %>%
  mutate(
    changed = reference_annual_behavior != era5_annual_behavior,
    reference_annual_behavior = factor(reference_annual_behavior, levels = behavior_levels),
    era5_annual_behavior = factor(era5_annual_behavior, levels = behavior_levels)
  )

class_levels <- behavior_levels
era5_annual_class_change_matrix <- expand_grid(
  reference_annual_behavior = factor(class_levels, levels = class_levels),
  era5_annual_behavior = factor(class_levels, levels = class_levels)
) %>%
  left_join(era5_annual_class_change_summary, by = c("reference_annual_behavior", "era5_annual_behavior")) %>%
  mutate(
    n_sites = replace_na(n_sites, 0L),
    changed = reference_annual_behavior != era5_annual_behavior,
    label = if_else(n_sites > 0, as.character(n_sites), "")
  )

site_standardized_flux <- read.csv(site_standardized_flux_file) %>%
  transmute(
    SITE_ID = as.character(SITE_ID),
    flux_30min_umolC_m2_s,
    flux_30min_sd_umolC_m2_s,
    flux_30min_se_umolC_m2_s,
    behavior_30min = factor(behavior_30min, levels = behavior_levels)
  )

daily_flux_summary <- read.csv(daily_flux_summary_file) %>%
  transmute(
    SITE_ID = as.character(SITE_ID),
    flux_daily_mgC_m2_day,
    flux_daily_sd_mgC_m2_day,
    flux_daily_se_mgC_m2_day,
    behavior_daily = factor(behavior_daily, levels = behavior_levels)
  )

scaled_annual_budget <- read.csv(scaled_annual_budget_file) %>%
  transmute(
    SITE_ID = as.character(SITE_ID),
    flux_scaled_annual_gC_m2_yr = flux_scaled_daily_annual_gC_m2_yr,
    flux_scaled_annual_sd_gC_m2_yr = flux_scaled_daily_annual_sd_gC_m2_yr,
    flux_scaled_annual_se_gC_m2_yr = flux_scaled_daily_annual_se_gC_m2_yr,
    scaled_annual_behavior = factor(behavior_annual_scaled, levels = behavior_levels)
  )

era5_site_order <- era5_mean_annual_budget %>%
  mutate(era5_annual_behavior = factor(era5_annual_behavior, levels = behavior_levels)) %>%
  arrange(era5_annual_behavior, mean_era5_gapfilled_annual_budget_gC_m2_yr, SITE_ID) %>%
  pull(SITE_ID)

era5_all_site_flux_magnitude_summary <- era5_mean_annual_budget %>%
  transmute(
    SITE_ID = as.character(SITE_ID),
    annual_behavior = factor(era5_annual_behavior, levels = behavior_levels),
    flux_era5_annual_gC_m2_yr = mean_era5_gapfilled_annual_budget_gC_m2_yr,
    flux_era5_annual_sd_gC_m2_yr = sd_era5_gapfilled_annual_budget_gC_m2_yr,
    flux_era5_annual_se_gC_m2_yr = sd_era5_gapfilled_annual_budget_gC_m2_yr / sqrt(n_years)
  ) %>%
  left_join(site_standardized_flux, by = "SITE_ID") %>%
  left_join(daily_flux_summary, by = "SITE_ID") %>%
  left_join(scaled_annual_budget, by = "SITE_ID") %>%
  left_join(site_metadata %>% dplyr::select(SITE_ID, EcoType), by = "SITE_ID") %>%
  transmute(
    SITE_ID,
    EcoType,
    annual_behavior,
    `30 min` = flux_30min_umolC_m2_s,
    Daily = flux_daily_mgC_m2_day,
    `Annual` = flux_scaled_annual_gC_m2_yr,
    `Annual ERA5` = flux_era5_annual_gC_m2_yr,
    sd_30min = flux_30min_sd_umolC_m2_s,
    sd_daily = flux_daily_sd_mgC_m2_day,
    sd_annual_scaled = flux_scaled_annual_sd_gC_m2_yr,
    sd_annual_era5 = flux_era5_annual_sd_gC_m2_yr,
    se_30min = flux_30min_se_umolC_m2_s,
    se_daily = flux_daily_se_mgC_m2_day,
    se_annual_scaled = flux_scaled_annual_se_gC_m2_yr,
    se_annual_era5 = flux_era5_annual_se_gC_m2_yr,
    behavior_30min,
    behavior_daily,
    scaled_annual_behavior
  ) %>%
  pivot_longer(
    cols = c(`30 min`, Daily, `Annual`, `Annual ERA5`),
    names_to = "scale",
    values_to = "flux_native"
  ) %>%
  mutate(
    flux_sd_native = case_when(
      scale == "30 min"      ~ sd_30min,
      scale == "Daily"       ~ sd_daily,
      scale == "Annual"      ~ sd_annual_scaled,
      scale == "Annual ERA5" ~ sd_annual_era5
    ),
    flux_se_native = case_when(
      scale == "30 min"      ~ se_30min,
      scale == "Daily"       ~ se_daily,
      scale == "Annual"      ~ se_annual_scaled,
      scale == "Annual ERA5" ~ se_annual_era5
    ),
    flux_unit = case_when(
      scale == "30 min"      ~ "umol C m-2 s-1",
      scale == "Daily"       ~ "mg C m-2 d-1",
      scale == "Annual"      ~ "g C m-2 yr-1",
      scale == "Annual ERA5" ~ "g C m-2 yr-1"
    ),
    flux_lower_native = flux_native - flux_sd_native,
    flux_upper_native = flux_native + flux_sd_native,
    scale = factor(scale, levels = c("30 min", "Daily", "Annual", "Annual ERA5")),
    SITE_ID_plot = factor(SITE_ID, levels = rev(era5_site_order)),
    scale_label = factor(
      paste0(scale, "\n", flux_unit),
      levels = paste0(
        c("30 min", "Daily", "Annual", "Annual ERA5"),
        "\n",
        c("umol C m-2 s-1", "mg C m-2 d-1", "g C m-2 yr-1", "g C m-2 yr-1")
      )
    )
  ) %>%
  filter(is.finite(flux_native), !is.na(annual_behavior))

era5_annual_site_map_data <- era5_mean_annual_budget %>%
  transmute(
    SITE_ID = as.character(SITE_ID),
    annual_behavior = factor(era5_annual_behavior, levels = behavior_levels),
    annual_flux_gC_m2_yr = mean_era5_gapfilled_annual_budget_gC_m2_yr,
    annual_flux_magnitude_gC_m2_yr = abs(mean_era5_gapfilled_annual_budget_gC_m2_yr)
  ) %>%
  left_join(site_metadata, by = "SITE_ID") %>%
  filter(
    is.finite(latitude),
    is.finite(longitude),
    is.finite(annual_flux_gC_m2_yr),
    !is.na(annual_behavior)
  )

era5_annual_method_flux_summary <- budget_comparison %>%
  transmute(
    SITE_ID = as.character(SITE_ID),
    annual_behavior = factor(era5_annual_behavior, levels = behavior_levels),
    `Scaled annual` = reference_annual_budget_gC_m2_yr,
    `ERA5 annual` = mean_era5_gapfilled_annual_budget_gC_m2_yr
  ) %>%
  pivot_longer(
    cols = c(`Scaled annual`, `ERA5 annual`),
    names_to = "annual_method",
    values_to = "annual_flux_gC_m2_yr"
  ) %>%
  mutate(annual_method = factor(annual_method, levels = c("Scaled annual", "ERA5 annual"))) %>%
  filter(is.finite(annual_flux_gC_m2_yr), !is.na(annual_behavior))

era5_site_diel_30min <- ch4_30min %>%
  left_join(
    era5_mean_annual_budget %>%
      transmute(SITE_ID = as.character(SITE_ID), annual_behavior = factor(era5_annual_behavior, levels = behavior_levels)),
    by = "SITE_ID"
  ) %>%
  left_join(site_metadata %>% dplyr::select(SITE_ID, EcoType_metadata = EcoType), by = "SITE_ID") %>%
  mutate(EcoType = coalesce(as.character(EcoType), EcoType_metadata)) %>%
  reframe(
    .by = c(SITE_ID, annual_behavior, EcoType, hour_num),
    n_30min = n(),
    flux_umolC_m2_s = mg_c_30min_to_umol_c_s(mean(CH4_mgC_30min, na.rm = TRUE)),
    source_probability = mean(CH4_mgC_30min > 0, na.rm = TRUE)
  ) %>%
  filter(!is.na(annual_behavior), is.finite(hour_num))

era5_diel_behavior_summary <- era5_site_diel_30min %>%
  reframe(
    .by = c(annual_behavior, hour_num),
    n_sites = n_distinct(SITE_ID),
    mean_flux_umolC_m2_s = mean(flux_umolC_m2_s, na.rm = TRUE),
    median_flux_umolC_m2_s = median(flux_umolC_m2_s, na.rm = TRUE),
    sd_flux_umolC_m2_s = sd(flux_umolC_m2_s, na.rm = TRUE),
    se_flux_umolC_m2_s = sd_flux_umolC_m2_s / sqrt(n_sites),
    mean_source_probability = mean(source_probability, na.rm = TRUE),
    sd_source_probability = sd(source_probability, na.rm = TRUE),
    se_source_probability = sd_source_probability / sqrt(n_sites)
  ) %>%
  mutate(
    se_flux_umolC_m2_s = replace_na(se_flux_umolC_m2_s, 0),
    se_source_probability = replace_na(se_source_probability, 0),
    annual_behavior = factor(annual_behavior, levels = behavior_levels)
  )

# Seasonal (monthly) source probability by behavior class
era5_site_seasonal_30min <- ch4_30min %>%
  left_join(
    era5_mean_annual_budget %>%
      transmute(SITE_ID = as.character(SITE_ID),
                annual_behavior = factor(era5_annual_behavior, levels = behavior_levels)),
    by = "SITE_ID"
  ) %>%
  filter(!is.na(annual_behavior), is.finite(month)) %>%
  reframe(
    .by = c(SITE_ID, annual_behavior, month),
    n_30min = n(),
    flux_umolC_m2_s = mg_c_30min_to_umol_c_s(mean(CH4_mgC_30min, na.rm = TRUE)),
    source_probability = mean(CH4_mgC_30min > 0, na.rm = TRUE)
  )

era5_seasonal_behavior_summary <- era5_site_seasonal_30min %>%
  reframe(
    .by = c(annual_behavior, month),
    n_sites = n_distinct(SITE_ID),
    mean_source_probability = mean(source_probability, na.rm = TRUE),
    sd_source_probability = sd(source_probability, na.rm = TRUE),
    se_source_probability = sd_source_probability / sqrt(n_sites)
  ) %>%
  mutate(
    se_source_probability = replace_na(se_source_probability, 0),
    annual_behavior = factor(annual_behavior, levels = behavior_levels)
  )

era5_behavior_site_counts <- era5_annual_site_map_data %>%
  transmute(
    SITE_ID,
    annual_behavior = factor(annual_behavior, levels = behavior_levels),
    EcoType = as.character(EcoType)
  ) %>%
  distinct()

era5_annual_behavior_counts <- era5_behavior_site_counts %>%
  count(annual_behavior, name = "n_sites") %>%
  complete(
    annual_behavior = factor(behavior_levels, levels = behavior_levels),
    fill = list(n_sites = 0L)
  )

era5_annual_behavior_ecotype_counts <- era5_behavior_site_counts %>%
  filter(!is.na(EcoType)) %>%
  count(annual_behavior, EcoType, name = "n_sites") %>%
  complete(
    annual_behavior = factor(behavior_levels, levels = behavior_levels),
    EcoType,
    fill = list(n_sites = 0L)
  )

write.csv(era5_annual_budget, "OUTPUT/NEON_ERA5_gapfilled_annual_budget_by_year.csv", row.names = FALSE)
write.csv(era5_mean_annual_budget, "OUTPUT/NEON_ERA5_gapfilled_mean_annual_budget.csv", row.names = FALSE)
write.csv(budget_comparison, "OUTPUT/NEON_ERA5_vs_model_standardized_budget_comparison.csv", row.names = FALSE)
write.csv(
  era5_all_site_flux_magnitude_summary,
  "OUTPUT/NEON_ERA5_all_site_flux_magnitude_summary.csv",
  row.names = FALSE
)
write.csv(
  era5_annual_site_map_data,
  "OUTPUT/NEON_ERA5_annual_site_map_data.csv",
  row.names = FALSE
)
write.csv(
  era5_annual_method_flux_summary,
  "OUTPUT/NEON_ERA5_scaled_vs_era5_annual_flux_by_site.csv",
  row.names = FALSE
)
write.csv(
  era5_diel_behavior_summary,
  "OUTPUT/NEON_ERA5_diel_behavior_summary.csv",
  row.names = FALSE
)
write.csv(
  era5_seasonal_behavior_summary,
  "OUTPUT/NEON_ERA5_seasonal_behavior_summary.csv",
  row.names = FALSE
)
write.csv(
  era5_annual_behavior_counts,
  "OUTPUT/NEON_ERA5_annual_behavior_site_counts.csv",
  row.names = FALSE
)
write.csv(
  era5_annual_behavior_ecotype_counts,
  "OUTPUT/NEON_ERA5_annual_behavior_ecotype_counts.csv",
  row.names = FALSE
)
write.csv(
  era5_annual_class_change_summary,
  "OUTPUT/NEON_ERA5_reference_annual_class_changes.csv",
  row.names = FALSE
)

comparison_cor <- suppressWarnings(cor(
  budget_comparison$mean_era5_gapfilled_annual_budget_gC_m2_yr,
  budget_comparison$model_standardized_annual_budget_gC_m2_yr,
  method = "spearman",
  use = "complete.obs"
))

comparison_rmse <- sqrt(mean(budget_comparison$budget_difference_era5_minus_model_standardized_gC_m2_yr^2, na.rm = TRUE))
comparison_bias <- mean(budget_comparison$budget_difference_era5_minus_model_standardized_gC_m2_yr, na.rm = TRUE)

# Budget comparison scatter: prominent 1:1 line, shaded sign-disagreement quadrants,
# label sites where methods disagree on sign.
budget_comparison <- budget_comparison %>%
  mutate(
    sign_agree = sign(mean_era5_gapfilled_annual_budget_gC_m2_yr) ==
                 sign(model_standardized_annual_budget_gC_m2_yr)
  )









top_difference_lines <- budget_comparison %>%
  slice_head(n = 10) %>%
  mutate(
    line = paste0(
      "- ", SITE_ID, ": ERA5 = ",
      signif(mean_era5_gapfilled_annual_budget_gC_m2_yr, 3),
      ", model-standardized = ",
      signif(model_standardized_annual_budget_gC_m2_yr, 3),
      ", difference = ",
      signif(budget_difference_era5_minus_model_standardized_gC_m2_yr, 3),
      " g C m-2 yr-1"
    )
  ) %>%
  pull(line)






era5_annual_behavior_counts <- budget_comparison %>%
  count(era5_annual_behavior, name = "n_sites") %>%
  mutate(line = paste0("- ", era5_annual_behavior, ": ", n_sites, " sites")) %>%
  pull(line)

era5_annual_class_change_lines <- era5_annual_class_change_summary %>%
  arrange(reference_annual_behavior, era5_annual_behavior) %>%
  mutate(
    line = paste0(
      "- Reference annual ", reference_annual_behavior, " -> ERA5 annual ", era5_annual_behavior,
      ": ", n_sites, " sites",
      if_else(changed, " changed", " unchanged")
    )
  ) %>%
  pull(line)

writeLines(
  c(
    "# NEON ERA5 Half-Hourly Gapfilled Annual Budgets",
    "",
    "## ERA5 Source",
    "- Hourly ERA5 point covariates were requested from the Open-Meteo Archive API for each NEON tower coordinate.",
    "- Variables: 2 m air temperature, 2 m relative humidity, 0-7 cm soil temperature, 0-7 cm and 7-28 cm soil moisture, precipitation (mm/hr), shortwave radiation.",
    "- Hourly values were linearly interpolated to 30-minute timestamps (precipitation: floor-hour value).",
    "",
    "## Model",
    paste0("- Training observations: ", nrow(training_rf), "."),
    paste0("- Site-level RF OOB R²: ", round(era5_gapfill_rf_site$r.squared, 3),
           " (mtry = ", best_mtry_site, "); OOB RMSE: ",
           signif(fit_metrics$rmse_site_mgC_m2_30min, 3), " mg C m-2 30 min-1."),
    paste0("- Population RF OOB R²: ", round(era5_gapfill_rf_pop$r.squared, 3),
           " (mtry = ", best_mtry_pop, "); OOB RMSE: ",
           signif(fit_metrics$rmse_pop_mgC_m2_30min, 3), " mg C m-2 30 min-1."),
    paste0("- 30-min OOB sign accuracy (site RF): ",
           round(100 * fit_metrics$sign_accuracy_site_30min, 1), "%."),
    paste0("- OOB annual budget sign accuracy (site RF): ",
           round(100 * oob_sign_summary$sign_accuracy_site, 1), "% of site-years."),
    paste0("- Sign-stable site-years (95% CI does not straddle zero): ",
           sum(sign_sensitivity$sign_stable, na.rm = TRUE), " of ", nrow(sign_sensitivity), "."),
    "- Response: total CH4 flux in mg C m-2 30 min-1.",
    paste0("- Predictors (", length(rf_predictors_pop), " population): ERA5 air/soil temperature, surface/deep soil moisture,",
           " precipitation, shortwave radiation, VPD; 7-day rolling mean Tair and cumulative precip;",
           " Tair×VSWC interaction; climatological site×time-of-day mean (cascade: month-hour → season-hour → hour → global);",
           " cyclic hour (sin/cos), DOY, month, season, ecosystem type.",
           " Site-level RF additionally includes SITE_ID."),
    "- Annual budgets retain observed half-hour fluxes where available and fill missing half-hours with best-model predictions.",
    "",
    "## Comparison To LUT (Balanced Site-Month-Hour) Annual Budget",
    paste0("- Sites compared: ", nrow(budget_comparison), "."),
    paste0("- Spearman correlation: ", signif(comparison_cor, 3), "."),
    paste0("- Mean ERA5-minus-LUT difference: ", signif(comparison_bias, 3), " g C m-2 yr-1."),
    paste0("- RMSE: ", signif(comparison_rmse, 3), " g C m-2 yr-1."),
    paste0("- Sign agreement: ", sum(budget_comparison$sign_agreement, na.rm = TRUE), " of ", nrow(budget_comparison), " sites."),
    "",
    "## Reference Annual vs ERA5 Annual-Budget Class",
    "- Reference annual classes are from lookup-filled daily annual budgets produced by `flow.30min.analysis.R`.",
    "- ERA5 annual class is based on the fraction of gapfilled site-years with a positive annual budget.",
    "- Weak-source: 100% of years positive; Weak-sink: 0% positive years; otherwise Fluctuating.",
    "### ERA5 Annual-Budget Behavior Counts",
    era5_annual_behavior_counts,
    "### Reference-vs-ERA5 Class Changes",
    era5_annual_class_change_lines,
    "",
    "## Largest Absolute Differences",
    top_difference_lines,
    "",
    "## Outputs",
    "- `OUTPUT/NEON_ERA5_hourly_site_covariates.csv.gz`",
    "- `OUTPUT/NEON_ERA5_30min_site_covariates.csv.gz`",
    "- `OUTPUT/NEON_ERA5_gapfilled_30min.csv.gz`",
    "- `OUTPUT/NEON_ERA5_gapfilled_annual_budget_by_year.csv`",
    "- `OUTPUT/NEON_ERA5_gapfilled_mean_annual_budget.csv`",
    "- `OUTPUT/NEON_ERA5_vs_model_standardized_budget_comparison.csv`",
    "- `OUTPUT/NEON_ERA5_all_site_flux_magnitude_summary.csv`",
    "- `OUTPUT/NEON_ERA5_annual_site_map_data.csv`",
    "- `OUTPUT/NEON_ERA5_scaled_vs_era5_annual_flux_by_site.csv`",
    "- `OUTPUT/NEON_ERA5_diel_behavior_summary.csv`",
    "- `OUTPUT/NEON_ERA5_annual_behavior_site_counts.csv`",
    "- `OUTPUT/NEON_ERA5_annual_behavior_ecotype_counts.csv`",
    "- `OUTPUT/NEON_ERA5_reference_annual_class_changes.csv`",
    "- `OUTPUT/NEON_ERA5_halfhour_gapfill_model_summary.txt`",
    "- `OUTPUT/NEON_ERA5_halfhour_gapfill_fit_metrics.csv`",
    "- `OUTPUT/NEON_ERA5_halfhour_gapfill_model_effects.csv`",
    "- `OUTPUT/NEON_ERA5_rf_oob_sign_accuracy_by_site.csv`",
    "- `OUTPUT/NEON_ERA5_rf_oob_annual_budget_sign.csv`",
    "- `OUTPUT/NEON_ERA5_rf_oob_sign_summary.csv`",
    "- `OUTPUT/NEON_ERA5_rf_sign_sensitivity.csv`",
    "- `OUTPUT/NEON_ERA5_rf_sign_sensitivity_by_site.csv`",
    "- `FIGURES/NEON_ERA5_halfhour_gapfill_model_fit_summary.png`",
    "- `FIGURES/NEON_ERA5_vs_model_standardized_budget_scatter.png`",
    "- `FIGURES/NEON_ERA5_vs_model_standardized_budget_differences.png`",
    "- `FIGURES/NEON_ERA5_reference_annual_class_changes.png`",
    "- `FIGURES/NEON_ERA5_all_site_category_flux_magnitudes.png`",
    "- `FIGURES/NEON_ERA5_annual_site_category_map.png`",
    "- `FIGURES/NEON_ERA5_flux_pattern_diel_behavior_panel.png`"
  ),
  "OUTPUT/NEON_ERA5_halfhour_gapfill_results.md"
)

message("Wrote ERA5 half-hourly gapfill and budget comparison outputs.")
