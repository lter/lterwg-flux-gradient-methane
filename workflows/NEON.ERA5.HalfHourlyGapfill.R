# ERA5-driven half-hourly total-flux CH4 gap filling for NEON sites.
#
# This workflow obtains hourly ERA5 point covariates for NEON tower sites,
# interpolates them to 30-minute timestamps, fits a total-flux GAM using ERA5
# temperature and soil volumetric water content, and compares the resulting
# annual budgets to the model-standardized 30-minute annual budget from
# NEON.30min.Gapfill.r and to the annual categories from flow.30min.analysis.R.

library(tidyverse)
library(ggplot2)
library(mgcv)
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
model_standardized_budget_file <- "OUTPUT/NON_30min_gapfill_annual_budgets.csv"
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

behavior_levels <- c("Consistent sink", "Fluctuating", "Consistent source")
# Color convention: blue = sink (uptake), grey = fluctuating, red = source (emission)
behavior_colors <- c(
  "Consistent sink"   = "#2166AC",
  "Fluctuating"       = "#4D4D4D",
  "Consistent source" = "#B2182B"
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
  query <- list(
    latitude = latitude,
    longitude = longitude,
    start_date = as.character(start_date),
    end_date = as.character(end_date),
    hourly = "temperature_2m,soil_moisture_0_to_7cm",
    timezone = "UTC",
    models = "era5"
  )
  url <- paste0(
    "https://archive-api.open-meteo.com/v1/archive?",
    paste(paste0(names(query), "=", utils::URLencode(as.character(query), reserved = TRUE)), collapse = "&")
  )

  parsed <- jsonlite::fromJSON(url)

  if (is.null(parsed$hourly$time) || is.null(parsed$hourly$temperature_2m) ||
      is.null(parsed$hourly$soil_moisture_0_to_7cm)) {
    stop("Open-Meteo response did not include expected ERA5 hourly variables for ", site)
  }

  tibble(
    SITE_ID = site,
    time_hour = as.POSIXct(parsed$hourly$time, format = "%Y-%m-%dT%H:%M", tz = "UTC"),
    ERA5_Tair_C = as.numeric(parsed$hourly$temperature_2m),
    ERA5_VSWC = as.numeric(parsed$hourly$soil_moisture_0_to_7cm),
    ERA5_latitude = parsed$latitude,
    ERA5_longitude = parsed$longitude,
    ERA5_elevation = parsed$elevation %||% NA_real_
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

if (!file.exists(era5_hourly_file)) {
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
} else {
  message("Using cached ERA5 hourly covariates: ", era5_hourly_file)
  era5_hourly <- data.table::fread(era5_hourly_file) %>%
    as_tibble() %>%
    mutate(time_hour = as.POSIXct(time_hour, tz = "UTC"))
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

if (!file.exists(era5_halfhour_file)) {
  message("Interpolating hourly ERA5 covariates to 30-minute site-year grid.")
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
        dplyr::select(SITE_ID, time_floor_hour = time_hour, ERA5_Tair_C_floor = ERA5_Tair_C, ERA5_VSWC_floor = ERA5_VSWC),
      by = c("SITE_ID", "time_floor_hour")
    ) %>%
    left_join(
      era5_hourly %>%
        dplyr::select(SITE_ID, time_ceiling_hour = time_hour, ERA5_Tair_C_ceiling = ERA5_Tair_C, ERA5_VSWC_ceiling = ERA5_VSWC),
      by = c("SITE_ID", "time_ceiling_hour")
    ) %>%
    mutate(
      ERA5_Tair_C = if_else(
        is.finite(ERA5_Tair_C_floor) & is.finite(ERA5_Tair_C_ceiling),
        ERA5_Tair_C_floor + hour_weight * (ERA5_Tair_C_ceiling - ERA5_Tair_C_floor),
        coalesce(ERA5_Tair_C_floor, ERA5_Tair_C_ceiling)
      ),
      ERA5_VSWC = if_else(
        is.finite(ERA5_VSWC_floor) & is.finite(ERA5_VSWC_ceiling),
        ERA5_VSWC_floor + hour_weight * (ERA5_VSWC_ceiling - ERA5_VSWC_floor),
        coalesce(ERA5_VSWC_floor, ERA5_VSWC_ceiling)
      )
    ) %>%
    dplyr::select(SITE_ID, Year, time.rounded, Date, month, doy, hour_num, season, EcoType, ERA5_Tair_C, ERA5_VSWC)

  data.table::fwrite(era5_halfhour, era5_halfhour_file)
} else {
  message("Using cached ERA5 half-hour covariates: ", era5_halfhour_file)
  era5_halfhour <- data.table::fread(era5_halfhour_file) %>%
    as_tibble() %>%
    mutate(
      time.rounded = as.POSIXct(time.rounded, tz = "UTC"),
      Date = as.Date(Date),
      season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
      EcoType = factor(EcoType)
    )
}

gapfill_data <- era5_halfhour %>%
  left_join(
    ch4_30min %>%
      dplyr::select(SITE_ID, time.rounded, CH4_mgC_30min),
    by = c("SITE_ID", "time.rounded")
  ) %>%
  mutate(
    observed_flux = is.finite(CH4_mgC_30min),
    sin_hour = sin(2 * pi * hour_num / 24),
    cos_hour = cos(2 * pi * hour_num / 24),
    SITE_ID = factor(SITE_ID),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    EcoType = factor(EcoType)
  )

training_data <- gapfill_data %>%
  filter(
    observed_flux,
    is.finite(CH4_mgC_30min),
    is.finite(ERA5_Tair_C),
    is.finite(ERA5_VSWC),
    is.finite(hour_num),
    is.finite(doy),
    !is.na(season),
    !is.na(EcoType)
  ) %>%
  droplevels()

if (nrow(training_data) < 500) {
  stop("Too few observations with ERA5 covariates to fit model: ", nrow(training_data))
}

site_training_counts <- training_data %>%
  count(SITE_ID, name = "n_training_obs")

era5_gapfill_model <- mgcv::bam(
  CH4_mgC_30min ~
    s(hour_num, bs = "cc", k = 12) +
    s(doy, bs = "cc", k = 20) +
    s(ERA5_Tair_C, k = 10) +
    s(ERA5_VSWC, k = 10) +
    ti(ERA5_Tair_C, ERA5_VSWC, k = c(5, 5)) +
    season +
    EcoType +
    s(SITE_ID, bs = "re"),
  data = training_data,
  family = gaussian(),
  method = "fREML",
  discrete = FALSE,
  knots = list(hour_num = c(0, 24), doy = c(0.5, 366.5))
)

capture.output(
  {
    cat("NEON ERA5-driven half-hourly total-flux gapfill model\n")
    cat("Training observations:", nrow(training_data), "\n")
    cat("ERA5 hourly file:", era5_hourly_file, "\n")
    cat("ERA5 half-hour file:", era5_halfhour_file, "\n\n")
    print(summary(era5_gapfill_model))
  },
  file = "OUTPUT/NEON_ERA5_halfhour_gapfill_model_summary.txt"
)

training_fit <- training_data %>%
  mutate(
    fitted_CH4_mgC_30min = as.numeric(predict(era5_gapfill_model, type = "response")),
    residual_CH4_mgC_30min = CH4_mgC_30min - fitted_CH4_mgC_30min
  )

fit_metrics <- training_fit %>%
  summarise(
    n_training = dplyr::n(),
    rmse_mgC_m2_30min = sqrt(mean(residual_CH4_mgC_30min^2, na.rm = TRUE)),
    mae_mgC_m2_30min = mean(abs(residual_CH4_mgC_30min), na.rm = TRUE),
    bias_mgC_m2_30min = mean(residual_CH4_mgC_30min, na.rm = TRUE),
    observed_sd_mgC_m2_30min = sd(CH4_mgC_30min, na.rm = TRUE),
    fitted_sd_mgC_m2_30min = sd(fitted_CH4_mgC_30min, na.rm = TRUE),
    correlation_observed_fitted = suppressWarnings(cor(CH4_mgC_30min, fitted_CH4_mgC_30min, use = "complete.obs"))
  )

write.csv(fit_metrics, "OUTPUT/NEON_ERA5_halfhour_gapfill_fit_metrics.csv", row.names = FALSE)

set.seed(20260519)
fit_plot_data <- training_fit %>%
  slice_sample(n = min(10000, nrow(training_fit)))

fit_axis_lim <- quantile(
  abs(c(fit_plot_data$CH4_mgC_30min, fit_plot_data$fitted_CH4_mgC_30min)),
  0.995,
  na.rm = TRUE
)

plot_fit_observed <- fit_plot_data %>%
  ggplot(aes(x = fitted_CH4_mgC_30min, y = CH4_mgC_30min)) +
  geom_bin2d(bins = 70) +
  geom_abline(slope = 1, intercept = 0, color = "white", linewidth = 0.9) +
  geom_hline(yintercept = 0, color = "grey75", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey75", linetype = "dashed", linewidth = 0.5) +
  coord_cartesian(xlim = c(-fit_axis_lim, fit_axis_lim), ylim = c(-fit_axis_lim, fit_axis_lim)) +
  scale_fill_viridis_c(option = "magma", trans = "sqrt", name = "Count") +
  labs(
    title = "A. Observed vs fitted",
    subtitle = paste0(
      "RMSE = ", signif(fit_metrics$rmse_mgC_m2_30min, 3),
      "; r = ", signif(fit_metrics$correlation_observed_fitted, 3)
    ),
    x = "Fitted CH4 flux (mg C m-2 30 min-1)",
    y = "Observed CH4 flux (mg C m-2 30 min-1)"
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

plot_fit_residuals <- fit_plot_data %>%
  ggplot(aes(x = fitted_CH4_mgC_30min, y = residual_CH4_mgC_30min)) +
  geom_bin2d(bins = 70) +
  geom_hline(yintercept = 0, color = "white", linewidth = 0.9) +
  coord_cartesian(
    xlim = c(-fit_axis_lim, fit_axis_lim),
    ylim = quantile(fit_plot_data$residual_CH4_mgC_30min, c(0.005, 0.995), na.rm = TRUE)
  ) +
  scale_fill_viridis_c(option = "magma", trans = "sqrt", name = "Count") +
  labs(
    title = "B. Residuals vs fitted",
    subtitle = paste0("Bias = ", signif(fit_metrics$bias_mgC_m2_30min, 3),
                      "; MAE = ", signif(fit_metrics$mae_mgC_m2_30min, 3)),
    x = "Fitted CH4 flux (mg C m-2 30 min-1)",
    y = "Residual (observed - fitted)"
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

reference_values <- training_data %>%
  summarise(
    hour_num = 12,
    doy = 182,
    ERA5_Tair_C = median(ERA5_Tair_C, na.rm = TRUE),
    ERA5_VSWC = median(ERA5_VSWC, na.rm = TRUE),
    season = names(sort(table(season), decreasing = TRUE))[1],
    EcoType = names(sort(table(EcoType), decreasing = TRUE))[1],
    SITE_ID = levels(SITE_ID)[1]
  )

make_effect_grid <- function(driver, values) {
  tibble(
    driver = driver,
    driver_value = values,
    hour_num = reference_values$hour_num,
    doy = reference_values$doy,
    ERA5_Tair_C = reference_values$ERA5_Tair_C,
    ERA5_VSWC = reference_values$ERA5_VSWC,
    season = reference_values$season,
    EcoType = reference_values$EcoType,
    SITE_ID = reference_values$SITE_ID
  ) %>%
    mutate(
      hour_num = if_else(driver == "hour_num", driver_value, hour_num),
      doy = if_else(driver == "doy", driver_value, doy),
      ERA5_Tair_C = if_else(driver == "ERA5_Tair_C", driver_value, ERA5_Tair_C),
      ERA5_VSWC = if_else(driver == "ERA5_VSWC", driver_value, ERA5_VSWC),
      season = factor(season, levels = levels(training_data$season)),
      EcoType = factor(EcoType, levels = levels(training_data$EcoType)),
      SITE_ID = factor(SITE_ID, levels = levels(training_data$SITE_ID))
    )
}

effect_grid <- bind_rows(
  make_effect_grid("hour_num", seq(0, 23.5, by = 0.5)),
  make_effect_grid("doy", seq(1, 366, length.out = 160)),
  make_effect_grid(
    "ERA5_Tair_C",
    seq(quantile(training_data$ERA5_Tair_C, 0.02, na.rm = TRUE),
        quantile(training_data$ERA5_Tair_C, 0.98, na.rm = TRUE),
        length.out = 160)
  ),
  make_effect_grid(
    "ERA5_VSWC",
    seq(quantile(training_data$ERA5_VSWC, 0.02, na.rm = TRUE),
        quantile(training_data$ERA5_VSWC, 0.98, na.rm = TRUE),
        length.out = 160)
  )
)

effect_pred <- predict(
  era5_gapfill_model,
  newdata = effect_grid,
  type = "response",
  se.fit = TRUE,
  exclude = "s(SITE_ID)"
)

effect_grid <- effect_grid %>%
  mutate(
    pred_CH4_mgC_30min = as.numeric(effect_pred$fit),
    se_CH4_mgC_30min = as.numeric(effect_pred$se.fit),
    lower_CH4_mgC_30min = pred_CH4_mgC_30min - 1.64 * se_CH4_mgC_30min,
    upper_CH4_mgC_30min = pred_CH4_mgC_30min + 1.64 * se_CH4_mgC_30min,
    driver_label = recode(
      driver,
      hour_num = "Hour of day",
      doy = "Day of year",
      ERA5_Tair_C = "ERA5 air temperature (C)",
      ERA5_VSWC = "ERA5 soil moisture (0-7 cm)"
    )
  )

write.csv(effect_grid, "OUTPUT/NEON_ERA5_halfhour_gapfill_model_effects.csv", row.names = FALSE)

plot_model_effects <- effect_grid %>%
  ggplot(aes(x = driver_value, y = pred_CH4_mgC_30min)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.5) +
  geom_ribbon(aes(ymin = lower_CH4_mgC_30min, ymax = upper_CH4_mgC_30min),
              fill = "grey35", alpha = 0.18) +
  geom_line(color = "grey10", linewidth = 1.0) +
  facet_wrap(~driver_label, scales = "free_x", ncol = 2) +
  labs(
    title = "C. Population-level model effects",
    subtitle = "Predictions hold other covariates at typical values and exclude the site random effect; ribbons are +/- 90% model SE.",
    x = NULL,
    y = "Predicted CH4 flux (mg C m-2 30 min-1)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8, color = "grey35"),
    strip.background = element_rect(fill = "grey94", color = "grey40"),
    strip.text = element_text(face = "bold", size = 8.5),
    panel.grid.minor = element_blank()
  )

plot_residual_distribution <- fit_plot_data %>%
  ggplot(aes(x = residual_CH4_mgC_30min)) +
  geom_vline(xintercept = 0, color = "grey35", linetype = "dashed") +
  geom_histogram(bins = 80, fill = "#5E81AC", color = "white", linewidth = 0.15) +
  coord_cartesian(xlim = quantile(fit_plot_data$residual_CH4_mgC_30min, c(0.005, 0.995), na.rm = TRUE)) +
  labs(
    title = "D. Residual distribution",
    x = "Residual (mg C m-2 30 min-1)",
    y = "Training observations"
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

plot_era5_model_fit_summary <- (plot_fit_observed | plot_fit_residuals) /
  (plot_model_effects | plot_residual_distribution) +
  plot_layout(heights = c(1, 1.25), guides = "collect") +
  plot_annotation(
    title = "ERA5 Half-Hourly Gapfill GAM Fit Summary",
    subtitle = paste0(
      "Training rows = ", fit_metrics$n_training,
      "; deviance explained = ", signif(summary(era5_gapfill_model)$dev.expl, 3),
      "; adjusted R2 = ", signif(summary(era5_gapfill_model)$r.sq, 3)
    ),
    caption = "Observed/fitted diagnostics include the site random effect. Model-effect panels exclude the site random effect to show population-level responses."
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10, color = "grey35"),
    plot.caption = element_text(size = 8, color = "grey35"),
    legend.position = "bottom"
  )

ggsave(
  "FIGURES/NEON_ERA5_halfhour_gapfill_model_fit_summary.png",
  plot_era5_model_fit_summary,
  width = 13,
  height = 10,
  units = "in",
  dpi = 300
)


model_site_levels <- levels(training_data$SITE_ID)

prediction_data <- gapfill_data %>%
  mutate(SITE_ID_original = as.character(SITE_ID)) %>%
  left_join(site_training_counts %>% mutate(SITE_ID = as.character(SITE_ID)), by = "SITE_ID") %>%
  mutate(
    n_training_obs = replace_na(n_training_obs, 0L),
    use_site_effect = SITE_ID_original %in% model_site_levels & n_training_obs >= min_training_obs_per_site,
    SITE_ID = factor(
      if_else(SITE_ID_original %in% model_site_levels, SITE_ID_original, model_site_levels[1]),
      levels = model_site_levels
    )
  )

prediction_population <- predict(
  era5_gapfill_model,
  newdata = prediction_data,
  type = "response",
  exclude = "s(SITE_ID)"
)

prediction_site <- rep(NA_real_, nrow(prediction_data))
site_effect_rows <- which(prediction_data$use_site_effect)

if (length(site_effect_rows) > 0) {
  prediction_site[site_effect_rows] <- predict(
    era5_gapfill_model,
    newdata = prediction_data[site_effect_rows, , drop = FALSE],
    type = "response"
  )
}

era5_gapfilled_30min <- prediction_data %>%
  mutate(
    pred_CH4_mgC_30min_population = as.numeric(prediction_population),
    pred_CH4_mgC_30min_site = as.numeric(prediction_site),
    pred_CH4_mgC_30min = if_else(use_site_effect, pred_CH4_mgC_30min_site, pred_CH4_mgC_30min_population),
    gapfilled_CH4_mgC_30min = if_else(observed_flux, CH4_mgC_30min, pred_CH4_mgC_30min),
    prediction_type = case_when(
      observed_flux ~ "observed",
      use_site_effect ~ "site_model_gapfill",
      TRUE ~ "population_model_gapfill"
    ),
    SITE_ID = SITE_ID_original
  ) %>%
  dplyr::select(
    SITE_ID, Year, time.rounded, Date, month, doy, hour_num, season, EcoType,
    ERA5_Tair_C, ERA5_VSWC, CH4_mgC_30min, observed_flux,
    pred_CH4_mgC_30min, gapfilled_CH4_mgC_30min, prediction_type, n_training_obs
  )

data.table::fwrite(era5_gapfilled_30min, "OUTPUT/NEON_ERA5_gapfilled_30min.csv.gz")

era5_annual_budget <- era5_gapfilled_30min %>%
  reframe(
    .by = c(SITE_ID, Year),
    n_30min = dplyr::n(),
    n_observed = sum(observed_flux, na.rm = TRUE),
    observed_coverage = mean(observed_flux, na.rm = TRUE),
    n_training_obs = max(n_training_obs, na.rm = TRUE),
    annual_budget_gC_m2_yr = sum(gapfilled_CH4_mgC_30min, na.rm = TRUE) / 1000,
    model_only_annual_budget_gC_m2_yr = sum(pred_CH4_mgC_30min, na.rm = TRUE) / 1000,
    observed_partial_gC_m2 = sum(CH4_mgC_30min[observed_flux], na.rm = TRUE) / 1000,
    gapfilled_partial_gC_m2 = sum(gapfilled_CH4_mgC_30min[!observed_flux], na.rm = TRUE) / 1000
  ) %>%
  arrange(SITE_ID, Year)

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
      era5_prop_source_years >= 0.75 ~ "Consistent source",
      era5_prop_source_years <= 0.25 ~ "Consistent sink",
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
    `Annual scaled` = flux_scaled_annual_gC_m2_yr,
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
    cols = c(`30 min`, Daily, `Annual scaled`, `Annual ERA5`),
    names_to = "scale",
    values_to = "flux_native"
  ) %>%
  mutate(
    flux_sd_native = case_when(
      scale == "30 min" ~ sd_30min,
      scale == "Daily" ~ sd_daily,
      scale == "Annual scaled" ~ sd_annual_scaled,
      scale == "Annual ERA5" ~ sd_annual_era5
    ),
    flux_se_native = case_when(
      scale == "30 min" ~ se_30min,
      scale == "Daily" ~ se_daily,
      scale == "Annual scaled" ~ se_annual_scaled,
      scale == "Annual ERA5" ~ se_annual_era5
    ),
    flux_unit = case_when(
      scale == "30 min" ~ "umol C m-2 s-1",
      scale == "Daily" ~ "mg C m-2 d-1",
      scale == "Annual scaled" ~ "g C m-2 yr-1",
      scale == "Annual ERA5" ~ "g C m-2 yr-1"
    ),
    flux_lower_native = flux_native - flux_sd_native,
    flux_upper_native = flux_native + flux_sd_native,
    scale = factor(scale, levels = c("30 min", "Daily", "Annual scaled", "Annual ERA5")),
    SITE_ID_plot = factor(SITE_ID, levels = rev(era5_site_order)),
    scale_label = factor(
      paste0(scale, "\n", flux_unit),
      levels = paste0(
        c("30 min", "Daily", "Annual scaled", "Annual ERA5"),
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

axis_lim <- max(
  abs(c(budget_comparison$model_standardized_annual_budget_gC_m2_yr,
        budget_comparison$mean_era5_gapfilled_annual_budget_gC_m2_yr,
        budget_comparison$model_standardized_lwr_gC_m2_yr,
        budget_comparison$model_standardized_upr_gC_m2_yr)),
  na.rm = TRUE
) * 1.12

plot_budget_comparison <- budget_comparison %>%
  ggplot(aes(x = model_standardized_annual_budget_gC_m2_yr,
             y = mean_era5_gapfilled_annual_budget_gC_m2_yr,
             color = reference_annual_behavior)) +
  # sign-disagreement quadrant shading
  annotate("rect", xmin = -axis_lim, xmax = 0, ymin = 0, ymax =  axis_lim,
           fill = "#ff7043", alpha = 0.08) +
  annotate("rect", xmin =  0, xmax =  axis_lim, ymin = -axis_lim, ymax = 0,
           fill = "#ff7043", alpha = 0.08) +
  annotate("text", x = -axis_lim * 0.55, y =  axis_lim * 0.75,
           label = "Model-standardized: sink\nERA5: source", size = 3.0,
           color = "#c62828", fontface = "italic", hjust = 0.5) +
  annotate("text", x =  axis_lim * 0.55, y = -axis_lim * 0.75,
           label = "Model-standardized: source\nERA5: sink", size = 3.0,
           color = "#c62828", fontface = "italic", hjust = 0.5) +
  geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.8) +
  # prominent 1:1 reference line
  geom_abline(slope = 1, intercept = 0, color = "grey20", linewidth = 1.2) +
  geom_errorbar(aes(xmin = model_standardized_lwr_gC_m2_yr, xmax = model_standardized_upr_gC_m2_yr),
                orientation = "y", alpha = 0.30, width = 0, linewidth = 0.8) +
  geom_point(aes(shape = sign_agree, size = sign_agree), alpha = 0.87) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "Sign agrees", `FALSE` = "Sign disagrees"),
                     name = NULL) +
  scale_size_manual(values = c(`TRUE` = 2.2, `FALSE` = 3.5), guide = "none") +
  # label only sites where sign disagrees
  ggrepel::geom_text_repel(
    data = budget_comparison %>% filter(!sign_agree),
    aes(label = SITE_ID), size = 2.7, max.overlaps = 30, show.legend = FALSE
  ) +
  scale_color_manual(values = behavior_colors, na.translate = FALSE) +
  coord_cartesian(xlim = c(-axis_lim, axis_lim), ylim = c(-axis_lim, axis_lim)) +
  labs(
    title    = "ERA5 Half-Hour Gapfilled Budget vs Model-Standardized 30-Minute Budget",
    subtitle = paste0("Spearman rho = ", signif(comparison_cor, 3),
                      " - RMSE = ", signif(comparison_rmse, 3),
                      " g C m-2 yr-1 - Sign agrees: ",
                      sum(budget_comparison$sign_agree, na.rm = TRUE), " of ",
                      nrow(budget_comparison), " sites"),
    x        = expression(paste("Model-standardized annual budget (g C ", m^-2, " yr"^-1, ")")),
    y        = expression(paste("ERA5 half-hour gapfilled annual budget (g C ", m^-2, " yr"^-1, ")")),
    color    = "Reference annual class",
    caption  = "Orange quadrants = methods disagree on sink/source sign. Bars = 95% model-standardized simulation CI."
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey35"),
    plot.caption  = element_text(size = 7.5, color = "grey40"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("FIGURES/NEON_ERA5_vs_model_standardized_budget_scatter.png",
       plot_budget_comparison, width = 9, height = 8, units = "in", dpi = 300)


plot_budget_difference <- budget_comparison %>%
  mutate(SITE_ID = fct_reorder(SITE_ID, budget_difference_era5_minus_model_standardized_gC_m2_yr)) %>%
  ggplot(aes(x = budget_difference_era5_minus_model_standardized_gC_m2_yr, y = SITE_ID, fill = reference_annual_behavior)) +
  geom_vline(xintercept = 0, color = "grey45", linetype = "dashed") +
  geom_col(width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = behavior_colors, na.translate = FALSE) +
  labs(
    title = "Difference Between ERA5 Half-Hour and Model-Standardized Annual Budgets",
    x = expression(paste("ERA5 gapfilled - model-standardized budget (g C ", m^-2, " yr"^-1, ")")),
    y = NULL,
    fill = "Reference annual class"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom", axis.text.y = element_text(size = 7))

ggsave("FIGURES/NEON_ERA5_vs_model_standardized_budget_differences.png", plot_budget_difference, width = 9, height = 9, units = "in", dpi = 300)


plot_era5_annual_class_changes <- era5_annual_class_change_matrix %>%
  mutate(
    tile_fill = case_when(
      !changed & n_sites > 0 ~ "agreement",
      changed  & n_sites > 0 ~ "changed",
      TRUE                    ~ "empty"
    )
  ) %>%
  ggplot(aes(x = era5_annual_behavior, y = reference_annual_behavior)) +
  geom_tile(aes(fill = tile_fill), color = "white", linewidth = 1.4) +
  geom_text(aes(label = label), fontface = "bold", size = 8, color = "grey15") +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_y_discrete(limits = rev(class_levels), drop = FALSE) +
  scale_fill_manual(
    values = c("agreement" = "#c8e6c9", "changed" = "#ffe0b2", "empty" = "grey97"),
    labels = c("agreement" = "Unchanged class", "changed" = "Class changed",
               "empty"     = "No sites"),
    name   = NULL,
    na.translate = FALSE
  ) +
  labs(
    title    = "Reference Annual Behavior vs ERA5 Annual-Budget Class",
    subtitle = "Rows: annual classes from lookup-filled daily budgets - Columns: ERA5 annual classes from source-year fraction\nAnnual source if >= 75% of years are positive; annual sink if <= 25%; otherwise fluctuating",
    x        = "ERA5 annual-budget class",
    y        = "Reference annual class"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "grey35"),
    legend.position  = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x  = element_text(face = "bold", angle = 25, hjust = 0),
    axis.text.y  = element_text(face = "bold")
  )

ggsave("FIGURES/NEON_ERA5_reference_annual_class_changes.png",
       plot_era5_annual_class_changes, width = 7.5, height = 6.5, units = "in", dpi = 300)

plot_era5_all_site_flux_magnitude <- era5_all_site_flux_magnitude_summary %>%
  ggplot(aes(x = flux_native, y = SITE_ID_plot, color = annual_behavior, shape = EcoType)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_errorbar(
    aes(xmin = flux_lower_native, xmax = flux_upper_native),
    orientation = "y",
    width = 0.18,
    linewidth = 0.45,
    color = "grey35",
    alpha = 0.65,
    na.rm = TRUE
  ) +
  geom_point(size = 3.35, alpha = 0.55, stroke = 0.55) +
  facet_grid(
    rows = vars(annual_behavior),
    cols = vars(scale_label),
    scales = "free",
    space = "free_y"
  ) +
  scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  theme_bw(base_size = 9.6) +
  labs(
    x = "Flux Magnitude",
    y = NULL,
    shape = "Ecosystem type",
    title = expression(paste("NEON CH"[4]," Flux Magnitudes"))
  ) +
  guides(
    color = "none",
    shape = guide_legend(override.aes = list(size = 3.7, alpha = 1))
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    strip.text.x = element_text(face = "bold", size = 10, lineheight = 0.95),
    strip.text.y = element_text(face = "bold", size = 8),
    panel.spacing.x = unit(0.55, "lines"),
    panel.spacing.y = unit(0.35, "lines")
  )

ggsave(
  "FIGURES/NEON_ERA5_all_site_category_flux_magnitudes.png",
  plot_era5_all_site_flux_magnitude,
  width = 12.5,
  height = 10.5,
  units = "in",
  dpi = 300
)


north_america_map <- rnaturalearth::ne_countries(
  continent = "North America",
  returnclass = "sf"
)

plot_era5_annual_site_map <- ggplot() +
  geom_sf(data = north_america_map, fill = "grey94", color = "white", linewidth = 0.25) +
  geom_point(
    data = era5_annual_site_map_data,
    aes(
      x = longitude,
      y = latitude,
      color = annual_behavior,
      shape = EcoType,
      size = annual_flux_magnitude_gC_m2_yr
    ),
    alpha = 0.55,
    stroke = 0.65
  ) +
  scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_size_continuous(
    range = c(2.4, 8),
    breaks = scales::breaks_pretty(n = 4),
    name = expression(paste("ERA5 Annual Flux (g C ", m^-2, " ", yr^-1, ")"))
  ) +
  coord_sf(xlim = c(-170, -60), ylim = c(15, 72), expand = FALSE) +
  theme_bw(base_size = 10) +
  labs(
    x = NULL,
    y = NULL,
    color = "ERA5-Annual Behavior class",
    shape = "Ecosystem type",
    title = "ERA5 Annual Flux Behavior Across Sites",
    subtitle = "Color shows ERA5 annual behavior class; symbol size shows absolute ERA5 annual flux magnitude."
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4.2, alpha = 1)),
    shape = guide_legend(override.aes = list(size = 4.2, alpha = 1))
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 8.8),
    legend.text = element_text(size = 8.2),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.2)
  )

ggsave(
  "FIGURES/NEON_ERA5_annual_site_category_map.png",
  plot_era5_annual_site_map,
  width = 10.5,
  height = 8,
  units = "in",
  dpi = 300
)

# BEHAVIOR PATTERN FIGURE # ####

annual_method_colors <- c(
  "Scaled annual" = "#7F7F7F",
  "ERA5 annual" = "#009E73"
)

plot_era5_annual_method_boxplots <- era5_annual_method_flux_summary %>% filter(annual_method ==  "ERA5 annual") %>% 
  ggplot(aes(x = annual_behavior, y = annual_flux_gC_m2_yr, fill = annual_behavior, col = annual_behavior)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.62,
    alpha=0.5,
    position = position_dodge(width = 0.72)
  ) +
  geom_point(
    aes(group = annual_method,col=annual_behavior),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72, seed = 20260522),
    size = 2.1,
    alpha = 0.55,
    stroke = 0.35
  ) +
  scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    x = "ERA5 annual class",
    y = expression(paste("Annual CH"[4], " flux (g C ", m^-2, " ", yr^-1, ")")),
    title = "A. Annual Flux Magnitude by Behavior Class"
  ) +
  guides(fill = guide_legend(nrow = 1, order = 1, override.aes = list(alpha = 1))) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "none",
    legend.box = "vertical"
    
  )

plot_era5_diel_flux <- era5_diel_behavior_summary %>%
  ggplot(aes(x = hour_num, y = mean_flux_umolC_m2_s, color = annual_behavior, fill = annual_behavior)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_ribbon(
    aes(
      ymin = mean_flux_umolC_m2_s - se_flux_umolC_m2_s,
      ymax = mean_flux_umolC_m2_s + se_flux_umolC_m2_s
    ),
    color = NA,
    alpha = 0.16
  ) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(0, 23.5)) +
  theme_bw(base_size = 10) +
  labs(
    x = "Hour of day",
    y = expression(paste("30-min CH"[4], " flux (", mu, "mol C ", m^-2, " ", s^-1, ")")),
    color = "Behavior Class",
    title = "C. Diel Flux Pattern"
  ) +
  guides(fill = "none", color = guide_legend(nrow = 1, order = 2, override.aes = list(linewidth = 1.2))) +
  theme(
    legend.position = "none",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 12)
  )

plot_era5_diel_source_probability <- era5_diel_behavior_summary %>%
  ggplot(aes(x = hour_num, y = mean_source_probability, color = annual_behavior, fill = annual_behavior)) +
  geom_ribbon(
    aes(
      ymin = pmax(0, mean_source_probability - se_source_probability),
      ymax = pmin(1, mean_source_probability + se_source_probability)
    ),
    color = NA,
    alpha = 0.16
  ) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(0, 23.5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_bw(base_size = 12) +
  labs(
    x = "Hour of day",
    y = "Probability of positive flux",
    color = "Behavior Class",
    title = "D. Diel source probability"
  ) +
  guides(fill = "none", color = guide_legend(nrow = 1, order = 2, override.aes = list(linewidth = 1.2))) +
  theme(
    legend.position = "none",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 12)
  )

plot_era5_annual_behavior_counts <- era5_annual_behavior_counts %>%
  ggplot(aes(x = n_sites, y = fct_rev(annual_behavior), fill = annual_behavior)) +
  geom_col(width = 0.68, color = "grey30", linewidth = 0.25) +
  geom_text(aes(label = n_sites), hjust = -0.25, size = 3.1) +
  scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16)), breaks = scales::breaks_width(5)) +
  theme_bw(base_size = 12) +
  labs(
    x = "Number of sites",
    y = NULL,
    title = "B. Sites per ERA5 annual class"
  ) +
  guides(fill = "none") +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 8.7)
  )


library(ggpubr)

plot_era5_flux_pattern_panel <- ggarrange( plot_era5_annual_method_boxplots,
                                           plot_era5_annual_behavior_counts,
                                           plot_era5_diel_flux , plot_era5_diel_source_probability,
                                           ncol=2, nrow=2 ) %>%  annotate_figure(
                top = text_grob("Flux Patterns By ERA5 Annual Behavior Class", 
                                color = "black", face = "bold", size = 14))

ggsave(
  "FIGURES/NEON_ERA5_flux_pattern_diel_behavior_panel.png",
  plot_era5_flux_pattern_panel,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
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


plot_era5_annual_behavior_ecotypes <- era5_annual_behavior_ecotype_counts %>%
  ggplot(aes(x = n_sites, y = fct_rev(annual_behavior), fill = EcoType)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = ecotype_colors, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04)), breaks = scales::breaks_width(5)) +
  theme_bw(base_size = 10) +
  labs(
    x = "Number of sites",
    y = NULL,
    fill = "Ecosystem type",
    title = "By ecosystem type"
  ) +
  guides(fill = guide_legend(nrow = 1, order = 3, override.aes = list(alpha = 1))) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8)
  )




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
    "- Variables: 2 m air temperature and 0-7 cm soil volumetric water content.",
    "- Hourly values were linearly interpolated to 30-minute timestamps.",
    "",
    "## Model",
    paste0("- Training observations with ERA5 covariates: ", nrow(training_data), "."),
    paste0("- Model fit RMSE: ", signif(fit_metrics$rmse_mgC_m2_30min, 3),
           " mg C m-2 30 min-1; observed-fitted correlation: ",
           signif(fit_metrics$correlation_observed_fitted, 3), "."),
    "- Response: total CH4 flux in mg C m-2 30 min-1.",
    "- Predictors: cyclic hour, cyclic day of year, ERA5 air temperature, ERA5 soil moisture, their tensor interaction, season, ecosystem type, and site random effect.",
    "- Annual budgets retain observed half-hour fluxes where available and fill missing half-hours with ERA5-driven model predictions.",
    "",
    "## Comparison To Model-Standardized 30-Minute Annual Budget",
    paste0("- Sites compared: ", nrow(budget_comparison), "."),
    paste0("- Spearman correlation: ", signif(comparison_cor, 3), "."),
    paste0("- Mean ERA5-minus-model-standardized difference: ", signif(comparison_bias, 3), " g C m-2 yr-1."),
    paste0("- RMSE: ", signif(comparison_rmse, 3), " g C m-2 yr-1."),
    paste0("- Sign agreement: ", sum(budget_comparison$sign_agreement, na.rm = TRUE), " of ", nrow(budget_comparison), " sites."),
    "",
    "## Reference Annual vs ERA5 Annual-Budget Class",
    "- Reference annual classes are from lookup-filled daily annual budgets produced by `flow.30min.analysis.R`.",
    "- ERA5 annual class is based on the fraction of gapfilled site-years with a positive annual budget.",
    "- Annual source: at least 75% positive years; annual sink: no more than 25% positive years; otherwise fluctuating.",
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
    "- `FIGURES/NEON_ERA5_halfhour_gapfill_model_fit_summary.png`",
    "- `FIGURES/NEON_ERA5_vs_model_standardized_budget_scatter.png`",
    "- `FIGURES/NEON_ERA5_vs_model_standardized_budget_differences.png`",
    "- `FIGURES/NEON_ERA5_reference_annual_class_changes.png`",
    "- `FIGURES/NEON_ERA5_all_site_category_flux_magnitudes.png`",
    "- `FIGURES/NEON_ERA5_annual_site_category_map.png`",
    "- `FIGURES/NEON_ERA5_flux_pattern_diel_behavior_panel.png`",
  ),
  "OUTPUT/NEON_ERA5_halfhour_gapfill_results.md"
)

message("Wrote ERA5 half-hourly gapfill and budget comparison outputs.")
