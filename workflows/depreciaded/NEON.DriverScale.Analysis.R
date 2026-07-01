# Driver analysis for NEON CH4 source probability and flux magnitude by scale.
#
# This workflow intentionally separates within-site temporal drivers from
# among-site/slow drivers. The 30-minute and daily tables use site-centered
# short-term drivers. The annual table uses ERA5 annual anomalies plus site-level
# attributes. Effect uncertainty and variable importance use site-block
# bootstrap weights rather than row-level resampling.

library(tidyverse)
library(ggplot2)
library(mgcv)
library(patchwork)

input_dir <- Sys.getenv(
  "LOCALDIR_CH4",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/Methane"
)

output_dir <- Sys.getenv("NEON_DRIVER_OUTPUT_DIR", unset = file.path(input_dir, "OUTPUT"))
figure_dir <- Sys.getenv("NEON_DRIVER_FIGURE_DIR", unset = file.path(input_dir, "FIGURES"))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

input_path <- function(filename) file.path(input_dir, "OUTPUT", filename)
output_path <- function(filename) file.path(output_dir, filename)
figure_path <- function(filename) file.path(figure_dir, filename)

set.seed(20260520)

n_boot <- as.integer(Sys.getenv("NEON_DRIVER_BOOT_N", unset = "20"))
min_abs_flux <- 1e-9
seconds_per_30min <- 30 * 60
ug_c_per_umol_c <- 12.011
source_threshold <- 1
behavior_levels <- c("Weak-sink", "Fluctuating", "Weak-source")
behavior_colors <- c(
  "Weak-sink" = "#2166AC",
  "Fluctuating" = "#4D4D4D",
  "Weak-source" = "#B2182B"
)
scale_levels <- c("30min", "daily", "annual")
scale_labels <- c("30min" = "30 min", "daily" = "Daily", "annual" = "Annual")

model_data_file <- input_path("30min_ch4_model_data.csv")
daily_flux_file <- input_path("NEON_scale_daily_flux_all_sites.csv")
annual_era5_file <- input_path("NEON_ERA5_gapfilled_annual_budget_by_year.csv")
era5_halfhour_file <- input_path("NEON_ERA5_30min_site_covariates.csv.gz")
site_attribute_file <- input_path("NEON_site_attribute_multivariate_matrix.csv")

required_files <- c(model_data_file, daily_flux_file, annual_era5_file, era5_halfhour_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required inputs: ", paste(missing_files, collapse = ", "))
}

mg_c_30min_to_umol_c_s <- function(x) {
  x * 1000 / ug_c_per_umol_c / seconds_per_30min
}

center_by_site <- function(data, variable) {
  site_mean_name <- paste0(variable, "_site_mean")
  within_name <- paste0(variable, "_within_site")

  data %>%
    group_by(SITE_ID) %>%
    mutate(
      "{site_mean_name}" := mean(.data[[variable]], na.rm = TRUE),
      "{within_name}" := .data[[variable]] - .data[[site_mean_name]]
    ) %>%
    ungroup()
}

finite_or_na <- function(x) {
  x[!is.finite(x)] <- NA_real_
  x
}

fill_numeric_median <- function(data, variables) {
  out <- data
  for (variable in intersect(variables, names(out))) {
    fill_value <- median(out[[variable]], na.rm = TRUE)
    if (!is.finite(fill_value)) fill_value <- 0
    out[[variable]][!is.finite(out[[variable]]) | is.na(out[[variable]])] <- fill_value
  }
  out
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  names(sort(table(x), decreasing = TRUE))[1]
}

safe_gam <- function(formula, data, family, weights = NULL) {
  fit_fun <- if (nrow(data) > 5000) mgcv::bam else mgcv::gam
  extra_args <- if (identical(fit_fun, mgcv::bam)) {
    list(discrete = TRUE, nthreads = 2)
  } else {
    list()
  }

  tryCatch(
    {
      base_args <- list(
        formula = formula,
        data = data,
        family = family,
        method = if (identical(fit_fun, mgcv::bam)) "fREML" else "REML",
        select = TRUE
      )
      base_args <- c(base_args, extra_args)
      if (is.null(weights)) {
        do.call(fit_fun, base_args)
      } else {
        do.call(fit_fun, c(base_args, list(weights = weights)))
      }
    },
    error = function(e) {
      message("Model fit failed: ", conditionMessage(e))
      NULL
    }
  )
}

model_metric <- function(model) {
  if (is.null(model)) return(NA_real_)
  summary(model)$dev.expl
}

source_family <- binomial(link = "logit")
magnitude_family <- gaussian()

read_model_data <- function() {
  read.csv(model_data_file) %>%
    mutate(
      SITE_ID = as.character(SITE_ID),
      Date = as.Date(Date),
      time.rounded = as.POSIXct(time.rounded, tz = "UTC"),
      Year = as.integer(format(Date, "%Y")),
      month = as.integer(month),
      doy = as.numeric(doy),
      hour_num = as.numeric(hour_num),
      season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
      EcoType = factor(EcoType),
      CH4_mgC_30min = as.numeric(CH4_mgC_30min),
      flux_30min_umolC_m2_s = mg_c_30min_to_umol_c_s(CH4_mgC_30min),
      source_state = as.integer(flux_30min_umolC_m2_s > 0),
      log_abs_flux = log10(abs(flux_30min_umolC_m2_s) + min_abs_flux),
      PAR = pmax(as.numeric(PAR), 0),
      log_PAR = log1p(PAR),
      Tair_C = as.numeric(Tair_C),
      VSWCMean = as.numeric(VSWCMean)
    )
}

ch4_30min_raw <- read_model_data()

if (!"CH4_gradient_mgC_30min" %in% names(ch4_30min_raw)) {
  ch4_30min_raw$CH4_gradient_mgC_30min <- NA_real_
}
if (!"CH4_storage_mgC_30min" %in% names(ch4_30min_raw)) {
  ch4_30min_raw$CH4_storage_mgC_30min <- NA_real_
}
if (!"storage_abs_fraction" %in% names(ch4_30min_raw)) {
  ch4_30min_raw$storage_abs_fraction <- NA_real_
}

# Source/sink behavior is derived here, not in flow.30min.analysis.R. Months are
# classified after equal-weighting observed half-hour-of-day bins within each
# site-month, so behavior is not driven by uneven sampling of particular hours.
monthly_observed_ch4 <- ch4_30min_raw %>%
  reframe(
    .by = c(SITE_ID, YearMon),
    n_30min = n(),
    observed_mean_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    observed_median_CH4_mgC_30min = median(CH4_mgC_30min, na.rm = TRUE),
    observed_mean_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    observed_mean_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE),
    observed_prop_source_30min = mean(CH4_mgC_30min > 0, na.rm = TRUE),
    observed_prop_source_gradient_30min = mean(CH4_gradient_mgC_30min > 0, na.rm = TRUE),
    mean_Tair_C = mean(Tair_C, na.rm = TRUE),
    mean_VSWC = mean(VSWCMean, na.rm = TRUE)
  )

halfhour_site_month_ch4 <- ch4_30min_raw %>%
  reframe(
    .by = c(SITE_ID, YearMon, hour_num),
    n_30min_bin = n(),
    mean_CH4_mgC_30min_bin = mean(CH4_mgC_30min, na.rm = TRUE),
    median_CH4_mgC_30min_bin = median(CH4_mgC_30min, na.rm = TRUE),
    mean_gradient_mgC_30min_bin = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    median_gradient_mgC_30min_bin = median(CH4_gradient_mgC_30min, na.rm = TRUE),
    mean_storage_mgC_30min_bin = mean(CH4_storage_mgC_30min, na.rm = TRUE),
    median_storage_mgC_30min_bin = median(CH4_storage_mgC_30min, na.rm = TRUE),
    prop_source_30min_bin = mean(CH4_mgC_30min > 0, na.rm = TRUE),
    prop_source_gradient_30min_bin = mean(CH4_gradient_mgC_30min > 0, na.rm = TRUE),
    median_storage_abs_fraction_bin = median(storage_abs_fraction, na.rm = TRUE)
  )

monthly_site_ch4 <- halfhour_site_month_ch4 %>%
  reframe(
    .by = c(SITE_ID, YearMon),
    n_halfhour_bins = n_distinct(hour_num),
    mean_CH4_mgC_30min = mean(mean_CH4_mgC_30min_bin, na.rm = TRUE),
    median_CH4_mgC_30min = median(median_CH4_mgC_30min_bin, na.rm = TRUE),
    mean_gradient_mgC_30min = mean(mean_gradient_mgC_30min_bin, na.rm = TRUE),
    median_gradient_mgC_30min = median(median_gradient_mgC_30min_bin, na.rm = TRUE),
    mean_storage_mgC_30min = mean(mean_storage_mgC_30min_bin, na.rm = TRUE),
    median_storage_mgC_30min = median(median_storage_mgC_30min_bin, na.rm = TRUE),
    prop_source_30min = mean(prop_source_30min_bin, na.rm = TRUE),
    prop_source_gradient_30min = mean(prop_source_gradient_30min_bin, na.rm = TRUE),
    median_storage_abs_fraction = median(median_storage_abs_fraction_bin, na.rm = TRUE)
  ) %>%
  left_join(monthly_observed_ch4, by = c("SITE_ID", "YearMon")) %>%
  mutate(
    source_month = mean_CH4_mgC_30min > 0,
    source_gradient_month = mean_gradient_mgC_30min > 0
  )

site_behavior_core <- monthly_site_ch4 %>%
  rename(
    monthly_mean_CH4_mgC_30min = mean_CH4_mgC_30min,
    monthly_mean_gradient_mgC_30min = mean_gradient_mgC_30min,
    monthly_mean_storage_mgC_30min = mean_storage_mgC_30min
  ) %>%
  arrange(SITE_ID, YearMon) %>%
  group_by(SITE_ID) %>%
  mutate(sign_change = source_month != lag(source_month)) %>%
  reframe(
    n_months = n(),
    prop_source_months = mean(source_month, na.rm = TRUE),
    prop_sink_months = 1 - prop_source_months,
    prop_source_gradient_months = mean(source_gradient_month, na.rm = TRUE),
    mean_CH4_mgC_30min = mean(monthly_mean_CH4_mgC_30min, na.rm = TRUE),
    median_CH4_mgC_30min = median(monthly_mean_CH4_mgC_30min, na.rm = TRUE),
    sd_monthly_CH4 = sd(monthly_mean_CH4_mgC_30min, na.rm = TRUE),
    mean_gradient_mgC_30min = mean(monthly_mean_gradient_mgC_30min, na.rm = TRUE),
    median_gradient_mgC_30min = median(monthly_mean_gradient_mgC_30min, na.rm = TRUE),
    mean_storage_mgC_30min = mean(monthly_mean_storage_mgC_30min, na.rm = TRUE),
    median_storage_mgC_30min = median(monthly_mean_storage_mgC_30min, na.rm = TRUE),
    median_storage_abs_fraction = median(median_storage_abs_fraction, na.rm = TRUE),
    sign_changes = sum(sign_change, na.rm = TRUE)
  ) %>%
  mutate(
    CH4_behavior = case_when(
      prop_source_months >= source_threshold ~ "Weak-source",
      prop_source_months <= 1 - source_threshold ~ "Weak-sink",
      TRUE ~ "Fluctuating"
    ),
    CH4_gradient_behavior = case_when(
      prop_source_gradient_months >= source_threshold ~ "Weak-source",
      prop_source_gradient_months <= 1 - source_threshold ~ "Weak-sink",
      TRUE ~ "Fluctuating"
    ),
    behavior_changed_from_gradient = CH4_behavior != CH4_gradient_behavior,
    CH4_behavior = factor(CH4_behavior, levels = behavior_levels),
    CH4_gradient_behavior = factor(CH4_gradient_behavior, levels = behavior_levels)
  )

site_attributes <- if (file.exists(site_attribute_file)) {
  read.csv(site_attribute_file) %>%
    mutate(SITE_ID = as.character(SITE_ID)) %>%
    dplyr::select(
      SITE_ID, any_of(c(
        "MAP", "MAT", "acidity", "carbonTot", "nitrogenTot", "sulfurTot",
        "clayTotal", "sandTotal", "bulkDensOvenDry", "LAI.mean",
        "canopyHeight_m", "EcoType", "CH4_behavior"
      ))
    )
} else {
  ch4_30min_raw %>%
    reframe(
      .by = SITE_ID,
      EcoType = mode_value(as.character(EcoType)),
      MAP = mean(MAP, na.rm = TRUE),
      MAT = mean(MAT, na.rm = TRUE),
      acidity = mean(acidity, na.rm = TRUE),
      carbonTot = mean(carbonTot, na.rm = TRUE),
      nitrogenTot = mean(nitrogenTot, na.rm = TRUE),
      sulfurTot = mean(sulfurTot, na.rm = TRUE),
      clayTotal = mean(clayTotal, na.rm = TRUE),
      sandTotal = mean(sandTotal, na.rm = TRUE),
      bulkDensOvenDry = mean(bulkDensOvenDry, na.rm = TRUE),
      LAI.mean = mean(LAI.mean, na.rm = TRUE),
      canopyHeight_m = mean(canopyHeight_m, na.rm = TRUE)
    )
}

site_numeric_drivers <- c(
  "MAP", "MAT", "acidity", "carbonTot", "nitrogenTot", "sulfurTot",
  "clayTotal", "sandTotal", "bulkDensOvenDry", "LAI.mean", "canopyHeight_m"
)

site_ecotype_lookup <- ch4_30min_raw %>%
  reframe(.by = SITE_ID, EcoType_from_30min = mode_value(as.character(EcoType)))

if (!"EcoType" %in% names(site_attributes)) {
  site_attributes <- site_attributes %>%
    left_join(site_ecotype_lookup, by = "SITE_ID") %>%
    rename(EcoType = EcoType_from_30min)
} else {
  site_attributes <- site_attributes %>%
    left_join(site_ecotype_lookup, by = "SITE_ID") %>%
    mutate(EcoType = coalesce(as.character(EcoType), EcoType_from_30min)) %>%
    dplyr::select(-EcoType_from_30min)
}

site_attributes <- site_attributes %>%
  dplyr::select(-any_of("CH4_behavior")) %>%
  mutate(EcoType = factor(EcoType)) %>%
  fill_numeric_median(site_numeric_drivers)

site_behavior <- site_behavior_core %>%
  left_join(site_attributes, by = "SITE_ID")

component_behavior_comparison <- site_behavior_core %>%
  transmute(
    SITE_ID,
    n_months,
    prop_source_total = prop_source_months,
    prop_source_gradient = prop_source_gradient_months,
    behavior_total = CH4_behavior,
    behavior_gradient = CH4_gradient_behavior,
    behavior_changed_from_gradient,
    mean_total_mgC_30min = mean_CH4_mgC_30min,
    mean_gradient_mgC_30min,
    mean_storage_mgC_30min,
    median_storage_abs_fraction
  )

ch4_30min_raw <- ch4_30min_raw %>%
  left_join(site_behavior_core %>% dplyr::select(SITE_ID, CH4_behavior), by = "SITE_ID") %>%
  mutate(
    CH4_behavior = factor(CH4_behavior, levels = behavior_levels),
    EcoType = factor(EcoType)
  )

# readr::write_csv(monthly_site_ch4, output_path("30min_monthly_site_ch4.csv"))
# readr::write_csv(halfhour_site_month_ch4, output_path("30min_halfhour_balancing_bins.csv"))
# readr::write_csv(site_behavior, output_path("30min_site_behavior.csv"))
# readr::write_csv(component_behavior_comparison, output_path("30min_total_vs_gradient_behavior_comparison.csv"))

temporal_30min_drivers <- c("Tair_C", "VSWCMean", "log_PAR")

driver_30min_data <- ch4_30min_raw %>%
  filter(
    is.finite(flux_30min_umolC_m2_s),
    is.finite(log_abs_flux),
    is.finite(Tair_C),
    is.finite(VSWCMean),
    is.finite(log_PAR),
    is.finite(hour_num),
    is.finite(doy),
    !is.na(season),
    !is.na(EcoType)
  ) %>%
  reduce(
    temporal_30min_drivers,
    .init = .,
    .f = center_by_site
  ) %>%
  dplyr::select(
    SITE_ID, time.rounded, Date, Year, month, season, doy, hour_num,
    EcoType, source_state, flux_30min_umolC_m2_s, log_abs_flux,
    Tair_C, Tair_C_within_site, Tair_C_site_mean,
    VSWCMean, VSWCMean_within_site, VSWCMean_site_mean,
    log_PAR, log_PAR_within_site, log_PAR_site_mean
  ) %>%
  mutate(SITE_ID = factor(SITE_ID), EcoType = factor(EcoType))

daily_flux <- read.csv(daily_flux_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    Date = as.Date(Date),
    daily_mgC_m2_day = as.numeric(daily_mgC_m2_day),
    source_state = as.integer(daily_mgC_m2_day > 0),
    log_abs_flux = log10(abs(daily_mgC_m2_day) + min_abs_flux)
  ) %>%
  dplyr::select(
    SITE_ID, Date, any_of(c(
      "n_30min", "daily_mgC_m2_day_observed", "daily_gradient_mgC_m2_day_observed",
      "daily_storage_mgC_m2_day_observed", "daily_mgC_m2_day",
      "daily_gradient_mgC_m2_day", "daily_storage_mgC_m2_day"
    )),
    source_state, log_abs_flux
  )

daily_scaled_flux <- ch4_30min_raw %>%
  reframe(
    .by = c(SITE_ID, Date),
    n_30min_scaled = n(),
    daily_scaled_mgC_m2_day = sum(CH4_mgC_30min, na.rm = TRUE) * 48 / n_30min_scaled
  ) %>%
  filter(n_30min_scaled > 0)

# Lookup-filled daily fluxes are a comparison product here, not the primary
# scale-based daily input. Observed half-hours are retained; missing slots are
# filled from progressively broader half-hour-of-day lookup means:
# site-month, site-season, site-bi-season, site-annual, then matching global
# fallbacks. Bi-seasonal means are warm (Spring/Summer) and cool
# (Autumn/Winter) groups.
half_hour_bins <- sort(unique(ch4_30min_raw$hour_num[is.finite(ch4_30min_raw$hour_num)]))

ch4_30min_lookup_source <- ch4_30min_raw %>%
  mutate(
    season_lookup = as.character(season),
    biseason_lookup = case_when(
      season_lookup %in% c("Spring", "Summer") ~ "Warm",
      season_lookup %in% c("Autumn", "Winter") ~ "Cool",
      TRUE ~ NA_character_
    )
  )

daily_observed_bins <- ch4_30min_lookup_source %>%
  reframe(
    .by = c(SITE_ID, Date, hour_num),
    observed_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE)
  )

site_month_hour_lookup <- ch4_30min_lookup_source %>%
  reframe(
    .by = c(SITE_ID, month, hour_num),
    lookup_site_month_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    n_site_month_hour_lookup = n()
  )

site_season_hour_lookup <- ch4_30min_lookup_source %>%
  reframe(
    .by = c(SITE_ID, season_lookup, hour_num),
    lookup_site_season_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    n_site_season_hour_lookup = n()
  )

site_biseason_hour_lookup <- ch4_30min_lookup_source %>%
  filter(!is.na(biseason_lookup)) %>%
  reframe(
    .by = c(SITE_ID, biseason_lookup, hour_num),
    lookup_site_biseason_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    n_site_biseason_hour_lookup = n()
  )

site_hour_lookup <- ch4_30min_lookup_source %>%
  reframe(
    .by = c(SITE_ID, hour_num),
    lookup_site_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    n_site_hour_lookup = n()
  )

global_month_hour_lookup <- ch4_30min_lookup_source %>%
  reframe(
    .by = c(month, hour_num),
    lookup_global_month_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE)
  )

global_season_hour_lookup <- ch4_30min_lookup_source %>%
  reframe(
    .by = c(season_lookup, hour_num),
    lookup_global_season_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE)
  )

global_biseason_hour_lookup <- ch4_30min_lookup_source %>%
  filter(!is.na(biseason_lookup)) %>%
  reframe(
    .by = c(biseason_lookup, hour_num),
    lookup_global_biseason_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE)
  )

global_hour_lookup <- ch4_30min_lookup_source %>%
  reframe(
    .by = hour_num,
    lookup_global_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE)
  )

site_lookup <- ch4_30min_lookup_source %>%
  reframe(
    .by = SITE_ID,
    lookup_site_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE)
  )

global_lookup <- mean(ch4_30min_lookup_source$CH4_mgC_30min, na.rm = TRUE)

daily_lookup_dates <- ch4_30min_lookup_source %>%
  reframe(
    .by = c(SITE_ID, Date),
    month = first(month),
    season_lookup = mode_value(season_lookup),
    biseason_lookup = mode_value(biseason_lookup)
  )

daily_lookup_grid <- daily_lookup_dates %>%
  tidyr::expand_grid(hour_num = half_hour_bins) %>%
  left_join(daily_observed_bins, by = c("SITE_ID", "Date", "hour_num")) %>%
  left_join(site_month_hour_lookup, by = c("SITE_ID", "month", "hour_num")) %>%
  left_join(site_season_hour_lookup, by = c("SITE_ID", "season_lookup", "hour_num")) %>%
  left_join(site_biseason_hour_lookup, by = c("SITE_ID", "biseason_lookup", "hour_num")) %>%
  left_join(site_hour_lookup, by = c("SITE_ID", "hour_num")) %>%
  left_join(global_month_hour_lookup, by = c("month", "hour_num")) %>%
  left_join(global_season_hour_lookup, by = c("season_lookup", "hour_num")) %>%
  left_join(global_biseason_hour_lookup, by = c("biseason_lookup", "hour_num")) %>%
  left_join(global_hour_lookup, by = "hour_num") %>%
  left_join(site_lookup, by = "SITE_ID") %>%
  mutate(
    fill_source = case_when(
      is.finite(observed_CH4_mgC_30min) ~ "observed",
      is.finite(lookup_site_month_hour_CH4_mgC_30min) ~ "site_month_hour_lookup",
      is.finite(lookup_site_season_hour_CH4_mgC_30min) ~ "site_season_hour_lookup",
      is.finite(lookup_site_biseason_hour_CH4_mgC_30min) ~ "site_biseason_hour_lookup",
      is.finite(lookup_site_hour_CH4_mgC_30min) ~ "site_annual_hour_lookup",
      is.finite(lookup_global_month_hour_CH4_mgC_30min) ~ "global_month_hour_lookup",
      is.finite(lookup_global_season_hour_CH4_mgC_30min) ~ "global_season_hour_lookup",
      is.finite(lookup_global_biseason_hour_CH4_mgC_30min) ~ "global_biseason_hour_lookup",
      is.finite(lookup_global_hour_CH4_mgC_30min) ~ "global_annual_hour_lookup",
      is.finite(lookup_site_CH4_mgC_30min) ~ "site_annual_mean_lookup",
      TRUE ~ "global_lookup"
    ),
    filled_CH4_mgC_30min = coalesce(
      observed_CH4_mgC_30min,
      lookup_site_month_hour_CH4_mgC_30min,
      lookup_site_season_hour_CH4_mgC_30min,
      lookup_site_biseason_hour_CH4_mgC_30min,
      lookup_site_hour_CH4_mgC_30min,
      lookup_global_month_hour_CH4_mgC_30min,
      lookup_global_season_hour_CH4_mgC_30min,
      lookup_global_biseason_hour_CH4_mgC_30min,
      lookup_global_hour_CH4_mgC_30min,
      lookup_site_CH4_mgC_30min,
      global_lookup
    ),
    was_observed = fill_source == "observed"
  )

daily_lookup_flux <- daily_lookup_grid %>%
  reframe(
    .by = c(SITE_ID, Date),
    n_observed_30min = sum(was_observed, na.rm = TRUE),
    n_lookup_filled_30min = sum(!was_observed, na.rm = TRUE),
    n_site_month_hour_lookup_filled_30min = sum(fill_source == "site_month_hour_lookup", na.rm = TRUE),
    n_site_season_hour_lookup_filled_30min = sum(fill_source == "site_season_hour_lookup", na.rm = TRUE),
    n_site_biseason_hour_lookup_filled_30min = sum(fill_source == "site_biseason_hour_lookup", na.rm = TRUE),
    n_site_annual_hour_lookup_filled_30min = sum(fill_source == "site_annual_hour_lookup", na.rm = TRUE),
    n_global_lookup_filled_30min = sum(str_detect(fill_source, "^global"), na.rm = TRUE),
    daily_lookup_mgC_m2_day = sum(filled_CH4_mgC_30min, na.rm = TRUE),
    daily_lookup_mean_mgC_m2_30min = mean(filled_CH4_mgC_30min, na.rm = TRUE)
  )

daily_scaled_vs_lookup <- daily_scaled_flux %>%
  left_join(daily_lookup_flux, by = c("SITE_ID", "Date")) %>%
  mutate(
    daily_lookup_minus_scaled_mgC_m2_day = daily_lookup_mgC_m2_day - daily_scaled_mgC_m2_day,
    daily_lookup_to_scaled_ratio = daily_lookup_mgC_m2_day / daily_scaled_mgC_m2_day
  )

# readr::write_csv(site_hour_lookup, output_path("NEON_driver_daily_site_halfhour_lookup.csv"))
# readr::write_csv(daily_lookup_grid, output_path("NEON_driver_daily_lookup_grid.csv"))
# readr::write_csv(daily_lookup_flux, output_path("NEON_driver_daily_lookup_flux_all_sites.csv"))
# readr::write_csv(daily_scaled_vs_lookup, output_path("NEON_driver_daily_scaled_vs_lookup_comparison.csv"))

daily_drivers <- ch4_30min_raw %>%
  filter(!is.na(Date)) %>%
  reframe(
    .by = c(SITE_ID, Date),
    Year = first(Year),
    month = as.integer(format(first(Date), "%m")),
    doy = as.numeric(format(first(Date), "%j")),
    season = mode_value(as.character(season)),
    EcoType = mode_value(as.character(EcoType)),
    mean_Tair_C = mean(Tair_C, na.rm = TRUE),
    mean_VSWC = mean(VSWCMean, na.rm = TRUE),
    mean_log_PAR = mean(log_PAR, na.rm = TRUE),
    max_PAR = max(PAR, na.rm = TRUE)
  )

temporal_daily_drivers <- c("mean_Tair_C", "mean_VSWC", "mean_log_PAR")

daily_driver_data <- daily_flux %>%
  left_join(daily_drivers, by = c("SITE_ID", "Date")) %>%
  filter(
    is.finite(daily_mgC_m2_day),
    is.finite(log_abs_flux),
    is.finite(mean_Tair_C),
    is.finite(mean_VSWC),
    is.finite(mean_log_PAR),
    is.finite(doy),
    !is.na(season),
    !is.na(EcoType)
  ) %>%
  reduce(
    temporal_daily_drivers,
    .init = .,
    .f = center_by_site
  ) %>%
  mutate(
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    SITE_ID = factor(SITE_ID),
    EcoType = factor(EcoType)
  )

era5_annual_drivers <- read.csv(era5_halfhour_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    Year = as.integer(Year),
    ERA5_Tair_C = as.numeric(ERA5_Tair_C),
    ERA5_VSWC = as.numeric(ERA5_VSWC)
  ) %>%
  reframe(
    .by = c(SITE_ID, Year),
    ERA5_Tair_C_annual = mean(ERA5_Tair_C, na.rm = TRUE),
    ERA5_Tair_C_sd_annual = sd(ERA5_Tair_C, na.rm = TRUE),
    ERA5_VSWC_annual = mean(ERA5_VSWC, na.rm = TRUE),
    ERA5_VSWC_sd_annual = sd(ERA5_VSWC, na.rm = TRUE)
  )

annual_era5_driver_data <- read.csv(annual_era5_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    Year = as.integer(Year),
    annual_budget_gC_m2_yr = as.numeric(annual_budget_gC_m2_yr),
    source_state = as.integer(annual_budget_gC_m2_yr > 0),
    log_abs_flux = log10(abs(annual_budget_gC_m2_yr) + min_abs_flux)
  ) %>%
  left_join(era5_annual_drivers, by = c("SITE_ID", "Year")) %>%
  left_join(site_attributes, by = "SITE_ID") %>%
  reduce(
    c("ERA5_Tair_C_annual", "ERA5_VSWC_annual"),
    .init = .,
    .f = center_by_site
  ) %>%
  fill_numeric_median(c(
    "ERA5_Tair_C_annual", "ERA5_Tair_C_sd_annual", "ERA5_VSWC_annual",
    "ERA5_VSWC_sd_annual", "ERA5_Tair_C_annual_within_site",
    "ERA5_VSWC_annual_within_site", "ERA5_Tair_C_annual_site_mean",
    "ERA5_VSWC_annual_site_mean", site_numeric_drivers
  )) %>%
  filter(is.finite(annual_budget_gC_m2_yr), is.finite(log_abs_flux)) %>%
  mutate(SITE_ID = factor(SITE_ID), EcoType = factor(EcoType))

daily_behavior <- daily_driver_data %>%
  reframe(
    .by = SITE_ID,
    n_days_daily = n(),
    prop_source_daily = mean(daily_mgC_m2_day > 0, na.rm = TRUE),
    mean_daily_mgC_m2_day = mean(daily_mgC_m2_day, na.rm = TRUE),
    median_daily_mgC_m2_day = median(daily_mgC_m2_day, na.rm = TRUE),
    sd_daily_mgC_m2_day = sd(daily_mgC_m2_day, na.rm = TRUE)
  ) %>%
  mutate(
    CH4_daily_behavior = case_when(
      prop_source_daily >= source_threshold ~ "Weak-source",
      prop_source_daily <= 1 - source_threshold ~ "Weak-sink",
      TRUE ~ "Fluctuating"
    ),
    CH4_daily_behavior = factor(CH4_daily_behavior, levels = behavior_levels)
  )

annual_behavior <- annual_era5_driver_data %>%
  reframe(
    .by = SITE_ID,
    n_years_annual = n(),
    prop_source_annual = mean(annual_budget_gC_m2_yr > 0, na.rm = TRUE),
    mean_annual_budget_gC_m2_yr = mean(annual_budget_gC_m2_yr, na.rm = TRUE),
    median_annual_budget_gC_m2_yr = median(annual_budget_gC_m2_yr, na.rm = TRUE),
    sd_annual_budget_gC_m2_yr = sd(annual_budget_gC_m2_yr, na.rm = TRUE)
  ) %>%
  mutate(
    CH4_annual_behavior = case_when(
      prop_source_annual >= source_threshold ~ "Weak-source",
      prop_source_annual <= 1 - source_threshold ~ "Weak-sink",
      TRUE ~ "Fluctuating"
    ),
    CH4_annual_behavior = factor(CH4_annual_behavior, levels = behavior_levels)
  )

source_sink_scale_summary <- site_behavior_core %>%
  dplyr::select(
    SITE_ID,
    n_months_30min = n_months,
    prop_source_30min = prop_source_months,
    flux_30min_mgC_m2_30min = mean_CH4_mgC_30min,
    behavior_30min = CH4_behavior
  ) %>%
  left_join(daily_behavior, by = "SITE_ID") %>%
  left_join(annual_behavior, by = "SITE_ID") %>%
  mutate(
    behavior_30min = factor(behavior_30min, levels = behavior_levels),
    behavior_30min_to_daily_changed = behavior_30min != CH4_daily_behavior,
    behavior_daily_to_annual_changed = CH4_daily_behavior != CH4_annual_behavior,
    behavior_30min_to_annual_changed = behavior_30min != CH4_annual_behavior
  )

source_sink_long <- source_sink_scale_summary %>%
  transmute(
    SITE_ID,
    `30 min` = behavior_30min,
    Daily = CH4_daily_behavior,
    Annual = CH4_annual_behavior,
    prop_source_30min,
    prop_source_daily,
    prop_source_annual
  ) %>%
  pivot_longer(cols = c(`30 min`, Daily, Annual), names_to = "scale", values_to = "behavior") %>%
  mutate(
    prop_source = case_when(
      scale == "30 min" ~ prop_source_30min,
      scale == "Daily" ~ prop_source_daily,
      scale == "Annual" ~ prop_source_annual
    ),
    scale = factor(scale, levels = c("30 min", "Daily", "Annual")),
    behavior = factor(behavior, levels = behavior_levels)
  )

source_sink_counts <- source_sink_long %>%
  count(scale, behavior, name = "n_sites") %>%
  complete(
    scale,
    behavior = factor(behavior_levels, levels = behavior_levels),
    fill = list(n_sites = 0L)
  )

class_change_30min_daily <- source_sink_scale_summary %>%
  count(behavior_30min, CH4_daily_behavior, name = "n_sites") %>%
  complete(
    behavior_30min = factor(behavior_levels, levels = behavior_levels),
    CH4_daily_behavior = factor(behavior_levels, levels = behavior_levels),
    fill = list(n_sites = 0L)
  )

class_change_daily_annual <- source_sink_scale_summary %>%
  count(CH4_daily_behavior, CH4_annual_behavior, name = "n_sites") %>%
  complete(
    CH4_daily_behavior = factor(behavior_levels, levels = behavior_levels),
    CH4_annual_behavior = factor(behavior_levels, levels = behavior_levels),
    fill = list(n_sites = 0L)
  )

class_change_30min_annual <- source_sink_scale_summary %>%
  count(behavior_30min, CH4_annual_behavior, name = "n_sites") %>%
  complete(
    behavior_30min = factor(behavior_levels, levels = behavior_levels),
    CH4_annual_behavior = factor(behavior_levels, levels = behavior_levels),
    fill = list(n_sites = 0L)
  )

# readr::write_csv(daily_behavior, output_path("NEON_driver_scale_daily_behavior.csv"))
# readr::write_csv(annual_behavior, output_path("NEON_driver_scale_annual_behavior.csv"))
# readr::write_csv(source_sink_scale_summary, output_path("NEON_driver_scale_source_sink_summary.csv"))
# readr::write_csv(source_sink_counts, output_path("NEON_driver_scale_source_sink_counts.csv"))
# readr::write_csv(class_change_30min_daily, output_path("NEON_driver_scale_class_change_30min_to_daily.csv"))
# readr::write_csv(class_change_daily_annual, output_path("NEON_driver_scale_class_change_daily_to_annual.csv"))
# readr::write_csv(class_change_30min_annual, output_path("NEON_driver_scale_class_change_30min_to_annual.csv"))

# readr::write_csv(driver_30min_data, output_path("NEON_driver_scale_30min_driver_data.csv"))
# readr::write_csv(daily_driver_data, output_path("NEON_driver_scale_daily_driver_data.csv"))
# readr::write_csv(annual_era5_driver_data, output_path("NEON_driver_scale_annual_era5_driver_data.csv"))

scale_specs <- list(
  `30min` = list(
    data = driver_30min_data,
    response_flux = "flux_30min_umolC_m2_s",
    source_formula = source_state ~
      s(Tair_C_within_site, k = 6) +
      s(VSWCMean_within_site, k = 6) +
      s(log_PAR_within_site, k = 6) +
      s(hour_num, bs = "cc", k = 12) +
      s(doy, bs = "cc", k = 20) +
      season + EcoType + s(SITE_ID, bs = "re"),
    magnitude_formula = log_abs_flux ~
      s(Tair_C_within_site, k = 6) +
      s(VSWCMean_within_site, k = 6) +
      s(log_PAR_within_site, k = 6) +
      s(hour_num, bs = "cc", k = 12) +
      s(doy, bs = "cc", k = 20) +
      season + EcoType + s(SITE_ID, bs = "re"),
    driver_terms = list(
      Tair_within = "s(Tair_C_within_site)",
      VSWC_within = "s(VSWCMean_within_site)",
      PAR_within = "s(log_PAR_within_site)",
      hour_of_day = "s(hour_num)",
      day_of_year = "s(doy)",
      season = "season",
      EcoType = "EcoType"
    ),
    partial_drivers = c("Tair_C_within_site", "VSWCMean_within_site", "log_PAR_within_site", "hour_num", "doy")
  ),
  daily = list(
    data = daily_driver_data,
    response_flux = "daily_mgC_m2_day",
    source_formula = source_state ~
      s(mean_Tair_C_within_site, k = 6) +
      s(mean_VSWC_within_site, k = 6) +
      s(mean_log_PAR_within_site, k = 6) +
      s(doy, bs = "cc", k = 20) +
      season + EcoType + s(SITE_ID, bs = "re"),
    magnitude_formula = log_abs_flux ~
      s(mean_Tair_C_within_site, k = 6) +
      s(mean_VSWC_within_site, k = 6) +
      s(mean_log_PAR_within_site, k = 6) +
      s(doy, bs = "cc", k = 20) +
      season + EcoType + s(SITE_ID, bs = "re"),
    driver_terms = list(
      Tair_within = "s(mean_Tair_C_within_site)",
      VSWC_within = "s(mean_VSWC_within_site)",
      PAR_within = "s(mean_log_PAR_within_site)",
      day_of_year = "s(doy)",
      season = "season",
      EcoType = "EcoType"
    ),
    partial_drivers = c("mean_Tair_C_within_site", "mean_VSWC_within_site", "mean_log_PAR_within_site", "doy")
  ),
  annual = list(
    data = annual_era5_driver_data,
    response_flux = "annual_budget_gC_m2_yr",
    source_formula = source_state ~
      ERA5_Tair_C_annual_within_site +
      ERA5_VSWC_annual_within_site +
      ERA5_Tair_C_annual_site_mean +
      ERA5_VSWC_annual_site_mean +
      MAP + acidity + carbonTot + sulfurTot + EcoType,
    magnitude_formula = log_abs_flux ~
      ERA5_Tair_C_annual_within_site +
      ERA5_VSWC_annual_within_site +
      ERA5_Tair_C_annual_site_mean +
      ERA5_VSWC_annual_site_mean +
      MAP + acidity + carbonTot + sulfurTot + EcoType,
    driver_terms = list(
      ERA5_Tair_within = "ERA5_Tair_C_annual_within_site",
      ERA5_VSWC_within = "ERA5_VSWC_annual_within_site",
      ERA5_Tair_site_mean = "ERA5_Tair_C_annual_site_mean",
      ERA5_VSWC_site_mean = "ERA5_VSWC_annual_site_mean",
      MAP = "MAP",
      acidity = "acidity",
      carbonTot = "carbonTot",
      sulfurTot = "sulfurTot",
      EcoType = "EcoType"
    ),
    partial_drivers = c(
      "ERA5_Tair_C_annual_within_site", "ERA5_VSWC_annual_within_site",
      "ERA5_Tair_C_annual_site_mean", "ERA5_VSWC_annual_site_mean",
      "MAP", "acidity", "carbonTot", "sulfurTot"
    )
  )
)

fit_scale_models <- function(spec) {
  data <- spec$data
  list(
    source = safe_gam(spec$source_formula, data, source_family),
    magnitude = safe_gam(spec$magnitude_formula, data, magnitude_family)
  )
}

scale_models <- purrr::map(scale_specs, fit_scale_models)

model_summary <- purrr::imap_dfr(scale_specs, function(spec, scale_name) {
  models <- scale_models[[scale_name]]
  tibble(
    scale = scale_name,
    response = c("source_probability", "magnitude_log_abs"),
    n_rows = nrow(spec$data),
    n_sites = n_distinct(spec$data$SITE_ID),
    deviance_explained = c(model_metric(models$source), model_metric(models$magnitude)),
    aic = c(if (is.null(models$source)) NA_real_ else AIC(models$source),
      if (is.null(models$magnitude)) NA_real_ else AIC(models$magnitude))
  )
}) %>%
  mutate(scale = factor(scale, levels = scale_levels, labels = scale_labels))

# readr::write_csv(model_summary, output_path("NEON_driver_scale_model_summary.csv"))

site_boot_weights <- function(data) {
  sites <- levels(factor(data$SITE_ID))
  sampled <- sample(sites, length(sites), replace = TRUE)
  counts <- table(factor(sampled, levels = sites))
  as.numeric(counts[as.character(data$SITE_ID)])
}

importance_from_model <- function(model, driver_terms) {
  if (is.null(model)) {
    return(tibble(
      driver = names(driver_terms),
      full_deviance_explained = NA_real_,
      importance_statistic = NA_real_
    ))
  }

  model_summary <- summary(model)
  smooth_table <- as.data.frame(model_summary$s.table)
  smooth_table$term <- rownames(smooth_table)
  param_table <- as.data.frame(model_summary$p.table)
  param_table$term <- rownames(param_table)

  smooth_stat_col <- intersect(c("Chi.sq", "F"), names(smooth_table))[1]
  param_stat_col <- intersect(c("z value", "t value"), names(param_table))[1]

  purrr::imap_dfr(driver_terms, function(term, driver) {
    if (startsWith(term, "s(")) {
      idx <- startsWith(smooth_table$term, term)
      statistic <- if (any(idx) && !is.na(smooth_stat_col)) {
        sum(abs(smooth_table[[smooth_stat_col]][idx]), na.rm = TRUE)
      } else {
        NA_real_
      }
    } else {
      idx <- startsWith(param_table$term, term)
      statistic <- if (any(idx) && !is.na(param_stat_col)) {
        sum(abs(param_table[[param_stat_col]][idx]), na.rm = TRUE)
      } else {
        NA_real_
      }
    }

    tibble(
      driver = driver,
      full_deviance_explained = model_metric(model),
      importance_statistic = statistic
    )
  })
}

importance_once <- function(spec, response_name, weights = NULL) {
  formula <- if (response_name == "source_probability") spec$source_formula else spec$magnitude_formula
  family <- if (response_name == "source_probability") source_family else magnitude_family
  data <- spec$data

  if (response_name == "source_probability" && length(unique(data$source_state[weights %||% rep(1, nrow(data)) > 0])) < 2) {
    return(tibble())
  }

  full_model <- safe_gam(formula, data, family, weights = weights)
  importance_from_model(full_model, spec$driver_terms)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

importance_base <- purrr::imap_dfr(scale_specs, function(spec, scale_name) {
  bind_rows(
    importance_once(spec, "source_probability") %>% mutate(response = "source_probability"),
    importance_once(spec, "magnitude_log_abs") %>% mutate(response = "magnitude_log_abs")
  ) %>%
    mutate(scale = scale_name, bootstrap = 0L)
})

bootstrap_importance <- purrr::imap_dfr(scale_specs, function(spec, scale_name) {
  purrr::map_dfr(seq_len(n_boot), function(i) {
    weights <- site_boot_weights(spec$data)
    bind_rows(
      importance_once(spec, "source_probability", weights = weights) %>% mutate(response = "source_probability"),
      importance_once(spec, "magnitude_log_abs", weights = weights) %>% mutate(response = "magnitude_log_abs")
    ) %>%
      mutate(scale = scale_name, bootstrap = i)
  })
})

importance_all <- bind_rows(importance_base, bootstrap_importance)
# readr::write_csv(importance_all, output_path("NEON_driver_scale_variable_importance_bootstrap.csv"))

importance_summary <- importance_all %>%
  filter(bootstrap > 0) %>%
  reframe(
    .by = c(scale, response, driver),
    importance_median = median(importance_statistic, na.rm = TRUE),
    importance_lwr = quantile(importance_statistic, 0.025, na.rm = TRUE),
    importance_upr = quantile(importance_statistic, 0.975, na.rm = TRUE),
    supported = is.finite(importance_lwr) & importance_lwr > 0 & importance_median > 0.001
  ) %>%
  filter(is.finite(importance_median)) %>%
  mutate(scale = factor(scale, levels = scale_levels, labels = scale_labels)) %>%
  arrange(response, scale, desc(importance_median))

# readr::write_csv(importance_summary, output_path("NEON_driver_scale_variable_importance_summary.csv"))

make_reference_row <- function(data) {
  reference <- data[1, , drop = FALSE]
  for (nm in names(reference)) {
    if (is.numeric(data[[nm]]) || is.integer(data[[nm]])) {
      reference[[nm]] <- median(data[[nm]], na.rm = TRUE)
    } else if (is.factor(data[[nm]])) {
      reference[[nm]] <- factor(mode_value(as.character(data[[nm]])), levels = levels(data[[nm]]))
    } else {
      reference[[nm]] <- mode_value(data[[nm]])
    }
  }
  reference$SITE_ID <- factor(levels(factor(data$SITE_ID))[1], levels = levels(factor(data$SITE_ID)))
  reference
}

partial_effects_one <- function(spec, scale_name, response_name, model, driver) {
  if (is.null(model) || !driver %in% names(spec$data)) return(tibble())
  data <- spec$data
  reference <- make_reference_row(data)
  driver_values <- quantile(data[[driver]], probs = seq(0.05, 0.95, length.out = 60), na.rm = TRUE)
  driver_values <- unique(as.numeric(driver_values))
  if (length(driver_values) < 2) return(tibble())

  newdata <- reference[rep(1, length(driver_values)), , drop = FALSE]
  newdata[[driver]] <- driver_values

  random_terms <- grep("SITE_ID", names(model$smooth), value = TRUE)
  pred <- tryCatch(
    predict(model, newdata = newdata, type = "link", se.fit = TRUE, exclude = "s(SITE_ID)"),
    error = function(e) predict(model, newdata = newdata, type = "link", se.fit = TRUE)
  )

  fit_link <- as.numeric(pred$fit)
  se_link <- as.numeric(pred$se.fit)
  if (response_name == "source_probability") {
    fit <- plogis(fit_link)
    lwr <- plogis(fit_link - 1.96 * se_link)
    upr <- plogis(fit_link + 1.96 * se_link)
  } else {
    fit <- fit_link
    lwr <- fit_link - 1.96 * se_link
    upr <- fit_link + 1.96 * se_link
  }

  tibble(
    scale = scale_name,
    response = response_name,
    driver = driver,
    driver_value = driver_values,
    fit = fit,
    lwr = lwr,
    upr = upr
  )
}

partial_effects <- purrr::imap_dfr(scale_specs, function(spec, scale_name) {
  models <- scale_models[[scale_name]]
  purrr::map_dfr(spec$partial_drivers, function(driver) {
    bind_rows(
      partial_effects_one(spec, scale_name, "source_probability", models$source, driver),
      partial_effects_one(spec, scale_name, "magnitude_log_abs", models$magnitude, driver)
    )
  })
}) %>%
  mutate(scale = factor(scale, levels = scale_levels, labels = scale_labels))

# readr::write_csv(partial_effects, output_path("NEON_driver_scale_partial_effects.csv"))

driver_scale_summary <- importance_summary %>%
  mutate(
    scale_raw = names(scale_labels)[match(as.character(scale), scale_labels)],
    driver_family = case_when(
      str_detect(driver, "within") ~ "within_site_temporal",
      str_detect(driver, "site_mean|MAP|MAT|acidity|carbonTot|sulfurTot|EcoType") ~ "among_site_or_slow",
      driver %in% c("hour_of_day", "day_of_year", "season") ~ "seasonal_or_diel",
      TRUE ~ "other"
    )
  ) %>%
  reframe(
    .by = c(response, driver),
    acts_30min = any(scale_raw == "30min" & supported, na.rm = TRUE),
    acts_daily = any(scale_raw == "daily" & supported, na.rm = TRUE),
    acts_annual = any(scale_raw == "annual" & supported, na.rm = TRUE),
    acts_at_scales = paste(c("30min", "daily", "annual")[c(acts_30min, acts_daily, acts_annual)], collapse = ", "),
    max_importance_median = max(importance_median, na.rm = TRUE),
    driver_family = first(driver_family)
  ) %>%
  mutate(acts_at_scales = if_else(acts_at_scales == "", "none", acts_at_scales)) %>%
  arrange(response, desc(max_importance_median))

# readr::write_csv(driver_scale_summary, output_path("NEON_driver_scale_cross_scale_summary.csv"))

plot_importance <- importance_summary %>%
  mutate(
    driver = fct_reorder(driver, importance_median),
    response = recode(
      response,
      source_probability = "Source probability",
      magnitude_log_abs = "Flux magnitude"
    )
  ) %>%
  ggplot(aes(x = importance_median, y = driver, color = supported)) +
  geom_errorbar(aes(xmin = importance_lwr, xmax = importance_upr), orientation = "y", width = 0.18, alpha = 0.7) +
  geom_point(size = 2.2) +
  facet_grid(response ~ scale, scales = "free") +
  scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey45")) +
  theme_bw(base_size = 10.5) +
  labs(
    x = "Variable importance statistic\n(smooth F/Chi-square or summed parametric |t/z|)",
    y = NULL,
    color = "Bootstrap-supported",
    title = "Driver Importance By Temporal Scale",
    subtitle = paste0("Intervals are site-block bootstrap 95% intervals; n_boot = ", n_boot, ".")
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ggsave(
#   figure_path("NEON_driver_scale_variable_importance.png"),
#   plot = plot_importance,
#   width = 12,
#   height = 8,
#   units = "in",
#   dpi = 300
# )
# ggsave(
#   figure_path("NEON_driver_scale_variable_importance.pdf"),
#   plot = plot_importance,
#   width = 12,
#   height = 8,
#   units = "in"
# )

plot_partials <- partial_effects %>%
  mutate(
    response = recode(
      response,
      source_probability = "Source probability",
      magnitude_log_abs = "Flux magnitude"
    )
  ) %>%
  ggplot(aes(x = driver_value, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey75", alpha = 0.55) +
  geom_line(color = "#2166AC", linewidth = 0.7) +
  facet_grid(response + scale ~ driver, scales = "free") +
  theme_bw(base_size = 9) +
  labs(
    x = "Driver value",
    y = "Partial effect",
    title = "Source-Probability And Magnitude Partial Effects By Scale",
    subtitle = "Predictions hold other drivers at typical values and exclude site random effects."
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 7),
    panel.grid.minor = element_blank()
  )

# ggsave(
#   figure_path("NEON_driver_scale_partial_effects.png"),
#   plot = plot_partials,
#   width = 16,
#   height = 11,
#   units = "in",
#   dpi = 300
# )
# ggsave(
#   figure_path("NEON_driver_scale_partial_effects.pdf"),
#   plot = plot_partials,
#   width = 16,
#   height = 11,
#   units = "in"
# )

plot_source_sink_counts <- source_sink_counts %>%
  ggplot(aes(x = scale, y = n_sites, fill = behavior)) +
  geom_col(width = 0.68) +
  scale_fill_manual(values = behavior_colors, drop = FALSE) +
  theme_bw(base_size = 10.5) +
  labs(
    x = NULL,
    y = "Number of sites",
    fill = "Source/sink class",
    title = "CH4 Source/Sink Class Counts By Temporal Scale"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ggsave(
#   figure_path("NEON_driver_scale_source_sink_counts.png"),
#   plot = plot_source_sink_counts,
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 300
# )
# ggsave(
#   figure_path("NEON_driver_scale_source_sink_counts.pdf"),
#   plot = plot_source_sink_counts,
#   width = 7,
#   height = 5,
#   units = "in"
# )

plot_transition_matrix <- function(data, x_col, y_col, title, x_lab, y_lab) {
  data %>%
    mutate(
      changed = .data[[x_col]] != .data[[y_col]],
      label = if_else(n_sites > 0, as.character(n_sites), "")
    ) %>%
    ggplot(aes(x = .data[[y_col]], y = .data[[x_col]], fill = changed)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label), fontface = "bold", size = 4) +
    scale_fill_manual(
      values = c("FALSE" = "#D9E8F5", "TRUE" = "#F3B6A8"),
      labels = c("Same class", "Changed class"),
      drop = FALSE
    ) +
    coord_equal() +
    theme_bw(base_size = 10.5) +
    labs(x = x_lab, y = y_lab, fill = NULL, title = title) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
}

plot_source_sink_transitions <- (
  plot_transition_matrix(
    class_change_30min_daily,
    "behavior_30min",
    "CH4_daily_behavior",
    "30-Minute To Daily",
    "Daily class",
    "30-minute class"
  ) |
    plot_transition_matrix(
      class_change_daily_annual,
      "CH4_daily_behavior",
      "CH4_annual_behavior",
      "Daily To Annual",
      "Annual class",
      "Daily class"
    ) |
    plot_transition_matrix(
      class_change_30min_annual,
      "behavior_30min",
      "CH4_annual_behavior",
      "30-Minute To Annual",
      "Annual class",
      "30-minute class"
    )
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ggsave(
#   figure_path("NEON_driver_scale_source_sink_transitions.png"),
#   plot = plot_source_sink_transitions,
#   width = 11,
#   height = 4.5,
#   units = "in",
#   dpi = 300
# )
# ggsave(
#   figure_path("NEON_driver_scale_source_sink_transitions.pdf"),
#   plot = plot_source_sink_transitions,
#   width = 11,
#   height = 4.5,
#   units = "in"
# )

scale_counts <- bind_rows(
  tibble(scale = "30min", n_rows = nrow(driver_30min_data), n_sites = n_distinct(driver_30min_data$SITE_ID)),
  tibble(scale = "daily", n_rows = nrow(daily_driver_data), n_sites = n_distinct(daily_driver_data$SITE_ID)),
  tibble(scale = "annual", n_rows = nrow(annual_era5_driver_data), n_sites = n_distinct(annual_era5_driver_data$SITE_ID))
)

summary_lines <- c(
  "# NEON Driver Scale Analysis",
  "",
  "## Design",
  "- Source/sink classes are derived in this script from balanced 30-minute months, sampling-adjusted daily fluxes, and ERA5 annual budgets.",
  "- Lookup-filled daily fluxes are generated here only as a comparison to the scaled daily flux product.",
  "- Source probability and flux magnitude are modeled separately at 30-minute, daily, and annual scales.",
  "- 30-minute and daily temporal drivers are centered within site before modeling. This estimates within-site temporal effects separately from among-site differences.",
  "- Annual models use ERA5 annual anomalies within site plus ERA5/site mean and soil/ecosystem attributes. They are intentionally simpler linear models, with site-block bootstrap handling repeated-site dependence.",
  "- Variable importance uncertainty uses site-block bootstrap weights, not row-level resampling.",
  paste0("- Bootstrap replicates: ", n_boot, "."),
  "- For final inference, rerun with a larger bootstrap count, for example `NEON_DRIVER_BOOT_N=100 Rscript workflows/NEON.DriverScale.Analysis.R`.",
  "",
  "## Model-Ready Tables",
  paste0("- 30-minute: ", scale_counts$n_rows[scale_counts$scale == "30min"], " rows, ", scale_counts$n_sites[scale_counts$scale == "30min"], " sites."),
  paste0("- Daily: ", scale_counts$n_rows[scale_counts$scale == "daily"], " rows, ", scale_counts$n_sites[scale_counts$scale == "daily"], " sites."),
  paste0("- Annual ERA5: ", scale_counts$n_rows[scale_counts$scale == "annual"], " rows, ", scale_counts$n_sites[scale_counts$scale == "annual"], " sites."),
  "",
  "## Outputs",
  "- `OUTPUT/30min_site_behavior.csv`",
  "- `OUTPUT/30min_halfhour_balancing_bins.csv`",
  "- `OUTPUT/30min_total_vs_gradient_behavior_comparison.csv`",
  "- `OUTPUT/NEON_driver_scale_source_sink_summary.csv`",
  "- `OUTPUT/NEON_driver_scale_source_sink_counts.csv`",
  "- `OUTPUT/NEON_driver_scale_class_change_30min_to_daily.csv`",
  "- `OUTPUT/NEON_driver_scale_class_change_daily_to_annual.csv`",
  "- `OUTPUT/NEON_driver_scale_class_change_30min_to_annual.csv`",
  "- `OUTPUT/NEON_driver_daily_site_halfhour_lookup.csv`",
  "- `OUTPUT/NEON_driver_daily_lookup_grid.csv`",
  "- `OUTPUT/NEON_driver_daily_lookup_flux_all_sites.csv`",
  "- `OUTPUT/NEON_driver_daily_scaled_vs_lookup_comparison.csv`",
  "- `OUTPUT/NEON_driver_scale_30min_driver_data.csv`",
  "- `OUTPUT/NEON_driver_scale_daily_driver_data.csv`",
  "- `OUTPUT/NEON_driver_scale_annual_era5_driver_data.csv`",
  "- `OUTPUT/NEON_driver_scale_model_summary.csv`",
  "- `OUTPUT/NEON_driver_scale_variable_importance_bootstrap.csv`",
  "- `OUTPUT/NEON_driver_scale_variable_importance_summary.csv`",
  "- `OUTPUT/NEON_driver_scale_partial_effects.csv`",
  "- `OUTPUT/NEON_driver_scale_cross_scale_summary.csv`",
  "- `FIGURES/NEON_driver_scale_source_sink_counts.png`",
  "- `FIGURES/NEON_driver_scale_source_sink_transitions.png`",
  "- `FIGURES/NEON_driver_scale_variable_importance.png`",
  "- `FIGURES/NEON_driver_scale_partial_effects.png`"
)

# writeLines(summary_lines, output_path("NEON_driver_scale_analysis_results.md"))

message("NEON driver scale analysis complete (file outputs disabled).")
