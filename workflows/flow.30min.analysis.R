# NEON CH4 native-scale flux products.
#
# This script is intentionally limited to flux-rate and budget products:
#   1. standardized mean 30-minute total flux from balanced site-month-hour bins,
#   2. sampling-adjusted daily total fluxes,
#   3. annual budgets from gap-filled daily sums.
#
# Source/sink classification, source-probability models, magnitude models, and
# driver analyses are handled in NEON.DriverScale.Analysis.R.

library(tidyverse)
library(ggplot2)
library(patchwork)

localdir.ch4 <- Sys.getenv(
  "LOCALDIR_CH4",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/Methane"
)
localdir <- Sys.getenv(
  "LOCALDIR_FLUXGRADIENT",
  unset = "/Volumes/MaloneLab/Research/FluxGradient"
)

if (!dir.exists(localdir.ch4)) {
  stop("CH4 data directory does not exist: ", localdir.ch4)
}

setwd(localdir.ch4)
dir.create("OUTPUT", showWarnings = FALSE, recursive = TRUE)
dir.create("FIGURES", showWarnings = FALSE, recursive = TRUE)

ch4_input_file <- file.path(localdir.ch4, "SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata")
soil_input_file <- file.path(localdir.ch4, "Soildata_YearMon.Rdata")
metadata_file <- file.path(localdir, "Ameriflux_NEON field-sites.csv")
soil_biogeo_file <- file.path(localdir, "Soil_Biogeochem_RootBiomass.csv")
canopy_file <- file.path(localdir, "canopy_commbined.csv")
era5_mean_annual_file <- "OUTPUT/NEON_ERA5_gapfilled_mean_annual_budget.csv"
era5_annual_by_year_file <- "OUTPUT/NEON_ERA5_gapfilled_annual_budget_by_year.csv"

required_input_files <- c(ch4_input_file, soil_input_file, metadata_file)
missing_input_files <- required_input_files[!file.exists(required_input_files)]
if (length(missing_input_files) > 0) {
  stop("Missing required inputs: ", paste(missing_input_files, collapse = ", "))
}

load(ch4_input_file)
load(soil_input_file)

ch4_option <- "flux_total"
gradient_option <- "FG_ENSEMBLE_RSHP"
storage_option <- "storage_flux_filled"
flux_to_mgC_30min <- 2 * 0.0000288872 * 1000
seconds_per_30min <- 30 * 60
ug_c_per_umol_c <- 12.011
scale_levels <- c("30 min", "Daily", "Annual scaled")
source_threshold <- 0.75
behavior_levels <- c("Consistent sink", "Fluctuating", "Consistent source")
behavior_colors <- c(
  "Consistent sink" = "#2166AC",
  "Fluctuating" = "#4D4D4D",
  "Consistent source" = "#B2182B"
)
ecotype_colors <- c(
  "Cropland" = "#E69F00",
  "Forest" = "#009E73",
  "Grassland" = "#F0E442",
  "Shrubland" = "#CC79A7",
  "Wetland" = "#56B4E9"
)

mg_c_30min_to_umol_c_s <- function(x) {
  x * 1000 / ug_c_per_umol_c / seconds_per_30min
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

site_metadata <- read.csv(metadata_file) %>%
  mutate(
    EcoType = case_when(
      Vegetation.Abbreviation..IGBP. %in% c("ENF", "DBF", "MF", "EBF", "SAV") ~ "Forest",
      Vegetation.Abbreviation..IGBP. == "WET" ~ "Wetland",
      Vegetation.Abbreviation..IGBP. == "GRA" ~ "Grassland",
      Vegetation.Abbreviation..IGBP. %in% c("CVM", "CRO") ~ "Cropland",
      Vegetation.Abbreviation..IGBP. == "OSH" ~ "Shrubland",
      TRUE ~ NA_character_
    )
  ) %>%
  rename(
    SITE_ID = Site_Id.NEON,
    MAP = Mean.Average.Precipitation..mm.,
    MAT = Mean.Average.Tempurature..degrees.C.
  )

soil_biogeo <- if (file.exists(soil_biogeo_file)) {
  read.csv(soil_biogeo_file) %>% rename(SITE_ID = siteID)
} else {
  tibble(SITE_ID = character())
}

canopy_site <- if (file.exists(canopy_file)) {
  read.csv(canopy_file) %>%
    rename(SITE_ID = Site) %>%
    reframe(
      .by = SITE_ID,
      canopyHeight_m = mean(canopyHeight_m, na.rm = TRUE),
      LAI.mean = mean(LAI.mean, na.rm = TRUE),
      CHM.mean = mean(CHM.mean, na.rm = TRUE)
    )
} else {
  tibble(SITE_ID = character())
}

site_attributes <- site_metadata %>%
  dplyr::select(
    SITE_ID, EcoType, MAP, MAT, Vegetation.Abbreviation..IGBP.,
    Climate.Class.Abbreviation..Koeppen., Latitude..degrees.,
    Longitude..degrees.
  ) %>%
  left_join(soil_biogeo, by = "SITE_ID") %>%
  left_join(canopy_site, by = "SITE_ID")

ch4_30min <- purrr::imap_dfr(SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage, function(site_df, site) {
  site_df %>% mutate(SITE_ID = site)
}) %>%
  mutate(
    time.rounded = as.POSIXct(time.rounded, tz = "UTC"),
    YearMon = format(time.rounded, "%Y-%m"),
    Date = as.Date(time.rounded),
    Year = as.integer(format(time.rounded, "%Y")),
    month = as.integer(format(time.rounded, "%m")),
    doy = as.numeric(format(time.rounded, "%j")),
    hour_num = as.numeric(format(time.rounded, "%H")) + as.numeric(format(time.rounded, "%M")) / 60,
    hour_factor = factor(as.numeric(format(time.rounded, "%H")), levels = 0:23),
    sin_hour = sin(2 * pi * hour_num / 24),
    cos_hour = cos(2 * pi * hour_num / 24),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    CH4_raw = .data[[ch4_option]],
    CH4_total_mgC_30min = CH4_raw * flux_to_mgC_30min,
    CH4_gradient_mgC_30min = .data[[gradient_option]] * flux_to_mgC_30min,
    CH4_storage_mgC_30min = .data[[storage_option]] * flux_to_mgC_30min,
    CH4_mgC_30min = CH4_total_mgC_30min,
    source_30min = CH4_mgC_30min > 0,
    source_gradient_30min = CH4_gradient_mgC_30min > 0,
    storage_fraction = if_else(
      is.finite(CH4_total_mgC_30min) & abs(CH4_total_mgC_30min) > .Machine$double.eps,
      CH4_storage_mgC_30min / CH4_total_mgC_30min,
      NA_real_
    ),
    storage_abs_fraction = if_else(
      is.finite(CH4_total_mgC_30min) & abs(CH4_total_mgC_30min) > .Machine$double.eps,
      abs(CH4_storage_mgC_30min) / abs(CH4_total_mgC_30min),
      NA_real_
    ),
    PAR = pmax(PAR, 0),
    log_PAR = log1p(PAR)
  ) %>%
  left_join(Site_SoilData %>% rename(SITE_ID = Site), by = c("SITE_ID", "YearMon")) %>%
  left_join(site_attributes, by = "SITE_ID") %>%
  filter(
    count > 0,
    is.finite(CH4_mgC_30min),
    is.finite(hour_num),
    is.finite(Tair_C),
    is.finite(VSWCMean),
    is.finite(log_PAR),
    !is.na(season),
    !is.na(EcoType),
    !is.na(SITE_ID)
  ) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    EcoType = as.factor(EcoType)
  )

readr::write_csv(ch4_30min, "OUTPUT/30min_ch4_model_data.csv")

# Standardized mean 30-minute flux: each observed half-hour-of-day bin receives
# equal weight inside each site-month before site-level means are calculated.
monthly_observed_ch4 <- ch4_30min %>%
  reframe(
    .by = c(SITE_ID, YearMon),
    n_30min = n(),
    observed_mean_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    observed_median_CH4_mgC_30min = median(CH4_mgC_30min, na.rm = TRUE),
    observed_mean_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    observed_mean_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE),
    observed_prop_positive_30min = mean(CH4_mgC_30min > 0, na.rm = TRUE),
    mean_Tair_C = mean(Tair_C, na.rm = TRUE),
    mean_VSWC = mean(VSWCMean, na.rm = TRUE),
    mean_log_PAR = mean(log_PAR, na.rm = TRUE)
  )

halfhour_site_month_ch4 <- ch4_30min %>%
  reframe(
    .by = c(SITE_ID, YearMon, hour_num),
    n_30min_bin = n(),
    mean_CH4_mgC_30min_bin = mean(CH4_mgC_30min, na.rm = TRUE),
    median_CH4_mgC_30min_bin = median(CH4_mgC_30min, na.rm = TRUE),
    mean_gradient_mgC_30min_bin = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    median_gradient_mgC_30min_bin = median(CH4_gradient_mgC_30min, na.rm = TRUE),
    mean_storage_mgC_30min_bin = mean(CH4_storage_mgC_30min, na.rm = TRUE),
    median_storage_mgC_30min_bin = median(CH4_storage_mgC_30min, na.rm = TRUE),
    prop_positive_30min_bin = mean(CH4_mgC_30min > 0, na.rm = TRUE),
    prop_positive_gradient_30min_bin = mean(CH4_gradient_mgC_30min > 0, na.rm = TRUE),
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
    prop_positive_30min = mean(prop_positive_30min_bin, na.rm = TRUE),
    prop_positive_gradient_30min = mean(prop_positive_gradient_30min_bin, na.rm = TRUE),
    median_storage_abs_fraction = median(median_storage_abs_fraction_bin, na.rm = TRUE)
  ) %>%
  left_join(monthly_observed_ch4, by = c("SITE_ID", "YearMon"))

site_standardized_30min_flux <- monthly_site_ch4 %>%
  rename(
    monthly_mean_CH4_mgC_30min = mean_CH4_mgC_30min,
    monthly_mean_gradient_mgC_30min = mean_gradient_mgC_30min,
    monthly_mean_storage_mgC_30min = mean_storage_mgC_30min
  ) %>%
  reframe(
    .by = SITE_ID,
    n_months_30min = n(),
    mean_CH4_mgC_30min = mean(monthly_mean_CH4_mgC_30min, na.rm = TRUE),
    median_CH4_mgC_30min = median(monthly_mean_CH4_mgC_30min, na.rm = TRUE),
    sd_monthly_CH4_mgC_30min = sd(monthly_mean_CH4_mgC_30min, na.rm = TRUE),
    se_monthly_CH4_mgC_30min = sd_monthly_CH4_mgC_30min / sqrt(n_months_30min),
    mean_gradient_mgC_30min = mean(monthly_mean_gradient_mgC_30min, na.rm = TRUE),
    mean_storage_mgC_30min = mean(monthly_mean_storage_mgC_30min, na.rm = TRUE),
    median_storage_abs_fraction = median(median_storage_abs_fraction, na.rm = TRUE),
    prop_positive_months = mean(monthly_mean_CH4_mgC_30min > 0, na.rm = TRUE)
  ) %>%
  mutate(
    flux_30min_umolC_m2_s = mg_c_30min_to_umol_c_s(mean_CH4_mgC_30min),
    flux_30min_sd_umolC_m2_s = mg_c_30min_to_umol_c_s(sd_monthly_CH4_mgC_30min),
    flux_30min_se_umolC_m2_s = mg_c_30min_to_umol_c_s(se_monthly_CH4_mgC_30min),
    behavior_30min = case_when(
      prop_positive_months >= source_threshold ~ "Consistent source",
      prop_positive_months <= 1 - source_threshold ~ "Consistent sink",
      TRUE ~ "Fluctuating"
    ),
    behavior_30min = factor(behavior_30min, levels = behavior_levels)
  ) %>%
  left_join(site_attributes, by = "SITE_ID") %>%
  arrange(SITE_ID)

half_hour_bins <- sort(unique(ch4_30min$hour_num[is.finite(ch4_30min$hour_num)]))

ch4_daily_lookup_source <- ch4_30min %>%
  mutate(
    season_lookup = as.character(season),
    biseason_lookup = case_when(
      season_lookup %in% c("Spring", "Summer") ~ "Warm",
      season_lookup %in% c("Autumn", "Winter") ~ "Cool",
      TRUE ~ NA_character_
    )
  )

daily_date_attributes <- ch4_daily_lookup_source %>%
  reframe(
    .by = c(SITE_ID, Date),
    Year = first(Year),
    month = first(month),
    doy = first(doy),
    season = mode_value(season_lookup),
    biseason_lookup = mode_value(biseason_lookup),
    mean_Tair_C = mean(Tair_C, na.rm = TRUE),
    mean_VSWC = mean(VSWCMean, na.rm = TRUE),
    mean_log_PAR = mean(log_PAR, na.rm = TRUE),
    EcoType = first(as.character(EcoType))
  )

daily_observed_bins <- ch4_daily_lookup_source %>%
  reframe(
    .by = c(SITE_ID, Date, hour_num),
    n_observed_in_bin = n(),
    observed_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    observed_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    observed_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE)
  )

site_month_hour_lookup <- ch4_daily_lookup_source %>%
  reframe(
    .by = c(SITE_ID, month, hour_num),
    lookup_site_month_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    lookup_site_month_hour_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    lookup_site_month_hour_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE),
    n_site_month_hour_lookup = n()
  )

site_season_hour_lookup <- ch4_daily_lookup_source %>%
  reframe(
    .by = c(SITE_ID, season_lookup, hour_num),
    lookup_site_season_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    lookup_site_season_hour_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    lookup_site_season_hour_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE),
    n_site_season_hour_lookup = n()
  )

site_biseason_hour_lookup <- ch4_daily_lookup_source %>%
  filter(!is.na(biseason_lookup)) %>%
  reframe(
    .by = c(SITE_ID, biseason_lookup, hour_num),
    lookup_site_biseason_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    lookup_site_biseason_hour_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    lookup_site_biseason_hour_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE),
    n_site_biseason_hour_lookup = n()
  )

site_hour_lookup <- ch4_daily_lookup_source %>%
  reframe(
    .by = c(SITE_ID, hour_num),
    lookup_site_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    lookup_site_hour_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    lookup_site_hour_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE),
    n_site_hour_lookup = n()
  )

global_month_hour_lookup <- ch4_daily_lookup_source %>%
  reframe(
    .by = c(month, hour_num),
    lookup_global_month_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    lookup_global_month_hour_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    lookup_global_month_hour_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE)
  )

global_season_hour_lookup <- ch4_daily_lookup_source %>%
  reframe(
    .by = c(season_lookup, hour_num),
    lookup_global_season_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    lookup_global_season_hour_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    lookup_global_season_hour_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE)
  )

global_biseason_hour_lookup <- ch4_daily_lookup_source %>%
  filter(!is.na(biseason_lookup)) %>%
  reframe(
    .by = c(biseason_lookup, hour_num),
    lookup_global_biseason_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    lookup_global_biseason_hour_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    lookup_global_biseason_hour_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE)
  )

global_hour_lookup <- ch4_daily_lookup_source %>%
  reframe(
    .by = hour_num,
    lookup_global_hour_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    lookup_global_hour_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    lookup_global_hour_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE)
  )

site_lookup <- ch4_daily_lookup_source %>%
  reframe(
    .by = SITE_ID,
    lookup_site_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    lookup_site_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    lookup_site_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE)
  )

global_lookup <- ch4_daily_lookup_source %>%
  summarise(
    lookup_global_CH4_mgC_30min = mean(CH4_mgC_30min, na.rm = TRUE),
    lookup_global_gradient_mgC_30min = mean(CH4_gradient_mgC_30min, na.rm = TRUE),
    lookup_global_storage_mgC_30min = mean(CH4_storage_mgC_30min, na.rm = TRUE)
  )

daily_lookup_grid <- daily_date_attributes %>%
  dplyr::select(SITE_ID, Date, month, season_lookup = season, biseason_lookup) %>%
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
  bind_cols(global_lookup[rep(1, nrow(.)), ]) %>%
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
      lookup_global_CH4_mgC_30min
    ),
    filled_gradient_mgC_30min = coalesce(
      observed_gradient_mgC_30min,
      lookup_site_month_hour_gradient_mgC_30min,
      lookup_site_season_hour_gradient_mgC_30min,
      lookup_site_biseason_hour_gradient_mgC_30min,
      lookup_site_hour_gradient_mgC_30min,
      lookup_global_month_hour_gradient_mgC_30min,
      lookup_global_season_hour_gradient_mgC_30min,
      lookup_global_biseason_hour_gradient_mgC_30min,
      lookup_global_hour_gradient_mgC_30min,
      lookup_site_gradient_mgC_30min,
      lookup_global_gradient_mgC_30min
    ),
    filled_storage_mgC_30min = coalesce(
      observed_storage_mgC_30min,
      lookup_site_month_hour_storage_mgC_30min,
      lookup_site_season_hour_storage_mgC_30min,
      lookup_site_biseason_hour_storage_mgC_30min,
      lookup_site_hour_storage_mgC_30min,
      lookup_global_month_hour_storage_mgC_30min,
      lookup_global_season_hour_storage_mgC_30min,
      lookup_global_biseason_hour_storage_mgC_30min,
      lookup_global_hour_storage_mgC_30min,
      lookup_site_storage_mgC_30min,
      lookup_global_storage_mgC_30min
    ),
    was_observed = fill_source == "observed"
  )

daily_flux <- daily_lookup_grid %>%
  reframe(
    .by = c(SITE_ID, Date),
    n_30min = sum(was_observed, na.rm = TRUE),
    n_observed_30min = sum(was_observed, na.rm = TRUE),
    n_lookup_filled_30min = sum(!was_observed, na.rm = TRUE),
    n_site_month_hour_lookup_filled_30min = sum(fill_source == "site_month_hour_lookup", na.rm = TRUE),
    n_site_season_hour_lookup_filled_30min = sum(fill_source == "site_season_hour_lookup", na.rm = TRUE),
    n_site_biseason_hour_lookup_filled_30min = sum(fill_source == "site_biseason_hour_lookup", na.rm = TRUE),
    n_site_annual_hour_lookup_filled_30min = sum(fill_source == "site_annual_hour_lookup", na.rm = TRUE),
    n_global_lookup_filled_30min = sum(str_detect(fill_source, "^global"), na.rm = TRUE),
    daily_mgC_m2_day_observed = sum(observed_CH4_mgC_30min, na.rm = TRUE),
    daily_gradient_mgC_m2_day_observed = sum(observed_gradient_mgC_30min, na.rm = TRUE),
    daily_storage_mgC_m2_day_observed = sum(observed_storage_mgC_30min, na.rm = TRUE),
    daily_mgC_m2_day = sum(filled_CH4_mgC_30min, na.rm = TRUE),
    daily_gradient_mgC_m2_day = sum(filled_gradient_mgC_30min, na.rm = TRUE),
    daily_storage_mgC_m2_day = sum(filled_storage_mgC_30min, na.rm = TRUE),
    daily_mean_mgC_m2_30min = mean(filled_CH4_mgC_30min, na.rm = TRUE),
    daily_observed_mean_mgC_m2_30min = mean(observed_CH4_mgC_30min, na.rm = TRUE)
  ) %>%
  left_join(daily_date_attributes, by = c("SITE_ID", "Date")) %>%
  mutate(
    daily_fill_method = "observed_plus_month_season_biseason_annual_lookup",
    daily_observed_coverage = n_observed_30min / length(half_hour_bins)
  ) %>%
  arrange(SITE_ID, Date)

fill_source_levels <- c(
  "observed",
  "site_month_hour_lookup",
  "site_season_hour_lookup",
  "site_biseason_hour_lookup",
  "site_annual_hour_lookup",
  "global_month_hour_lookup",
  "global_season_hour_lookup",
  "global_biseason_hour_lookup",
  "global_annual_hour_lookup",
  "site_annual_mean_lookup",
  "global_lookup"
)

fill_source_labels <- c(
  observed = "Observed",
  site_month_hour_lookup = "Site-month-hour",
  site_season_hour_lookup = "Site-season-hour",
  site_biseason_hour_lookup = "Site-biseason-hour",
  site_annual_hour_lookup = "Site-annual-hour",
  global_month_hour_lookup = "All-site month-hour",
  global_season_hour_lookup = "All-site season-hour",
  global_biseason_hour_lookup = "All-site biseason-hour",
  global_annual_hour_lookup = "All-site annual-hour",
  site_annual_mean_lookup = "Site annual mean",
  global_lookup = "All-site mean"
)

fill_source_colors <- c(
  observed = "#222222",
  site_month_hour_lookup = "#1B9E77",
  site_season_hour_lookup = "#66A61E",
  site_biseason_hour_lookup = "#A6D854",
  site_annual_hour_lookup = "#B2DF8A",
  global_month_hour_lookup = "#7570B3",
  global_season_hour_lookup = "#A676C8",
  global_biseason_hour_lookup = "#CAB2D6",
  global_annual_hour_lookup = "#8DA0CB",
  site_annual_mean_lookup = "#E6AB02",
  global_lookup = "#D95F02"
)

fill_source_lookup <- tibble(
  fill_source = factor(fill_source_levels, levels = fill_source_levels),
  fill_source_label = unname(fill_source_labels)
)

daily_fill_source_summary <- daily_lookup_grid %>%
  mutate(
    fill_source = factor(fill_source, levels = fill_source_levels)
  ) %>%
  count(fill_source, name = "n_halfhour_slots") %>%
  right_join(fill_source_lookup, by = "fill_source") %>%
  mutate(
    n_halfhour_slots = replace_na(n_halfhour_slots, 0L),
    prop_halfhour_slots = n_halfhour_slots / sum(n_halfhour_slots),
    pct_halfhour_slots = 100 * prop_halfhour_slots
  )

daily_fill_source_by_site <- daily_lookup_grid %>%
  mutate(
    fill_source = factor(fill_source, levels = fill_source_levels)
  ) %>%
  count(SITE_ID, fill_source, name = "n_halfhour_slots") %>%
  left_join(fill_source_lookup, by = "fill_source") %>%
  group_by(SITE_ID) %>%
  mutate(
    total_halfhour_slots = sum(n_halfhour_slots),
    prop_halfhour_slots = n_halfhour_slots / total_halfhour_slots,
    pct_halfhour_slots = 100 * prop_halfhour_slots
  ) %>%
  ungroup()

daily_summary <- daily_flux %>%
  reframe(
    .by = SITE_ID,
    n_days_daily = n(),
    median_n_30min_per_day = median(n_observed_30min, na.rm = TRUE),
    min_n_30min_per_day = min(n_observed_30min, na.rm = TRUE),
    median_lookup_filled_30min_per_day = median(n_lookup_filled_30min, na.rm = TRUE),
    mean_daily_observed_coverage = mean(daily_observed_coverage, na.rm = TRUE),
    prop_positive_daily = mean(daily_mgC_m2_day > 0, na.rm = TRUE),
    flux_daily_mgC_m2_day = mean(daily_mgC_m2_day, na.rm = TRUE),
    flux_daily_sd_mgC_m2_day = sd(daily_mgC_m2_day, na.rm = TRUE),
    flux_daily_se_mgC_m2_day = flux_daily_sd_mgC_m2_day / sqrt(n_days_daily),
    flux_daily_median_mgC_m2_day = median(daily_mgC_m2_day, na.rm = TRUE),
    flux_daily_iqr_mgC_m2_day = IQR(daily_mgC_m2_day, na.rm = TRUE),
    flux_daily_mean_mgC_m2_30min = mean(daily_mean_mgC_m2_30min, na.rm = TRUE)
  ) %>%
  mutate(
    behavior_daily = case_when(
      prop_positive_daily >= source_threshold ~ "Consistent source",
      prop_positive_daily <= 1 - source_threshold ~ "Consistent sink",
      TRUE ~ "Fluctuating"
    ),
    behavior_daily = factor(behavior_daily, levels = behavior_levels)
  ) %>%
  left_join(site_attributes, by = "SITE_ID") %>%
  arrange(SITE_ID)

scaled_daily_annual_budget_by_year <- daily_flux %>%
  reframe(
    .by = c(SITE_ID, Year),
    n_observed_days = n(),
    mean_daily_observed_coverage = mean(daily_observed_coverage, na.rm = TRUE),
    annual_budget_gC_m2_yr = sum(daily_mgC_m2_day, na.rm = TRUE) / 1000,
    annual_gradient_budget_gC_m2_yr = sum(daily_gradient_mgC_m2_day, na.rm = TRUE) / 1000,
    annual_storage_budget_gC_m2_yr = sum(daily_storage_mgC_m2_day, na.rm = TRUE) / 1000
  ) %>%
  left_join(site_attributes, by = "SITE_ID") %>%
  arrange(SITE_ID, Year)

scaled_daily_annual_budget_summary <- scaled_daily_annual_budget_by_year %>%
  reframe(
    .by = SITE_ID,
    n_years_scaled_daily = n(),
    prop_positive_scaled_daily_annual = mean(annual_budget_gC_m2_yr > 0, na.rm = TRUE),
    flux_scaled_daily_annual_gC_m2_yr = mean(annual_budget_gC_m2_yr, na.rm = TRUE),
    flux_scaled_daily_annual_sd_gC_m2_yr = sd(annual_budget_gC_m2_yr, na.rm = TRUE),
    flux_scaled_daily_annual_se_gC_m2_yr = flux_scaled_daily_annual_sd_gC_m2_yr / sqrt(n_years_scaled_daily),
    mean_observed_days_per_year = mean(n_observed_days, na.rm = TRUE),
    mean_daily_observed_coverage = mean(mean_daily_observed_coverage, na.rm = TRUE)
  ) %>%
  mutate(
    behavior_annual_scaled = case_when(
      prop_positive_scaled_daily_annual >= source_threshold ~ "Consistent source",
      prop_positive_scaled_daily_annual <= 1 - source_threshold ~ "Consistent sink",
      TRUE ~ "Fluctuating"
    ),
    behavior_annual_scaled = factor(behavior_annual_scaled, levels = behavior_levels)
  ) %>%
  left_join(site_attributes, by = "SITE_ID") %>%
  arrange(SITE_ID)

era5_annual_budget_summary <- if (file.exists(era5_mean_annual_file)) {
  read.csv(era5_mean_annual_file) %>%
    mutate(SITE_ID = as.character(SITE_ID)) %>%
    transmute(
      SITE_ID,
      n_years_era5 = n_years,
      prop_positive_era5_annual = era5_prop_source_years,
      flux_era5_annual_gC_m2_yr = mean_era5_gapfilled_annual_budget_gC_m2_yr,
      flux_era5_annual_sd_gC_m2_yr = sd_era5_gapfilled_annual_budget_gC_m2_yr,
      flux_era5_annual_se_gC_m2_yr = sd_era5_gapfilled_annual_budget_gC_m2_yr / sqrt(n_years),
      era5_annual_behavior = factor(era5_annual_behavior, levels = behavior_levels),
      mean_observed_coverage,
      n_training_obs
    ) %>%
    left_join(site_attributes, by = "SITE_ID") %>%
    arrange(SITE_ID)
} else if (file.exists(era5_annual_by_year_file)) {
  read.csv(era5_annual_by_year_file) %>%
    mutate(SITE_ID = as.character(SITE_ID)) %>%
    reframe(
      .by = SITE_ID,
      n_years_era5 = n(),
      prop_positive_era5_annual = mean(annual_budget_gC_m2_yr > 0, na.rm = TRUE),
      flux_era5_annual_gC_m2_yr = mean(annual_budget_gC_m2_yr, na.rm = TRUE),
      flux_era5_annual_sd_gC_m2_yr = sd(annual_budget_gC_m2_yr, na.rm = TRUE),
      flux_era5_annual_se_gC_m2_yr = flux_era5_annual_sd_gC_m2_yr / sqrt(n_years_era5),
      era5_annual_behavior = case_when(
        prop_positive_era5_annual >= source_threshold ~ "Consistent source",
        prop_positive_era5_annual <= 1 - source_threshold ~ "Consistent sink",
        TRUE ~ "Fluctuating"
      ),
      mean_observed_coverage = mean(observed_coverage, na.rm = TRUE),
      n_training_obs = max(n_training_obs, na.rm = TRUE)
    ) %>%
    mutate(era5_annual_behavior = factor(era5_annual_behavior, levels = behavior_levels)) %>%
    left_join(site_attributes, by = "SITE_ID") %>%
    arrange(SITE_ID)
} else {
  warning(
    "ERA5 annual budget files not found. Run NEON.ERA5.HalfHourlyGapfill.R ",
    "to populate annual budget outputs."
  )
  tibble(
    SITE_ID = character(),
    n_years_era5 = integer(),
    prop_positive_era5_annual = numeric(),
    flux_era5_annual_gC_m2_yr = numeric(),
    flux_era5_annual_sd_gC_m2_yr = numeric(),
    flux_era5_annual_se_gC_m2_yr = numeric(),
    era5_annual_behavior = factor(levels = behavior_levels)
  )
}

annual_budget_summary <- scaled_daily_annual_budget_summary

scale_site_summary <- site_standardized_30min_flux %>%
  dplyr::select(
    SITE_ID, EcoType, n_months_30min, prop_positive_months,
    flux_30min_umolC_m2_s, flux_30min_sd_umolC_m2_s, flux_30min_se_umolC_m2_s,
    behavior_30min
  ) %>%
  left_join(
    daily_summary %>%
      dplyr::select(
        SITE_ID, n_days_daily, prop_positive_daily,
        flux_daily_mgC_m2_day, flux_daily_sd_mgC_m2_day,
        flux_daily_se_mgC_m2_day, flux_daily_median_mgC_m2_day,
        behavior_daily
      ),
    by = "SITE_ID"
  ) %>%
  left_join(
    annual_budget_summary %>%
      dplyr::select(
        SITE_ID, n_years_scaled_daily, prop_positive_scaled_daily_annual,
        flux_scaled_daily_annual_gC_m2_yr, flux_scaled_daily_annual_sd_gC_m2_yr,
        flux_scaled_daily_annual_se_gC_m2_yr, behavior_annual_scaled
      ),
    by = "SITE_ID"
  ) %>%
  arrange(SITE_ID)

scale_long_summary <- scale_site_summary %>%
  transmute(
    SITE_ID,
    EcoType,
    `30 min` = flux_30min_umolC_m2_s,
    Daily = flux_daily_mgC_m2_day,
    `Annual scaled` = flux_scaled_daily_annual_gC_m2_yr,
    sd_30min = flux_30min_sd_umolC_m2_s,
    sd_daily = flux_daily_sd_mgC_m2_day,
    sd_annual = flux_scaled_daily_annual_sd_gC_m2_yr,
    se_30min = flux_30min_se_umolC_m2_s,
    se_daily = flux_daily_se_mgC_m2_day,
    se_annual = flux_scaled_daily_annual_se_gC_m2_yr,
    behavior_30min,
    behavior_daily,
    behavior_annual_scaled
  ) %>%
  pivot_longer(cols = all_of(scale_levels), names_to = "scale", values_to = "flux_native") %>%
  mutate(
    flux_sd_native = case_when(
      scale == "30 min" ~ sd_30min,
      scale == "Daily" ~ sd_daily,
      scale == "Annual scaled" ~ sd_annual
    ),
    flux_se_native = case_when(
      scale == "30 min" ~ se_30min,
      scale == "Daily" ~ se_daily,
      scale == "Annual scaled" ~ se_annual
    ),
    flux_unit = case_when(
      scale == "30 min" ~ "umol C m-2 s-1",
      scale == "Daily" ~ "mg C m-2 d-1",
      scale == "Annual scaled" ~ "g C m-2 yr-1"
    ),
    behavior = case_when(
      scale == "30 min" ~ as.character(behavior_annual_scaled),
      scale == "Daily" ~ as.character(behavior_annual_scaled),
      scale == "Annual scaled" ~ as.character(behavior_annual_scaled)
    ),
    behavior = factor(behavior, levels = behavior_levels),
    flux_lower_native = flux_native - flux_sd_native,
    flux_upper_native = flux_native + flux_sd_native,
    scale = factor(scale, levels = scale_levels)
  )

site_diel_30min <- ch4_30min %>%
  left_join(
    annual_budget_summary %>%
      dplyr::select(SITE_ID, behavior_annual_scaled),
    by = "SITE_ID"
  ) %>%
  reframe(
    .by = c(SITE_ID, behavior_annual_scaled, EcoType, hour_num),
    n_30min = n(),
    flux_umolC_m2_s = mg_c_30min_to_umol_c_s(mean(CH4_mgC_30min, na.rm = TRUE)),
    gradient_umolC_m2_s = mg_c_30min_to_umol_c_s(mean(CH4_gradient_mgC_30min, na.rm = TRUE)),
    storage_umolC_m2_s = mg_c_30min_to_umol_c_s(mean(CH4_storage_mgC_30min, na.rm = TRUE)),
    source_probability = mean(CH4_mgC_30min > 0, na.rm = TRUE)
  ) %>%
  filter(!is.na(behavior_annual_scaled), is.finite(hour_num))

diel_behavior_summary <- site_diel_30min %>%
  reframe(
    .by = c(behavior_annual_scaled, hour_num),
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
    behavior_annual_scaled = factor(behavior_annual_scaled, levels = behavior_levels)
  )

diel_component_summary <- site_diel_30min %>%
  dplyr::select(SITE_ID, behavior_annual_scaled, hour_num, gradient_umolC_m2_s, storage_umolC_m2_s) %>%
  pivot_longer(
    cols = c(gradient_umolC_m2_s, storage_umolC_m2_s),
    names_to = "component",
    values_to = "component_flux_umolC_m2_s"
  ) %>%
  mutate(
    component = recode(
      component,
      gradient_umolC_m2_s = "Gradient",
      storage_umolC_m2_s = "Storage"
    )
  ) %>%
  reframe(
    .by = c(behavior_annual_scaled, component, hour_num),
    n_sites = n_distinct(SITE_ID),
    mean_component_flux_umolC_m2_s = mean(component_flux_umolC_m2_s, na.rm = TRUE),
    sd_component_flux_umolC_m2_s = sd(component_flux_umolC_m2_s, na.rm = TRUE),
    se_component_flux_umolC_m2_s = sd_component_flux_umolC_m2_s / sqrt(n_sites)
  ) %>%
  mutate(
    se_component_flux_umolC_m2_s = replace_na(se_component_flux_umolC_m2_s, 0),
    behavior_annual_scaled = factor(behavior_annual_scaled, levels = behavior_levels),
    component = factor(component, levels = c("Gradient", "Storage"))
  )

annual_behavior_site_counts <- annual_budget_summary %>%
  transmute(
    SITE_ID,
    behavior = factor(behavior_annual_scaled, levels = behavior_levels),
    EcoType = as.character(EcoType)
  ) %>%
  distinct()

annual_behavior_counts <- annual_behavior_site_counts %>%
  count(behavior, name = "n_sites") %>%
  complete(
    behavior = factor(behavior_levels, levels = behavior_levels),
    fill = list(n_sites = 0L)
  )

annual_behavior_ecotype_counts <- annual_behavior_site_counts %>%
  filter(!is.na(EcoType)) %>%
  count(behavior, EcoType, name = "n_sites") %>%
  complete(
    behavior = factor(behavior_levels, levels = behavior_levels),
    EcoType,
    fill = list(n_sites = 0L)
  )

annual_site_order <- scale_site_summary %>%
  mutate(behavior_annual_scaled = factor(behavior_annual_scaled, levels = behavior_levels)) %>%
  arrange(behavior_annual_scaled, flux_scaled_daily_annual_gC_m2_yr, SITE_ID) %>%
  pull(SITE_ID)

all_site_flux_magnitude_summary <- scale_long_summary %>%
  mutate(
    annual_behavior = factor(behavior_annual_scaled, levels = behavior_levels),
    SITE_ID_plot = factor(SITE_ID, levels = rev(annual_site_order)),
    scale_label = factor(
      paste0(scale, "\n", flux_unit),
      levels = paste0(scale_levels, "\n", c("umol C m-2 s-1", "mg C m-2 d-1", "g C m-2 yr-1"))
    )
  ) %>%
  filter(is.finite(flux_native), !is.na(annual_behavior))

annual_site_map_data <- annual_budget_summary %>%
  transmute(
    SITE_ID,
    EcoType,
    annual_behavior = factor(behavior_annual_scaled, levels = behavior_levels),
    annual_flux_gC_m2_yr = flux_scaled_daily_annual_gC_m2_yr,
    annual_flux_magnitude_gC_m2_yr = abs(flux_scaled_daily_annual_gC_m2_yr),
    latitude = Latitude..degrees.,
    longitude = Longitude..degrees.
  ) %>%
  filter(
    is.finite(latitude),
    is.finite(longitude),
    is.finite(annual_flux_gC_m2_yr),
    !is.na(annual_behavior)
  )

readr::write_csv(monthly_site_ch4, "OUTPUT/30min_monthly_site_ch4.csv")
readr::write_csv(halfhour_site_month_ch4, "OUTPUT/30min_halfhour_balancing_bins.csv")
readr::write_csv(site_standardized_30min_flux, "OUTPUT/30min_site_standardized_flux.csv")
readr::write_csv(site_month_hour_lookup, "OUTPUT/NEON_daily_site_month_halfhour_lookup.csv")
readr::write_csv(daily_lookup_grid, "OUTPUT/NEON_scale_daily_flux_lookup_grid.csv")
readr::write_csv(daily_flux, "OUTPUT/NEON_scale_daily_flux_all_sites.csv")
readr::write_csv(daily_fill_source_summary, "OUTPUT/NEON_daily_fill_source_summary.csv")
readr::write_csv(daily_fill_source_by_site, "OUTPUT/NEON_daily_fill_source_by_site.csv")
readr::write_csv(daily_summary, "OUTPUT/NEON_scale_daily_flux_summary.csv")
readr::write_csv(scaled_daily_annual_budget_by_year, "OUTPUT/NEON_scale_scaled_daily_annual_budget_by_year.csv")
readr::write_csv(scaled_daily_annual_budget_summary, "OUTPUT/NEON_scale_scaled_daily_annual_budget_summary.csv")
readr::write_csv(era5_annual_budget_summary, "OUTPUT/NEON_scale_ERA5_annual_budget_summary.csv")
readr::write_csv(annual_budget_summary, "OUTPUT/NEON_scale_annual_budget_summary.csv")
readr::write_csv(scale_site_summary, "OUTPUT/NEON_scale_site_flux_budget_summary.csv")
readr::write_csv(scale_long_summary, "OUTPUT/NEON_scale_long_flux_budget_summary.csv")
readr::write_csv(site_diel_30min, "OUTPUT/NEON_site_diel_30min_by_behavior.csv")
readr::write_csv(diel_behavior_summary, "OUTPUT/NEON_diel_30min_behavior_summary.csv")
readr::write_csv(diel_component_summary, "OUTPUT/NEON_diel_30min_component_behavior_summary.csv")
readr::write_csv(annual_behavior_counts, "OUTPUT/NEON_annual_behavior_site_counts.csv")
readr::write_csv(annual_behavior_ecotype_counts, "OUTPUT/NEON_annual_behavior_ecotype_counts.csv")
readr::write_csv(all_site_flux_magnitude_summary, "OUTPUT/NEON_all_site_flux_magnitude_summary.csv")
readr::write_csv(annual_site_map_data, "OUTPUT/NEON_annual_site_map_data.csv")
readr::write_csv(
  scale_long_summary %>%
    count(scale, behavior, name = "n_sites") %>%
    complete(
      scale,
      behavior = factor(behavior_levels, levels = behavior_levels),
      fill = list(n_sites = 0L)
    ),
  "OUTPUT/NEON_scale_behavior_counts.csv"
)

# Backward-compatible names for plotting/table consumers.
readr::write_csv(scale_site_summary, "OUTPUT/NEON_scale_site_comparison.csv")
readr::write_csv(scale_long_summary, "OUTPUT/NEON_scale_long_comparison.csv")

plot_scale_flux <- scale_long_summary %>%
  filter(is.finite(flux_native)) %>%
  mutate(SITE_ID_plot = fct_reorder(SITE_ID, flux_native, .na_rm = TRUE)) %>%
  ggplot(aes(x = flux_native, y = SITE_ID_plot, color = behavior)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_errorbar(
    aes(xmin = flux_lower_native, xmax = flux_upper_native),
    orientation = "y",
    width = 0.22,
    linewidth = 0.55,
    color = "grey30",
    alpha = 0.75,
    na.rm = TRUE
  ) +
  geom_point(size = 2.2, alpha = 0.9) +
  facet_wrap(~ scale + flux_unit, scales = "free_x", nrow = 1) +
  scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  theme_bw(base_size = 10.5) +
  labs(
    x = "Flux in native units",
    y = NULL,
    color = "Behavior class",
    title = "NEON CH4 Native-Scale Flux Products",
    subtitle = "30-min values are standardized by balanced site-month-hour bins; daily values use month/season/biseason/annual lookup fills. Bars show +/- 1 SD."
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 8.5)
  )

ggsave("FIGURES/NEON_scale_native_flux_by_site.png", plot_scale_flux, width = 12, height = 8, units = "in", dpi = 300)
ggsave("FIGURES/NEON_scale_native_flux_by_site.pdf", plot_scale_flux, width = 12, height = 8, units = "in")

plot_flux_by_behavior <- scale_long_summary %>%
  filter(is.finite(flux_native), !is.na(behavior)) %>%
  ggplot(aes(x = behavior, y = flux_native, color = behavior)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_boxplot(outlier.shape = NA, width = 0.55, color = "grey35", fill = "grey92") +
  geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.75) +
  facet_wrap(~ scale + flux_unit, scales = "free_y", nrow = 1) +
  scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  theme_bw(base_size = 10) +
  labs(
    x = NULL,
    y = "Flux in native units",
    title = "A. Native-scale flux distributions"
  ) +
  guides(color = "none") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 25, hjust = 1),
    strip.text = element_text(face = "bold", size = 7.8)
  )

plot_diel_flux <- diel_behavior_summary %>%
  ggplot(aes(x = hour_num, y = mean_flux_umolC_m2_s, color = behavior_annual_scaled, fill = behavior_annual_scaled)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_ribbon(
    aes(ymin = mean_flux_umolC_m2_s - se_flux_umolC_m2_s, ymax = mean_flux_umolC_m2_s + se_flux_umolC_m2_s),
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
    color = "Annual behavior class",
    title = "B. Mean diel flux pattern"
  ) +
  guides(fill = "none") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 10)
  )

plot_diel_source_probability <- diel_behavior_summary %>%
  ggplot(aes(x = hour_num, y = mean_source_probability, color = behavior_annual_scaled, fill = behavior_annual_scaled)) +
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
  theme_bw(base_size = 10) +
  labs(
    x = "Hour of day",
    y = "Probability of positive flux",
    color = "Annual behavior class",
    title = "C. Diel source probability"
  ) +
  guides(fill = "none") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 10)
  )

plot_annual_behavior_counts <- annual_behavior_counts %>%
  ggplot(aes(x = n_sites, y = fct_rev(behavior), fill = behavior)) +
  geom_col(width = 0.68, color = "grey30", linewidth = 0.25) +
  geom_text(aes(label = n_sites), hjust = -0.25, size = 3.1) +
  scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16)), breaks = scales::breaks_width(5)) +
  theme_bw(base_size = 10) +
  labs(
    x = "Number of sites",
    y = NULL,
    title = "D. Sites per annual class"
  ) +
  guides(fill = "none") +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8.7)
  )

plot_annual_behavior_ecotypes <- annual_behavior_ecotype_counts %>%
  ggplot(aes(x = n_sites, y = fct_rev(behavior), fill = EcoType)) +
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
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8)
  )

plot_site_count_panel <- plot_annual_behavior_counts | plot_annual_behavior_ecotypes

plot_flux_pattern_diel_panel <- (
  plot_flux_by_behavior / ((plot_diel_flux | plot_diel_source_probability) / plot_site_count_panel)
) +
  plot_layout(heights = c(0.95, 1.7), guides = "collect") +
  plot_annotation(
    title = "NEON CH4 Flux Patterns By Behavior Class",
    subtitle = "Categories are defined from annual fluxes; finer-scale panels show 30-min and daily behavior within those annual classes. Diel ribbons show +/- 1 SE among sites."
  ) &
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("FIGURES/NEON_flux_pattern_diel_behavior_panel.png", plot_flux_pattern_diel_panel, width = 13, height = 10, units = "in", dpi = 300)
ggsave("FIGURES/NEON_flux_pattern_diel_behavior_panel.pdf", plot_flux_pattern_diel_panel, width = 13, height = 10, units = "in")

plot_all_site_flux_magnitude <- all_site_flux_magnitude_summary %>%
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
    x = "Flux magnitude in native units",
    y = NULL,
    shape = "Ecosystem type",
    title = "NEON CH4 Site Categories And Flux Magnitudes",
    subtitle = "Rows define annual flux categories; points show each site's 30-min, daily, and annual flux magnitude. Bars show +/- 1 SD."
  ) +
  guides(
    color = "none",
    shape = guide_legend(override.aes = list(size = 3.7, alpha = 1))
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 8.8),
    legend.text = element_text(size = 8.2),
    axis.text.y = element_text(size = 7.1),
    axis.text.x = element_text(size = 7.7),
    strip.text.x = element_text(face = "bold", size = 7.9, lineheight = 0.95),
    strip.text.y = element_text(face = "bold", size = 8.1),
    panel.spacing.x = unit(0.55, "lines"),
    panel.spacing.y = unit(0.35, "lines")
  )

ggsave("FIGURES/NEON_all_site_category_flux_magnitudes.png", plot_all_site_flux_magnitude, width = 12.5, height = 10.5, units = "in", dpi = 300)
ggsave("FIGURES/NEON_all_site_category_flux_magnitudes.pdf", plot_all_site_flux_magnitude, width = 12.5, height = 10.5, units = "in")

north_america_map <- rnaturalearth::ne_countries(
  continent = "North America",
  returnclass = "sf"
)

plot_annual_site_map <- ggplot() +
  geom_sf(data = north_america_map, fill = "grey94", color = "white", linewidth = 0.25) +
  geom_point(
    data = annual_site_map_data,
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
    name = expression(paste("|Annual flux| (g C ", m^-2, " ", yr^-1, ")"))
  ) +
  coord_sf(xlim = c(-170, -60), ylim = c(15, 72), expand = FALSE) +
  theme_bw(base_size = 10) +
  labs(
    x = NULL,
    y = NULL,
    color = "Annual behavior class",
    shape = "Ecosystem type",
    title = "NEON CH4 Annual Flux Categories Across Sites",
    subtitle = "Color shows annual behavior class; symbol size shows absolute annual flux magnitude."
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

ggsave("FIGURES/NEON_annual_site_category_map.png", plot_annual_site_map, width = 10.5, height = 8, units = "in", dpi = 300)
ggsave("FIGURES/NEON_annual_site_category_map.pdf", plot_annual_site_map, width = 10.5, height = 8, units = "in")

plot_daily_fill_source_frequency <- daily_fill_source_summary %>%
  mutate(fill_source_label = fct_reorder(fill_source_label, pct_halfhour_slots)) %>%
  ggplot(aes(x = pct_halfhour_slots, y = fill_source_label, fill = fill_source)) +
  geom_col(width = 0.72, color = "grey30", linewidth = 0.2) +
  geom_text(
    aes(label = if_else(pct_halfhour_slots >= 0.1, sprintf("%.1f%%", pct_halfhour_slots), "<0.1%")),
    hjust = -0.12,
    size = 3
  ) +
  scale_fill_manual(values = fill_source_colors, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.16))
  ) +
  theme_bw(base_size = 10) +
  labs(
    x = "Half-hour slots used in daily flux calculation",
    y = NULL,
    title = "Frequency Of Data Sources Used For Daily CH4 Gap Filling",
    subtitle = "Observed half-hour fluxes are used first; missing slots are filled by increasingly broad lookup means."
  ) +
  guides(fill = "none") +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8.8)
  )

plot_daily_fill_source_by_site <- daily_fill_source_by_site %>%
  left_join(
    annual_budget_summary %>%
      transmute(SITE_ID, annual_behavior = factor(behavior_annual_scaled, levels = behavior_levels)),
    by = "SITE_ID"
  ) %>%
  mutate(
    SITE_ID_plot = fct_reorder(SITE_ID, prop_halfhour_slots, .fun = function(x) sum(x, na.rm = TRUE)),
    fill_source_label = factor(fill_source_label, levels = unname(fill_source_labels))
  ) %>%
  ggplot(aes(x = SITE_ID_plot, y = pct_halfhour_slots, fill = fill_source)) +
  geom_col(width = 0.78, color = NA) +
  facet_grid(~ annual_behavior, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = fill_source_colors,
    breaks = fill_source_levels,
    labels = fill_source_labels,
    drop = FALSE,
    na.translate = FALSE
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_size = 9.4) +
  labs(
    x = NULL,
    y = "Half-hour slots",
    fill = "Data source",
    title = "Daily Gap-Fill Data Sources By Site",
    subtitle = "Bars show the fraction of site-date half-hour slots that were observed or filled from each lookup level."
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 7.7),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.8),
    strip.text = element_text(face = "bold", size = 8.2)
  )

plot_daily_fill_source_panel <- plot_daily_fill_source_frequency / plot_daily_fill_source_by_site +
  plot_layout(heights = c(0.85, 1.15))

ggsave("FIGURES/NEON_daily_gapfill_source_frequency.png", plot_daily_fill_source_panel, width = 12, height = 10, units = "in", dpi = 300)
ggsave("FIGURES/NEON_daily_gapfill_source_frequency.pdf", plot_daily_fill_source_panel, width = 12, height = 10, units = "in")

plot_daily_distribution <- daily_flux %>%
  ggplot(aes(x = daily_mgC_m2_day, y = fct_reorder(SITE_ID, daily_mgC_m2_day, median, .na_rm = TRUE))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_boxplot(outlier.alpha = 0.12, width = 0.55, fill = "grey90", color = "grey30") +
  theme_bw(base_size = 10.5) +
  labs(
    x = "Gap-filled daily CH4 flux (mg C m-2 d-1)",
    y = NULL,
    title = "Daily Total CH4 Flux By Site",
    subtitle = "Missing half-hour slots are filled from site month-hour, season-hour, biseason-hour, then annual-hour means."
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave("FIGURES/NEON_daily_flux_by_site.png", plot_daily_distribution, width = 8, height = 9, units = "in", dpi = 300)
ggsave("FIGURES/NEON_daily_flux_by_site.pdf", plot_daily_distribution, width = 8, height = 9, units = "in")

summary_lines <- c(
  "# NEON CH4 Native-Scale Flux Products",
  "",
  "## Scope",
  "- This script creates standardized mean 30-minute fluxes, gap-filled daily fluxes, and annual budgets from those daily fluxes.",
  "- Source/sink classification and driver models are produced by `NEON.DriverScale.Analysis.R`.",
  "",
  "## Outputs",
  "- `OUTPUT/30min_ch4_model_data.csv`",
  "- `OUTPUT/30min_halfhour_balancing_bins.csv`",
  "- `OUTPUT/30min_monthly_site_ch4.csv`",
  "- `OUTPUT/30min_site_standardized_flux.csv`",
  "- `OUTPUT/NEON_daily_site_month_halfhour_lookup.csv`",
  "- `OUTPUT/NEON_scale_daily_flux_lookup_grid.csv`",
  "- `OUTPUT/NEON_scale_daily_flux_all_sites.csv`",
  "- `OUTPUT/NEON_scale_daily_flux_summary.csv`",
  "- `OUTPUT/NEON_scale_scaled_daily_annual_budget_by_year.csv`",
  "- `OUTPUT/NEON_scale_scaled_daily_annual_budget_summary.csv`",
  "- `OUTPUT/NEON_scale_ERA5_annual_budget_summary.csv`",
  "- `OUTPUT/NEON_scale_annual_budget_summary.csv`",
  "- `OUTPUT/NEON_scale_site_flux_budget_summary.csv`",
  "- `OUTPUT/NEON_scale_long_flux_budget_summary.csv`",
  "- `OUTPUT/NEON_scale_behavior_counts.csv`",
  "- `FIGURES/NEON_scale_native_flux_by_site.png`",
  "- `FIGURES/NEON_daily_flux_by_site.png`"
)

writeLines(summary_lines, "OUTPUT/30min_flux_products_results.md")

message("Wrote NEON native-scale flux products.")
