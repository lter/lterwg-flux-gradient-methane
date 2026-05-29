# ERA5-based upland CH4 two-stage upscaling workflow.
#
# Goal:
#   1. classify non-wetland/non-inundated upland NEON sites as weak sinks,
#      weak sources, or fluctuating using ERA5 gapfilled flux rates;
#   2. evaluate which temporal resolution is most appropriate for that
#      classification;
#   3. estimate annual flux rates conditional on ecosystem type and exchange
#      class; and
#   4. demonstrate an area-weighted annual upscaling compared with the Global
#      Methane Budget soil sink.
#
# Global budget reference used here:
#   Global Methane Budget 2000-2020 (ESSD, 2025) reports soil uptake of
#   35 [35, 36] Tg CH4 yr-1. The script treats this as -35 Tg CH4 yr-1
#   because the NEON sign convention is positive = emission, negative = uptake.

library(tidyverse)
library(ggplot2)
library(data.table)

localdir.ch4 <- Sys.getenv(
  "LOCALDIR_CH4",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/Methane"
)

if (!dir.exists(localdir.ch4)) {
  stop("CH4 data directory does not exist: ", localdir.ch4)
}

setwd(localdir.ch4)
dir.create("OUTPUT", showWarnings = FALSE, recursive = TRUE)
dir.create("FIGURES", showWarnings = FALSE, recursive = TRUE)

era5_30min_file <- "OUTPUT/NEON_ERA5_gapfilled_30min.csv.gz"
era5_annual_file <- "OUTPUT/NEON_ERA5_gapfilled_annual_budget_by_year.csv"
site_behavior_file <- "OUTPUT/30min_site_behavior.csv"
area_assumption_file <- "OUTPUT/ERA_Upscaling_ecosystem_area_assumptions.csv"

required_files <- c(era5_30min_file, era5_annual_file, site_behavior_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required ERA5/native files: ", paste(missing_files, collapse = ", "))
}

behavior_levels <- c("Weak sink", "Fluctuating", "Weak source")
source_threshold <- 0.75
sink_threshold <- 0.25
gmb_soil_sink_tg_ch4_yr <- -35
gmb_soil_sink_low_tg_ch4_yr <- -36
gmb_soil_sink_high_tg_ch4_yr <- -35
source_probability_threshold_default <- 0.80

exchange_colors <- c(
  "Weak sink" = "#2166AC",
  "Fluctuating" = "#4D4D4D",
  "Weak source" = "#B2182B"
)

classify_exchange <- function(prop_positive) {
  case_when(
    is.na(prop_positive) ~ NA_character_,
    prop_positive >= source_threshold ~ "Weak source",
    prop_positive <= sink_threshold ~ "Weak sink",
    TRUE ~ "Fluctuating"
  )
}

calc_entropy <- function(p) {
  p <- p[p > 0 & is.finite(p)]
  if (length(p) == 0) return(NA_real_)
  -sum(p * log(p))
}

annualize_30min_gC <- function(mean_mgC_30min) {
  mean_mgC_30min * 48 * 365 / 1000
}

annualize_daily_gC <- function(mean_mgC_day) {
  mean_mgC_day * 365 / 1000
}

gC_m2_yr_to_tg_ch4 <- function(flux_gC_m2_yr, area_mha) {
  flux_gC_m2_yr * area_mha * 0.0133333333333333
}

interval_se <- function(low, high) {
  se <- (high - low) / (2 * 1.96)
  ifelse(is.finite(se) & se > 0, se, NA_real_)
}

site_behavior <- read.csv(site_behavior_file) %>%
  mutate(SITE_ID = as.character(SITE_ID))

upland_sites <- site_behavior %>%
  filter(
    !is.na(EcoType),
    !str_detect(EcoType, regex("wetland|inundat|flood|marsh|swamp|bog|fen|lake|rice", ignore_case = TRUE))
  ) %>%
  distinct(SITE_ID, EcoType, MAP, MAT)

if (nrow(upland_sites) == 0) {
  stop("No non-wetland/non-inundated sites available after filtering.")
}

era5_30min <- data.table::fread(era5_30min_file) %>%
  as_tibble() %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    Date = as.Date(Date),
    Year = as.integer(Year),
    gapfilled_CH4_mgC_30min = as.numeric(gapfilled_CH4_mgC_30min),
    pred_CH4_mgC_30min = as.numeric(pred_CH4_mgC_30min),
    ERA5_Tair_C = as.numeric(ERA5_Tair_C),
    ERA5_VSWC = as.numeric(ERA5_VSWC)
  ) %>%
  inner_join(upland_sites, by = c("SITE_ID", "EcoType")) %>%
  filter(is.finite(gapfilled_CH4_mgC_30min), is.finite(Year))

era5_annual <- read.csv(era5_annual_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    Year = as.integer(Year),
    annual_budget_gC_m2_yr = as.numeric(annual_budget_gC_m2_yr),
    observed_coverage = as.numeric(observed_coverage)
  ) %>%
  inner_join(upland_sites, by = "SITE_ID") %>%
  filter(is.finite(annual_budget_gC_m2_yr), is.finite(Year))

daily_era5 <- era5_30min %>%
  reframe(
    .by = c(SITE_ID, EcoType, Date, Year),
    daily_mgC_m2_day = sum(gapfilled_CH4_mgC_30min, na.rm = TRUE),
    mean_ERA5_Tair_C = mean(ERA5_Tair_C, na.rm = TRUE),
    mean_ERA5_VSWC = mean(ERA5_VSWC, na.rm = TRUE),
    n_30min = dplyr::n(),
    prop_observed_30min = mean(observed_flux, na.rm = TRUE)
  )

annual_drivers <- era5_30min %>%
  reframe(
    .by = c(SITE_ID, EcoType, Year),
    mean_ERA5_Tair_C = mean(ERA5_Tair_C, na.rm = TRUE),
    mean_ERA5_VSWC = mean(ERA5_VSWC, na.rm = TRUE),
    n_30min_grid = dplyr::n(),
    observed_coverage_from_30min = mean(observed_flux, na.rm = TRUE)
  )

site_scale_summary <- bind_rows(
  era5_30min %>%
    reframe(
      .by = c(SITE_ID, EcoType),
      temporal_resolution = "Half-hour",
      n_units = dplyr::n(),
      prop_positive = mean(gapfilled_CH4_mgC_30min > 0, na.rm = TRUE),
      annualized_flux_gC_m2_yr = annualize_30min_gC(mean(gapfilled_CH4_mgC_30min, na.rm = TRUE))
    ),
  daily_era5 %>%
    reframe(
      .by = c(SITE_ID, EcoType),
      temporal_resolution = "Daily",
      n_units = dplyr::n(),
      prop_positive = mean(daily_mgC_m2_day > 0, na.rm = TRUE),
      annualized_flux_gC_m2_yr = annualize_daily_gC(mean(daily_mgC_m2_day, na.rm = TRUE))
    ),
  era5_annual %>%
    reframe(
      .by = c(SITE_ID, EcoType),
      temporal_resolution = "Annual",
      n_units = dplyr::n(),
      prop_positive = mean(annual_budget_gC_m2_yr > 0, na.rm = TRUE),
      annualized_flux_gC_m2_yr = mean(annual_budget_gC_m2_yr, na.rm = TRUE)
    )
) %>%
  mutate(
    temporal_resolution = factor(temporal_resolution, levels = c("Half-hour", "Daily", "Annual")),
    exchange_class = factor(classify_exchange(prop_positive), levels = behavior_levels)
  )

annual_reference <- site_scale_summary %>%
  filter(temporal_resolution == "Annual") %>%
  dplyr::select(SITE_ID, annual_exchange_class = exchange_class, annual_flux_gC_m2_yr = annualized_flux_gC_m2_yr)

temporal_resolution_diagnostics <- site_scale_summary %>%
  left_join(annual_reference, by = "SITE_ID") %>%
  reframe(
    .by = temporal_resolution,
    n_sites = n_distinct(SITE_ID),
    n_weak_sink = sum(exchange_class == "Weak sink", na.rm = TRUE),
    n_fluctuating = sum(exchange_class == "Fluctuating", na.rm = TRUE),
    n_weak_source = sum(exchange_class == "Weak source", na.rm = TRUE),
    class_entropy = calc_entropy(prop.table(table(exchange_class))),
    agreement_with_annual_class = mean(exchange_class == annual_exchange_class, na.rm = TRUE),
    annual_flux_spearman = suppressWarnings(cor(annualized_flux_gC_m2_yr, annual_flux_gC_m2_yr, method = "spearman", use = "complete.obs")),
    annual_flux_rmse_gC_m2_yr = sqrt(mean((annualized_flux_gC_m2_yr - annual_flux_gC_m2_yr)^2, na.rm = TRUE))
  ) %>%
  mutate(
    # Primary target is annual net exchange, so agreement with annual class and
    # annual flux RMSE dominate the resolution score.
    resolution_score = agreement_with_annual_class +
      scales::rescale(annual_flux_spearman, to = c(0, 0.25), from = c(-1, 1)) -
      scales::rescale(annual_flux_rmse_gC_m2_yr, to = c(0, 0.25)),
    selected_for_stage1 = temporal_resolution == temporal_resolution[which.max(resolution_score)]
  ) %>%
  arrange(desc(resolution_score))

selected_resolution <- as.character(temporal_resolution_diagnostics$temporal_resolution[1])

stage1_site_data <- site_scale_summary %>%
  filter(as.character(temporal_resolution) == selected_resolution) %>%
  left_join(upland_sites, by = c("SITE_ID", "EcoType")) %>%
  mutate(exchange_class = factor(exchange_class, levels = behavior_levels)) %>%
  filter(!is.na(exchange_class))

class_model_data <- stage1_site_data %>%
  mutate(EcoType = factor(EcoType)) %>%
  filter(!is.na(EcoType))

# Stage 1 is intentionally a smoothed empirical probability model rather than a
# separable multinomial regression. Some ecosystem types have very few NEON
# sites, so Laplace smoothing keeps every class possible while preserving the
# observed class structure.
stage1_alpha <- 0.5

class_probability_by_ecotype <- class_model_data %>%
  count(EcoType, exchange_class, name = "n_class_sites") %>%
  complete(EcoType, exchange_class = factor(behavior_levels, levels = behavior_levels), fill = list(n_class_sites = 0L)) %>%
  group_by(EcoType) %>%
  mutate(
    n_sites = sum(n_class_sites),
    smoothed_probability = (n_class_sites + stage1_alpha) / (n_sites + stage1_alpha * length(behavior_levels))
  ) %>%
  ungroup() %>%
  dplyr::select(EcoType, exchange_class, n_class_sites, n_sites, smoothed_probability) %>%
  pivot_wider(
    names_from = exchange_class,
    values_from = c(n_class_sites, smoothed_probability),
    names_glue = "{.value}_{make.names(exchange_class)}"
  ) %>%
  transmute(
    EcoType,
    n_sites,
    observed_weak_sink = n_class_sites_Weak.sink / n_sites,
    observed_fluctuating = n_class_sites_Fluctuating / n_sites,
    observed_weak_source = n_class_sites_Weak.source / n_sites,
    prob_weak_sink = smoothed_probability_Weak.sink,
    prob_fluctuating = smoothed_probability_Fluctuating,
    prob_weak_source = smoothed_probability_Weak.source
  )

stage1_class_counts <- class_model_data %>%
  count(EcoType, exchange_class, name = "n_class_sites") %>%
  complete(EcoType, exchange_class = factor(behavior_levels, levels = behavior_levels), fill = list(n_class_sites = 0L))

stage1_prior_scenarios <- tribble(
  ~stage1_prior_scenario, ~prior_weak_sink, ~prior_fluctuating, ~prior_weak_source, ~prior_note,
  "Weak empirical smoothing", 0.5, 0.5, 0.5,
  "Current weak Laplace smoothing; probabilities mostly follow observed NEON annual class counts.",
  "Equal class prior", 2.0, 2.0, 2.0,
  "Balanced prior with equal pseudo-counts for weak sink, fluctuating, and weak source classes.",
  "Sink-favoring prior", 4.0, 1.0, 1.0,
  "Conservative upland prior that assumes weak sink unless weak-source evidence is strong."
) %>%
  pivot_longer(
    cols = c(prior_weak_sink, prior_fluctuating, prior_weak_source),
    names_to = "prior_class",
    values_to = "prior_count"
  ) %>%
  mutate(
    exchange_class = recode(
      prior_class,
      prior_weak_sink = "Weak sink",
      prior_fluctuating = "Fluctuating",
      prior_weak_source = "Weak source"
    ),
    exchange_class = factor(exchange_class, levels = behavior_levels)
  ) %>%
  dplyr::select(stage1_prior_scenario, prior_note, exchange_class, prior_count)

class_probability_prior_scenarios <- stage1_class_counts %>%
  crossing(stage1_prior_scenario = unique(stage1_prior_scenarios$stage1_prior_scenario)) %>%
  left_join(stage1_prior_scenarios, by = c("stage1_prior_scenario", "exchange_class")) %>%
  group_by(stage1_prior_scenario, EcoType) %>%
  mutate(
    n_sites = sum(n_class_sites),
    prior_total = sum(prior_count),
    smoothed_probability = (n_class_sites + prior_count) / (n_sites + prior_total),
    prior_note = first(prior_note)
  ) %>%
  ungroup() %>%
  dplyr::select(stage1_prior_scenario, prior_note, EcoType, exchange_class, n_class_sites, n_sites, prior_count, smoothed_probability) %>%
  pivot_wider(
    names_from = exchange_class,
    values_from = c(n_class_sites, prior_count, smoothed_probability),
    names_glue = "{.value}_{make.names(exchange_class)}"
  ) %>%
  transmute(
    stage1_prior_scenario,
    prior_note,
    EcoType,
    n_sites,
    observed_weak_sink = n_class_sites_Weak.sink / n_sites,
    observed_fluctuating = n_class_sites_Fluctuating / n_sites,
    observed_weak_source = n_class_sites_Weak.source / n_sites,
    prior_weak_sink = prior_count_Weak.sink,
    prior_fluctuating = prior_count_Fluctuating,
    prior_weak_source = prior_count_Weak.source,
    prob_weak_sink = smoothed_probability_Weak.sink,
    prob_fluctuating = smoothed_probability_Fluctuating,
    prob_weak_source = smoothed_probability_Weak.source
  )

stage1_site_probabilities <- class_model_data %>%
  left_join(class_probability_by_ecotype, by = "EcoType")

annual_flux_model_data <- era5_annual %>%
  left_join(annual_drivers, by = c("SITE_ID", "EcoType", "Year")) %>%
  left_join(stage1_site_data %>% dplyr::select(SITE_ID, exchange_class), by = "SITE_ID") %>%
  mutate(
    EcoType = factor(EcoType),
    exchange_class = factor(exchange_class, levels = behavior_levels)
  ) %>%
  filter(!is.na(exchange_class), is.finite(annual_budget_gC_m2_yr))

stage2_flux_model <- lm(
  annual_budget_gC_m2_yr ~ exchange_class + EcoType + mean_ERA5_Tair_C + mean_ERA5_VSWC,
  data = annual_flux_model_data
)

stage2_rate_by_ecotype_class <- annual_flux_model_data %>%
  reframe(
    .by = c(EcoType, exchange_class),
    n_site_years = dplyr::n(),
    n_sites = n_distinct(SITE_ID),
    observed_mean_gC_m2_yr = mean(annual_budget_gC_m2_yr, na.rm = TRUE),
    observed_median_gC_m2_yr = median(annual_budget_gC_m2_yr, na.rm = TRUE),
    observed_sd_gC_m2_yr = sd(annual_budget_gC_m2_yr, na.rm = TRUE),
    mean_ERA5_Tair_C = mean(mean_ERA5_Tair_C, na.rm = TRUE),
    mean_ERA5_VSWC = mean(mean_ERA5_VSWC, na.rm = TRUE)
  ) %>%
  complete(EcoType, exchange_class = factor(behavior_levels, levels = behavior_levels)) %>%
  group_by(EcoType) %>%
  mutate(
    mean_ERA5_Tair_C = replace_na(mean_ERA5_Tair_C, mean(mean_ERA5_Tair_C, na.rm = TRUE)),
    mean_ERA5_VSWC = replace_na(mean_ERA5_VSWC, mean(mean_ERA5_VSWC, na.rm = TRUE))
  ) %>%
  ungroup()

prediction_grid <- stage2_rate_by_ecotype_class %>%
  mutate(
    EcoType = factor(EcoType, levels = levels(annual_flux_model_data$EcoType)),
    exchange_class = factor(exchange_class, levels = behavior_levels),
    mean_ERA5_Tair_C = replace_na(mean_ERA5_Tair_C, mean(annual_flux_model_data$mean_ERA5_Tair_C, na.rm = TRUE)),
    mean_ERA5_VSWC = replace_na(mean_ERA5_VSWC, mean(annual_flux_model_data$mean_ERA5_VSWC, na.rm = TRUE))
  )

flux_predictions <- predict(stage2_flux_model, newdata = prediction_grid, se.fit = TRUE)

stage2_rate_by_ecotype_class <- prediction_grid %>%
  mutate(
    model_rate_gC_m2_yr = as.numeric(flux_predictions$fit),
    model_rate_se_gC_m2_yr = as.numeric(flux_predictions$se.fit)
  ) %>%
  dplyr::select(
    EcoType, exchange_class, n_site_years, n_sites,
    observed_mean_gC_m2_yr, observed_median_gC_m2_yr, observed_sd_gC_m2_yr,
    model_rate_gC_m2_yr, model_rate_se_gC_m2_yr
  )

chamber_process_sink_reference <- tibble(
  EcoType = c("Cropland", "Forest", "Grassland", "Shrubland"),
  external_sink_rate_gC_m2_yr = c(-0.0410625, -0.254625, -0.294, -0.2597625),
  external_sink_se_gC_m2_yr = c(0.04329719387755102, 0.10163817827832189, 0.08002260658003012, 0.0674676560064131),
  n_sink_references = c(1L, 10L, 3L, 4L),
  sink_reference_sources = c(
    "Soil chamber literature",
    "Process-based model; Soil chamber literature",
    "Process-based model; Soil chamber literature",
    "Process-based model; Soil chamber literature"
  ),
  sink_reference_classes = c(
    "Cropland soil chambers",
    "Boreal deciduous forest; Boreal evergreen forest; Mixed forest; Rural forest soil chambers; Temperate broadleaf forest; Temperate deciduous forest; Temperate needleleaf forest; Temperate/subarctic forest soil chambers; Tropical deciduous forest; Tropical evergreen forest",
    "Grassland soil chambers; Grassland/steppe; Savanna",
    "Dense shrubland; Desert soil chambers; Open shrubland; Polar desert/rock/ice"
  )
)

literature_rate_references <- chamber_process_sink_reference %>%
  right_join(
    tibble(EcoType = unique(as.character(stage2_rate_by_ecotype_class$EcoType))),
    by = "EcoType"
  ) %>%
  arrange(EcoType)

stage2_rate_scenarios <- bind_rows(
    stage2_rate_by_ecotype_class %>%
      mutate(
        rate_scenario = "NEON ERA5 rates",
        calibrated_rate_gC_m2_yr = model_rate_gC_m2_yr,
        calibrated_rate_se_gC_m2_yr = model_rate_se_gC_m2_yr,
        calibration_note = "NEON ERA5 stage-2 annual flux model."
      ),
    stage2_rate_by_ecotype_class %>%
      left_join(literature_rate_references, by = "EcoType") %>%
      mutate(
        rate_scenario = "Chamber/process sink + NEON source",
        calibrated_rate_gC_m2_yr = if_else(
          exchange_class == "Weak sink" & is.finite(external_sink_rate_gC_m2_yr),
          external_sink_rate_gC_m2_yr,
          model_rate_gC_m2_yr
        ),
        calibrated_rate_se_gC_m2_yr = if_else(
          exchange_class == "Weak sink" & is.finite(external_sink_se_gC_m2_yr),
          external_sink_se_gC_m2_yr,
          model_rate_se_gC_m2_yr
        ),
        calibration_note = if_else(
          exchange_class == "Weak sink",
          "Weak-sink rate replaced with ecosystem-matched soil chamber/process-model uptake benchmark; weak-source and fluctuating rates retained from NEON ERA5.",
          "NEON ERA5 stage-2 annual flux model retained for non-sink class."
        )
      )
) %>%
  dplyr::select(
    rate_scenario, EcoType, exchange_class, n_site_years, n_sites,
    observed_mean_gC_m2_yr, observed_median_gC_m2_yr, observed_sd_gC_m2_yr,
    model_rate_gC_m2_yr, model_rate_se_gC_m2_yr,
    calibrated_rate_gC_m2_yr, calibrated_rate_se_gC_m2_yr,
    external_sink_rate_gC_m2_yr, external_sink_se_gC_m2_yr,
    calibration_note
  )

if (!file.exists(area_assumption_file)) {
  area_assumptions <- tibble(
    EcoType = c("Cropland", "Forest", "Grassland", "Shrubland"),
    area_mha = c(1550, 4058, 5200, 1200),
    area_note = "Demonstration placeholder: replace with analysis-specific non-wetland, non-inundated area estimates before final inference."
  )
  write.csv(area_assumptions, area_assumption_file, row.names = FALSE)
} else {
  area_assumptions <- read.csv(area_assumption_file) %>%
    mutate(EcoType = as.character(EcoType), area_mha = as.numeric(area_mha))
}

flux_rates_and_global_area <- stage2_rate_scenarios %>%
  left_join(area_assumptions, by = "EcoType") %>%
  arrange(EcoType, rate_scenario, exchange_class) %>%
  transmute(
    EcoType,
    global_area_mha = area_mha,
    area_note,
    rate_scenario,
    exchange_class,
    flux_rate_gC_m2_yr = calibrated_rate_gC_m2_yr,
    flux_rate_se_gC_m2_yr = calibrated_rate_se_gC_m2_yr,
    neon_model_rate_gC_m2_yr = model_rate_gC_m2_yr,
    chamber_process_sink_rate_gC_m2_yr = external_sink_rate_gC_m2_yr,
    n_site_years,
    n_sites,
    calibration_note
  )

upscaling_components <- area_assumptions %>%
  filter(EcoType %in% unique(stage1_site_data$EcoType), is.finite(area_mha), area_mha > 0) %>%
  left_join(class_probability_by_ecotype, by = "EcoType") %>%
  pivot_longer(
    cols = c(prob_weak_sink, prob_fluctuating, prob_weak_source),
    names_to = "probability_name",
    values_to = "class_probability"
  ) %>%
  mutate(
    exchange_class = recode(
      probability_name,
      prob_weak_sink = "Weak sink",
      prob_fluctuating = "Fluctuating",
      prob_weak_source = "Weak source"
    ),
    exchange_class = factor(exchange_class, levels = behavior_levels)
  ) %>%
  left_join(stage2_rate_by_ecotype_class, by = c("EcoType", "exchange_class")) %>%
  mutate(
    area_weighted_rate_gC_m2_yr = class_probability * model_rate_gC_m2_yr,
    upscaled_tg_ch4_yr = gC_m2_yr_to_tg_ch4(area_weighted_rate_gC_m2_yr, area_mha),
    upscaled_tg_ch4_yr_se_component = abs(gC_m2_yr_to_tg_ch4(class_probability * model_rate_se_gC_m2_yr, area_mha))
  )

upscaled_global_estimate <- upscaling_components %>%
  summarise(
    selected_temporal_resolution = selected_resolution,
    n_neon_upland_sites = n_distinct(stage1_site_data$SITE_ID),
    total_area_mha = sum(unique(area_assumptions$area_mha[area_assumptions$EcoType %in% unique(stage1_site_data$EcoType)]), na.rm = TRUE),
    upscaled_net_exchange_tg_ch4_yr = sum(upscaled_tg_ch4_yr, na.rm = TRUE),
    approximate_model_se_tg_ch4_yr = sqrt(sum(upscaled_tg_ch4_yr_se_component^2, na.rm = TRUE)),
    global_methane_budget_soil_sink_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_global_soil_sink_magnitude = 100 * upscaled_net_exchange_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr),
    interpretation = case_when(
      upscaled_net_exchange_tg_ch4_yr < 0 ~ "Net upland uptake under NEON ERA5 two-stage framework",
      upscaled_net_exchange_tg_ch4_yr > 0 ~ "Net upland emission under NEON ERA5 two-stage framework",
      TRUE ~ "Near-zero net upland exchange under NEON ERA5 two-stage framework"
    )
  )

source_threshold_sensitivity <- crossing(
  source_probability_threshold = seq(0.50, 0.95, by = 0.05),
  class_probability_by_ecotype
) %>%
  mutate(
    selected_exchange_class = if_else(
      prob_weak_source >= source_probability_threshold,
      "Weak source",
      "Weak sink"
    ),
    selected_exchange_class = factor(selected_exchange_class, levels = behavior_levels)
  ) %>%
  left_join(
    area_assumptions %>%
      filter(EcoType %in% unique(stage1_site_data$EcoType), is.finite(area_mha), area_mha > 0),
    by = "EcoType"
  ) %>%
  left_join(
    stage2_rate_by_ecotype_class %>%
      dplyr::select(EcoType, exchange_class, model_rate_gC_m2_yr, model_rate_se_gC_m2_yr),
    by = c("EcoType", "selected_exchange_class" = "exchange_class")
  ) %>%
  mutate(
    upscaled_tg_ch4_yr = gC_m2_yr_to_tg_ch4(model_rate_gC_m2_yr, area_mha),
    upscaled_tg_ch4_yr_se_component = abs(gC_m2_yr_to_tg_ch4(model_rate_se_gC_m2_yr, area_mha))
  ) %>%
  group_by(source_probability_threshold) %>%
  summarise(
    n_ecotypes_as_weak_source = sum(selected_exchange_class == "Weak source", na.rm = TRUE),
    n_ecotypes_as_weak_sink = sum(selected_exchange_class == "Weak sink", na.rm = TRUE),
    upscaled_net_exchange_tg_ch4_yr = sum(upscaled_tg_ch4_yr, na.rm = TRUE),
    approximate_model_se_tg_ch4_yr = sqrt(sum(upscaled_tg_ch4_yr_se_component^2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    global_methane_budget_soil_sink_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_global_soil_sink_magnitude = 100 * upscaled_net_exchange_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr)
  )

source_threshold_default_components <- class_probability_by_ecotype %>%
  mutate(
    source_probability_threshold = source_probability_threshold_default,
    selected_exchange_class = if_else(
      prob_weak_source >= source_probability_threshold,
      "Weak source",
      "Weak sink"
    ),
    selected_exchange_class = factor(selected_exchange_class, levels = behavior_levels)
  ) %>%
  left_join(
    area_assumptions %>%
      filter(EcoType %in% unique(stage1_site_data$EcoType), is.finite(area_mha), area_mha > 0),
    by = "EcoType"
  ) %>%
  left_join(
    stage2_rate_by_ecotype_class %>%
      dplyr::select(EcoType, exchange_class, model_rate_gC_m2_yr, model_rate_se_gC_m2_yr),
    by = c("EcoType", "selected_exchange_class" = "exchange_class")
  ) %>%
  mutate(
    upscaled_tg_ch4_yr = gC_m2_yr_to_tg_ch4(model_rate_gC_m2_yr, area_mha),
    upscaled_tg_ch4_yr_se_component = abs(gC_m2_yr_to_tg_ch4(model_rate_se_gC_m2_yr, area_mha))
  )

source_threshold_default_estimate <- source_threshold_default_components %>%
  summarise(
    selected_temporal_resolution = selected_resolution,
    source_probability_threshold = source_probability_threshold_default,
    n_neon_upland_sites = n_distinct(stage1_site_data$SITE_ID),
    total_area_mha = sum(area_mha, na.rm = TRUE),
    upscaled_net_exchange_tg_ch4_yr = sum(upscaled_tg_ch4_yr, na.rm = TRUE),
    approximate_model_se_tg_ch4_yr = sqrt(sum(upscaled_tg_ch4_yr_se_component^2, na.rm = TRUE)),
    global_methane_budget_soil_sink_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_global_soil_sink_magnitude = 100 * upscaled_net_exchange_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr),
    interpretation = case_when(
      upscaled_net_exchange_tg_ch4_yr < 0 ~ "Net upland uptake under sink-default source-threshold framework",
      upscaled_net_exchange_tg_ch4_yr > 0 ~ "Net upland emission under sink-default source-threshold framework",
      TRUE ~ "Near-zero net upland exchange under sink-default source-threshold framework"
    )
  )

source_threshold_sensitivity_rate_scenarios <- crossing(
  source_probability_threshold = seq(0.50, 0.95, by = 0.05),
  class_probability_by_ecotype
) %>%
  mutate(
    selected_exchange_class = if_else(
      prob_weak_source >= source_probability_threshold,
      "Weak source",
      "Weak sink"
    ),
    selected_exchange_class = factor(selected_exchange_class, levels = behavior_levels)
  ) %>%
  left_join(
    area_assumptions %>%
      filter(EcoType %in% unique(stage1_site_data$EcoType), is.finite(area_mha), area_mha > 0),
    by = "EcoType"
  ) %>%
  left_join(
    stage2_rate_scenarios %>%
      dplyr::select(
        rate_scenario, EcoType, exchange_class,
        calibrated_rate_gC_m2_yr, calibrated_rate_se_gC_m2_yr
      ),
    by = c("EcoType", "selected_exchange_class" = "exchange_class"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    upscaled_tg_ch4_yr = gC_m2_yr_to_tg_ch4(calibrated_rate_gC_m2_yr, area_mha),
    upscaled_tg_ch4_yr_se_component = abs(gC_m2_yr_to_tg_ch4(calibrated_rate_se_gC_m2_yr, area_mha))
  ) %>%
  group_by(rate_scenario, source_probability_threshold) %>%
  summarise(
    n_ecotypes_as_weak_source = sum(selected_exchange_class == "Weak source", na.rm = TRUE),
    n_ecotypes_as_weak_sink = sum(selected_exchange_class == "Weak sink", na.rm = TRUE),
    upscaled_net_exchange_tg_ch4_yr = sum(upscaled_tg_ch4_yr, na.rm = TRUE),
    approximate_model_se_tg_ch4_yr = sqrt(sum(upscaled_tg_ch4_yr_se_component^2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    selected_temporal_resolution = selected_resolution,
    n_neon_upland_sites = n_distinct(stage1_site_data$SITE_ID),
    global_methane_budget_soil_sink_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_global_soil_sink_magnitude = 100 * upscaled_net_exchange_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr)
  )

source_threshold_default_components_rate_scenarios <- class_probability_by_ecotype %>%
  mutate(
    source_probability_threshold = source_probability_threshold_default,
    selected_exchange_class = if_else(
      prob_weak_source >= source_probability_threshold,
      "Weak source",
      "Weak sink"
    ),
    selected_exchange_class = factor(selected_exchange_class, levels = behavior_levels)
  ) %>%
  left_join(
    area_assumptions %>%
      filter(EcoType %in% unique(stage1_site_data$EcoType), is.finite(area_mha), area_mha > 0),
    by = "EcoType"
  ) %>%
  left_join(
    stage2_rate_scenarios %>%
      dplyr::select(
        rate_scenario, EcoType, exchange_class,
        calibrated_rate_gC_m2_yr, calibrated_rate_se_gC_m2_yr,
        model_rate_gC_m2_yr, external_sink_rate_gC_m2_yr,
        calibration_note
      ),
    by = c("EcoType", "selected_exchange_class" = "exchange_class"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    upscaled_tg_ch4_yr = gC_m2_yr_to_tg_ch4(calibrated_rate_gC_m2_yr, area_mha),
    upscaled_tg_ch4_yr_se_component = abs(gC_m2_yr_to_tg_ch4(calibrated_rate_se_gC_m2_yr, area_mha))
  )

source_threshold_default_estimate_rate_scenarios <- source_threshold_default_components_rate_scenarios %>%
  group_by(rate_scenario) %>%
  summarise(
    selected_temporal_resolution = selected_resolution,
    source_probability_threshold = source_probability_threshold_default,
    n_neon_upland_sites = n_distinct(stage1_site_data$SITE_ID),
    total_area_mha = sum(area_mha, na.rm = TRUE),
    upscaled_net_exchange_tg_ch4_yr = sum(upscaled_tg_ch4_yr, na.rm = TRUE),
    approximate_model_se_tg_ch4_yr = sqrt(sum(upscaled_tg_ch4_yr_se_component^2, na.rm = TRUE)),
    global_methane_budget_soil_sink_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_global_soil_sink_magnitude = 100 * upscaled_net_exchange_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr),
    .groups = "drop"
  ) %>%
  mutate(
    interpretation = case_when(
      upscaled_net_exchange_tg_ch4_yr < 0 ~ "Net upland uptake under sink-default source-threshold framework",
      upscaled_net_exchange_tg_ch4_yr > 0 ~ "Net upland emission under sink-default source-threshold framework",
      TRUE ~ "Near-zero net upland exchange under sink-default source-threshold framework"
    )
  )

stage1_prior_source_threshold_sensitivity <- crossing(
  source_probability_threshold = seq(0.50, 0.95, by = 0.05),
  class_probability_prior_scenarios
) %>%
  mutate(
    selected_exchange_class = if_else(
      prob_weak_source >= source_probability_threshold,
      "Weak source",
      "Weak sink"
    ),
    selected_exchange_class = factor(selected_exchange_class, levels = behavior_levels)
  ) %>%
  left_join(
    area_assumptions %>%
      filter(EcoType %in% unique(stage1_site_data$EcoType), is.finite(area_mha), area_mha > 0),
    by = "EcoType"
  ) %>%
  left_join(
    stage2_rate_scenarios %>%
      dplyr::select(
        rate_scenario, EcoType, exchange_class,
        calibrated_rate_gC_m2_yr, calibrated_rate_se_gC_m2_yr
      ),
    by = c("EcoType", "selected_exchange_class" = "exchange_class"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    upscaled_tg_ch4_yr = gC_m2_yr_to_tg_ch4(calibrated_rate_gC_m2_yr, area_mha),
    upscaled_tg_ch4_yr_se_component = abs(gC_m2_yr_to_tg_ch4(calibrated_rate_se_gC_m2_yr, area_mha))
  ) %>%
  group_by(stage1_prior_scenario, prior_note, rate_scenario, source_probability_threshold) %>%
  summarise(
    n_ecotypes_as_weak_source = sum(selected_exchange_class == "Weak source", na.rm = TRUE),
    n_ecotypes_as_weak_sink = sum(selected_exchange_class == "Weak sink", na.rm = TRUE),
    upscaled_net_exchange_tg_ch4_yr = sum(upscaled_tg_ch4_yr, na.rm = TRUE),
    approximate_model_se_tg_ch4_yr = sqrt(sum(upscaled_tg_ch4_yr_se_component^2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    selected_temporal_resolution = selected_resolution,
    n_neon_upland_sites = n_distinct(stage1_site_data$SITE_ID),
    global_methane_budget_soil_sink_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_global_soil_sink_magnitude = 100 * upscaled_net_exchange_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr)
  )

stage1_prior_default_estimate <- stage1_prior_source_threshold_sensitivity %>%
  filter(source_probability_threshold == source_probability_threshold_default) %>%
  arrange(stage1_prior_scenario, rate_scenario)

terrestrial_sink_comparison <- stage1_prior_default_estimate %>%
  transmute(
    comparison_type = "ERA5 upscaling estimate",
    stage1_prior_scenario,
    rate_scenario,
    source_probability_threshold,
    exchange_tg_ch4_yr = upscaled_net_exchange_tg_ch4_yr,
    approximate_se_tg_ch4_yr = approximate_model_se_tg_ch4_yr,
    terrestrial_sink_reference_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_terrestrial_sink_magnitude = percent_of_global_soil_sink_magnitude
  ) %>%
  bind_rows(
    tibble(
      comparison_type = "Global Methane Budget terrestrial soil sink",
      stage1_prior_scenario = "Reference",
      rate_scenario = "Global Methane Budget",
      source_probability_threshold = NA_real_,
      exchange_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
      approximate_se_tg_ch4_yr = NA_real_,
      terrestrial_sink_reference_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
      percent_of_terrestrial_sink_magnitude = -100
    )
  )

readr::write_csv(site_scale_summary, "OUTPUT/ERA_Upscaling_temporal_site_classifications.csv")
readr::write_csv(temporal_resolution_diagnostics, "OUTPUT/ERA_Upscaling_temporal_resolution_diagnostics.csv")
readr::write_csv(stage1_site_probabilities, "OUTPUT/ERA_Upscaling_stage1_site_class_probabilities.csv")
readr::write_csv(class_probability_by_ecotype, "OUTPUT/ERA_Upscaling_stage1_class_probability_by_ecotype.csv")
readr::write_csv(class_probability_prior_scenarios, "OUTPUT/ERA_Upscaling_stage1_class_probability_prior_scenarios.csv")
readr::write_csv(annual_flux_model_data, "OUTPUT/ERA_Upscaling_stage2_annual_flux_model_data.csv")
readr::write_csv(stage2_rate_by_ecotype_class, "OUTPUT/ERA_Upscaling_stage2_flux_rates_by_ecotype_class.csv")
readr::write_csv(literature_rate_references, "OUTPUT/ERA_Upscaling_literature_rate_references.csv")
readr::write_csv(stage2_rate_scenarios, "OUTPUT/ERA_Upscaling_stage2_flux_rates_rate_scenarios.csv")
readr::write_csv(flux_rates_and_global_area, "OUTPUT/ERA_Upscaling_flux_rates_and_global_area.csv")
readr::write_csv(upscaling_components, "OUTPUT/ERA_Upscaling_global_components.csv")
readr::write_csv(upscaled_global_estimate, "OUTPUT/ERA_Upscaling_global_net_exchange_estimate.csv")
readr::write_csv(source_threshold_sensitivity, "OUTPUT/ERA_Upscaling_source_threshold_sensitivity.csv")
readr::write_csv(source_threshold_default_components, "OUTPUT/ERA_Upscaling_source_threshold_default_components.csv")
readr::write_csv(source_threshold_default_estimate, "OUTPUT/ERA_Upscaling_source_threshold_default_estimate.csv")
readr::write_csv(source_threshold_sensitivity_rate_scenarios, "OUTPUT/ERA_Upscaling_source_threshold_sensitivity_rate_scenarios.csv")
readr::write_csv(source_threshold_default_components_rate_scenarios, "OUTPUT/ERA_Upscaling_source_threshold_default_components_rate_scenarios.csv")
readr::write_csv(source_threshold_default_estimate_rate_scenarios, "OUTPUT/ERA_Upscaling_source_threshold_default_estimate_rate_scenarios.csv")
readr::write_csv(stage1_prior_source_threshold_sensitivity, "OUTPUT/ERA_Upscaling_stage1_prior_source_threshold_sensitivity.csv")
readr::write_csv(stage1_prior_default_estimate, "OUTPUT/ERA_Upscaling_stage1_prior_default_estimate.csv")
readr::write_csv(terrestrial_sink_comparison, "OUTPUT/ERA_Upscaling_terrestrial_sink_comparison.csv")

capture.output(
  {
    cat("ERA-Upscaling two-stage framework\n\n")
    cat("Selected temporal resolution:", selected_resolution, "\n\n")
    cat("Stage 1 class probability model:\n")
    cat("Smoothed empirical probabilities by ecosystem type with Laplace alpha =", stage1_alpha, "\n")
    print(class_probability_by_ecotype)
    cat("\nStage 2 annual flux model:\n")
    print(summary(stage2_flux_model))
    cat("\nGlobal demonstration estimate:\n")
    print(upscaled_global_estimate)
    cat("\nSink-default source-threshold estimate:\n")
    print(source_threshold_default_estimate)
    cat("\nChamber/process weak-sink rate references:\n")
    print(literature_rate_references)
    cat("\nSink-default source-threshold estimates by rate scenario:\n")
    print(source_threshold_default_estimate_rate_scenarios)
    cat("\nStage 1 prior sensitivity at default source threshold:\n")
    print(stage1_prior_default_estimate)
  },
  file = "OUTPUT/ERA_Upscaling_model_summary.txt"
)

plot_temporal_resolution <- temporal_resolution_diagnostics %>%
  pivot_longer(
    cols = c(n_weak_sink, n_fluctuating, n_weak_source),
    names_to = "class_count",
    values_to = "n_class_sites"
  ) %>%
  mutate(
    exchange_class = recode(
      class_count,
      n_weak_sink = "Weak sink",
      n_fluctuating = "Fluctuating",
      n_weak_source = "Weak source"
    ),
    exchange_class = factor(exchange_class, levels = behavior_levels)
  ) %>%
  ggplot(aes(x = temporal_resolution, y = n_class_sites, fill = exchange_class)) +
  geom_col(width = 0.68) +
  geom_point(
    data = temporal_resolution_diagnostics,
    aes(x = temporal_resolution, y = n_sites + 1.3),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    color = "grey20",
    size = 3
  ) +
  geom_text(
    data = temporal_resolution_diagnostics,
    aes(x = temporal_resolution, y = n_sites + 2.3, label = paste0("score=", signif(resolution_score, 2))),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  scale_fill_manual(values = exchange_colors, drop = FALSE) +
  labs(
    title = "ERA5 Temporal Resolution Diagnostics For Upland Exchange Classification",
    subtitle = paste("Selected resolution:", selected_resolution),
    x = NULL,
    y = "NEON upland sites",
    fill = "Exchange class"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  "FIGURES/ERA_Upscaling_temporal_resolution_diagnostics.png",
  plot_temporal_resolution,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300
)

plot_upscaled_components <- upscaling_components %>%
  mutate(component_label = paste(EcoType, exchange_class, sep = " - ")) %>%
  ggplot(aes(x = reorder(component_label, upscaled_tg_ch4_yr), y = upscaled_tg_ch4_yr, fill = exchange_class)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = exchange_colors, drop = FALSE) +
  labs(
    title = "ERA5 Two-Stage Upland CH4 Upscaling Components",
    subtitle = "Area assumptions are demonstration placeholders written to OUTPUT/ERA_Upscaling_ecosystem_area_assumptions.csv",
    x = NULL,
    y = "Upscaled exchange (Tg CH4 yr-1)",
    fill = "Exchange class"
  ) +
  theme_bw(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  "FIGURES/ERA_Upscaling_global_components.png",
  plot_upscaled_components,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)

plot_budget_comparison <- bind_rows(
  upscaled_global_estimate %>%
    transmute(
      source = "NEON ERA5 two-stage upland estimate",
      exchange_tg_ch4_yr = upscaled_net_exchange_tg_ch4_yr
    ),
  tibble(
    source = "Global Methane Budget soil uptake",
    exchange_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr
  )
) %>%
  mutate(source = factor(source, levels = source)) %>%
  ggplot(aes(x = source, y = exchange_tg_ch4_yr, fill = source)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_errorbar(
    data = upscaled_global_estimate,
    aes(
      x = "NEON ERA5 two-stage upland estimate",
      ymin = upscaled_net_exchange_tg_ch4_yr - approximate_model_se_tg_ch4_yr,
      ymax = upscaled_net_exchange_tg_ch4_yr + approximate_model_se_tg_ch4_yr
    ),
    inherit.aes = FALSE,
    width = 0.12
  ) +
  annotate(
    "errorbar",
    x = 2,
    ymin = gmb_soil_sink_low_tg_ch4_yr,
    ymax = gmb_soil_sink_high_tg_ch4_yr,
    width = 0.12
  ) +
  coord_flip() +
  labs(
    title = "ERA5 Upland Upscaling Compared With Global Methane Budget Soil Sink",
    subtitle = "Global Methane Budget 2000-2020 reports soil uptake of 35 [35, 36] Tg CH4 yr-1; negative values are uptake.",
    x = NULL,
    y = "Net exchange (Tg CH4 yr-1)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  "FIGURES/ERA_Upscaling_global_budget_comparison.png",
  plot_budget_comparison,
  width = 9,
  height = 4.8,
  units = "in",
  dpi = 300
)

plot_source_threshold_sensitivity <- source_threshold_sensitivity %>%
  ggplot(aes(x = source_probability_threshold, y = upscaled_net_exchange_tg_ch4_yr)) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = gmb_soil_sink_low_tg_ch4_yr,
    ymax = gmb_soil_sink_high_tg_ch4_yr,
    fill = "#2166AC",
    alpha = 0.12
  ) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_ribbon(
    aes(
      ymin = upscaled_net_exchange_tg_ch4_yr - approximate_model_se_tg_ch4_yr,
      ymax = upscaled_net_exchange_tg_ch4_yr + approximate_model_se_tg_ch4_yr
    ),
    fill = "grey45",
    alpha = 0.18
  ) +
  geom_line(linewidth = 1.0, color = "grey15") +
  geom_point(aes(fill = factor(n_ecotypes_as_weak_source)), shape = 21, size = 3.1, color = "grey15") +
  scale_fill_brewer(palette = "Reds", name = "Ecosystems\nas source") +
  labs(
    title = "Sink-Default Source-Threshold Sensitivity",
    subtitle = "Ecosystem types are weak sinks unless smoothed weak-source probability exceeds the threshold.",
    x = "Weak-source probability threshold",
    y = "Upscaled exchange (Tg CH4 yr-1)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  "FIGURES/ERA_Upscaling_source_threshold_sensitivity.png",
  plot_source_threshold_sensitivity,
  width = 9,
  height = 5.5,
  units = "in",
  dpi = 300
)

scenario_colors <- c(
  "NEON ERA5 rates" = "#4D4D4D",
  "Chamber/process sink + NEON source" = "#1B9E77"
)

plot_literature_rate_threshold_sensitivity <- source_threshold_sensitivity_rate_scenarios %>%
  mutate(rate_scenario = factor(rate_scenario, levels = names(scenario_colors))) %>%
  ggplot(aes(x = source_probability_threshold, y = upscaled_net_exchange_tg_ch4_yr, color = rate_scenario)) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = gmb_soil_sink_low_tg_ch4_yr,
    ymax = gmb_soil_sink_high_tg_ch4_yr,
    fill = "#2166AC",
    alpha = 0.12
  ) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_line(linewidth = 1.0) +
  geom_point(aes(shape = factor(n_ecotypes_as_weak_source)), size = 2.8) +
  scale_color_manual(values = scenario_colors, name = "Rate scenario") +
  scale_shape_manual(values = c("0" = 16, "1" = 17, "2" = 15, "3" = 18, "4" = 8), name = "Ecosystems\nas source") +
  labs(
    title = "ERA5 Upland Upscaling With Chamber/Process Sink Calibration",
    subtitle = "Weak-sink rates use ecosystem-matched chamber/process benchmarks; weak-source rates remain NEON ERA5 based.",
    x = "Weak-source probability threshold",
    y = "Upscaled exchange (Tg CH4 yr-1)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  "FIGURES/ERA_Upscaling_literature_calibrated_source_threshold_sensitivity.png",
  plot_literature_rate_threshold_sensitivity,
  width = 10,
  height = 5.8,
  units = "in",
  dpi = 300
)

ggsave(
  "FIGURES/ERA_Upscaling_chamber_process_source_threshold_sensitivity.png",
  plot_literature_rate_threshold_sensitivity,
  width = 10,
  height = 5.8,
  units = "in",
  dpi = 300
)

plot_default_rate_scenario_data <- bind_rows(
  source_threshold_default_estimate_rate_scenarios %>%
    transmute(
      source = rate_scenario,
      exchange_tg_ch4_yr = upscaled_net_exchange_tg_ch4_yr,
      approximate_se_tg_ch4_yr = approximate_model_se_tg_ch4_yr
    ),
  tibble(
    source = "Global Methane Budget soil uptake",
    exchange_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    approximate_se_tg_ch4_yr = NA_real_
  )
) %>%
  mutate(source = factor(source, levels = source))

plot_default_rate_scenario_comparison <- plot_default_rate_scenario_data %>%
  ggplot(aes(x = source, y = exchange_tg_ch4_yr, fill = source)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_errorbar(
    aes(
      ymin = exchange_tg_ch4_yr - approximate_se_tg_ch4_yr,
      ymax = exchange_tg_ch4_yr + approximate_se_tg_ch4_yr
    ),
    width = 0.12,
    na.rm = TRUE
  ) +
  annotate(
    "errorbar",
    x = which(levels(plot_default_rate_scenario_data$source) == "Global Methane Budget soil uptake"),
    ymin = gmb_soil_sink_low_tg_ch4_yr,
    ymax = gmb_soil_sink_high_tg_ch4_yr,
    width = 0.12
  ) +
  coord_flip() +
  labs(
    title = "Default Threshold Upscaling With Chamber/Process Sink Calibration",
    subtitle = paste0("Weak-source probability threshold = ", source_probability_threshold_default),
    x = NULL,
    y = "Net exchange (Tg CH4 yr-1)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  "FIGURES/ERA_Upscaling_literature_calibrated_budget_comparison.png",
  plot_default_rate_scenario_comparison,
  width = 9.5,
  height = 5.2,
  units = "in",
  dpi = 300
)

ggsave(
  "FIGURES/ERA_Upscaling_chamber_process_budget_comparison.png",
  plot_default_rate_scenario_comparison,
  width = 9.5,
  height = 5.2,
  units = "in",
  dpi = 300
)

prior_colors <- c(
  "Weak empirical smoothing" = "#4D4D4D",
  "Equal class prior" = "#7570B3",
  "Sink-favoring prior" = "#1B9E77"
)

plot_stage1_prior_probabilities <- class_probability_prior_scenarios %>%
  dplyr::select(stage1_prior_scenario, EcoType, prob_weak_sink, prob_fluctuating, prob_weak_source) %>%
  pivot_longer(
    cols = starts_with("prob_"),
    names_to = "probability_name",
    values_to = "class_probability"
  ) %>%
  mutate(
    exchange_class = recode(
      probability_name,
      prob_weak_sink = "Weak sink",
      prob_fluctuating = "Fluctuating",
      prob_weak_source = "Weak source"
    ),
    exchange_class = factor(exchange_class, levels = behavior_levels),
    stage1_prior_scenario = factor(stage1_prior_scenario, levels = names(prior_colors))
  ) %>%
  ggplot(aes(x = EcoType, y = class_probability, fill = exchange_class)) +
  geom_col(width = 0.72) +
  facet_wrap(~stage1_prior_scenario, ncol = 1) +
  scale_fill_manual(values = exchange_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Stage 1 Class Probabilities Under Alternative Priors",
    subtitle = "Equal and sink-favoring priors test whether source probabilities are driven by NEON class imbalance.",
    x = NULL,
    y = "Class probability",
    fill = "Exchange class"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  "FIGURES/ERA_Upscaling_stage1_prior_class_probabilities.png",
  plot_stage1_prior_probabilities,
  width = 9,
  height = 8,
  units = "in",
  dpi = 300
)

plot_stage1_prior_threshold_sensitivity <- stage1_prior_source_threshold_sensitivity %>%
  mutate(
    stage1_prior_scenario = factor(stage1_prior_scenario, levels = names(prior_colors)),
    rate_scenario = factor(rate_scenario, levels = names(scenario_colors))
  ) %>%
  ggplot(aes(x = source_probability_threshold, y = upscaled_net_exchange_tg_ch4_yr, color = stage1_prior_scenario)) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = gmb_soil_sink_low_tg_ch4_yr,
    ymax = gmb_soil_sink_high_tg_ch4_yr,
    fill = "#2166AC",
    alpha = 0.12
  ) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_line(linewidth = 1.0) +
  geom_point(aes(shape = factor(n_ecotypes_as_weak_source)), size = 2.6) +
  facet_wrap(~rate_scenario, ncol = 1) +
  scale_color_manual(values = prior_colors, name = "Stage 1 prior") +
  scale_shape_manual(values = c("0" = 16, "1" = 17, "2" = 15, "3" = 18, "4" = 8), name = "Ecosystems\nas source") +
  labs(
    title = "Stage 1 Prior Sensitivity For Sink-Default Upscaling",
    subtitle = "Prior choices change whether ecosystem weak-source probabilities clear the source threshold.",
    x = "Weak-source probability threshold",
    y = "Upscaled exchange (Tg CH4 yr-1)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  "FIGURES/ERA_Upscaling_stage1_prior_source_threshold_sensitivity.png",
  plot_stage1_prior_threshold_sensitivity,
  width = 10,
  height = 7.5,
  units = "in",
  dpi = 300
)

plot_terrestrial_sink_comparison <- stage1_prior_default_estimate %>%
  mutate(
    stage1_prior_scenario = factor(stage1_prior_scenario, levels = names(prior_colors)),
    rate_scenario = factor(rate_scenario, levels = names(scenario_colors))
  ) %>%
  ggplot(aes(x = stage1_prior_scenario, y = upscaled_net_exchange_tg_ch4_yr, fill = stage1_prior_scenario)) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = gmb_soil_sink_low_tg_ch4_yr,
    ymax = gmb_soil_sink_high_tg_ch4_yr,
    fill = "#2166AC",
    alpha = 0.12
  ) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_col(width = 0.68, color = "grey20") +
  geom_errorbar(
    aes(
      ymin = upscaled_net_exchange_tg_ch4_yr - approximate_model_se_tg_ch4_yr,
      ymax = upscaled_net_exchange_tg_ch4_yr + approximate_model_se_tg_ch4_yr
    ),
    width = 0.16
  ) +
  facet_wrap(~rate_scenario, ncol = 1) +
  scale_fill_manual(values = prior_colors, guide = "none") +
  labs(
    title = "ERA5 Upland Estimates Compared With Terrestrial Soil CH4 Sink",
    subtitle = paste0(
      "Default weak-source threshold = ", source_probability_threshold_default,
      "; blue band/line show Global Methane Budget soil uptake of 35 [35, 36] Tg CH4 yr-1."
    ),
    x = NULL,
    y = "Net exchange (Tg CH4 yr-1)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

ggsave(
  "FIGURES/ERA_Upscaling_terrestrial_sink_comparison.png",
  plot_terrestrial_sink_comparison,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)

message("Wrote ERA-Upscaling outputs using ", n_distinct(stage1_site_data$SITE_ID), " non-wetland NEON upland sites.")
