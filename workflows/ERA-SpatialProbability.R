# ERA5 condition-based upland CH4 source-probability model.
#
# Purpose:
#   1. estimate P(weak source) for non-wetland upland sites from ecosystem type
#      and environmental conditions rather than ecosystem class counts alone;
#   2. validate that model with leave-one-site-out prediction before any global
#      spatial extrapolation is attempted; and
#   3. define the gridded data contract needed for spatial estimates from
#      2000-2025.
#
# This script does not download global gridded data. It writes a manifest of the
# required spatial inputs and stops before global prediction unless those inputs
# already exist locally.

library(tidyverse)
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
rate_file <- "OUTPUT/ERA_Upscaling_stage2_flux_rates_rate_scenarios.csv"
area_file <- "OUTPUT/ERA_Upscaling_ecosystem_area_assumptions.csv"

required_files <- c(era5_30min_file, era5_annual_file, site_behavior_file, rate_file, area_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

source_threshold <- 1
sink_threshold <- 0
source_probability_threshold_default <- 0.80
gmb_soil_sink_tg_ch4_yr <- -35

behavior_levels <- c("Weak sink", "Fluctuating", "Weak source")

classify_exchange <- function(prop_positive) {
  case_when(
    is.na(prop_positive) ~ NA_character_,
    prop_positive >= source_threshold ~ "Weak source",
    prop_positive <= sink_threshold ~ "Weak sink",
    TRUE ~ "Fluctuating"
  )
}

gC_m2_yr_to_tg_ch4 <- function(flux_gC_m2_yr, area_mha) {
  flux_gC_m2_yr * area_mha * 0.0133333333333333
}

logit_clip <- function(p, eps = 1e-6) {
  pmin(pmax(p, eps), 1 - eps)
}

calc_auc <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  if (length(unique(obs)) < 2) return(NA_real_)
  ranks <- rank(pred, ties.method = "average")
  n_pos <- sum(obs == 1)
  n_neg <- sum(obs == 0)
  (sum(ranks[obs == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

safe_predict_glm <- function(model, newdata) {
  pred <- tryCatch(
    predict(model, newdata = newdata, type = "response"),
    error = function(e) NA_real_
  )
  as.numeric(pred)
}

site_behavior <- read.csv(site_behavior_file) %>%
  mutate(SITE_ID = as.character(SITE_ID))

upland_sites <- site_behavior %>%
  filter(
    !is.na(EcoType),
    !str_detect(EcoType, regex("wetland|inundat|flood|marsh|swamp|bog|fen|lake|rice", ignore_case = TRUE))
  ) %>%
  distinct(
    SITE_ID, EcoType, MAP, MAT, acidity, carbonTot, nitrogenTot,
    sulfurTot, clayTotal, sandTotal, bulkDensOvenDry, LAI.mean, canopyHeight_m
  )

era5_30min <- data.table::fread(era5_30min_file) %>%
  as_tibble() %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    Year = as.integer(Year),
    ERA5_Tair_C = as.numeric(ERA5_Tair_C),
    ERA5_VSWC = as.numeric(ERA5_VSWC),
    gapfilled_CH4_mgC_30min = as.numeric(gapfilled_CH4_mgC_30min)
  ) %>%
  inner_join(upland_sites %>% dplyr::select(SITE_ID, EcoType), by = c("SITE_ID", "EcoType")) %>%
  filter(is.finite(Year), is.finite(ERA5_Tair_C), is.finite(ERA5_VSWC))

era5_annual <- read.csv(era5_annual_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    Year = as.integer(Year),
    annual_budget_gC_m2_yr = as.numeric(annual_budget_gC_m2_yr)
  ) %>%
  inner_join(upland_sites, by = "SITE_ID") %>%
  filter(is.finite(Year), is.finite(annual_budget_gC_m2_yr))

site_annual_class <- era5_annual %>%
  reframe(
    .by = c(SITE_ID, EcoType),
    n_years = dplyr::n(),
    prop_positive_annual = mean(annual_budget_gC_m2_yr > 0, na.rm = TRUE),
    mean_annual_budget_gC_m2_yr = mean(annual_budget_gC_m2_yr, na.rm = TRUE)
  ) %>%
  mutate(
    exchange_class = factor(classify_exchange(prop_positive_annual), levels = behavior_levels),
    weak_source = as.integer(exchange_class == "Weak source")
  )

site_era5_conditions <- era5_30min %>%
  reframe(
    .by = c(SITE_ID, EcoType),
    mean_ERA5_Tair_C = mean(ERA5_Tair_C, na.rm = TRUE),
    sd_ERA5_Tair_C = sd(ERA5_Tair_C, na.rm = TRUE),
    mean_ERA5_VSWC = mean(ERA5_VSWC, na.rm = TRUE),
    sd_ERA5_VSWC = sd(ERA5_VSWC, na.rm = TRUE)
  )

condition_model_data <- site_annual_class %>%
  left_join(site_era5_conditions, by = c("SITE_ID", "EcoType")) %>%
  left_join(upland_sites, by = c("SITE_ID", "EcoType")) %>%
  mutate(
    EcoType = factor(EcoType),
    MAP = as.numeric(MAP),
    MAT = as.numeric(MAT),
    acidity = as.numeric(acidity),
    carbonTot = as.numeric(carbonTot),
    clayTotal = as.numeric(clayTotal),
    sandTotal = as.numeric(sandTotal)
  ) %>%
  filter(
    !is.na(exchange_class),
    is.finite(weak_source),
    is.finite(mean_ERA5_Tair_C),
    is.finite(mean_ERA5_VSWC),
    is.finite(MAP),
    is.finite(MAT)
  )

if (nrow(condition_model_data) < 10 || length(unique(condition_model_data$weak_source)) < 2) {
  stop("Not enough non-wetland upland sites or class variation to fit the condition model.")
}

condition_model_formula <- weak_source ~ EcoType + scale(mean_ERA5_Tair_C) +
  scale(mean_ERA5_VSWC) + scale(MAP) + scale(MAT)

add_class_balance_weights <- function(dat) {
  n_total <- nrow(dat)
  n_class <- dat %>%
    count(weak_source, name = "n_class")

  dat %>%
    left_join(n_class, by = "weak_source") %>%
    mutate(class_balance_weight = n_total / (length(unique(weak_source)) * n_class)) %>%
    dplyr::select(-n_class)
}

fit_condition_model <- function(dat, balanced = FALSE) {
  dat <- if (balanced) add_class_balance_weights(dat) else mutate(dat, class_balance_weight = 1)
  glm(
    condition_model_formula,
    data = dat,
    family = if (balanced) quasibinomial() else binomial(),
    weights = class_balance_weight
  )
}

fit_glm_model <- function(dat, formula, balanced = FALSE) {
  dat <- if (balanced) add_class_balance_weights(dat) else mutate(dat, class_balance_weight = 1)
  glm(
    formula,
    data = dat,
    family = if (balanced) quasibinomial() else binomial(),
    weights = class_balance_weight
  )
}

condition_model_scenarios <- tibble(
  model_scenario = "Class-balanced condition model",
  balanced = TRUE
)

condition_models <- condition_model_scenarios %>%
  mutate(
    model = map(balanced, ~ fit_condition_model(condition_model_data, balanced = .x))
  )

condition_model_predictions <- condition_models %>%
  transmute(
    model_scenario,
    pred_source_probability = map(model, ~ as.numeric(predict(.x, newdata = condition_model_data, type = "response")))
  ) %>%
  unnest(pred_source_probability) %>%
  group_by(model_scenario) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

condition_model_data_with_id <- condition_model_data %>%
  mutate(row_id = row_number())

loo_predictions <- map_dfr(seq_len(nrow(condition_model_data)), function(i) {
  train_base <- condition_model_data[-i, , drop = FALSE]
  test <- condition_model_data[i, , drop = FALSE]

  condition_model_scenarios %>%
    mutate(
      loo_model = map(balanced, ~ tryCatch(
        fit_condition_model(train_base, balanced = .x),
        error = function(e) NULL
      )),
      loo_pred_source_probability = map_dbl(loo_model, ~ if (is.null(.x)) NA_real_ else safe_predict_glm(.x, test))
    ) %>%
    transmute(
      model_scenario,
      SITE_ID = test$SITE_ID,
      EcoType = as.character(test$EcoType),
      observed_exchange_class = as.character(test$exchange_class),
      observed_weak_source = test$weak_source,
      loo_pred_source_probability
    )
})

validation_metrics <- loo_predictions %>%
  group_by(model_scenario) %>%
  summarise(
    n_sites = dplyr::n(),
    n_weak_source = sum(observed_weak_source == 1, na.rm = TRUE),
    n_not_weak_source = sum(observed_weak_source == 0, na.rm = TRUE),
    auc = calc_auc(observed_weak_source, loo_pred_source_probability),
    brier_score = mean((observed_weak_source - loo_pred_source_probability)^2, na.rm = TRUE),
    log_loss = -mean(
      observed_weak_source * log(logit_clip(loo_pred_source_probability)) +
        (1 - observed_weak_source) * log(1 - logit_clip(loo_pred_source_probability)),
      na.rm = TRUE
    ),
    accuracy_threshold_0_5 = mean((loo_pred_source_probability >= 0.5) == observed_weak_source, na.rm = TRUE),
    sensitivity_threshold_0_5 = sum(loo_pred_source_probability >= 0.5 & observed_weak_source == 1, na.rm = TRUE) /
      sum(observed_weak_source == 1, na.rm = TRUE),
    specificity_threshold_0_5 = sum(loo_pred_source_probability < 0.5 & observed_weak_source == 0, na.rm = TRUE) /
      sum(observed_weak_source == 0, na.rm = TRUE),
    n_predicted_source_at_0_8 = sum(loo_pred_source_probability >= source_probability_threshold_default, na.rm = TRUE),
    .groups = "drop"
  )

condition_site_predictions <- condition_model_data_with_id %>%
  dplyr::select(
    row_id, SITE_ID, EcoType, exchange_class, weak_source, prop_positive_annual,
    mean_ERA5_Tair_C, mean_ERA5_VSWC,
    MAP, MAT, mean_annual_budget_gC_m2_yr
  ) %>%
  left_join(condition_model_predictions, by = "row_id") %>%
  left_join(loo_predictions, by = c("model_scenario", "SITE_ID", "EcoType"))

monthly_condition_model_data <- era5_30min %>%
  reframe(
    .by = c(SITE_ID, EcoType, Year, month, season),
    monthly_budget_mgC_m2 = sum(gapfilled_CH4_mgC_30min, na.rm = TRUE),
    mean_ERA5_Tair_C = mean(ERA5_Tair_C, na.rm = TRUE),
    mean_ERA5_VSWC = mean(ERA5_VSWC, na.rm = TRUE),
    n_30min = dplyr::n()
  ) %>%
  left_join(upland_sites, by = c("SITE_ID", "EcoType")) %>%
  mutate(
    EcoType = factor(EcoType),
    month = factor(month),
    season = factor(season),
    weak_source = as.integer(monthly_budget_mgC_m2 > 0),
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

monthly_condition_model_formula <- weak_source ~ EcoType + season +
  scale(mean_ERA5_Tair_C) + scale(mean_ERA5_VSWC) + scale(MAP) + scale(MAT)

monthly_condition_models <- condition_model_scenarios %>%
  mutate(
    model_scenario = paste("Monthly", model_scenario),
    model = map(balanced, ~ fit_glm_model(monthly_condition_model_data, monthly_condition_model_formula, balanced = .x))
  )

monthly_condition_model_predictions <- monthly_condition_models %>%
  transmute(
    model_scenario,
    pred_source_probability = map(model, ~ as.numeric(predict(.x, newdata = monthly_condition_model_data, type = "response")))
  ) %>%
  unnest(pred_source_probability) %>%
  group_by(model_scenario) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

monthly_condition_model_data_with_id <- monthly_condition_model_data %>%
  mutate(row_id = row_number())

monthly_loo_predictions <- map_dfr(unique(monthly_condition_model_data$SITE_ID), function(heldout_site) {
  train_base <- monthly_condition_model_data %>% filter(SITE_ID != heldout_site)
  test <- monthly_condition_model_data %>% filter(SITE_ID == heldout_site)

  map2_dfr(condition_model_scenarios$model_scenario, condition_model_scenarios$balanced, function(model_name, balanced_model) {
    model_i <- tryCatch(
      fit_glm_model(train_base, monthly_condition_model_formula, balanced = balanced_model),
      error = function(e) NULL
    )

    test %>%
      transmute(
        model_scenario = paste("Monthly", model_name),
        SITE_ID,
        EcoType = as.character(EcoType),
        Year,
        month = as.integer(as.character(month)),
        season = as.character(season),
        observed_weak_source = weak_source,
        monthly_budget_mgC_m2,
        loo_pred_source_probability = if (is.null(model_i)) NA_real_ else safe_predict_glm(model_i, test)
      )
  })
})

monthly_validation_metrics <- monthly_loo_predictions %>%
  group_by(model_scenario) %>%
  summarise(
    n_site_months = dplyr::n(),
    n_sites = n_distinct(SITE_ID),
    n_source_months = sum(observed_weak_source == 1, na.rm = TRUE),
    n_sink_months = sum(observed_weak_source == 0, na.rm = TRUE),
    auc = calc_auc(observed_weak_source, loo_pred_source_probability),
    brier_score = mean((observed_weak_source - loo_pred_source_probability)^2, na.rm = TRUE),
    log_loss = -mean(
      observed_weak_source * log(logit_clip(loo_pred_source_probability)) +
        (1 - observed_weak_source) * log(1 - logit_clip(loo_pred_source_probability)),
      na.rm = TRUE
    ),
    accuracy_threshold_0_5 = mean((loo_pred_source_probability >= 0.5) == observed_weak_source, na.rm = TRUE),
    sensitivity_threshold_0_5 = sum(loo_pred_source_probability >= 0.5 & observed_weak_source == 1, na.rm = TRUE) /
      sum(observed_weak_source == 1, na.rm = TRUE),
    specificity_threshold_0_5 = sum(loo_pred_source_probability < 0.5 & observed_weak_source == 0, na.rm = TRUE) /
      sum(observed_weak_source == 0, na.rm = TRUE),
    n_predicted_source_at_0_8 = sum(loo_pred_source_probability >= source_probability_threshold_default, na.rm = TRUE),
    .groups = "drop"
  )

monthly_condition_predictions <- monthly_condition_model_data_with_id %>%
  mutate(month = as.integer(as.character(month))) %>%
  left_join(monthly_condition_model_predictions, by = "row_id") %>%
  left_join(
    monthly_loo_predictions %>%
      dplyr::select(model_scenario, SITE_ID, Year, month, loo_pred_source_probability),
    by = c("model_scenario", "SITE_ID", "Year", "month")
  )

rate_scenarios <- read.csv(rate_file) %>%
  mutate(
    EcoType = as.character(EcoType),
    exchange_class = as.character(exchange_class),
    calibrated_rate_gC_m2_yr = as.numeric(calibrated_rate_gC_m2_yr),
    calibrated_rate_se_gC_m2_yr = as.numeric(calibrated_rate_se_gC_m2_yr)
  )

area_assumptions <- read.csv(area_file) %>%
  mutate(EcoType = as.character(EcoType), area_mha = as.numeric(area_mha))

ecotype_mean_conditions <- condition_site_predictions %>%
  reframe(
    .by = c(model_scenario, EcoType),
    condition_pred_source_probability = mean(pred_source_probability, na.rm = TRUE),
    n_sites = dplyr::n(),
    mean_ERA5_Tair_C = mean(mean_ERA5_Tair_C, na.rm = TRUE),
    mean_ERA5_VSWC = mean(mean_ERA5_VSWC, na.rm = TRUE),
    MAP = mean(MAP, na.rm = TRUE),
    MAT = mean(MAT, na.rm = TRUE)
  )

condition_area_weighted_demo <- ecotype_mean_conditions %>%
  mutate(
    selected_exchange_class = if_else(
      condition_pred_source_probability >= source_probability_threshold_default,
      "Weak source",
      "Weak sink"
    )
  ) %>%
  left_join(area_assumptions, by = "EcoType") %>%
  left_join(
    rate_scenarios %>%
      dplyr::select(rate_scenario, EcoType, exchange_class, calibrated_rate_gC_m2_yr, calibrated_rate_se_gC_m2_yr),
    by = c("EcoType", "selected_exchange_class" = "exchange_class"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    upscaled_tg_ch4_yr = gC_m2_yr_to_tg_ch4(calibrated_rate_gC_m2_yr, area_mha),
    upscaled_tg_ch4_yr_se_component = abs(gC_m2_yr_to_tg_ch4(calibrated_rate_se_gC_m2_yr, area_mha))
  )

condition_area_weighted_estimate <- condition_area_weighted_demo %>%
  group_by(model_scenario, rate_scenario) %>%
  summarise(
    source_probability_threshold = source_probability_threshold_default,
    n_ecotypes_as_weak_source = sum(selected_exchange_class == "Weak source", na.rm = TRUE),
    n_ecotypes_as_weak_sink = sum(selected_exchange_class == "Weak sink", na.rm = TRUE),
    upscaled_net_exchange_tg_ch4_yr = sum(upscaled_tg_ch4_yr, na.rm = TRUE),
    approximate_model_se_tg_ch4_yr = sqrt(sum(upscaled_tg_ch4_yr_se_component^2, na.rm = TRUE)),
    global_methane_budget_soil_sink_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_global_soil_sink_magnitude = 100 * upscaled_net_exchange_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr),
    .groups = "drop"
  )

monthly_ecotype_mean_conditions <- monthly_condition_predictions %>%
  reframe(
    .by = c(model_scenario, EcoType),
    condition_pred_source_probability = mean(pred_source_probability, na.rm = TRUE),
    n_site_months = dplyr::n(),
    n_sites = n_distinct(SITE_ID),
    mean_ERA5_Tair_C = mean(mean_ERA5_Tair_C, na.rm = TRUE),
    mean_ERA5_VSWC = mean(mean_ERA5_VSWC, na.rm = TRUE),
    MAP = mean(MAP, na.rm = TRUE),
    MAT = mean(MAT, na.rm = TRUE)
  )

monthly_area_weighted_demo <- monthly_ecotype_mean_conditions %>%
  mutate(
    selected_exchange_class = if_else(
      condition_pred_source_probability >= source_probability_threshold_default,
      "Weak source",
      "Weak sink"
    )
  ) %>%
  left_join(area_assumptions, by = "EcoType") %>%
  left_join(
    rate_scenarios %>%
      dplyr::select(rate_scenario, EcoType, exchange_class, calibrated_rate_gC_m2_yr, calibrated_rate_se_gC_m2_yr),
    by = c("EcoType", "selected_exchange_class" = "exchange_class"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    upscaled_tg_ch4_yr = gC_m2_yr_to_tg_ch4(calibrated_rate_gC_m2_yr, area_mha),
    upscaled_tg_ch4_yr_se_component = abs(gC_m2_yr_to_tg_ch4(calibrated_rate_se_gC_m2_yr, area_mha))
  )

monthly_area_weighted_estimate <- monthly_area_weighted_demo %>%
  group_by(model_scenario, rate_scenario) %>%
  summarise(
    source_probability_threshold = source_probability_threshold_default,
    n_ecotypes_as_weak_source = sum(selected_exchange_class == "Weak source", na.rm = TRUE),
    n_ecotypes_as_weak_sink = sum(selected_exchange_class == "Weak sink", na.rm = TRUE),
    upscaled_net_exchange_tg_ch4_yr = sum(upscaled_tg_ch4_yr, na.rm = TRUE),
    approximate_model_se_tg_ch4_yr = sqrt(sum(upscaled_tg_ch4_yr_se_component^2, na.rm = TRUE)),
    global_methane_budget_soil_sink_tg_ch4_yr = gmb_soil_sink_tg_ch4_yr,
    percent_of_global_soil_sink_magnitude = 100 * upscaled_net_exchange_tg_ch4_yr / abs(gmb_soil_sink_tg_ch4_yr),
    .groups = "drop"
  )

spatial_input_manifest <- tribble(
  ~input_name, ~required_fields, ~temporal_coverage, ~role, ~local_expected_path,
  "Non-wetland upland land-cover grid",
  "cell_id, lon, lat, year, EcoType, cell_area_m2, non_wetland_upland_mask",
  "annual 2000-2025 or static land cover with annual masks",
  "Defines which grid cells are eligible and maps each cell to Cropland, Forest, Grassland, or Shrubland.",
  "SPATIAL/landcover_nonwetland_upland_2000_2025.parquet",
  "Wetland/inundation exclusion mask",
  "cell_id, year, wetland_or_inundated",
  "annual 2000-2025",
  "Excludes wetlands, lakes, rice, inundated/flooded areas, bogs, fens, marshes, and swamps before prediction.",
  "SPATIAL/wetland_inundation_mask_2000_2025.parquet",
  "ERA5 annual condition grid",
  "cell_id, year, mean_ERA5_Tair_C, mean_ERA5_VSWC",
  "annual 2000-2025",
  "Provides the condition covariates used by the source-probability model.",
  "SPATIAL/era5_annual_conditions_2000_2025.parquet",
  "Climate normals or annual climate covariates",
  "cell_id, year, MAP, MAT",
  "annual or climatological 2000-2025",
  "Provides MAP and MAT terms used in the source-probability model.",
  "SPATIAL/climate_covariates_2000_2025.parquet"
)

global_spatial_ready <- all(file.exists(spatial_input_manifest$local_expected_path))

if (global_spatial_ready) {
  stop("Spatial input files are present, but gridded prediction has not been implemented in this first validation script.")
}

readr::write_csv(condition_model_data, "OUTPUT/ERA_SpatialProbability_condition_model_training_data.csv")
readr::write_csv(condition_site_predictions, "OUTPUT/ERA_SpatialProbability_site_predictions.csv")
readr::write_csv(validation_metrics, "OUTPUT/ERA_SpatialProbability_validation_metrics.csv")
readr::write_csv(condition_area_weighted_demo, "OUTPUT/ERA_SpatialProbability_area_weighted_demo_components.csv")
readr::write_csv(condition_area_weighted_estimate, "OUTPUT/ERA_SpatialProbability_area_weighted_demo_estimate.csv")
readr::write_csv(monthly_condition_model_data, "OUTPUT/ERA_SpatialProbability_monthly_condition_model_training_data.csv")
readr::write_csv(monthly_condition_predictions, "OUTPUT/ERA_SpatialProbability_monthly_predictions.csv")
readr::write_csv(monthly_validation_metrics, "OUTPUT/ERA_SpatialProbability_monthly_validation_metrics.csv")
readr::write_csv(monthly_area_weighted_demo, "OUTPUT/ERA_SpatialProbability_monthly_area_weighted_demo_components.csv")
readr::write_csv(monthly_area_weighted_estimate, "OUTPUT/ERA_SpatialProbability_monthly_area_weighted_demo_estimate.csv")
readr::write_csv(spatial_input_manifest, "OUTPUT/ERA_SpatialProbability_global_spatial_input_manifest.csv")

capture.output(
  {
    cat("ERA spatial probability model\n\n")
    cat("Condition model formula:\n")
    print(condition_model_formula)
    cat("\nModel summaries:\n")
    walk2(condition_models$model_scenario, condition_models$model, function(model_name, model_object) {
      cat("\n", model_name, "\n", sep = "")
      print(summary(model_object))
    })
    cat("\nLeave-one-site-out validation metrics:\n")
    print(validation_metrics)
    cat("\nMonthly leave-one-site-out validation metrics:\n")
    print(monthly_validation_metrics)
    cat("\nArea-weighted demonstration estimate using ecosystem-mean conditions:\n")
    print(condition_area_weighted_estimate)
    cat("\nMonthly area-weighted demonstration estimate using ecosystem-mean monthly conditions:\n")
    print(monthly_area_weighted_estimate)
    cat("\nSpatial global prediction status:\n")
    cat("Required gridded spatial inputs present:", global_spatial_ready, "\n")
    print(spatial_input_manifest)
  },
  file = "OUTPUT/ERA_SpatialProbability_model_summary.txt"
)

plot_validation <- condition_site_predictions %>%
  mutate(exchange_class = factor(exchange_class, levels = behavior_levels)) %>%
  ggplot(aes(x = exchange_class, y = loo_pred_source_probability, fill = exchange_class)) +
  geom_hline(yintercept = source_probability_threshold_default, linetype = "dashed", color = "grey25") +
  geom_boxplot(width = 0.55, alpha = 0.65, outlier.shape = NA) +
  geom_jitter(aes(shape = EcoType), width = 0.12, height = 0, size = 2.3, color = "grey15") +
  facet_wrap(~model_scenario, ncol = 1) +
  scale_fill_manual(values = c("Weak sink" = "#2166AC", "Fluctuating" = "#4D4D4D", "Weak source" = "#B2182B"), drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Condition-Based Source Probability Model: Leave-One-Site-Out Predictions",
    subtitle = paste0("Dashed line = ", source_probability_threshold_default, " weak-source threshold"),
    x = "Observed annual ERA5 behavior class",
    y = "Predicted P(weak source)",
    fill = "Observed class",
    shape = "Ecosystem"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  "FIGURES/ERA_SpatialProbability_validation.png",
  plot_validation,
  width = 9,
  height = 8,
  units = "in",
  dpi = 300
)

plot_demo_estimate <- condition_area_weighted_estimate %>%
  ggplot(aes(x = rate_scenario, y = upscaled_net_exchange_tg_ch4_yr, fill = rate_scenario)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -35, fill = "#2166AC", alpha = 0.12) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_errorbar(
    aes(
      ymin = upscaled_net_exchange_tg_ch4_yr - approximate_model_se_tg_ch4_yr,
      ymax = upscaled_net_exchange_tg_ch4_yr + approximate_model_se_tg_ch4_yr
    ),
    width = 0.12
  ) +
  facet_wrap(~model_scenario, ncol = 1) +
  coord_flip() +
  labs(
    title = "Condition-Based Probability Model: Area-Weighted Demonstration",
    subtitle = "Uses ecosystem-mean NEON conditions and placeholder global areas; not a spatial global estimate.",
    x = NULL,
    y = "Net exchange (Tg CH4 yr-1)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  "FIGURES/ERA_SpatialProbability_area_weighted_demo.png",
  plot_demo_estimate,
  width = 9,
  height = 4.8,
  units = "in",
  dpi = 300
)

plot_monthly_validation <- monthly_loo_predictions %>%
  mutate(
    observed_month_class = if_else(observed_weak_source == 1, "Source month", "Sink month"),
    observed_month_class = factor(observed_month_class, levels = c("Sink month", "Source month"))
  ) %>%
  ggplot(aes(x = observed_month_class, y = loo_pred_source_probability, fill = observed_month_class)) +
  geom_hline(yintercept = source_probability_threshold_default, linetype = "dashed", color = "grey25") +
  geom_boxplot(width = 0.55, alpha = 0.65, outlier.shape = NA) +
  geom_jitter(aes(shape = EcoType), width = 0.12, height = 0, size = 1.5, color = "grey15", alpha = 0.65) +
  facet_wrap(~model_scenario, ncol = 1) +
  scale_fill_manual(values = c("Sink month" = "#2166AC", "Source month" = "#B2182B"), drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Monthly Condition-Based Source Probability Model: Leave-One-Site-Out Predictions",
    subtitle = paste0("Site-blocked validation; dashed line = ", source_probability_threshold_default, " weak-source threshold"),
    x = "Observed monthly ERA5 source state",
    y = "Predicted P(source month)",
    fill = "Observed month",
    shape = "Ecosystem"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  "FIGURES/ERA_SpatialProbability_monthly_validation.png",
  plot_monthly_validation,
  width = 9,
  height = 8,
  units = "in",
  dpi = 300
)

plot_monthly_demo_estimate <- monthly_area_weighted_estimate %>%
  ggplot(aes(x = rate_scenario, y = upscaled_net_exchange_tg_ch4_yr, fill = rate_scenario)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -35, fill = "#2166AC", alpha = 0.12) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_errorbar(
    aes(
      ymin = upscaled_net_exchange_tg_ch4_yr - approximate_model_se_tg_ch4_yr,
      ymax = upscaled_net_exchange_tg_ch4_yr + approximate_model_se_tg_ch4_yr
    ),
    width = 0.12
  ) +
  facet_wrap(~model_scenario, ncol = 1) +
  coord_flip() +
  labs(
    title = "Monthly Condition-Based Probability Model: Area-Weighted Demonstration",
    subtitle = "Uses ecosystem-mean monthly predictions and placeholder global areas; not a spatial global estimate.",
    x = NULL,
    y = "Net exchange (Tg CH4 yr-1)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(
  "FIGURES/ERA_SpatialProbability_monthly_area_weighted_demo.png",
  plot_monthly_demo_estimate,
  width = 9,
  height = 4.8,
  units = "in",
  dpi = 300
)

message(
  "Wrote ERA-SpatialProbability outputs. Global spatial inputs present: ",
  global_spatial_ready
)
