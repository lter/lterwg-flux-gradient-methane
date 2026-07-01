# Model-standardized 30-minute CH4 gap filling and annual budgets.
#
# This script uses the total-flux 30-minute GAM from flow.30min.analysis.R to
# infer site behavior and annual CH4 budgets on a balanced design:
# every site gets every month and every half-hour-of-day slot. Covariates are
# filled from observed site/month/hour medians where available, then broader
# site/season/global medians when a cell was unobserved.

library(tidyverse)
library(ggplot2)
library(mgcv)

localdir.ch4 <- "/Volumes/MaloneLab/Research/FluxGradient/Methane"
setwd(localdir.ch4)

dir.create("OUTPUT", showWarnings = FALSE)
dir.create("FIGURES", showWarnings = FALSE)

model_file <- "OUTPUT/30min_ch4_models.Rdata"
model_data_file <- "OUTPUT/30min_ch4_model_data.csv"
site_behavior_file <- "OUTPUT/30min_site_behavior.csv"

if (!file.exists(model_file)) {
  stop("Missing OUTPUT/30min_ch4_models.Rdata. Run flow.30min.analysis.R first.")
}

if (!file.exists(model_data_file)) {
  stop("Missing OUTPUT/30min_ch4_model_data.csv. Run flow.30min.analysis.R first.")
}

if (!file.exists(site_behavior_file)) {
  stop("Missing OUTPUT/30min_site_behavior.csv. Run flow.30min.analysis.R first.")
}

load(model_file)

required_objects <- c("flux_model")
missing_objects <- setdiff(required_objects, ls())
if (length(missing_objects) > 0) {
  stop("Model file is missing required objects: ", paste(missing_objects, collapse = ", "))
}

ch4_30min <- read.csv(model_data_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    month = as.integer(format(as.Date(Date), "%m")),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    EcoType = as.character(EcoType)
  ) %>%
  filter(
    !is.na(SITE_ID),
    !is.na(month),
    !is.na(season),
    is.finite(hour_num),
    is.finite(Tair_C),
    is.finite(VSWCMean),
    is.finite(log_PAR)
  )

site_behavior <- read.csv(site_behavior_file) %>%
  mutate(SITE_ID = as.character(SITE_ID))

model_levels <- function(model, variable, fallback) {
  if (!is.null(model$xlevels[[variable]])) {
    return(model$xlevels[[variable]])
  }
  sort(unique(fallback))
}

site_levels <- model_levels(flux_model, "SITE_ID", ch4_30min$SITE_ID)
season_levels <- model_levels(flux_model, "season", as.character(ch4_30min$season))
ecotype_levels <- model_levels(flux_model, "EcoType", ch4_30min$EcoType)

month_lookup <- tibble(
  month = 1:12,
  month_name = month.abb,
  days_in_month = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31),
  season = factor(
    c("Winter", "Winter", "Spring", "Spring", "Spring", "Summer",
      "Summer", "Summer", "Autumn", "Autumn", "Autumn", "Winter"),
    levels = c("Winter", "Spring", "Summer", "Autumn")
  )
)

covariates <- c("Tair_C", "VSWCMean", "log_PAR")
hour_grid <- seq(0, 23.5, by = 0.5)

site_lookup <- ch4_30min %>%
  reframe(
    .by = SITE_ID,
    EcoType = names(sort(table(EcoType), decreasing = TRUE))[1]
  )

summarise_covariates <- function(data, group_vars, suffix) {
  data %>%
    reframe(
      .by = all_of(group_vars),
      Tair_C = median(Tair_C, na.rm = TRUE),
      VSWCMean = median(VSWCMean, na.rm = TRUE),
      log_PAR = median(log_PAR, na.rm = TRUE)
    ) %>%
    rename_with(~ paste0(.x, suffix), all_of(covariates))
}

site_month_hour <- summarise_covariates(ch4_30min, c("SITE_ID", "month", "hour_num"), "_site_month_hour")
site_month <- summarise_covariates(ch4_30min, c("SITE_ID", "month"), "_site_month")
site_season_hour <- summarise_covariates(ch4_30min, c("SITE_ID", "season", "hour_num"), "_site_season_hour")
site_season <- summarise_covariates(ch4_30min, c("SITE_ID", "season"), "_site_season")
site_all <- summarise_covariates(ch4_30min, "SITE_ID", "_site")
global_month_hour <- summarise_covariates(ch4_30min, c("month", "hour_num"), "_global_month_hour")
global_month <- summarise_covariates(ch4_30min, "month", "_global_month")
global_all <- ch4_30min %>%
  summarise(
    Tair_C_global = median(Tair_C, na.rm = TRUE),
    VSWCMean_global = median(VSWCMean, na.rm = TRUE),
    log_PAR_global = median(log_PAR, na.rm = TRUE)
  )

prediction_grid <- expand_grid(
  SITE_ID = sort(unique(ch4_30min$SITE_ID)),
  month = 1:12,
  hour_num = hour_grid
) %>%
  left_join(month_lookup, by = "month") %>%
  left_join(site_lookup, by = "SITE_ID") %>%
  left_join(site_month_hour, by = c("SITE_ID", "month", "hour_num")) %>%
  left_join(site_month, by = c("SITE_ID", "month")) %>%
  left_join(site_season_hour, by = c("SITE_ID", "season", "hour_num")) %>%
  left_join(site_season, by = c("SITE_ID", "season")) %>%
  left_join(site_all, by = "SITE_ID") %>%
  left_join(global_month_hour, by = c("month", "hour_num")) %>%
  left_join(global_month, by = "month") %>%
  mutate(
    Tair_C = coalesce(
      Tair_C_site_month_hour, Tair_C_site_month, Tair_C_site_season_hour,
      Tair_C_site_season, Tair_C_site, Tair_C_global_month_hour,
      Tair_C_global_month, global_all$Tair_C_global
    ),
    VSWCMean = coalesce(
      VSWCMean_site_month_hour, VSWCMean_site_month, VSWCMean_site_season_hour,
      VSWCMean_site_season, VSWCMean_site, VSWCMean_global_month_hour,
      VSWCMean_global_month, global_all$VSWCMean_global
    ),
    log_PAR = coalesce(
      log_PAR_site_month_hour, log_PAR_site_month, log_PAR_site_season_hour,
      log_PAR_site_season, log_PAR_site, log_PAR_global_month_hour,
      log_PAR_global_month, global_all$log_PAR_global
    ),
    covariate_source = case_when(
      if_all(ends_with("_site_month_hour"), ~ !is.na(.x)) ~ "site_month_hour",
      if_all(ends_with("_site_month"), ~ !is.na(.x)) ~ "site_month",
      if_all(ends_with("_site_season_hour"), ~ !is.na(.x)) ~ "site_season_hour",
      if_all(ends_with("_site_season"), ~ !is.na(.x)) ~ "site_season",
      if_all(ends_with("_site"), ~ !is.na(.x)) ~ "site",
      if_all(ends_with("_global_month_hour"), ~ !is.na(.x)) ~ "global_month_hour",
      if_all(ends_with("_global_month"), ~ !is.na(.x)) ~ "global_month",
      TRUE ~ "global"
    ),
    sin_hour = sin(2 * pi * hour_num / 24),
    cos_hour = cos(2 * pi * hour_num / 24),
    SITE_ID = factor(SITE_ID, levels = site_levels),
    season = factor(season, levels = season_levels),
    EcoType = factor(EcoType, levels = ecotype_levels)
  ) %>%
  dplyr::select(
    SITE_ID, month, month_name, season, days_in_month, hour_num,
    sin_hour, cos_hour, Tair_C, VSWCMean, log_PAR, EcoType, covariate_source
  )

prediction_grid$pred_flux_mgC_m2_30min <- predict(
  flux_model,
  newdata = prediction_grid,
  type = "response"
)

draw_coefficients <- function(model, n_sims = 1000) {
  beta <- coef(model)
  covariance <- vcov(model, unconditional = TRUE)
  covariance <- (covariance + t(covariance)) / 2
  eig <- eigen(covariance, symmetric = TRUE)
  values <- pmax(eig$values, 0)
  z <- matrix(rnorm(length(beta) * n_sims), nrow = length(beta), ncol = n_sims)
  sweep(eig$vectors %*% (sqrt(values) * z), 1, beta, "+")
}

set.seed(20260514)
n_sims <- 1000
model_matrix <- predict(flux_model, newdata = prediction_grid, type = "lpmatrix")
beta_draws <- draw_coefficients(flux_model, n_sims = n_sims)
pred_draws <- model_matrix %*% beta_draws

row_weight <- prediction_grid$days_in_month
month_group <- paste(prediction_grid$SITE_ID, prediction_grid$month, sep = "__")
annual_group <- as.character(prediction_grid$SITE_ID)

monthly_budget_draws <- rowsum(pred_draws * row_weight, group = month_group, reorder = FALSE)
annual_budget_draws <- rowsum(pred_draws * row_weight, group = annual_group, reorder = FALSE)

monthly_keys <- tibble(month_group = rownames(monthly_budget_draws)) %>%
  separate(month_group, into = c("SITE_ID", "month"), sep = "__", convert = TRUE) %>%
  left_join(month_lookup, by = "month")

monthly_budgets <- as_tibble(monthly_budget_draws, .name_repair = ~ paste0("sim_", seq_along(.x))) %>%
  mutate(SITE_ID = monthly_keys$SITE_ID, month = monthly_keys$month, month_name = monthly_keys$month_name, season = monthly_keys$season) %>%
  pivot_longer(starts_with("sim_"), names_to = "simulation", values_to = "budget_mgC_m2_month") %>%
  reframe(
    .by = c(SITE_ID, month, month_name, season),
    budget_mean_mgC_m2_month = mean(budget_mgC_m2_month),
    budget_median_mgC_m2_month = median(budget_mgC_m2_month),
    budget_lwr_mgC_m2_month = quantile(budget_mgC_m2_month, 0.025),
    budget_upr_mgC_m2_month = quantile(budget_mgC_m2_month, 0.975),
    prob_source_month = mean(budget_mgC_m2_month > 0)
  )

annual_budgets <- as_tibble(annual_budget_draws, .name_repair = ~ paste0("sim_", seq_along(.x))) %>%
  mutate(SITE_ID = rownames(annual_budget_draws)) %>%
  pivot_longer(starts_with("sim_"), names_to = "simulation", values_to = "budget_mgC_m2_yr") %>%
  reframe(
    .by = SITE_ID,
    annual_budget_mean_mgC_m2_yr = mean(budget_mgC_m2_yr),
    annual_budget_median_mgC_m2_yr = median(budget_mgC_m2_yr),
    annual_budget_lwr_mgC_m2_yr = quantile(budget_mgC_m2_yr, 0.025),
    annual_budget_upr_mgC_m2_yr = quantile(budget_mgC_m2_yr, 0.975),
    annual_budget_sd_mgC_m2_yr = sd(budget_mgC_m2_yr),
    prob_annual_source = mean(budget_mgC_m2_yr > 0)
  ) %>%
  mutate(
    annual_budget_mean_gC_m2_yr = annual_budget_mean_mgC_m2_yr / 1000,
    annual_budget_median_gC_m2_yr = annual_budget_median_mgC_m2_yr / 1000,
    annual_budget_lwr_gC_m2_yr = annual_budget_lwr_mgC_m2_yr / 1000,
    annual_budget_upr_gC_m2_yr = annual_budget_upr_mgC_m2_yr / 1000,
    annual_budget_sd_gC_m2_yr = annual_budget_sd_mgC_m2_yr / 1000
  )

standardized_behavior <- monthly_budgets %>%
  reframe(
    .by = SITE_ID,
    standardized_source_months = sum(budget_mean_mgC_m2_month > 0, na.rm = TRUE),
    standardized_prop_source_months = mean(budget_mean_mgC_m2_month > 0, na.rm = TRUE),
    standardized_behavior = case_when(
      standardized_prop_source_months >= 0.75 ~ "Consistent source",
      standardized_prop_source_months <= 0.25 ~ "Consistent sink",
      TRUE ~ "Fluctuating"
    )
  )

annual_budgets <- annual_budgets %>%
  left_join(standardized_behavior, by = "SITE_ID") %>%
  left_join(site_behavior %>% dplyr::select(SITE_ID, observed_behavior = CH4_behavior), by = "SITE_ID") %>%
  mutate(
    annual_behavior = case_when(
      prob_annual_source >= 0.75 ~ "Consistent source",
      prob_annual_source <= 0.25 ~ "Consistent sink",
      TRUE ~ "Fluctuating"
    ),
    standardized_behavior = factor(standardized_behavior, levels = c("Consistent sink", "Fluctuating", "Consistent source")),
    annual_behavior = factor(annual_behavior, levels = c("Consistent sink", "Fluctuating", "Consistent source")),
    observed_behavior = factor(observed_behavior, levels = c("Consistent sink", "Fluctuating", "Consistent source"))
  ) %>%
  arrange(standardized_behavior, annual_budget_mean_gC_m2_yr)

prediction_grid_out <- prediction_grid %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    season = as.character(season),
    EcoType = as.character(EcoType)
  )

write.csv(prediction_grid_out, "OUTPUT/NON_30min_gapfill_prediction_grid.csv", row.names = FALSE)
write.csv(monthly_budgets, "OUTPUT/NON_30min_gapfill_monthly_budgets.csv", row.names = FALSE)
write.csv(annual_budgets, "OUTPUT/NEON_30min_gapfill_annual_budgets.csv", row.names = FALSE)

# Color convention: blue = sink (uptake), grey = fluctuating, red = source (emission)
budget_colors <- c(
  "Consistent sink"   = "#2166AC",
  "Fluctuating"       = "#4D4D4D",
  "Consistent source" = "#B2182B"
)

class_change_summary <- annual_budgets %>%
  count(observed_behavior, standardized_behavior, name = "n_sites") %>%
  mutate(
    changed = observed_behavior != standardized_behavior,
    observed_behavior = factor(observed_behavior, levels = levels(annual_budgets$observed_behavior)),
    standardized_behavior = factor(standardized_behavior, levels = levels(annual_budgets$standardized_behavior))
  )

write.csv(class_change_summary, "OUTPUT/NON_30min_gapfill_class_changes.csv", row.names = FALSE)

annual_class_change_summary <- annual_budgets %>%
  count(observed_behavior, annual_behavior, name = "n_sites") %>%
  mutate(
    changed = observed_behavior != annual_behavior,
    observed_behavior = factor(observed_behavior, levels = levels(annual_budgets$observed_behavior)),
    annual_behavior = factor(annual_behavior, levels = levels(annual_budgets$annual_behavior))
  )

write.csv(
  annual_class_change_summary,
  "OUTPUT/NON_30min_gapfill_observed_vs_annual_class_changes.csv",
  row.names = FALSE
)

class_levels <- c("Consistent sink", "Fluctuating", "Consistent source")
class_change_matrix <- expand_grid(
  observed_behavior = factor(class_levels, levels = class_levels),
  standardized_behavior = factor(class_levels, levels = class_levels)
) %>%
  left_join(class_change_summary, by = c("observed_behavior", "standardized_behavior")) %>%
  mutate(
    n_sites = replace_na(n_sites, 0L),
    changed = observed_behavior != standardized_behavior,
    label = if_else(n_sites > 0, as.character(n_sites), "")
  )

annual_class_change_matrix <- expand_grid(
  observed_behavior = factor(class_levels, levels = class_levels),
  annual_behavior = factor(class_levels, levels = class_levels)
) %>%
  left_join(annual_class_change_summary, by = c("observed_behavior", "annual_behavior")) %>%
  mutate(
    n_sites = replace_na(n_sites, 0L),
    changed = observed_behavior != annual_behavior,
    label = if_else(n_sites > 0, as.character(n_sites), "")
  )

behavior_counts <- annual_budgets %>%
  count(standardized_behavior, name = "n_sites") %>%
  mutate(line = paste0("- ", standardized_behavior, ": ", n_sites, " sites")) %>%
  pull(line)

annual_behavior_counts <- annual_budgets %>%
  count(annual_behavior, name = "n_sites") %>%
  mutate(line = paste0("- ", annual_behavior, ": ", n_sites, " sites")) %>%
  pull(line)

annual_class_change_lines <- annual_class_change_summary %>%
  arrange(observed_behavior, annual_behavior) %>%
  mutate(
    line = paste0(
      "- Observed ", observed_behavior, " -> annual ", annual_behavior,
      ": ", n_sites, " sites",
      if_else(changed, " changed", " unchanged")
    )
  ) %>%
  pull(line)

top_budget_lines <- annual_budgets %>%
  arrange(desc(abs(annual_budget_mean_gC_m2_yr))) %>%
  slice_head(n = 8) %>%
  mutate(
    line = paste0(
      "- ", SITE_ID, " (", standardized_behavior, "): ",
      signif(annual_budget_mean_gC_m2_yr, 3), " g C m-2 yr-1 [",
      signif(annual_budget_lwr_gC_m2_yr, 3), ", ",
      signif(annual_budget_upr_gC_m2_yr, 3), "], P(source) = ",
      signif(prob_annual_source, 3)
    )
  ) %>%
  pull(line)

writeLines(
  c(
    "# NON 30-minute Model-Standardized Gapfill",
    "",
    "## Design",
    "- Grid: every site x 12 months x 48 half-hour-of-day slots.",
    "- Covariates: observed site/month/hour medians where available, falling back to broader site/season/global medians.",
    "- Annual budget: predicted mg C m-2 per 30 min summed across half-hours and weighted by days per month.",
    "- Uncertainty: 1,000 simulations from the fitted GAM coefficient covariance matrix; intervals describe model-estimation uncertainty in the expected budget.",
    "",
    "## Standardized Behavior Counts",
    behavior_counts,
    "",
    "## Annual-Budget Behavior Counts",
    annual_behavior_counts,
    "",
    "## Observed 30-Minute vs Annual-Budget Class Changes",
    annual_class_change_lines,
    "",
    "## Largest Absolute Annual Budgets",
    top_budget_lines,
    "",
    "## Outputs",
    "- `OUTPUT/NON_30min_gapfill_prediction_grid.csv`",
    "- `OUTPUT/NON_30min_gapfill_monthly_budgets.csv`",
    "- `OUTPUT/NEON_30min_gapfill_annual_budgets.csv`",
    "- `OUTPUT/NON_30min_gapfill_class_changes.csv`",
    "- `OUTPUT/NON_30min_gapfill_observed_vs_annual_class_changes.csv`"
  ),
  "OUTPUT/NON_30min_gapfill_results.md"
)

message("Wrote model-standardized 30-minute gapfill annual budgets.")
