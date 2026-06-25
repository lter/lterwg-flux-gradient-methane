# Supplementary model-result and driver/site-attribute plots for the
# 30-minute NEON total-flux CH4 workflow.

library(tidyverse)
library(ggplot2)
library(mgcv)
library(patchwork)
library(ggrepel)

localdir.ch4 <- "/Volumes/MaloneLab/Research/FluxGradient/Methane"
setwd(localdir.ch4)

dir.create("FIGURES", showWarnings = FALSE)
dir.create("OUTPUT", showWarnings = FALSE)

required_files <- c(
  "OUTPUT/30min_ch4_models.Rdata",
  "OUTPUT/30min_ch4_model_data.csv",
  "OUTPUT/30min_hourly_predictions.csv",
  "OUTPUT/30min_driver_predictions.csv",
  "OUTPUT/30min_site_behavior.csv",
  "OUTPUT/NEON_strong_sink_driver_comparison.csv",
  "OUTPUT/NEON_site_driver_values_for_sink_comparison.csv",
  "OUTPUT/NEON_30min_gapfill_annual_budgets.csv",
  "OUTPUT/NON_30min_gapfill_monthly_budgets.csv",
  "OUTPUT/NEON_consistency_magnitude_predictor_scores.csv",
  "OUTPUT/NEON_consistency_magnitude_tree_importance.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

load("OUTPUT/30min_ch4_models.Rdata")

behavior_levels <- c("Weak-sink", "Fluctuating", "Weak-source")
behavior_colors <- c(
  "Weak-sink" = "#2166AC",
  "Fluctuating" = "#4D4D4D",
  "Weak-source" = "#B2182B"
)
driver_group_colors <- c(
  "Moisture" = "#1B9E77",
  "Soil texture/physical" = "#D95F02",
  "Soil chemistry/depth" = "#7570B3",
  "Climate" = "#E7298A",
  "Vegetation/canopy" = "#66A61E"
)
season_colors <- c(
  "Winter" = "#5E81AC",
  "Spring" = "#4C9F70",
  "Summer" = "#D99000",
  "Autumn" = "#A65628"
)

site_behavior <- read.csv("OUTPUT/30min_site_behavior.csv") %>%
  mutate(CH4_behavior = factor(CH4_behavior, levels = behavior_levels))

ch4_30min <- read.csv("OUTPUT/30min_ch4_model_data.csv") %>%
  mutate(
    SITE_ID = factor(SITE_ID),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    EcoType = factor(EcoType),
    source_30min = as.logical(source_30min)
  ) %>%
  left_join(
    site_behavior %>% mutate(SITE_ID = as.character(SITE_ID)) %>% dplyr::select(SITE_ID, CH4_behavior),
    by = "SITE_ID"
  )

hour_predictions <- read.csv("OUTPUT/30min_hourly_predictions.csv") %>%
  mutate(season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")))

driver_predictions <- read.csv("OUTPUT/30min_driver_predictions.csv")

driver_comparison <- read.csv("OUTPUT/NEON_strong_sink_driver_comparison.csv")
site_drivers <- read.csv("OUTPUT/NEON_site_driver_values_for_sink_comparison.csv") %>%
  mutate(
    CH4_behavior = factor(CH4_behavior, levels = behavior_levels),
    behavior_comparison = factor(behavior_comparison)
  )

annual_budgets <- read.csv("OUTPUT/NEON_30min_gapfill_annual_budgets.csv") %>%
  mutate(
    standardized_behavior = factor(standardized_behavior, levels = behavior_levels),
    observed_behavior = factor(observed_behavior, levels = behavior_levels)
  )

monthly_budgets <- read.csv("OUTPUT/NON_30min_gapfill_monthly_budgets.csv") %>%
  mutate(
    month = as.integer(month),
    month_name = factor(month_name, levels = month.abb),
    monthly_budget_gC_m2 = budget_mean_mgC_m2_month / 1000
  )

predictor_scores <- read.csv("OUTPUT/NEON_consistency_magnitude_predictor_scores.csv")
tree_importance <- read.csv("OUTPUT/NEON_consistency_magnitude_tree_importance.csv")

set.seed(20260514)
diagnostic_data <- ch4_30min %>%
  mutate(
    flux_fitted = as.numeric(predict(flux_model, newdata = ch4_30min, type = "response")),
    flux_residual = CH4_mgC_30min - flux_fitted,
    source_pred = as.numeric(predict(source_model, newdata = ch4_30min, type = "response"))
  )

diagnostic_sample <- diagnostic_data %>%
  slice_sample(n = min(50000, nrow(diagnostic_data)))

plot_observed_fitted <- diagnostic_sample %>%
  ggplot(aes(x = flux_fitted, y = CH4_mgC_30min)) +
  geom_bin2d(bins = 80) +
  geom_abline(slope = 1, intercept = 0, color = "white", linetype = "dashed", linewidth = 0.7) +
  scale_fill_viridis_c(trans = "log10") +
  labs(
    title = "A. Observed vs fitted total flux",
    x = expression(paste("Fitted CH"[4], " flux (mg C ", m^-2, " 30 min"^-1, ")")),
    y = expression(paste("Observed CH"[4], " flux (mg C ", m^-2, " 30 min"^-1, ")")),
    fill = "Count"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"))

plot_residual_fitted <- diagnostic_sample %>%
  ggplot(aes(x = flux_fitted, y = flux_residual)) +
  geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
  geom_bin2d(bins = 80) +
  scale_fill_viridis_c(trans = "log10") +
  labs(
    title = "B. Residuals vs fitted",
    x = expression(paste("Fitted CH"[4], " flux")),
    y = "Residual",
    fill = "Count"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"))

qq_sample <- diagnostic_sample %>%
  arrange(flux_residual) %>%
  mutate(
    sample_quantile = flux_residual,
    theoretical_quantile = qnorm((row_number() - 0.5) / dplyr::n())
  )

plot_residual_qq <- qq_sample %>%
  ggplot(aes(x = theoretical_quantile, y = sample_quantile)) +
  geom_point(alpha = 0.18, size = 0.45) +
  geom_abline(
    intercept = median(qq_sample$sample_quantile, na.rm = TRUE),
    slope = IQR(qq_sample$sample_quantile, na.rm = TRUE) / IQR(qq_sample$theoretical_quantile, na.rm = TRUE),
    color = behavior_colors[["Weak-sink"]],
    linetype = "dashed"
  ) +
  labs(
    title = "C. Flux-model residual quantiles",
    x = "Theoretical normal quantile",
    y = "Residual quantile"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"))

plot_residual_season <- diagnostic_sample %>%
  ggplot(aes(x = season, y = flux_residual, color = season)) +
  geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  coord_cartesian(ylim = quantile(diagnostic_sample$flux_residual, c(0.01, 0.99), na.rm = TRUE)) +
  labs(
    title = "D. Residuals by season",
    x = NULL,
    y = "Residual"
  ) +
  theme_bw(base_size = 10.5) +
  scale_color_manual(values = season_colors, na.translate = FALSE) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

model_diagnostics_figure <- (plot_observed_fitted + plot_residual_fitted) /
  (plot_residual_qq + plot_residual_season) +
  plot_annotation(
    title = "Supplementary Figure S1. Thirty-minute total-flux GAM diagnostics",
    subtitle = "Diagnostics are plotted from a random sample of up to 50,000 observations.",
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  )

ggsave("FIGURES/SUPP_30min_flux_model_diagnostics.png", model_diagnostics_figure, width = 11, height = 9, units = "in", dpi = 300)
ggsave("FIGURES/SUPP_30min_flux_model_diagnostics.pdf", model_diagnostics_figure, width = 11, height = 9, units = "in")

calibration_data <- diagnostic_data %>%
  mutate(prob_bin = ntile(source_pred, 12)) %>%
  reframe(
    .by = prob_bin,
    mean_pred = mean(source_pred, na.rm = TRUE),
    observed_source = mean(source_30min, na.rm = TRUE),
    n = dplyr::n(),
    se = sqrt(observed_source * (1 - observed_source) / n),
    lwr = pmax(0, observed_source - 1.96 * se),
    upr = pmin(1, observed_source + 1.96 * se)
  )

plot_source_calibration <- calibration_data %>%
  ggplot(aes(x = mean_pred, y = observed_source)) +
  geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.01, alpha = 0.8) +
  geom_point(aes(size = n), color = behavior_colors[["Weak-source"]], alpha = 0.85) +
  scale_size(range = c(2, 6)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "A. Source-state calibration",
    x = "Mean predicted source probability",
    y = "Observed source frequency",
    size = "Observations"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

plot_source_by_behavior <- diagnostic_data %>%
  mutate(pred_bin = cut(source_pred, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
  reframe(
    .by = c(CH4_behavior, pred_bin),
    mean_pred = mean(source_pred, na.rm = TRUE),
    observed_source = mean(source_30min, na.rm = TRUE),
    n = dplyr::n()
  ) %>%
  filter(!is.na(CH4_behavior), n >= 20) %>%
  ggplot(aes(x = mean_pred, y = observed_source, color = CH4_behavior)) +
  geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  geom_point(aes(size = n), alpha = 0.85) +
  scale_color_manual(values = behavior_colors, na.translate = FALSE) +
  scale_size(range = c(1.5, 5)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "B. Calibration by behavior class",
    x = "Mean predicted source probability",
    y = "Observed source frequency",
    color = "Observed class",
    size = "Observations"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

plot_source_prob_distribution <- diagnostic_sample %>%
  ggplot(aes(x = source_pred, fill = CH4_behavior)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = behavior_colors, na.translate = FALSE) +
  labs(
    title = "C. Distribution of predicted source probabilities",
    x = "Predicted source probability",
    y = "Observations",
    fill = "Observed class"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

source_model_figure <- (plot_source_calibration + plot_source_by_behavior) / plot_source_prob_distribution +
  plot_layout(heights = c(1, 0.85)) +
  plot_annotation(
    title = "Supplementary Figure S2. Source-state model calibration",
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  )

ggsave("FIGURES/SUPP_30min_source_model_calibration.png", source_model_figure, width = 11, height = 8.5, units = "in", dpi = 300)
ggsave("FIGURES/SUPP_30min_source_model_calibration.pdf", source_model_figure, width = 11, height = 8.5, units = "in")

plot_hour_flux <- hour_predictions %>%
  ggplot(aes(x = hour_num, y = pred_flux, color = season)) +
  geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
  geom_line(linewidth = 0.95) +
  labs(
    title = "A. Diel total-flux prediction",
    x = "Hour of day",
    y = expression(paste("Predicted CH"[4], " flux (mg C ", m^-2, " 30 min"^-1, ")")),
    color = "Season"
  ) +
  theme_bw(base_size = 10.5) +
  scale_color_manual(values = season_colors, na.translate = FALSE) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

plot_hour_source <- hour_predictions %>%
  ggplot(aes(x = hour_num, y = pred_source_prob, color = season)) +
  geom_hline(yintercept = 0.5, color = "grey45", linetype = "dotted") +
  geom_line(linewidth = 0.95) +
  labs(
    title = "B. Diel source-probability prediction",
    x = "Hour of day",
    y = "Predicted source probability",
    color = "Season"
  ) +
  theme_bw(base_size = 10.5) +
  scale_color_manual(values = season_colors, na.translate = FALSE) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

driver_prediction_long <- driver_predictions %>%
  pivot_longer(c(pred_flux, pred_source_prob), names_to = "response", values_to = "prediction") %>%
  mutate(
    response = recode(response, pred_flux = "Total CH4 flux", pred_source_prob = "Source probability"),
    driver = recode(driver, Tair_C = "Air temperature", VSWCMean = "Soil water content", log_PAR = "log(PAR + 1)")
  )

plot_driver_response <- driver_prediction_long %>%
  ggplot(aes(x = driver_value, y = prediction)) +
  geom_hline(data = tibble(response = "Total CH4 flux", y = 0), aes(yintercept = y), inherit.aes = FALSE, color = "grey45", linetype = "dashed") +
  geom_line(linewidth = 0.95, color = "black") +
  facet_grid(response ~ driver, scales = "free") +
  labs(
    title = "C. Marginal model predictions",
    x = "Driver value",
    y = "Prediction"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"))

model_prediction_figure <- (plot_hour_flux + plot_hour_source) / plot_driver_response +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Supplementary Figure S3. Model-predicted diel and environmental responses",
    subtitle = "Population-level predictions exclude the site random effect.",
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  )

ggsave("FIGURES/SUPP_30min_model_predictions.png", model_prediction_figure, width = 12, height = 9, units = "in", dpi = 300)
ggsave("FIGURES/SUPP_30min_model_predictions.pdf", model_prediction_figure, width = 12, height = 9, units = "in")

driver_effect_plot <- driver_comparison %>%
  mutate(
    label = stringr::str_wrap(label, width = 22),
    label = factor(label, levels = rev(label[order(cliffs_delta)])),
    evidence = case_when(
      is.na(p_value) ~ "not tested",
      p_value < 0.05 ~ "p < 0.05",
      p_value < 0.10 ~ "p < 0.10",
      TRUE ~ "weak"
    ),
    evidence = factor(evidence, levels = c("p < 0.05", "p < 0.10", "weak", "not tested"))
  ) %>%
  ggplot(aes(x = cliffs_delta, y = label, color = driver_group, shape = evidence)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = cliffs_delta, yend = label), alpha = 0.7) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("p < 0.05" = 16, "p < 0.10" = 17, "weak" = 1, "not tested" = 4)) +
  scale_color_manual(values = driver_group_colors, na.translate = FALSE) +
  labs(
    title = "A. Site-attribute contrasts between behavior classes",
    subtitle = "Negative values are lower in consistent-source sites than fluctuating sites.",
    x = "Cliff's delta",
    y = NULL,
    color = "Attribute group",
    shape = "Evidence"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

top_driver_vars <- driver_comparison %>%
  slice_max(abs(cliffs_delta), n = min(12, nrow(driver_comparison)), with_ties = FALSE) %>%
  pull(variable)

top_driver_labels <- driver_comparison %>%
  filter(variable %in% top_driver_vars) %>%
  dplyr::select(variable, label, driver_group)

driver_distribution_plot <- site_drivers %>%
  dplyr::select(SITE_ID, behavior_comparison, all_of(top_driver_vars)) %>%
  pivot_longer(all_of(top_driver_vars), names_to = "variable", values_to = "value") %>%
  filter(is.finite(value), !is.na(behavior_comparison)) %>%
  group_by(variable) %>%
  mutate(value_z = as.numeric(scale(value))) %>%
  ungroup() %>%
  left_join(top_driver_labels, by = "variable") %>%
  mutate(label = factor(stringr::str_wrap(label, 18), levels = stringr::str_wrap(top_driver_labels$label, 18))) %>%
  ggplot(aes(x = behavior_comparison, y = value_z, color = behavior_comparison)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.15, width = 0.55) +
  geom_jitter(width = 0.13, alpha = 0.75, size = 1.7) +
  facet_wrap(~label, scales = "free_y", ncol = 4) +
  scale_color_manual(values = behavior_colors, na.translate = FALSE) +
  labs(
    title = "B. Top site-attribute distributions",
    x = NULL,
    y = "Standardized site value",
    color = "Behavior contrast"
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
    legend.position = "bottom"
  )

driver_category_figure <- driver_effect_plot / driver_distribution_plot +
  plot_layout(heights = c(1.1, 1.25)) +
  plot_annotation(
    title = "Supplementary Figure S4. Conditions associated with observed behavior categories",
    subtitle = "Site-level comparisons describe categories; they are not causal tests of half-hour flux drivers.",
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  )

ggsave("FIGURES/SUPP_total_flux_behavior_driver_contrasts.png", driver_category_figure, width = 13, height = 13, units = "in", dpi = 300)
ggsave("FIGURES/SUPP_total_flux_behavior_driver_contrasts.pdf", driver_category_figure, width = 13, height = 13, units = "in")

top_predictor_scores <- predictor_scores %>%
  slice_head(n = min(20, nrow(predictor_scores))) %>%
  mutate(
    variable = stringr::str_wrap(variable, 24),
    variable = factor(variable, levels = rev(variable))
  )

plot_predictor_scores <- top_predictor_scores %>%
  dplyr::select(variable, budget_spearman, source_month_spearman) %>%
  pivot_longer(c(budget_spearman, source_month_spearman), names_to = "association", values_to = "rho") %>%
  mutate(association = recode(association, budget_spearman = "Annual budget", source_month_spearman = "Source-month fraction")) %>%
  ggplot(aes(x = rho, y = variable, color = association)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  geom_point(size = 2.7, position = position_dodge(width = 0.45)) +
  labs(
    title = "A. Ranked univariate predictor associations",
    x = "Spearman rho",
    y = NULL,
    color = "Response"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

plot_tree_importance <- tree_importance %>%
  mutate(
    variable = stringr::str_wrap(variable, 24),
    variable = fct_reorder(variable, relative_importance)
  ) %>%
  ggplot(aes(x = relative_importance, y = variable, fill = model)) +
  geom_col(width = 0.72) +
  facet_wrap(~model, scales = "free_y") +
  labs(
    title = "B. Decision-tree variable importance",
    x = "Relative importance within tree",
    y = NULL,
    fill = "Tree"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

plot_budget_predictor <- annual_budgets %>%
  left_join(site_drivers %>% dplyr::select(SITE_ID, MAP, MAT, siltTotal, sandTotal), by = "SITE_ID") %>%
  pivot_longer(c(MAP, MAT, siltTotal, sandTotal), names_to = "variable", values_to = "value") %>%
  mutate(variable = recode(variable, MAP = "Mean annual precipitation", MAT = "Mean annual temperature", siltTotal = "Silt", sandTotal = "Sand")) %>%
  filter(is.finite(value)) %>%
  ggplot(aes(x = value, y = annual_budget_mean_gC_m2_yr, color = standardized_behavior)) +
  geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.45, linetype = "dotted") +
  facet_wrap(~variable, scales = "free_x", ncol = 2) +
  scale_color_manual(values = behavior_colors, na.translate = FALSE) +
  labs(
    title = "C. Annual budget versus selected site attributes",
    x = NULL,
    y = expression(paste("Annual CH"[4], " budget (g C ", m^-2, " yr"^-1, ")")),
    color = "Standardized behavior"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

predictor_summary_figure <- (plot_predictor_scores + plot_tree_importance) / plot_budget_predictor +
  plot_layout(heights = c(1.1, 1.1)) +
  plot_annotation(
    title = "Supplementary Figure S5. Candidate predictors of consistency and magnitude",
    subtitle = "Predictor screens and trees are exploratory summaries of site-attribute space.",
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  )

ggsave("FIGURES/SUPP_consistency_magnitude_predictors.png", predictor_summary_figure, width = 13, height = 11, units = "in", dpi = 300)
ggsave("FIGURES/SUPP_consistency_magnitude_predictors.pdf", predictor_summary_figure, width = 13, height = 11, units = "in")

monthly_behavior <- monthly_budgets %>%
  left_join(annual_budgets %>% dplyr::select(SITE_ID, standardized_behavior, annual_budget_mean_gC_m2_yr), by = "SITE_ID") %>%
  mutate(
    SITE_ID = fct_reorder(SITE_ID, annual_budget_mean_gC_m2_yr),
    standardized_behavior = factor(standardized_behavior, levels = behavior_levels)
  )

plot_monthly_budget_uncertainty <- monthly_behavior %>%
  ggplot(aes(x = month, y = budget_mean_mgC_m2_month / 1000, color = standardized_behavior, fill = standardized_behavior)) +
  geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
  geom_ribbon(
    aes(ymin = budget_lwr_mgC_m2_month / 1000, ymax = budget_upr_mgC_m2_month / 1000),
    alpha = 0.09,
    color = NA
  ) +
  geom_line(aes(group = SITE_ID), alpha = 0.35, linewidth = 0.45) +
  stat_summary(aes(group = standardized_behavior), fun = mean, geom = "line", linewidth = 1.1) +
  scale_color_manual(values = behavior_colors, na.translate = FALSE) +
  scale_fill_manual(values = behavior_colors, na.translate = FALSE) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title = "A. Monthly standardized budgets by site",
    x = NULL,
    y = expression(paste("Monthly CH"[4], " budget (g C ", m^-2, ")")),
    color = "Standardized behavior",
    fill = "Standardized behavior"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

plot_annual_uncertainty_width <- annual_budgets %>%
  mutate(SITE_ID = fct_reorder(SITE_ID, annual_budget_upr_gC_m2_yr - annual_budget_lwr_gC_m2_yr)) %>%
  ggplot(aes(x = annual_budget_upr_gC_m2_yr - annual_budget_lwr_gC_m2_yr, y = SITE_ID, color = standardized_behavior)) +
  geom_point(size = 2.2) +
  facet_grid(standardized_behavior ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
  scale_color_manual(values = behavior_colors, na.translate = FALSE) +
  labs(
    title = "B. Annual budget uncertainty width",
    x = expression(paste("95% interval width (g C ", m^-2, " yr"^-1, ")")),
    y = NULL,
    color = "Standardized behavior"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom", axis.text.y = element_text(size = 7))

budget_uncertainty_figure <- plot_monthly_budget_uncertainty / plot_annual_uncertainty_width +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Supplementary Figure S6. Model-standardized budget uncertainty",
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  )

ggsave("FIGURES/SUPP_gapfill_budget_uncertainty.png", budget_uncertainty_figure, width = 12, height = 12, units = "in", dpi = 300)
ggsave("FIGURES/SUPP_gapfill_budget_uncertainty.pdf", budget_uncertainty_figure, width = 12, height = 12, units = "in")

figure_index <- c(
  "# Supplementary Model and Driver Figures",
  "",
  "## Shared palettes",
  "- Behavior classes: consistent sink `#2166AC`, fluctuating `#4D4D4D`, consistent source `#B2182B`.",
  "- Driver groups: moisture `#1B9E77`, soil texture/physical `#D95F02`, soil chemistry/depth `#7570B3`, climate `#E7298A`, vegetation/canopy `#66A61E`.",
  "- Seasons: winter `#5E81AC`, spring `#4C9F70`, summer `#D99000`, autumn `#A65628`.",
  "",
  "- `FIGURES/SUPP_30min_flux_model_diagnostics.png`: observed versus fitted total flux, residual diagnostics, and seasonal residual structure for the Gaussian 30-minute GAM.",
  "- `FIGURES/SUPP_30min_source_model_calibration.png`: source-state calibration and predicted source-probability distributions.",
  "- `FIGURES/SUPP_30min_model_predictions.png`: population-level diel predictions and marginal air-temperature, soil-moisture, and PAR predictions.",
  "- `FIGURES/SUPP_total_flux_behavior_driver_contrasts.png`: site-level condition contrasts between observed consistent-source and fluctuating categories.",
  "- `FIGURES/SUPP_consistency_magnitude_predictors.png`: predictor screens, tree importance, and annual budget relationships with selected site attributes.",
  "- `FIGURES/SUPP_gapfill_budget_uncertainty.png`: model-standardized monthly budget trajectories and annual uncertainty widths.",
  "",
  "PDF versions with the same base names were also written to `FIGURES/`."
)

writeLines(figure_index, "OUTPUT/NEON_supplementary_model_driver_figures.md")

message("Wrote supplementary model and driver plots.")
