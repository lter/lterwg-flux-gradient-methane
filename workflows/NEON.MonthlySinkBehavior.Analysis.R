# Monthly sink behavior analysis for NEON 30-minute total-flux CH4.
#
# This script asks what distinguishes sink months from source months. It works
# at the site-month level, where sink behavior is much more common than at the
# site-level "consistent sink" classification.

library(tidyverse)
library(ggplot2)
library(mgcv)
library(patchwork)

localdir.ch4 <- "/Volumes/MaloneLab/Research/FluxGradient/Methane"
setwd(localdir.ch4)

dir.create("OUTPUT", showWarnings = FALSE)
dir.create("FIGURES", showWarnings = FALSE)

required_files <- c(
  "OUTPUT/30min_monthly_site_ch4.csv",
  "OUTPUT/30min_site_behavior.csv",
  "OUTPUT/NON_30min_gapfill_monthly_budgets.csv",
  "OUTPUT/NON_30min_gapfill_prediction_grid.csv",
  "OUTPUT/NEON_30min_gapfill_annual_budgets.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

behavior_levels <- c("Weak-sink", "Fluctuating", "Weak-source")
behavior_colors <- c(
  "Weak-sink" = "#B2182B",
  "Fluctuating" = "#4D4D4D",
  "Weak-source" = "#2166AC"
)
sink_source_colors <- c("Sink month" = "#B2182B", "Source month" = "#2166AC")
season_colors <- c(
  "Winter" = "#5E81AC",
  "Spring" = "#4C9F70",
  "Summer" = "#D99000",
  "Autumn" = "#A65628"
)

month_lookup <- tibble(
  month = 1:12,
  month_name = factor(month.abb, levels = month.abb),
  season = factor(
    c("Winter", "Winter", "Spring", "Spring", "Spring", "Summer",
      "Summer", "Summer", "Autumn", "Autumn", "Autumn", "Winter"),
    levels = c("Winter", "Spring", "Summer", "Autumn")
  )
)

site_behavior <- read.csv("OUTPUT/30min_site_behavior.csv") %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    CH4_behavior = factor(CH4_behavior, levels = behavior_levels)
  )

monthly_observed <- read.csv("OUTPUT/30min_monthly_site_ch4.csv") %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    month = as.integer(format(as.Date(paste0(YearMon, "-01")), "%m")),
    sink_month = !source_month,
    sink_state = factor(if_else(sink_month, "Sink month", "Source month"), levels = c("Sink month", "Source month"))
  ) %>%
  left_join(month_lookup, by = "month") %>%
  left_join(site_behavior %>% dplyr::select(SITE_ID, CH4_behavior, MAP, MAT, sandTotal, siltTotal, clayTotal, dryMass, acidity, ctonRatio, estimatedOC, canopyHeight_m, LAI.mean, CHM.mean), by = "SITE_ID") %>%
  group_by(SITE_ID) %>%
  mutate(
    Tair_site_mean = mean(mean_Tair_C, na.rm = TRUE),
    VSWC_site_mean = mean(mean_VSWC, na.rm = TRUE),
    Tair_anom = mean_Tair_C - Tair_site_mean,
    VSWC_anom = mean_VSWC - VSWC_site_mean
  ) %>%
  ungroup()

standardized_monthly <- read.csv("OUTPUT/NON_30min_gapfill_monthly_budgets.csv") %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    month = as.integer(month),
    month_name = factor(month_name, levels = month.abb),
    sink_month = budget_mean_mgC_m2_month <= 0,
    sink_state = factor(if_else(sink_month, "Sink month", "Source month"), levels = c("Sink month", "Source month")),
    monthly_budget_gC_m2 = budget_mean_mgC_m2_month / 1000
  )

standardized_annual <- read.csv("OUTPUT/NEON_30min_gapfill_annual_budgets.csv") %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    standardized_behavior = factor(standardized_behavior, levels = behavior_levels)
  )

prediction_grid_monthly <- read.csv("OUTPUT/NON_30min_gapfill_prediction_grid.csv") %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    month = as.integer(month)
  ) %>%
  reframe(
    .by = c(SITE_ID, month),
    standardized_Tair_C = mean(Tair_C, na.rm = TRUE),
    standardized_VSWC = mean(VSWCMean, na.rm = TRUE),
    standardized_log_PAR = mean(log_PAR, na.rm = TRUE)
  ) %>%
  group_by(SITE_ID) %>%
  mutate(
    standardized_Tair_anom = standardized_Tair_C - mean(standardized_Tair_C, na.rm = TRUE),
    standardized_VSWC_anom = standardized_VSWC - mean(standardized_VSWC, na.rm = TRUE),
    standardized_log_PAR_anom = standardized_log_PAR - mean(standardized_log_PAR, na.rm = TRUE)
  ) %>%
  ungroup()

standardized_monthly <- standardized_monthly %>%
  left_join(standardized_annual %>% dplyr::select(SITE_ID, standardized_behavior), by = "SITE_ID") %>%
  left_join(prediction_grid_monthly, by = c("SITE_ID", "month"))

monthly_timing <- bind_rows(
  monthly_observed %>%
    reframe(
      .by = c(month, month_name, season),
      n_site_months = dplyr::n(),
      n_sink_months = sum(sink_month, na.rm = TRUE),
      sink_fraction = mean(sink_month, na.rm = TRUE)
    ) %>%
    mutate(source = "Observed balanced monthly flux"),
  standardized_monthly %>%
    reframe(
      .by = c(month, month_name, season),
      n_site_months = dplyr::n(),
      n_sink_months = sum(sink_month, na.rm = TRUE),
      sink_fraction = mean(sink_month, na.rm = TRUE)
    ) %>%
    mutate(source = "Model-standardized monthly budget")
)

monthly_site_sink_summary <- monthly_observed %>%
  reframe(
    .by = c(SITE_ID, CH4_behavior),
    n_months = dplyr::n(),
    n_sink_months = sum(sink_month, na.rm = TRUE),
    observed_sink_fraction = mean(sink_month, na.rm = TRUE),
    mean_sink_month_flux = mean(mean_CH4_mgC_30min[sink_month], na.rm = TRUE),
    mean_source_month_flux = mean(mean_CH4_mgC_30min[!sink_month], na.rm = TRUE)
  ) %>%
  left_join(
    standardized_monthly %>%
      reframe(
        .by = SITE_ID,
        standardized_sink_months = sum(sink_month, na.rm = TRUE),
        standardized_sink_fraction = mean(sink_month, na.rm = TRUE)
      ),
    by = "SITE_ID"
  ) %>%
  arrange(desc(observed_sink_fraction))

compare_month_variable <- function(data, variable, label, source_name) {
  sink_values <- data[[variable]][data$sink_month]
  source_values <- data[[variable]][!data$sink_month]
  if (sum(is.finite(sink_values)) < 2 || sum(is.finite(source_values)) < 2) {
    return(tibble())
  }
  tibble(
    source = source_name,
    variable = variable,
    label = label,
    n_sink = sum(is.finite(sink_values)),
    n_source = sum(is.finite(source_values)),
    median_sink = median(sink_values, na.rm = TRUE),
    median_source = median(source_values, na.rm = TRUE),
    median_difference_sink_minus_source = median_sink - median_source,
    p_value = suppressWarnings(wilcox.test(sink_values, source_values, exact = FALSE)$p.value)
  )
}

observed_contrasts <- tribble(
  ~variable, ~label,
  "mean_Tair_C", "Monthly air temperature",
  "Tair_anom", "Within-site temperature anomaly",
  "mean_VSWC", "Monthly VSWC",
  "VSWC_anom", "Within-site VSWC anomaly",
  "mean_gradient_mgC_30min", "Gradient component",
  "mean_storage_mgC_30min", "Storage component",
  "median_storage_abs_fraction", "Storage absolute fraction",
  "n_halfhour_bins", "Observed half-hour bins"
) %>%
  pmap_dfr(~ compare_month_variable(monthly_observed, ..1, ..2, "Observed balanced monthly flux")) %>%
  mutate(p_adj_bh = p.adjust(p_value, method = "BH"))

standardized_contrasts <- tribble(
  ~variable, ~label,
  "standardized_Tair_C", "Monthly air temperature",
  "standardized_Tair_anom", "Within-site temperature anomaly",
  "standardized_VSWC", "Monthly VSWC",
  "standardized_VSWC_anom", "Within-site VSWC anomaly",
  "standardized_log_PAR", "Monthly log(PAR + 1)",
  "standardized_log_PAR_anom", "Within-site log(PAR + 1) anomaly",
  "monthly_budget_gC_m2", "Monthly budget"
) %>%
  pmap_dfr(~ compare_month_variable(standardized_monthly, ..1, ..2, "Model-standardized monthly budget")) %>%
  mutate(p_adj_bh = p.adjust(p_value, method = "BH"))

monthly_sink_contrasts <- bind_rows(observed_contrasts, standardized_contrasts)

monthly_sink_model_data <- monthly_observed %>%
  filter(
    is.finite(Tair_anom),
    is.finite(VSWC_anom),
    is.finite(n_halfhour_bins),
    !is.na(season),
    !is.na(SITE_ID)
  ) %>%
  mutate(
    SITE_ID = factor(SITE_ID),
    season = factor(season, levels = levels(month_lookup$season))
  )

monthly_sink_model <- mgcv::bam(
  sink_month ~
    season +
    s(Tair_anom, k = 6) +
    s(VSWC_anom, k = 6) +
    s(n_halfhour_bins, k = 6) +
    s(SITE_ID, bs = "re"),
  data = monthly_sink_model_data,
  family = binomial(),
  method = "fREML",
  discrete = FALSE
)

reference_values <- monthly_sink_model_data %>%
  summarise(
    Tair_anom = median(Tair_anom, na.rm = TRUE),
    VSWC_anom = median(VSWC_anom, na.rm = TRUE),
    n_halfhour_bins = median(n_halfhour_bins, na.rm = TRUE),
    SITE_ID = names(sort(table(SITE_ID), decreasing = TRUE))[1]
  )

effect_grid <- bind_rows(
  tibble(
    driver = "Within-site temperature anomaly",
    driver_value = seq(
      quantile(monthly_sink_model_data$Tair_anom, 0.02, na.rm = TRUE),
      quantile(monthly_sink_model_data$Tair_anom, 0.98, na.rm = TRUE),
      length.out = 100
    ),
    Tair_anom = driver_value,
    VSWC_anom = reference_values$VSWC_anom,
    n_halfhour_bins = reference_values$n_halfhour_bins
  ),
  tibble(
    driver = "Within-site VSWC anomaly",
    driver_value = seq(
      quantile(monthly_sink_model_data$VSWC_anom, 0.02, na.rm = TRUE),
      quantile(monthly_sink_model_data$VSWC_anom, 0.98, na.rm = TRUE),
      length.out = 100
    ),
    Tair_anom = reference_values$Tair_anom,
    VSWC_anom = driver_value,
    n_halfhour_bins = reference_values$n_halfhour_bins
  ),
  tibble(
    driver = "Observed half-hour bins",
    driver_value = seq(
      quantile(monthly_sink_model_data$n_halfhour_bins, 0.02, na.rm = TRUE),
      quantile(monthly_sink_model_data$n_halfhour_bins, 0.98, na.rm = TRUE),
      length.out = 100
    ),
    Tair_anom = reference_values$Tair_anom,
    VSWC_anom = reference_values$VSWC_anom,
    n_halfhour_bins = driver_value
  )
) %>%
  expand_grid(season = levels(monthly_sink_model_data$season)) %>%
  mutate(SITE_ID = factor(reference_values$SITE_ID, levels = levels(monthly_sink_model_data$SITE_ID)))

effect_predictions <- predict(
  monthly_sink_model,
  newdata = effect_grid,
  type = "link",
  se.fit = TRUE,
  exclude = "s(SITE_ID)"
)

effect_grid <- effect_grid %>%
  mutate(
    fit_link = as.numeric(effect_predictions$fit),
    se_link = as.numeric(effect_predictions$se.fit),
    sink_probability = plogis(fit_link),
    sink_probability_lwr = plogis(fit_link - 1.96 * se_link),
    sink_probability_upr = plogis(fit_link + 1.96 * se_link)
  )

write.csv(monthly_sink_contrasts, "OUTPUT/NEON_monthly_sink_condition_contrasts.csv", row.names = FALSE)
write.csv(monthly_timing, "OUTPUT/NEON_monthly_sink_timing.csv", row.names = FALSE)
write.csv(monthly_site_sink_summary, "OUTPUT/NEON_monthly_sink_site_summary.csv", row.names = FALSE)
write.csv(effect_grid, "OUTPUT/NEON_monthly_sink_model_effects.csv", row.names = FALSE)

save(monthly_sink_model, file = "OUTPUT/NEON_monthly_sink_model.Rdata")

capture.output(
  {
    cat("Monthly sink-state model\n")
    print(summary(monthly_sink_model))
    cat("\nObserved monthly sink counts\n")
    print(monthly_observed %>% count(sink_month))
    cat("\nModel-standardized monthly sink counts\n")
    print(standardized_monthly %>% count(sink_month))
    cat("\nObserved contrasts\n")
    print(observed_contrasts)
    cat("\nModel-standardized contrasts\n")
    print(standardized_contrasts)
  },
  file = "OUTPUT/NEON_monthly_sink_model_summary.txt"
)

plot_sink_timing <- monthly_timing %>%
  ggplot(aes(x = month, y = sink_fraction, color = source, group = source)) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(size = n_sink_months), alpha = 0.85) +
  facet_wrap(~source, ncol = 1) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Monthly Sink Behavior Is Seasonally Structured",
    subtitle = "Sink months are site-months with negative balanced monthly total flux or negative model-standardized monthly budget.",
    x = NULL,
    y = "Fraction of site-months that are sinks",
    size = "Sink months",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

ggsave("FIGURES/NEON_monthly_sink_timing.png", plot_sink_timing, width = 9, height = 7, units = "in", dpi = 300)
ggsave("FIGURES/NEON_monthly_sink_timing.pdf", plot_sink_timing, width = 9, height = 7, units = "in")

condition_plot_data <- monthly_observed %>%
  dplyr::select(
    sink_state,
    `Within-site temperature anomaly` = Tair_anom,
    `Within-site VSWC anomaly` = VSWC_anom,
    `Gradient component` = mean_gradient_mgC_30min,
    `Storage component` = mean_storage_mgC_30min,
    `Observed half-hour bins` = n_halfhour_bins
  ) %>%
  pivot_longer(-sink_state, names_to = "condition", values_to = "value") %>%
  filter(is.finite(value))

plot_sink_conditions <- condition_plot_data %>%
  ggplot(aes(x = sink_state, y = value, color = sink_state)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.15, width = 0.55) +
  geom_jitter(width = 0.13, alpha = 0.45, size = 1.2) +
  facet_wrap(~condition, scales = "free_y", ncol = 3) +
  scale_color_manual(values = sink_source_colors, na.translate = FALSE) +
  labs(
    title = "Conditions Associated with Observed Monthly Sink States",
    x = NULL,
    y = "Monthly condition or flux component",
    color = NULL
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave("FIGURES/NEON_monthly_sink_conditions.png", plot_sink_conditions, width = 11, height = 7, units = "in", dpi = 300)
ggsave("FIGURES/NEON_monthly_sink_conditions.pdf", plot_sink_conditions, width = 11, height = 7, units = "in")

plot_sink_model_effects <- effect_grid %>%
  ggplot(aes(x = driver_value, y = sink_probability, color = season, fill = season)) +
  geom_ribbon(aes(ymin = sink_probability_lwr, ymax = sink_probability_upr), alpha = 0.13, color = NA) +
  geom_line(linewidth = 0.85) +
  facet_wrap(~driver, scales = "free_x", ncol = 3) +
  scale_color_manual(values = season_colors, na.translate = FALSE) +
  scale_fill_manual(values = season_colors, na.translate = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Monthly Sink-State Model Effects",
    subtitle = "Predictions exclude the site random effect and show expected sink probability for a reference site.",
    x = NULL,
    y = "Predicted monthly sink probability",
    color = "Season",
    fill = "Season"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave("FIGURES/NEON_monthly_sink_model_effects.png", plot_sink_model_effects, width = 12, height = 5.8, units = "in", dpi = 300)
ggsave("FIGURES/NEON_monthly_sink_model_effects.pdf", plot_sink_model_effects, width = 12, height = 5.8, units = "in")

plot_sink_site_summary <- monthly_site_sink_summary %>%
  mutate(SITE_ID = fct_reorder(SITE_ID, observed_sink_fraction)) %>%
  ggplot(aes(x = observed_sink_fraction, y = SITE_ID, color = CH4_behavior)) +
  geom_point(aes(size = n_months), alpha = 0.9) +
  geom_point(aes(x = standardized_sink_fraction), shape = 1, size = 2.4, color = "black", stroke = 0.7) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = behavior_colors, na.translate = FALSE) +
  labs(
    title = "Monthly Sink Frequency by Site",
    subtitle = "Filled points show observed balanced months; open circles show model-standardized months.",
    x = "Fraction of months classified as sinks",
    y = NULL,
    color = "Observed site class",
    size = "Observed months"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom", axis.text.y = element_text(size = 7))

ggsave("FIGURES/NEON_monthly_sink_frequency_by_site.png", plot_sink_site_summary, width = 9, height = 9, units = "in", dpi = 300)
ggsave("FIGURES/NEON_monthly_sink_frequency_by_site.pdf", plot_sink_site_summary, width = 9, height = 9, units = "in")

observed_sink_count <- sum(monthly_observed$sink_month, na.rm = TRUE)
observed_month_count <- nrow(monthly_observed)
standardized_sink_count <- sum(standardized_monthly$sink_month, na.rm = TRUE)
standardized_month_count <- nrow(standardized_monthly)

season_lines <- monthly_timing %>%
  group_by(source, season) %>%
  summarise(
    sink_fraction = mean(sink_fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(line = paste0("- ", source, " / ", season, ": ", scales::percent(sink_fraction, accuracy = 1))) %>%
  pull(line)

top_observed_contrast_lines <- observed_contrasts %>%
  arrange(p_value) %>%
  slice_head(n = 6) %>%
  mutate(
    line = paste0(
      "- ", label, ": sink median = ", signif(median_sink, 3),
      ", source median = ", signif(median_source, 3),
      ", difference = ", signif(median_difference_sink_minus_source, 3),
      ", p = ", signif(p_value, 3)
    )
  ) %>%
  pull(line)

top_site_lines <- monthly_site_sink_summary %>%
  arrange(desc(observed_sink_fraction)) %>%
  slice_head(n = 8) %>%
  mutate(
    line = paste0(
      "- ", SITE_ID, " (", CH4_behavior, "): observed sink fraction = ",
      scales::percent(observed_sink_fraction, accuracy = 1),
      ", standardized sink fraction = ",
      scales::percent(standardized_sink_fraction, accuracy = 1)
    )
  ) %>%
  pull(line)

writeLines(
  c(
    "# NEON Monthly Sink Behavior Analysis",
    "",
    "## Purpose",
    "This analysis asks what distinguishes monthly sink states from monthly source states. It is intentionally monthly rather than annual or site-class based, because sink behavior occurs in many site-months even though few sites are consistent sinks.",
    "",
    "## Sink Month Counts",
    paste0("- Observed balanced monthly data: ", observed_sink_count, " sink months out of ", observed_month_count, " site-months (", scales::percent(observed_sink_count / observed_month_count, accuracy = 0.1), ")."),
    paste0("- Model-standardized monthly budgets: ", standardized_sink_count, " sink months out of ", standardized_month_count, " site-months (", scales::percent(standardized_sink_count / standardized_month_count, accuracy = 0.1), ")."),
    "",
    "## Seasonal Sink Fractions",
    season_lines,
    "",
    "## Strongest Observed Monthly Contrasts",
    top_observed_contrast_lines,
    "",
    "## Sites with Highest Observed Monthly Sink Frequency",
    top_site_lines,
    "",
    "## Interpretation",
    "- Monthly sink behavior is common even though consistent sink sites are rare.",
    "- Sink months are mostly explained by the monthly gradient component becoming negative; storage is secondary.",
    "- Observed sink months are most frequent in summer, and within a site they tend to occur in warmer-than-average months.",
    "- Mean VSWC alone is not a strong separator of sink versus source months in the current monthly summaries.",
    "- The number of observed half-hour bins differs between sink and source months, so observed monthly sink inference should be interpreted alongside the model-standardized monthly budgets.",
    "- Standardized sink months remain concentrated in a subset of months and sites, supporting the interpretation that monthly sink behavior is real but often seasonal or near-threshold rather than site-permanent.",
    "",
    "## Outputs",
    "- `OUTPUT/NEON_monthly_sink_condition_contrasts.csv`",
    "- `OUTPUT/NEON_monthly_sink_timing.csv`",
    "- `OUTPUT/NEON_monthly_sink_site_summary.csv`",
    "- `OUTPUT/NEON_monthly_sink_model_effects.csv`",
    "- `OUTPUT/NEON_monthly_sink_model.Rdata`",
    "- `OUTPUT/NEON_monthly_sink_model_summary.txt`",
    "- `FIGURES/NEON_monthly_sink_timing.png`",
    "- `FIGURES/NEON_monthly_sink_conditions.png`",
    "- `FIGURES/NEON_monthly_sink_model_effects.png`",
    "- `FIGURES/NEON_monthly_sink_frequency_by_site.png`"
  ),
  "OUTPUT/NEON_monthly_sink_behavior_results.md"
)

message("Wrote monthly sink behavior analysis outputs.")
