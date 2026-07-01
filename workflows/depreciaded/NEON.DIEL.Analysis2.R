# NEON CH4 DIEL analysis, part 2:
# Two-part 30-minute models plus apparent Q10 estimates for source and sink regimes.

library(tidyverse)
library(ggplot2)
library(mgcv)
library(patchwork)
library(ggrepel)

localdir.ch4 <- '/Volumes/MaloneLab/Research/FluxGradient/Methane'
load( fs::path(localdir.ch4 ,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata")) )

setwd( localdir.ch4 )
dir.create("OUTPUT", showWarnings = FALSE)
dir.create("FIGURES", showWarnings = FALSE)

model_data_file <- "OUTPUT/30min_ch4_model_data.csv"
site_behavior_file <- "OUTPUT/30min_site_behavior.csv"

if (!file.exists(model_data_file)) {
  stop("Missing OUTPUT/30min_ch4_model_data.csv. Run flow.30min.analysis.R first.")
}

if (!file.exists(site_behavior_file)) {
  stop("Missing OUTPUT/30min_site_behavior.csv. Run flow.30min.analysis.R first.")
}

min_q10_n <- 100
min_q10_temp_range <- 5
tref_c <- 10

site_behavior <- read.csv(site_behavior_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    CH4_behavior = factor(CH4_behavior, levels = c("Weak-sink", "Fluctuating", "Weak-source"))
  )

required_model_cols <- c(
  "time.rounded", "Date", "SITE_ID", "CH4_mgC_30min", "source_30min",
  "hour_num", "sin_hour", "cos_hour", "Tair_C", "VSWCMean", "log_PAR",
  "season", "EcoType"
)

ch4_30min_raw <- read.csv(model_data_file)
missing_model_cols <- setdiff(required_model_cols, names(ch4_30min_raw))

if (length(missing_model_cols) > 0) {
  stop(
    "Missing required columns in OUTPUT/30min_ch4_model_data.csv: ",
    paste(missing_model_cols, collapse = ", "),
    ". Re-run flow.30min.analysis.R and check its output schema."
  )
}

if (!"CH4_behavior" %in% names(ch4_30min_raw)) {
  ch4_30min_raw <- ch4_30min_raw %>%
    mutate(SITE_ID = as.character(SITE_ID)) %>%
    left_join(
      site_behavior %>% dplyr::select(SITE_ID, CH4_behavior),
      by = "SITE_ID"
    )
}

ch4_30min <- ch4_30min_raw %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    time.rounded = as.POSIXct(time.rounded),
    Date = as.Date(Date),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    SITE_ID = factor(SITE_ID),
    EcoType = factor(EcoType),
    CH4_behavior = factor(CH4_behavior, levels = c("Weak-sink", "Fluctuating", "Weak-source")),
    source_30min = as.logical(source_30min),
    positive_CH4 = CH4_mgC_30min > 0,
    sink_CH4 = CH4_mgC_30min < 0,
    log_positive_CH4 = log(if_else(positive_CH4, CH4_mgC_30min, NA_real_)),
    log_sink_uptake = log(if_else(sink_CH4, abs(CH4_mgC_30min), NA_real_))
  ) %>%
  filter(
    is.finite(CH4_mgC_30min),
    is.finite(Tair_C),
    is.finite(VSWCMean),
    is.finite(log_PAR),
    !is.na(season),
    !is.na(EcoType),
    !is.na(CH4_behavior)
  )

# Two-part models:
# 1. source-state probability
# 2. positive CH4 emission magnitude
# 3. negative CH4 uptake magnitude
source_state_model_data <- ch4_30min %>% droplevels()
positive_flux_model_data <- ch4_30min %>%
  filter(positive_CH4) %>%
  droplevels()
sink_uptake_model_data <- ch4_30min %>%
  filter(sink_CH4) %>%
  droplevels()

source_state_model <- mgcv::bam(
  source_30min ~
    s(hour_num, k = 12) +
    sin_hour +
    cos_hour +
    s(Tair_C, k = 8) +
    s(VSWCMean, k = 8) +
    s(log_PAR, k = 8) +
    ti(Tair_C, VSWCMean, k = c(5, 5)) +
    season +
    EcoType +
    CH4_behavior +
    s(SITE_ID, bs = "re"),
  data = source_state_model_data,
  family = binomial(),
  method = "fREML",
  discrete = FALSE
)

positive_flux_model <- mgcv::bam(
  log_positive_CH4 ~
    s(hour_num, k = 12) +
    sin_hour +
    cos_hour +
    s(Tair_C, k = 8) +
    s(VSWCMean, k = 8) +
    s(log_PAR, k = 8) +
    ti(Tair_C, VSWCMean, k = c(5, 5)) +
    season +
    EcoType +
    CH4_behavior +
    s(SITE_ID, bs = "re"),
  data = positive_flux_model_data,
  family = gaussian(),
  method = "fREML",
  discrete = FALSE
)

sink_uptake_model <- mgcv::bam(
  log_sink_uptake ~
    s(hour_num, k = 12) +
    sin_hour +
    cos_hour +
    s(Tair_C, k = 8) +
    s(VSWCMean, k = 8) +
    s(log_PAR, k = 8) +
    ti(Tair_C, VSWCMean, k = c(5, 5)) +
    season +
    EcoType +
    CH4_behavior +
    s(SITE_ID, bs = "re"),
  data = sink_uptake_model_data,
  family = gaussian(),
  method = "fREML",
  discrete = FALSE
)

save(
  source_state_model,
  positive_flux_model,
  sink_uptake_model,
  file = "OUTPUT/NEON_DIEL_Analysis2_two_part_models.Rdata"
)

capture.output(
  {
    cat("Source-state model\n")
    print(summary(source_state_model))
    cat("\nPositive CH4 emission magnitude model\n")
    print(summary(positive_flux_model))
    cat("\nSink uptake magnitude model\n")
    print(summary(sink_uptake_model))
  },
  file = "OUTPUT/NEON_DIEL_Analysis2_model_summaries.txt"
)

# Apparent Q10 from log-linear temperature response within site-season.
fit_apparent_q10 <- function(df, response_col, regime_label) {
  df %>%
    group_by(SITE_ID, season, CH4_behavior, EcoType) %>%
    group_modify(function(group_df, key) {
      response <- group_df[[response_col]]
      model_df <- group_df %>%
        mutate(response = response) %>%
        filter(is.finite(response), is.finite(Tair_C))
      
      n_obs <- nrow(model_df)
      temp_range <- diff(range(model_df$Tair_C, na.rm = TRUE))
      
      if (n_obs < min_q10_n || !is.finite(temp_range) || temp_range < min_q10_temp_range) {
        return(tibble(
          regime = regime_label,
          n_obs = n_obs,
          temp_range = temp_range,
          slope = NA_real_,
          slope_p = NA_real_,
          r2 = NA_real_,
          apparent_Q10 = NA_real_,
          Rref_mgC_30min = NA_real_
        ))
      }
      
      fit <- lm(response ~ Tair_C, data = model_df)
      coef_table <- summary(fit)$coefficients
      slope <- coef_table["Tair_C", "Estimate"]
      intercept <- coef_table["(Intercept)", "Estimate"]
      
      tibble(
        regime = regime_label,
        n_obs = n_obs,
        temp_range = temp_range,
        slope = slope,
        slope_p = coef_table["Tair_C", "Pr(>|t|)"],
        r2 = summary(fit)$r.squared,
        apparent_Q10 = exp(10 * slope),
        Rref_mgC_30min = exp(intercept + slope * tref_c)
      )
    }) %>%
    ungroup()
}

source_q10 <- fit_apparent_q10(
  ch4_30min %>% filter(positive_CH4),
  "log_positive_CH4",
  "Positive emission"
)

sink_q10 <- fit_apparent_q10(
  ch4_30min %>% filter(sink_CH4),
  "log_sink_uptake",
  "Sink uptake magnitude"
)

apparent_q10 <- bind_rows(source_q10, sink_q10) %>%
  mutate(
    apparent_Q10_flag = case_when(
      is.na(apparent_Q10) ~ "not estimated",
      apparent_Q10 < 0.5 | apparent_Q10 > 10 ~ "extreme",
      TRUE ~ "usable"
    )
  )

write.csv(apparent_q10, "OUTPUT/NEON_DIEL_Analysis2_apparent_Q10_site_season.csv", row.names = FALSE)

q10_behavior_summary <- apparent_q10 %>%
  filter(apparent_Q10_flag == "usable") %>%
  reframe(
    .by = c(regime, CH4_behavior),
    n = n(),
    median_Q10 = median(apparent_Q10, na.rm = TRUE),
    q25_Q10 = quantile(apparent_Q10, 0.25, na.rm = TRUE),
    q75_Q10 = quantile(apparent_Q10, 0.75, na.rm = TRUE),
    median_Rref = median(Rref_mgC_30min, na.rm = TRUE)
  )

q10_tests <- apparent_q10 %>%
  filter(apparent_Q10_flag == "usable") %>%
  group_by(regime) %>%
  group_modify(function(df, key) {
    if (n_distinct(df$CH4_behavior) < 2 || nrow(df) < 8) {
      return(tibble(p.value = NA_real_))
    }
    tibble(p.value = kruskal.test(apparent_Q10 ~ CH4_behavior, data = df)$p.value)
  }) %>%
  ungroup()

write.csv(q10_behavior_summary, "OUTPUT/NEON_DIEL_Analysis2_Q10_behavior_summary.csv", row.names = FALSE)
write.csv(q10_tests, "OUTPUT/NEON_DIEL_Analysis2_Q10_behavior_tests.csv", row.names = FALSE)

# Prediction grids for two-part model effects.
reference_values <- ch4_30min %>%
  summarise(
    Tair_C = median(Tair_C, na.rm = TRUE),
    VSWCMean = median(VSWCMean, na.rm = TRUE),
    log_PAR = median(log_PAR, na.rm = TRUE),
    season = names(sort(table(season), decreasing = TRUE))[1],
    EcoType = names(sort(table(EcoType), decreasing = TRUE))[1],
    CH4_behavior = names(sort(table(CH4_behavior), decreasing = TRUE))[1],
    SITE_ID = names(sort(table(SITE_ID), decreasing = TRUE))[1]
  )

hour_grid <- expand_grid(
  hour_num = seq(0, 23.5, by = 0.5),
  CH4_behavior = levels(ch4_30min$CH4_behavior)
) %>%
  mutate(
    sin_hour = sin(2 * pi * hour_num / 24),
    cos_hour = cos(2 * pi * hour_num / 24),
    Tair_C = reference_values$Tair_C,
    VSWCMean = reference_values$VSWCMean,
    log_PAR = reference_values$log_PAR,
    season = reference_values$season,
    EcoType = reference_values$EcoType,
    SITE_ID = reference_values$SITE_ID
  )

driver_grid <- bind_rows(
  tibble(
    driver = "Tair_C",
    driver_value = seq(quantile(ch4_30min$Tair_C, 0.02, na.rm = TRUE),
                       quantile(ch4_30min$Tair_C, 0.98, na.rm = TRUE),
                       length.out = 100),
    Tair_C = driver_value,
    VSWCMean = reference_values$VSWCMean,
    log_PAR = reference_values$log_PAR
  ),
  tibble(
    driver = "VSWCMean",
    driver_value = seq(quantile(ch4_30min$VSWCMean, 0.02, na.rm = TRUE),
                       quantile(ch4_30min$VSWCMean, 0.98, na.rm = TRUE),
                       length.out = 100),
    Tair_C = reference_values$Tair_C,
    VSWCMean = driver_value,
    log_PAR = reference_values$log_PAR
  ),
  tibble(
    driver = "log_PAR",
    driver_value = seq(quantile(ch4_30min$log_PAR, 0.02, na.rm = TRUE),
                       quantile(ch4_30min$log_PAR, 0.98, na.rm = TRUE),
                       length.out = 100),
    Tair_C = reference_values$Tair_C,
    VSWCMean = reference_values$VSWCMean,
    log_PAR = driver_value
  )
) %>%
  mutate(
    hour_num = 12,
    sin_hour = sin(2 * pi * hour_num / 24),
    cos_hour = cos(2 * pi * hour_num / 24),
    season = reference_values$season,
    EcoType = reference_values$EcoType,
    CH4_behavior = reference_values$CH4_behavior,
    SITE_ID = reference_values$SITE_ID
  )

get_factor_levels <- function(model_data) {
  list(
    season = levels(model_data$season),
    EcoType = levels(model_data$EcoType),
    CH4_behavior = levels(model_data$CH4_behavior),
    SITE_ID = levels(model_data$SITE_ID)
  )
}

predict_supported <- function(model, newdata, factor_levels, transform = identity) {
  supported <- as.character(newdata$season) %in% factor_levels$season &
    as.character(newdata$EcoType) %in% factor_levels$EcoType &
    as.character(newdata$CH4_behavior) %in% factor_levels$CH4_behavior

  prediction <- rep(NA_real_, nrow(newdata))

  if (any(supported)) {
    prediction_data <- newdata[supported, , drop = FALSE] %>%
      mutate(
        season = factor(as.character(season), levels = factor_levels$season),
        EcoType = factor(as.character(EcoType), levels = factor_levels$EcoType),
        CH4_behavior = factor(as.character(CH4_behavior), levels = factor_levels$CH4_behavior),
        SITE_ID = factor(
          if_else(
            as.character(SITE_ID) %in% factor_levels$SITE_ID,
            as.character(SITE_ID),
            factor_levels$SITE_ID[1]
          ),
          levels = factor_levels$SITE_ID
        )
      )

    prediction[supported] <- transform(predict(
      model,
      newdata = prediction_data,
      type = "response",
      exclude = "s(SITE_ID)"
    ))
  }

  prediction
}

source_factor_levels <- get_factor_levels(source_state_model_data)
positive_factor_levels <- get_factor_levels(positive_flux_model_data)
sink_factor_levels <- get_factor_levels(sink_uptake_model_data)

hour_grid$source_prob <- predict_supported(source_state_model, hour_grid, source_factor_levels)
hour_grid$positive_flux <- predict_supported(positive_flux_model, hour_grid, positive_factor_levels, exp)
hour_grid$sink_uptake <- predict_supported(sink_uptake_model, hour_grid, sink_factor_levels, exp)

driver_grid$source_prob <- predict_supported(source_state_model, driver_grid, source_factor_levels)
driver_grid$positive_flux <- predict_supported(positive_flux_model, driver_grid, positive_factor_levels, exp)
driver_grid$sink_uptake <- predict_supported(sink_uptake_model, driver_grid, sink_factor_levels, exp)

write.csv(hour_grid, "OUTPUT/NEON_DIEL_Analysis2_hourly_predictions.csv", row.names = FALSE)
write.csv(driver_grid, "OUTPUT/NEON_DIEL_Analysis2_driver_predictions.csv", row.names = FALSE)

# Figures
plot_q10_behavior <- apparent_q10 %>%
  filter(apparent_Q10_flag == "usable") %>%
  ggplot(aes(x = CH4_behavior, y = apparent_Q10, color = CH4_behavior)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35) +
  geom_jitter(width = 0.15, alpha = 0.65, size = 1.8) +
  facet_wrap(~regime + season, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = c("Weak-sink" = "red3", "Fluctuating" = "grey35", "Weak-source" = "blue4"), na.translate = FALSE) +
  labs(x = "", y = "Apparent Q10", color = "Site behavior") +
  theme(
    strip.background = element_rect(fill = "transparent", color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave("FIGURES/NEON_DIEL_Analysis2_Q10_by_behavior.png", plot = plot_q10_behavior, width = 11, height = 7, units = "in")

plot_q10_vswc <- apparent_q10 %>%
  filter(apparent_Q10_flag == "usable") %>%
  left_join(
    ch4_30min %>%
      reframe(.by = c(SITE_ID, season), VSWCMean_site_season = mean(VSWCMean, na.rm = TRUE)),
    by = c("SITE_ID", "season")
  ) %>%
  ggplot(aes(x = VSWCMean_site_season, y = apparent_Q10, color = CH4_behavior)) +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~regime, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = c("Weak-sink" = "red3", "Fluctuating" = "grey35", "Weak-source" = "blue4"), na.translate = FALSE) +
  labs(x = "Mean site-season VSWC", y = "Apparent Q10", color = "Site behavior")

ggsave("FIGURES/NEON_DIEL_Analysis2_Q10_vs_VSWC.png", plot = plot_q10_vswc, width = 8, height = 5, units = "in")

plot_two_part_hour <- hour_grid %>%
  pivot_longer(cols = c(source_prob, positive_flux, sink_uptake), names_to = "response", values_to = "prediction") %>%
  mutate(
    response = recode(
      response,
      source_prob = "Source probability",
      positive_flux = "Positive CH4 magnitude",
      sink_uptake = "Sink uptake magnitude"
    )
  ) %>%
  ggplot(aes(x = hour_num, y = prediction, color = CH4_behavior)) +
  geom_line(linewidth = 1) +
  facet_wrap(~response, scales = "free_y", ncol = 1) +
  theme_bw() +
  scale_color_manual(values = c("Weak-sink" = "red3", "Fluctuating" = "grey35", "Weak-source" = "blue4"), na.translate = FALSE) +
  labs(x = "Hour", y = "Model prediction", color = "Site behavior")

ggsave("FIGURES/NEON_DIEL_Analysis2_two_part_hourly_predictions.png", plot = plot_two_part_hour, width = 8, height = 8.5, units = "in")

plot_two_part_drivers <- driver_grid %>%
  pivot_longer(cols = c(source_prob, positive_flux, sink_uptake), names_to = "response", values_to = "prediction") %>%
  mutate(
    response = recode(
      response,
      source_prob = "Source probability",
      positive_flux = "Positive CH4 magnitude",
      sink_uptake = "Sink uptake magnitude"
    )
  ) %>%
  ggplot(aes(x = driver_value, y = prediction)) +
  geom_line(linewidth = 1, color = "black") +
  facet_grid(response ~ driver, scales = "free") +
  theme_bw() +
  labs(x = "Driver value", y = "Model prediction") +
  theme(strip.background = element_rect(fill = "transparent", color = "black"))

ggsave("FIGURES/NEON_DIEL_Analysis2_two_part_driver_effects.png", plot = plot_two_part_drivers, width = 10, height = 7, units = "in")

plot_q10_map_data <- apparent_q10 %>%
  filter(regime == "Positive emission", apparent_Q10_flag == "usable") %>%
  reframe(
    .by = SITE_ID,
    median_source_Q10 = median(apparent_Q10, na.rm = TRUE),
    n_q10 = n()
  ) %>%
  left_join(site_behavior %>% dplyr::select(SITE_ID, CH4_behavior), by = "SITE_ID")

plot_q10_map <- site_behavior %>%
  dplyr::select(SITE_ID, Latitude..degrees., Longitude..degrees.) %>%
  left_join(plot_q10_map_data, by = "SITE_ID") %>%
  filter(!is.na(median_source_Q10)) %>%
  ggplot(aes(x = Longitude..degrees., y = Latitude..degrees.)) +
  geom_point(aes(color = median_source_Q10, shape = CH4_behavior, size = n_q10), alpha = 0.85) +
  theme_bw() +
  scale_color_viridis_c(option = "plasma") +
  labs(x = "Longitude", y = "Latitude", color = "Median source Q10", shape = "Site behavior", size = "Q10 windows")

ggsave("FIGURES/NEON_DIEL_Analysis2_source_Q10_map.png", plot = plot_q10_map, width = 7, height = 5, units = "in")

q10_counts <- apparent_q10 %>%
  count(regime, apparent_Q10_flag)

model_lines <- c(
  paste0("- Source-state deviance explained: ", signif(summary(source_state_model)$dev.expl, 3)),
  paste0("- Positive-emission magnitude adjusted R2: ", signif(summary(positive_flux_model)$r.sq, 3),
         "; deviance explained: ", signif(summary(positive_flux_model)$dev.expl, 3)),
  paste0("- Sink-uptake magnitude adjusted R2: ", signif(summary(sink_uptake_model)$r.sq, 3),
         "; deviance explained: ", signif(summary(sink_uptake_model)$dev.expl, 3))
)

q10_count_lines <- q10_counts %>%
  mutate(line = paste0("- ", regime, " / ", apparent_Q10_flag, ": ", n, " site-season windows")) %>%
  pull(line)

q10_summary_lines <- q10_behavior_summary %>%
  mutate(
    line = paste0(
      "- ", regime, " / ", CH4_behavior, ": median Q10 = ", signif(median_Q10, 3),
      " [", signif(q25_Q10, 3), ", ", signif(q75_Q10, 3), "], n = ", n
    )
  ) %>%
  pull(line)

q10_test_lines <- q10_tests %>%
  mutate(line = paste0("- ", regime, ": Kruskal-Wallis p = ", signif(p.value, 3))) %>%
  pull(line)

writeLines(
  c(
    "# NEON DIEL Analysis 2 Results",
    "",
    "This analysis separates three pieces of CH4 behavior: source-state probability, positive CH4 emission magnitude, and sink uptake magnitude.",
    "",
    "## Two-Part Model Summary",
    model_lines,
    "",
    "## Apparent Q10 Windows",
    paste0("Q10 was estimated only when a site-season regime had at least ", min_q10_n,
           " observations and at least ", min_q10_temp_range, " C of temperature range."),
    q10_count_lines,
    "",
    "## Apparent Q10 by Site Behavior",
    q10_summary_lines,
    "",
    "## Apparent Q10 Behavior Tests",
    q10_test_lines,
    "",
    "## Interpretation Notes",
    "- Positive-emission Q10 is the most defensible Q10 estimate here. Sink uptake Q10 is reported separately because uptake is a different process.",
    "- These are apparent Q10 values from log-linear temperature slopes, not causal temperature sensitivities. Moisture, season, and site context still matter.",
    "- Use the two-part model figures to interpret switching into source mode separately from emission/uptake magnitude.",
    "",
    "## Outputs",
    "- `OUTPUT/NEON_DIEL_Analysis2_model_summaries.txt`",
    "- `OUTPUT/NEON_DIEL_Analysis2_apparent_Q10_site_season.csv`",
    "- `OUTPUT/NEON_DIEL_Analysis2_Q10_behavior_summary.csv`",
    "- `FIGURES/NEON_DIEL_Analysis2_Q10_by_behavior.png`",
    "- `FIGURES/NEON_DIEL_Analysis2_Q10_vs_VSWC.png`",
    "- `FIGURES/NEON_DIEL_Analysis2_two_part_hourly_predictions.png`",
    "- `FIGURES/NEON_DIEL_Analysis2_two_part_driver_effects.png`",
    "- `FIGURES/NEON_DIEL_Analysis2_source_Q10_map.png`"
  ),
  "OUTPUT/NEON_DIEL_Analysis2_results.md"
)

