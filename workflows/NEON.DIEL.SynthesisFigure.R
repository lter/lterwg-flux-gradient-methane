# Synthesis figure for NEON CH4 half-hour driver analysis.
# Combines observed flux patterns, source/sink site classes, model driver effects,
# and apparent Q10 summaries from NEON.DIEL.Analysis2.R.

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)



model_data_file <- "OUTPUT/30min_ch4_model_data.csv"
site_behavior_file <- "OUTPUT/30min_site_behavior.csv"
driver_prediction_file <- "OUTPUT/NEON_DIEL_Analysis2_driver_predictions.csv"
q10_summary_file <- "OUTPUT/NEON_DIEL_Analysis2_Q10_behavior_summary.csv"

required_files <- c(
  model_data_file,
  site_behavior_file,
  driver_prediction_file,
  q10_summary_file
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required analysis output(s): ",
    paste(missing_files, collapse = ", "),
    ". Run flow.30min.analysis.R and NEON.DIEL.Analysis2.R first."
  )
}

required_model_cols <- c(
  "SITE_ID", "CH4_mgC_30min", "source_30min", "hour_num"
)
required_behavior_cols <- c("SITE_ID", "CH4_behavior")
required_driver_prediction_cols <- c(
  "driver", "driver_value", "source_prob", "positive_flux", "sink_uptake"
)
required_q10_summary_cols <- c(
  "regime", "CH4_behavior", "median_Q10", "q25_Q10", "q75_Q10"
)

ch4_30min_raw <- read.csv(model_data_file)
site_behavior_raw <- read.csv(site_behavior_file)
driver_predictions_raw <- read.csv(driver_prediction_file)
q10_summary_raw <- read.csv(q10_summary_file)

check_required_cols <- function(data, required_cols, file_label) {
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", file_label, ": ",
      paste(missing_cols, collapse = ", "),
      ". Re-run flow.30min.analysis.R and NEON.DIEL.Analysis2.R."
    )
  }
}

check_required_cols(ch4_30min_raw, required_model_cols, model_data_file)
check_required_cols(site_behavior_raw, required_behavior_cols, site_behavior_file)
check_required_cols(driver_predictions_raw, required_driver_prediction_cols, driver_prediction_file)
check_required_cols(q10_summary_raw, required_q10_summary_cols, q10_summary_file)

behavior_levels <- c("Consistent sink", "Fluctuating", "Consistent source")
behavior_colors <- c(
  "Consistent sink" = "red3",
  "Fluctuating" = "grey35",
  "Consistent source" = "blue4"
)

driver_labels <- c(
  "log_PAR" = "Light (log PAR)",
  "Tair_C" = "Air temperature (C)",
  "VSWCMean" = "Soil moisture (VSWC)"
)

response_labels <- c(
  "source_prob" = "Chance of emitting",
  "positive_flux" = "Emission magnitude",
  "sink_uptake" = "Uptake magnitude"
)

site_behavior <- site_behavior_raw %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    CH4_behavior = factor(CH4_behavior, levels = behavior_levels)
  )

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
    CH4_behavior = factor(CH4_behavior, levels = behavior_levels),
    source_30min = as.logical(source_30min)
  ) %>%
  filter(
    is.finite(CH4_mgC_30min),
    is.finite(hour_num),
    !is.na(CH4_behavior)
  )

observed_hourly <- ch4_30min %>%
  mutate(hour_bin = floor(hour_num)) %>%
  reframe(
    .by = c(CH4_behavior, hour_bin),
    mean_CH4 = mean(CH4_mgC_30min, na.rm = TRUE),
    se_CH4 = sd(CH4_mgC_30min, na.rm = TRUE) / sqrt(dplyr::n()),
    source_probability = mean(source_30min, na.rm = TRUE),
    n = dplyr::n()
  )

behavior_counts <- site_behavior %>%
  filter(!is.na(CH4_behavior)) %>%
  count(CH4_behavior, name = "n_sites") %>%
  mutate(
    label = paste0(n_sites, " sites"),
    CH4_behavior = factor(CH4_behavior, levels = behavior_levels)
  )

driver_predictions <- driver_predictions_raw %>%
  mutate(
    driver = factor(driver, levels = names(driver_labels), labels = driver_labels),
    CH4_behavior = if ("CH4_behavior" %in% names(.)) {
      factor(CH4_behavior, levels = behavior_levels)
    } else {
      factor(NA_character_, levels = behavior_levels)
    }
  ) %>%
  filter(!is.na(driver)) %>%
  pivot_longer(
    cols = c(source_prob, positive_flux, sink_uptake),
    names_to = "response",
    values_to = "prediction"
  ) %>%
  mutate(
    response = factor(response, levels = names(response_labels), labels = response_labels)
  ) %>%
  filter(is.finite(driver_value), is.finite(prediction), !is.na(response))

q10_summary <- q10_summary_raw %>%
  mutate(
    CH4_behavior = factor(CH4_behavior, levels = behavior_levels),
    regime = factor(regime, levels = c("Positive emission", "Sink uptake magnitude"))
  ) %>%
  filter(
    !is.na(CH4_behavior),
    !is.na(regime),
    is.finite(median_Q10),
    is.finite(q25_Q10),
    is.finite(q75_Q10)
  )

plot_behavior <- behavior_counts %>%
  ggplot(aes(x = CH4_behavior, y = n_sites, fill = CH4_behavior)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(aes(label = label), vjust = -0.35, size = 3.5) +
  scale_fill_manual(values = behavior_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "A. Site behavior classes",
    x = NULL,
    y = "Number of sites"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 22, hjust = 1),
    plot.title = element_text(face = "bold")
  )

plot_observed_diel <- observed_hourly %>%
  ggplot(aes(x = hour_bin, y = mean_CH4, color = CH4_behavior, fill = CH4_behavior)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey30") +
  geom_ribbon(
    aes(ymin = mean_CH4 - 1.96 * se_CH4, ymax = mean_CH4 + 1.96 * se_CH4),
    alpha = 0.14,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = behavior_colors, name = "Site behavior") +
  scale_fill_manual(values = behavior_colors, guide = "none") +
  labs(
    title = "B. Observed 30-minute CH4 flux diel pattern",
    x = "Hour of day",
    y = expression("Mean CH"[4] * " flux (mg C 30 min"^-1 * ")")
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

plot_observed_source <- observed_hourly %>%
  ggplot(aes(x = hour_bin, y = source_probability, color = CH4_behavior)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = behavior_colors, name = "Site behavior") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "C. Observed probability of positive CH4 flux",
    x = "Hour of day",
    y = "Positive-flux observations"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

plot_driver_effects <- driver_predictions %>%
  ggplot(aes(x = driver_value, y = prediction)) +
  geom_line(linewidth = 0.95, color = "black") +
  facet_grid(response ~ driver, scales = "free", switch = "y") +
  labs(
    title = "D. Model-implied driver responses",
    x = "Driver gradient",
    y = "Prediction"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey94", color = "grey30"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

plot_q10 <- q10_summary %>%
  ggplot(aes(x = CH4_behavior, y = median_Q10, color = CH4_behavior)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_errorbar(aes(ymin = q25_Q10, ymax = q75_Q10), width = 0.12, linewidth = 0.75) +
  geom_point(size = 2.6) +
  facet_wrap(~regime, ncol = 1) +
  scale_color_manual(values = behavior_colors, guide = "none") +
  labs(
    title = "E. Apparent Q10 by behavior",
    x = NULL,
    y = "Median apparent Q10\n(IQR)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 22, hjust = 1),
    strip.background = element_rect(fill = "grey94", color = "grey30"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

synthesis_figure <- (
  (plot_behavior | plot_observed_diel | plot_observed_source) /
    (plot_driver_effects | plot_q10)
) +
  plot_layout(heights = c(0.95, 1.35), widths = c(1, 1.55, 1.15), guides = "collect") +
  plot_annotation(
    title = "NEON CH4 flux synthesis: site state, diel behavior, and environmental drivers",
    subtitle = paste(
      "Source/sink status is better explained as switching into positive-flux state plus separate emission and uptake magnitudes;",
      "apparent Q10 is modest and not a strong separator of site behavior."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10.5)
    )
  ) &
  theme(legend.position = "bottom")

ggsave(
  "FIGURES/NEON_DIEL_driver_flux_synthesis.png",
  plot = synthesis_figure,
  width = 15,
  height = 10,
  units = "in",
  dpi = 300
)

ggsave(
  "FIGURES/NEON_DIEL_driver_flux_synthesis.pdf",
  plot = synthesis_figure,
  width = 15,
  height = 10,
  units = "in"
)

message("Wrote FIGURES/NEON_DIEL_driver_flux_synthesis.png and .pdf")
