# Check whether total-flux gapfilling changes source/sink sign.

library(tidyverse)
library(ggplot2)

script_arg <- commandArgs(FALSE) %>%
  str_subset("^--file=") %>%
  str_remove("^--file=") %>%
  first()

script_dir <- if (is.na(script_arg)) {
  getwd()
} else {
  dirname(normalizePath(script_arg, mustWork = FALSE))
}

workflow_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
output_dir <- file.path(script_dir, "OUTPUT")
figure_dir <- file.path(script_dir, "FIGURES")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

annual_file <- file.path(workflow_dir, "OUTPUT", "NEON_total_flux_annual_budget.csv")
gapfilled_file <- file.path(workflow_dir, "OUTPUT", "NEON_total_flux_gapfilled_30min.csv.gz")

if (!file.exists(annual_file)) {
  stop("Missing annual budget file: ", annual_file)
}

if (!file.exists(gapfilled_file)) {
  stop("Missing gapfilled 30-minute file: ", gapfilled_file)
}

annual <- readr::read_csv(annual_file, show_col_types = FALSE)
gapfilled_30min <- readr::read_csv(gapfilled_file, show_col_types = FALSE)

required_annual <- c(
  "SITE_ID", "Year", "observed_coverage", "n_observed_days",
  "observed_scaled_annual_gC_m2_yr", "annual_budget_gC_m2_yr",
  "quality_flag", "interpretable_budget"
)

required_gapfilled <- c(
  "SITE_ID", "Year", "gapfilled", "observed_CH4_mgC_30min",
  "pred_CH4_mgC_30min"
)

missing_annual <- setdiff(required_annual, names(annual))
missing_gapfilled <- setdiff(required_gapfilled, names(gapfilled_30min))

if (length(missing_annual) > 0) {
  stop("Annual budget file is missing required columns: ", paste(missing_annual, collapse = ", "))
}

if (length(missing_gapfilled) > 0) {
  stop("Gapfilled 30-minute file is missing required columns: ", paste(missing_gapfilled, collapse = ", "))
}

interval_sign_summary <- gapfilled_30min %>%
  group_by(SITE_ID, Year) %>%
  summarise(
    n_obs = sum(!gapfilled & is.finite(observed_CH4_mgC_30min)),
    n_fill = sum(gapfilled & is.finite(pred_CH4_mgC_30min)),
    prop_obs_positive = mean(
      observed_CH4_mgC_30min[!gapfilled & is.finite(observed_CH4_mgC_30min)] > 0,
      na.rm = TRUE
    ),
    prop_fill_positive = mean(
      pred_CH4_mgC_30min[gapfilled & is.finite(pred_CH4_mgC_30min)] > 0,
      na.rm = TRUE
    ),
    mean_obs_mgC = mean(
      observed_CH4_mgC_30min[!gapfilled & is.finite(observed_CH4_mgC_30min)],
      na.rm = TRUE
    ),
    mean_fill_mgC = mean(
      pred_CH4_mgC_30min[gapfilled & is.finite(pred_CH4_mgC_30min)],
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    delta_prop_positive = prop_fill_positive - prop_obs_positive,
    delta_mean_mgC = mean_fill_mgC - mean_obs_mgC
  )

sign_bias_diagnostic <- annual %>%
  left_join(interval_sign_summary, by = c("SITE_ID", "Year")) %>%
  filter(
    is.finite(observed_scaled_annual_gC_m2_yr),
    is.finite(annual_budget_gC_m2_yr)
  ) %>%
  mutate(
    observed_sign = case_when(
      observed_scaled_annual_gC_m2_yr > 0 ~ "source",
      observed_scaled_annual_gC_m2_yr < 0 ~ "sink",
      TRUE ~ "zero"
    ),
    gapfilled_sign = case_when(
      annual_budget_gC_m2_yr > 0 ~ "source",
      annual_budget_gC_m2_yr < 0 ~ "sink",
      TRUE ~ "zero"
    ),
    sign_flip = observed_sign != gapfilled_sign,
    flip_direction = case_when(
      observed_sign == "sink" & gapfilled_sign == "source" ~ "observed sink -> gapfilled source",
      observed_sign == "source" & gapfilled_sign == "sink" ~ "observed source -> gapfilled sink",
      TRUE ~ "same sign"
    ),
    annual_budget_difference_gC_m2_yr = annual_budget_gC_m2_yr - observed_scaled_annual_gC_m2_yr
  )

sign_bias_summary <- sign_bias_diagnostic %>%
  summarise(
    .by = interpretable_budget,
    n_site_years = n(),
    n_observed_source = sum(observed_sign == "source"),
    n_gapfilled_source = sum(gapfilled_sign == "source"),
    n_sign_flips = sum(sign_flip),
    n_sink_to_source = sum(flip_direction == "observed sink -> gapfilled source"),
    n_source_to_sink = sum(flip_direction == "observed source -> gapfilled sink"),
    median_budget_difference_gC_m2_yr = median(annual_budget_difference_gC_m2_yr, na.rm = TRUE),
    median_observed_positive_interval_fraction = median(prop_obs_positive, na.rm = TRUE),
    median_gapfilled_positive_interval_fraction = median(prop_fill_positive, na.rm = TRUE),
    median_delta_positive_interval_fraction = median(delta_prop_positive, na.rm = TRUE),
    median_delta_mean_mgC_30min = median(delta_mean_mgC, na.rm = TRUE)
  ) %>%
  arrange(desc(interpretable_budget))

readr::write_csv(
  sign_bias_diagnostic,
  file.path(output_dir, "NEON_total_flux_gapfill_sign_bias_diagnostic.csv")
)

readr::write_csv(
  sign_bias_summary,
  file.path(output_dir, "NEON_total_flux_gapfill_sign_bias_summary.csv")
)

sign_flip_plot <- sign_bias_diagnostic %>%
  mutate(
    interpretable_class = if_else(interpretable_budget, "Interpretable", "Low coverage / extrapolated"),
    flip_direction = factor(
      flip_direction,
      levels = c(
        "same sign",
        "observed sink -> gapfilled source",
        "observed source -> gapfilled sink"
      )
    )
  ) %>%
  ggplot(aes(
    x = observed_scaled_annual_gC_m2_yr,
    y = annual_budget_gC_m2_yr,
    color = flip_direction,
    shape = interpretable_class
  )) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, color = "grey70") +
  geom_point(alpha = 0.8, size = 2.4) +
  scale_color_manual(
    values = c(
      "same sign" = "grey35",
      "observed sink -> gapfilled source" = "#D55E00",
      "observed source -> gapfilled sink" = "#0072B2"
    )
  ) +
  theme_bw(base_size = 11) +
  labs(
    x = "Observed-only scaled annual budget (g C m-2 yr-1)",
    y = "Gapfilled annual budget (g C m-2 yr-1)",
    color = "Sign comparison",
    shape = NULL,
    title = "Does gapfilling change annual CH4 source/sink sign?"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "NEON_total_flux_gapfill_sign_bias.png"),
  sign_flip_plot,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)

summary_lines <- sign_bias_summary %>%
  mutate(
    class = if_else(interpretable_budget, "Interpretable site-years", "All other observed site-years"),
    line = paste0(
      "- ", class, ": n = ", n_site_years,
      ", observed source = ", n_observed_source,
      ", gapfilled source = ", n_gapfilled_source,
      ", sign flips = ", n_sign_flips,
      " (sink->source ", n_sink_to_source,
      ", source->sink ", n_source_to_sink, ")",
      ", median budget shift = ", signif(median_budget_difference_gC_m2_yr, 3),
      " g C m-2 yr-1",
      ", median interval positive-fraction shift = ",
      signif(median_delta_positive_interval_fraction, 3)
    )
  ) %>%
  pull(line)

flip_lines <- sign_bias_diagnostic %>%
  filter(interpretable_budget, sign_flip) %>%
  arrange(SITE_ID, Year) %>%
  mutate(
    line = paste0(
      "- ", SITE_ID, " ", Year,
      ": observed = ", signif(observed_scaled_annual_gC_m2_yr, 3),
      ", gapfilled = ", signif(annual_budget_gC_m2_yr, 3),
      ", coverage = ", scales::percent(observed_coverage, accuracy = 0.1),
      ", ", flip_direction
    )
  ) %>%
  pull(line)

writeLines(
  c(
    "# NEON Total-Flux Gapfill Sign Bias Check",
    "",
    "## Interpretation",
    "Gapfilling slightly shifts the annual budgets in the source-positive direction. The shift is modest relative to the fact that most observed-only scaled annual budgets are already positive, but it can change sign in low-to-moderate coverage site-years.",
    "",
    "The sign conclusion should therefore be reported with coverage diagnostics: gapfilling is not the sole reason the budgets are source-like, but it does increase the number of source-classified site-years.",
    "",
    "## Summary",
    summary_lines,
    "",
    "## Interpretable Site-Years With Sign Flips",
    if (length(flip_lines) > 0) flip_lines else "- None",
    "",
    "## Outputs",
    "- `OUTPUT/NEON_total_flux_gapfill_sign_bias_diagnostic.csv`",
    "- `OUTPUT/NEON_total_flux_gapfill_sign_bias_summary.csv`",
    "- `FIGURES/NEON_total_flux_gapfill_sign_bias.png`"
  ),
  file.path(output_dir, "NEON_total_flux_gapfill_sign_bias.md")
)

message("Wrote NEON total-flux gapfill sign-bias diagnostics.")
