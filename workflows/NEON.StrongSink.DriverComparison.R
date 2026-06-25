# Compare NEON CH4 total-flux behavior classes across moisture, soil,
# climate, and vegetation/canopy drivers. Behavior classes are read from
# flow.30min.analysis.R outputs, where site-month source/sink state is based on
# equal-weighted site-month-hour flux bins.
#
# Historical versions focused on "strong sink" sites. After switching the
# workflow to total flux, there may be no Weak-sink sites. This script
# therefore keeps the old sink comparison when possible, and otherwise compares
# the behavior classes present in OUTPUT/30min_site_behavior.csv.


# This script compares the site-level environmental drivers of the CH4 behavior classes from flow.30min.analysis.R.
# 
# In plain English: it asks, “What site attributes differ between the sites classified as Weak-source and the sites classified as Fluctuating?”
# 
# There are currently no Weak-sink sites, so despite the script name, it is not comparing strong sinks right now. It falls back to:
#   
#   Weak-source vs Fluctuating
# What it does:
#   
#   Loads corrected outputs from the 30-minute analysis:
#   
#   OUTPUT/30min_ch4_model_data.csv
# OUTPUT/30min_site_behavior.csv
# Uses the corrected behavior classes:
#   
#   These classes are based on sampling-bias-adjusted monthly mean CH4 flux.
# The adjustment first averages by SITE_ID + YearMon + hour_num.
# Builds sampling-bias-adjusted moisture summaries:
#   
#   averages VSWCMean and VSWCVar within SITE_ID + YearMon + hour_num;
# then summarizes to monthly and site-level moisture metrics.
# Joins those moisture metrics to the site behavior table, which already contains:
#   
#   soil chemistry/depth variables;
# soil texture/physical variables;
# climate variables like MAP, MAT;
# canopy variables like LAI.mean, canopyHeight_m, CHM.mean.
# Reshapes all drivers into long format and compares each driver between behavior groups using:
#   
#   Wilcoxon rank-sum test;
# Cliff’s delta effect size;
# BH-adjusted p-values.
# Writes comparison tables and figures:
#   
#   OUTPUT/NEON_strong_sink_driver_comparison.csv
# OUTPUT/NEON_site_driver_values_for_sink_comparison.csv
# FIGURES/NEON_strong_sink_driver_comparison.png
# FIGURES/NEON_total_flux_driver_distributions_all.png
# The current strongest contrasts are that Weak-source sites tend to have lower dry mass, lower pH/acidity values, higher MAT, lower C:N, and lower estimated organic carbon than Fluctuating sites. These are exploratory contrasts, not causal tests.

library(tidyverse)
library(ggplot2)
library(patchwork)

localdir.ch4 <- "/Volumes/MaloneLab/Research/FluxGradient/Methane"
if (!dir.exists(localdir.ch4)) {
  stop("Missing methane output directory: ", localdir.ch4)
}
setwd(localdir.ch4)

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

required_model_cols <- c("SITE_ID", "YearMon", "hour_num", "VSWCMean", "VSWCVar")
required_behavior_cols <- c("SITE_ID", "CH4_behavior")

ch4_30min_raw <- read.csv(model_data_file)
missing_model_cols <- setdiff(required_model_cols, names(ch4_30min_raw))

if (length(missing_model_cols) > 0) {
  stop(
    "Missing required columns in OUTPUT/30min_ch4_model_data.csv: ",
    paste(missing_model_cols, collapse = ", "),
    ". Re-run flow.30min.analysis.R and check its output schema."
  )
}

site_behavior_raw <- read.csv(site_behavior_file)
missing_behavior_cols <- setdiff(required_behavior_cols, names(site_behavior_raw))

if (length(missing_behavior_cols) > 0) {
  stop(
    "Missing required columns in OUTPUT/30min_site_behavior.csv: ",
    paste(missing_behavior_cols, collapse = ", "),
    ". Re-run flow.30min.analysis.R and check its output schema."
  )
}

behavior_levels <- c("Weak-sink", "Fluctuating", "Weak-source")
# Color convention: blue = sink (uptake), grey = fluctuating, red = source (emission)
behavior_colors <- c(
  "Weak-sink"   = "#2166AC",
  "Fluctuating"       = "#4D4D4D",
  "Weak-source" = "#B2182B"
)
driver_group_colors <- c(
  "Moisture" = "#1B9E77",
  "Soil texture/physical" = "#D95F02",
  "Soil chemistry/depth" = "#7570B3",
  "Climate" = "#E7298A",
  "Vegetation/canopy" = "#66A61E"
)

site_behavior <- site_behavior_raw %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    CH4_behavior = factor(CH4_behavior, levels = behavior_levels)
  )

present_behaviors <- site_behavior %>%
  filter(!is.na(CH4_behavior)) %>%
  count(CH4_behavior, name = "n_sites") %>%
  filter(n_sites > 0) %>%
  arrange(CH4_behavior)

if (nrow(present_behaviors) < 2) {
  stop("Need at least two CH4 behavior classes to compare drivers.")
}

if ("Weak-sink" %in% as.character(present_behaviors$CH4_behavior)) {
  comparison_mode <- "strong_sink_vs_other"
  focus_behavior <- "Weak-sink"
  reference_behavior <- "Other behavior classes"
  comparison_title <- "Weak-sink vs other total-flux behavior classes"
} else {
  comparison_mode <- "available_behavior_contrast"
  focus_behavior <- as.character(present_behaviors$CH4_behavior[nrow(present_behaviors)])
  reference_behavior <- as.character(present_behaviors$CH4_behavior[1])
  comparison_title <- paste(focus_behavior, "vs", reference_behavior)
}

comparison_levels <- c(reference_behavior, focus_behavior)
comparison_colors <- behavior_colors[comparison_levels]
missing_comparison_colors <- is.na(comparison_colors)
comparison_colors[missing_comparison_colors] <- "#4D4D4D"

site_behavior <- site_behavior %>%
  mutate(
    behavior_comparison = case_when(
      comparison_mode == "strong_sink_vs_other" & CH4_behavior == focus_behavior ~ focus_behavior,
      comparison_mode == "strong_sink_vs_other" & !is.na(CH4_behavior) ~ reference_behavior,
      as.character(CH4_behavior) %in% comparison_levels ~ as.character(CH4_behavior),
      TRUE ~ NA_character_
    ),
    behavior_comparison = factor(behavior_comparison, levels = comparison_levels)
  )

ch4_30min <- ch4_30min_raw %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    YearMon = as.character(YearMon)
  ) %>%
  filter(!is.na(SITE_ID))

# Correct moisture summaries for sampling-time bias with the same structure used
# for flux behavior: first average repeated records within each
# site-month-hour bin, then average bins equally within months and sites.
moisture_site_month_hour <- ch4_30min %>%
  filter(is.finite(VSWCMean), is.finite(VSWCVar)) %>%
  reframe(
    .by = c(SITE_ID, YearMon, hour_num),
    VSWCMean_bin = mean(VSWCMean, na.rm = TRUE),
    VSWCVar_bin = mean(VSWCVar, na.rm = TRUE),
    n_30min_bin = dplyr::n()
  )

moisture_monthly <- moisture_site_month_hour %>%
  reframe(
    .by = c(SITE_ID, YearMon),
    VSWC_monthly_mean = mean(VSWCMean_bin, na.rm = TRUE),
    VSWC_monthly_var = mean(VSWCVar_bin, na.rm = TRUE),
    n_hour_bins_with_vswc = n_distinct(hour_num),
    n_30min_with_vswc = sum(n_30min_bin)
  )

moisture_site <- moisture_monthly %>%
  reframe(
    .by = SITE_ID,
    VSWCMean_site = mean(VSWC_monthly_mean, na.rm = TRUE),
    VSWCVar_site = mean(VSWC_monthly_var, na.rm = TRUE),
    VSWCMean_obs_sd = sd(VSWC_monthly_mean, na.rm = TRUE),
    VSWC_monthly_sd = sd(VSWC_monthly_mean, na.rm = TRUE),
    VSWC_monthly_range = diff(range(VSWC_monthly_mean, na.rm = TRUE)),
    n_months_with_vswc = dplyr::n(),
    n_hour_bins_with_vswc = sum(n_hour_bins_with_vswc, na.rm = TRUE),
    n_30min = sum(n_30min_with_vswc, na.rm = TRUE)
  )

soil_driver_variables <- site_behavior %>%
  dplyr::select(where(is.numeric)) %>%
  names() %>%
  intersect(c(
    "sulfurTot", "biogeoTopDepth", "biogeoBottomDepth", "carbonTot",
    "nitrogenTot", "ctonRatio", "estimatedOC", "acidity",
    "bulkDensOvenDry", "sandTotal", "siltTotal", "clayTotal", "dryMass"
  ))

driver_groups <- tibble(
  variable = c(
    "VSWCMean_site", "VSWCVar_site", "VSWC_monthly_sd", "VSWC_monthly_range",
    soil_driver_variables,
    "MAP", "MAT", "LAI.mean", "canopyHeight_m", "CHM.mean"
  )
) %>%
  mutate(
    driver_group = case_when(
      variable %in% c("VSWCMean_site", "VSWCVar_site", "VSWC_monthly_sd", "VSWC_monthly_range") ~ "Moisture",
      variable %in% c("sandTotal", "siltTotal", "clayTotal", "bulkDensOvenDry", "dryMass") ~ "Soil texture/physical",
      variable %in% soil_driver_variables ~ "Soil chemistry/depth",
      variable %in% c("MAP", "MAT") ~ "Climate",
      TRUE ~ "Vegetation/canopy"
    ),
    label = recode(
      variable,
      VSWCMean_site = "Mean VSWC",
      VSWCVar_site = "Mean VSWC variance",
      VSWC_monthly_sd = "Wetness seasonality\n(SD monthly VSWC)",
      VSWC_monthly_range = "Wetness seasonality\n(range monthly VSWC)",
      sulfurTot = "Sulfur total",
      biogeoTopDepth = "Biogeochem top depth",
      biogeoBottomDepth = "Biogeochem bottom depth",
      carbonTot = "Total C",
      nitrogenTot = "Total N",
      ctonRatio = "C:N",
      estimatedOC = "Estimated OC",
      acidity = "pH/acidity",
      bulkDensOvenDry = "Bulk density",
      sandTotal = "Sand",
      siltTotal = "Silt",
      clayTotal = "Clay",
      dryMass = "Dry mass",
      MAP = "MAP",
      MAT = "MAT",
      LAI.mean = "LAI",
      canopyHeight_m = "Canopy height",
      CHM.mean = "CHM",
      .default = variable
    )
  ) %>%
  distinct(variable, .keep_all = TRUE)

site_drivers <- site_behavior %>%
  left_join(moisture_site, by = "SITE_ID")

available_driver_groups <- driver_groups %>%
  filter(variable %in% names(site_drivers))

if (nrow(available_driver_groups) == 0) {
  stop("No requested driver columns are available after joining the flow.30min.analysis.R outputs.")
}

driver_long <- site_drivers %>%
  filter(!is.na(behavior_comparison)) %>%
  dplyr::select(SITE_ID, CH4_behavior, behavior_comparison, all_of(available_driver_groups$variable)) %>%
  pivot_longer(
    cols = all_of(available_driver_groups$variable),
    names_to = "variable",
    values_to = "value"
  ) %>%
  left_join(available_driver_groups, by = "variable") %>%
  filter(is.finite(value)) %>%
  group_by(variable) %>%
  mutate(value_z = as.numeric(scale(value))) %>%
  ungroup()

compare_one_variable <- function(df) {
  focus <- df$value[df$behavior_comparison == focus_behavior]
  reference <- df$value[df$behavior_comparison == reference_behavior]
  n_focus <- length(focus)
  n_reference <- length(reference)

  if (n_focus < 2 || n_reference < 2) {
    return(tibble(
      n_focus = n_focus,
      n_reference = n_reference,
      median_focus = median(focus, na.rm = TRUE),
      median_reference = median(reference, na.rm = TRUE),
      median_difference = NA_real_,
      cliffs_delta = NA_real_,
      p_value = NA_real_
    ))
  }

  ranks <- rank(c(focus, reference), ties.method = "average")
  rank_sum_focus <- sum(ranks[seq_along(focus)])
  u_focus <- rank_sum_focus - n_focus * (n_focus + 1) / 2
  cliffs_delta <- 2 * (u_focus / (n_focus * n_reference)) - 1

  tibble(
    n_focus = n_focus,
    n_reference = n_reference,
    median_focus = median(focus, na.rm = TRUE),
    median_reference = median(reference, na.rm = TRUE),
    median_difference = median_focus - median_reference,
    cliffs_delta = cliffs_delta,
    p_value = suppressWarnings(wilcox.test(value ~ behavior_comparison, data = df, exact = FALSE)$p.value)
  )
}

driver_comparison <- driver_long %>%
  group_by(variable, driver_group, label) %>%
  group_modify(~ compare_one_variable(.x)) %>%
  ungroup() %>%
  mutate(
    comparison_mode = comparison_mode,
    focus_behavior = focus_behavior,
    reference_behavior = reference_behavior,
    p_adj_bh = p.adjust(p_value, method = "BH"),
    direction = case_when(
      median_difference > 0 ~ paste("higher in", focus_behavior),
      median_difference < 0 ~ paste("lower in", focus_behavior),
      TRUE ~ "no difference"
    ),
    abs_delta = abs(cliffs_delta),
    evidence = case_when(
      is.na(p_value) ~ "not tested",
      p_value < 0.05 ~ "p < 0.05",
      p_value < 0.10 ~ "p < 0.10",
      TRUE ~ "weak"
    )
  ) %>%
  arrange(desc(abs_delta))

write.csv(site_drivers, "OUTPUT/NEON_site_driver_values_for_sink_comparison.csv", row.names = FALSE)
write.csv(driver_comparison, "OUTPUT/NEON_strong_sink_driver_comparison.csv", row.names = FALSE)

plot_effects_data <- driver_comparison %>%
  filter(is.finite(cliffs_delta)) %>%
  mutate(
    label    = factor(label, levels = rev(unique(label[order(cliffs_delta)]))),
    evidence = factor(evidence, levels = c("p < 0.05", "p < 0.10", "weak", "not tested"))
  )

if (nrow(plot_effects_data) == 0) {
  stop("No finite Cliff's delta values are available for plot_effects.")
}

effect_focus <- unique(plot_effects_data$focus_behavior)[1]
effect_reference <- unique(plot_effects_data$reference_behavior)[1]
effect_title <- paste(effect_focus, "vs", effect_reference)

plot_effects <- plot_effects_data %>%
  ggplot(aes(x = cliffs_delta, y = label, color = driver_group, shape = evidence)) +
  # effect-size interpretation bands (|δ|: <0.1 negligible, 0.1-0.3 small, 0.3-0.5 medium, >0.5 large)
  annotate("rect", xmin = -1.0, xmax = -0.5, ymin = -Inf, ymax = Inf, fill = "#c62828", alpha = 0.06) +
  annotate("rect", xmin = -0.5, xmax = -0.3, ymin = -Inf, ymax = Inf, fill = "#ef9a9a", alpha = 0.08) +
  annotate("rect", xmin = -0.3, xmax = -0.1, ymin = -Inf, ymax = Inf, fill = "#ffccbc", alpha = 0.09) +
  annotate("rect", xmin = -0.1, xmax =  0.1, ymin = -Inf, ymax = Inf, fill = "grey90",  alpha = 0.50) +
  annotate("rect", xmin =  0.1, xmax =  0.3, ymin = -Inf, ymax = Inf, fill = "#c8e6c9", alpha = 0.09) +
  annotate("rect", xmin =  0.3, xmax =  0.5, ymin = -Inf, ymax = Inf, fill = "#a5d6a7", alpha = 0.08) +
  annotate("rect", xmin =  0.5, xmax =  1.0, ymin = -Inf, ymax = Inf, fill = "#2e7d32", alpha = 0.06) +
  # band labels along top margin
  annotate("text", x = c(-0.75,-0.40,-0.20, 0, 0.20, 0.40, 0.75), y = Inf, vjust = -0.3,
           label = c("large","medium","small","negligible","small","medium","large"),
           size = 2.6, color = "grey50", fontface = "italic") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.8) +
  geom_segment(aes(x = 0, xend = cliffs_delta, yend = label), linewidth = 0.9, alpha = 0.78) +
  geom_point(size = 3.5) +
  scale_shape_manual(values = c("p < 0.05" = 18, "p < 0.10" = 17, "weak" = 16, "not tested" = 4)) +
  scale_color_manual(values = driver_group_colors, na.translate = FALSE) +
  coord_cartesian(xlim = c(-0.85, 0.85), clip = "off") +
  labs(
    title    = paste("A. Driver differences:", effect_title),
    subtitle = paste("Positive Cliff's delta = variable higher in", effect_focus,
                     "- Background bands = conventional effect-size thresholds"),
    x        = "Cliff's delta effect size",
    y        = NULL,
    color    = "Driver group",
    shape    = "Evidence"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey35"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.margin = ggplot2::margin(t = 18, r = 8, b = 6, l = 6)
  )

plot_all_distributions <- driver_long %>%
  mutate(label = factor(label, levels = available_driver_groups$label)) %>%
  ggplot(aes(x = behavior_comparison, y = value_z, color = behavior_comparison)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.12, width = 0.55) +
  geom_jitter(width = 0.13, alpha = 0.75, size = 1.8) +
  facet_wrap(~driver_group + label, scales = "free_y", ncol = 4) +
  scale_color_manual(values = comparison_colors, guide = "none") +
  labs(
    title = "All site-level driver distributions",
    x = NULL,
    y = "Standardized site value"
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey94", color = "grey35"),
    strip.text = element_text(face = "bold", size = 8),
    panel.spacing = unit(0.45, "lines")
  )

top_driver_variables <- driver_comparison$variable[seq_len(min(8, nrow(driver_comparison)))]

plot_by_behavior <- driver_long %>%
  filter(variable %in% top_driver_variables) %>%
  mutate(label = factor(label, levels = driver_comparison$label[seq_len(min(8, nrow(driver_comparison)))])) %>%
  ggplot(aes(x = CH4_behavior, y = value_z, color = CH4_behavior)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.12, width = 0.55) +
  geom_jitter(width = 0.13, alpha = 0.75, size = 1.8) +
  facet_wrap(~label, scales = "free_y", ncol = 4) +
  scale_color_manual(values = behavior_colors, guide = "none") +
  labs(
    title = "B. Top driver contrasts shown by original CH4 behavior class",
    x = NULL,
    y = "Standardized site value"
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 8.5),
    strip.background = element_rect(fill = "grey94", color = "grey35"),
    strip.text = element_text(face = "bold", size = 8.5),
    panel.spacing = unit(0.6, "lines")
  )

sink_driver_figure <- (plot_effects / plot_by_behavior) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(
    title = "NEON CH4 Total-Flux Behavior Driver Comparison",
    subtitle = comparison_title,
    caption = "Wilcoxon tests and Cliff's delta are exploratory because site counts are small and drivers covary.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      plot.caption = element_text(size = 9, color = "grey30")
    )
  )

ggsave(
  "FIGURES/NEON_strong_sink_driver_comparison.png",
  plot = sink_driver_figure,
  width = 14,
  height = 12,
  units = "in",
  dpi = 300
)

ggsave(
  "FIGURES/NEON_strong_sink_driver_comparison.pdf",
  plot = sink_driver_figure,
  width = 14,
  height = 12,
  units = "in"
)

ggsave(
  "FIGURES/NEON_total_flux_driver_distributions_all.png",
  plot = plot_all_distributions,
  width = 13,
  height = 16,
  units = "in",
  dpi = 300
)

ggsave(
  "FIGURES/NEON_total_flux_driver_distributions_all.pdf",
  plot = plot_all_distributions,
  width = 13,
  height = 16,
  units = "in"
)

top_lines <- driver_comparison %>%
  slice_head(n = 8) %>%
  mutate(
    line = paste0(
      "- ", label, " (", driver_group, "): ", focus_behavior, " median = ",
      signif(median_focus, 3), ", ", reference_behavior, " median = ",
      signif(median_reference, 3), ", Cliff's delta = ", signif(cliffs_delta, 3),
      ", p = ", signif(p_value, 3), ", BH p = ", signif(p_adj_bh, 3),
      " (", direction, ")"
    )
  ) %>%
  pull(line)

writeLines(
  c(
    "# NEON CH4 Total-Flux Behavior Driver Comparison",
    "",
    paste0("Comparison mode: `", comparison_mode, "`."),
    paste0("Focus group: `", focus_behavior, "`. Reference group: `", reference_behavior, "`."),
    "Behavior classes are from `OUTPUT/30min_site_behavior.csv`, where monthly source/sink state is based on equal-weighted site-month-hour flux bins.",
    "Moisture driver summaries in this script are also based on equal-weighted site-month-hour bins before site-level means, seasonal SDs, and ranges are calculated.",
    "",
    "## Strongest Contrasts",
    top_lines,
    "",
    "## Interpretation",
    paste0("- Positive Cliff's delta means a variable is higher in `", focus_behavior, "`; negative means lower."),
    "- Treat p-values as exploratory because site counts are small and many site attributes covary.",
    "- Under the current total-flux behavior classes, there are no `Weak-sink` sites; this script therefore compares available behavior classes unless sinks reappear in future outputs.",
    "",
    "## Outputs",
    "- `OUTPUT/NEON_strong_sink_driver_comparison.csv`",
    "- `OUTPUT/NEON_site_driver_values_for_sink_comparison.csv`",
    "- `FIGURES/NEON_strong_sink_driver_comparison.png`",
    "- `FIGURES/NEON_strong_sink_driver_comparison.pdf`",
    "- `FIGURES/NEON_total_flux_driver_distributions_all.png`",
    "- `FIGURES/NEON_total_flux_driver_distributions_all.pdf`"
  ),
  "OUTPUT/NEON_strong_sink_driver_comparison_results.md"
)

message("Wrote NEON total-flux behavior driver comparison outputs.")
