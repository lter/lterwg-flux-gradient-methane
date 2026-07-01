# ERA-SpatialUpscaling-Figures-RF.R
#
# Standalone figure script for the Random Forest spatial upscaling pipeline.
# Loads pre-computed outputs from ERA-SpatialUpscaling-Monthly-RF.R and
# regenerates all publication figures without re-fitting models.
#
# Produces:
#   Fig1_source_probability.png      — P(source) density, maps, source/sink area
#   FigS1_model_performance.png      — Calibration, skill metrics, variable importance
#   Fig2_magnitude_models.png        — Obs vs fitted, seasonal flux cycle
#   Fig3_budget.png                  — Annual time series, bar chart, GMB sensitivity
#   Fig4_spatial_maps.png            — Per-approach flux maps for 2025
#
# Prerequisites (written by ERA-SpatialUpscaling-Monthly-RF.R):
#   <output_dir>/monthly_cell_predictions_2000_2025.rds
#   <output_dir>/era5_template.tif
#   <output_dir>/model_parameters.csv
#   <output_dir>/annual_budget_three_approaches.csv
#   <output_dir>/annual_budget_long.csv
#   <output_dir>/budget_summary_three_approaches.csv
#   <output_dir>/comparison_budget_all_approaches.csv
#   <output_dir>/comparison_class_skill_GLM_vs_RF.csv
#   <output_dir>/probability_calibration_skill.csv
#   <output_dir>/magnitude_model_fitted_values.csv
#   <output_dir>/magnitude_model_skill.csv
#   <output_dir>/gmb_threshold_sensitivity.csv
#   <output_dir>/rf_class_variable_importance.csv

library(tidyverse)
library(terra)
library(cowplot)

# ─────────────────────────────────────────────────────────────────────────────
# Paths  (must match ERA-SpatialUpscaling-Monthly-RF.R)
# ─────────────────────────────────────────────────────────────────────────────

# When sourced from ERA-SpatialUpscaling-Monthly-RF.R, output_dir and
# figure_dir are already defined — use them directly.
if (!exists("output_dir")) {
  rf_dir     <- "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF"
  output_dir <- file.path(rf_dir, "OUTPUT")
}
if (!exists("figure_dir")) {
  figure_dir <- file.path(dirname(output_dir), "FIGURES")
}
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Load saved outputs  (skipped when already in memory from the main script)
# ─────────────────────────────────────────────────────────────────────────────

if (!exists("cell_preds")) {
  message("Loading monthly_cell_predictions_2000_2025.rds ...")
  cell_preds <- readRDS(file.path(output_dir, "monthly_cell_predictions_2000_2025.rds"))
} else {
  message("Using cell_preds already in memory.")
}

if (!exists("template") || !inherits(template, "SpatRaster"))
  template <- rast(file.path(output_dir, "era5_template.tif"))

# Scalar parameters: always re-read from CSV for consistency
params              <- read.csv(file.path(output_dir, "model_parameters.csv"))
P_star              <- params$P_star
binary_threshold    <- params$binary_threshold
gmb_soil_sink_tg_ch4_yr <- params$gmb_soil_sink_tg_ch4_yr

# Tabular outputs
recode_approach <- function(df) {
  mutate(df, across(any_of("approach"), ~recode(.,
    "GMB-Continuous"  = "All-Sink",
    "GMB-constrained" = "All-Sink",
    "GMB-Dichotomous" = "All-Sink"
  )))
}

# Always reload: normalize approach labels regardless of what's in memory
annual_budget     <- read.csv(file.path(output_dir, "annual_budget_three_approaches.csv")) %>% recode_approach()
annual_budget_long <- read.csv(file.path(output_dir, "annual_budget_long.csv")) %>% recode_approach()
budget_summary    <- read.csv(file.path(output_dir, "budget_summary_three_approaches.csv")) %>% recode_approach()
comparison_budget <- read.csv(file.path(output_dir, "comparison_budget_all_approaches.csv")) %>% recode_approach()
# Always reload: RF version has 'model' column + brier_null/brier_skill_score
comparison_class <- read.csv(file.path(output_dir, "comparison_class_skill_GLM_vs_RF.csv"))
# Always reload: RF version has a 'model' column (A/B) that the GLM version lacks
probability_calibration_skill <- read.csv(file.path(output_dir, "probability_calibration_skill.csv"))
# Always reload: normalize magnitude_model labels regardless of what's in memory
magnitude_fitted_values <- read.csv(file.path(output_dir, "magnitude_model_fitted_values.csv")) %>%
  mutate(magnitude_model = recode(magnitude_model,
                                  "Weak sink"   = "Weak-sink",
                                  "Weak source" = "Weak-source"))
magnitude_model_skill <- read.csv(file.path(output_dir, "magnitude_model_skill.csv")) %>%
  mutate(magnitude_model = recode(magnitude_model,
                                  "Weak sink"   = "Weak-sink",
                                  "Weak source" = "Weak-source"))
if (!exists("gmb_sensitivity"))
  gmb_sensitivity <- read.csv(file.path(output_dir, "gmb_threshold_sensitivity.csv"))
# Always reload: needs 'model' column (A — weighted / B — 1:1 balanced)
rf_class_importance <- read.csv(file.path(output_dir, "rf_class_variable_importance.csv"))

# ─────────────────────────────────────────────────────────────────────────────
# Helper
# ─────────────────────────────────────────────────────────────────────────────

gC_m2_yr_to_tg_ch4 <- function(flux_gC_m2_yr, area_mha) {
  # g C m⁻² × Mha × (16/12) × 1e10 m²/Mha × 1e-12 Tg/g
  flux_gC_m2_yr * area_mha * 0.013333333
}

# ─────────────────────────────────────────────────────────────────────────────
# Shared aesthetics
# ─────────────────────────────────────────────────────────────────────────────

# Colors chosen to avoid overlap with:
#   sink (#2166AC blue), source (#B2182B red),
#   Forest (#1B7837 green), Grassland (#D9B86C tan),
#   Shrubland (#C2A5CF lavender), Arid (#D95F02 orange)
approach_colors <- c(
  "Continuous"     = "#009688",  # material teal
  "Dichotomous"    = "#9C27B0",  # material purple (distinct from pale shrubland lavender)
  "All-Sink"       = "#E91E63"   # material pink (distinct from dark source red)
)

model_colors <- c("Continuous" = "#009688")

ecotype_colors <- c(
  "Forest"    = "#1B7837",
  "Grassland" = "#D9B86C",
  "Shrubland" = "#C2A5CF",
  "Arid"      = "#D95F02"
)

sink_col   <- "#2166AC"
source_col <- "#B2182B"
gmb_col    <- "#2166AC"

fig_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold", size = 12),
    plot.tag          = element_text(face = "bold", size = 12),
    legend.position   = "bottom",
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 10),
    strip.background  = element_rect(fill = "grey92"),
    strip.text        = element_text(size = 11),
    axis.title        = element_text(size = 11),
    axis.text         = element_text(size = 10)
  )

map_theme <- fig_theme +
  theme(panel.grid = element_blank(), legend.key.width = unit(1.2, "cm"))

flux_scale <- function(lim, name = "g C m⁻² yr⁻¹") {
  scale_fill_gradient2(
    low      = sink_col, mid = "#D4C08A", high = source_col,
    midpoint = 0, limits = c(-lim, lim),
    breaks   = pretty(c(-lim, lim), n = 5),
    name     = name, na.value = "transparent"
  )
}

prob_scale <- scale_fill_gradient2(
  low = sink_col, mid = "#D4C08A", high = source_col,
  midpoint = 0.5, limits = c(0, 1),
  breaks = c(0, 0.25, 0.5, 0.75, 1),
  name = "P(source)", na.value = "transparent"
)

make_flux_map <- function(data, flux_col, title_label, flux_abs_max) {
  data %>%
    ggplot(aes(x = x, y = y, fill = .data[[flux_col]])) +
    geom_tile(width = res(template)[1], height = res(template)[2]) +
    coord_equal(expand = FALSE) +
    flux_scale(flux_abs_max) +
    labs(title = title_label, x = "Longitude", y = "Latitude") +
    map_theme
}

# ─────────────────────────────────────────────────────────────────────────────
# Derived datasets
# ─────────────────────────────────────────────────────────────────────────────

cp <- cell_preds %>%
  mutate(EcoType_plot = if_else(is_arid == 1L, "Arid", as.character(EcoType)))

# Annual mean P(source) and majority class — 2025 (Model A)
annual_cell_A_2025 <- cp %>%
  filter(Year == 2025) %>%
  group_by(cell, x, y) %>%
  summarise(
    mean_source_prob_A    = mean(source_prob_A, na.rm = TRUE),
    annual_exchange_class = if_else(mean(source_prob_A, na.rm = TRUE) >= 0.5,
                                    "Weak-source", "Weak-sink"),
    .groups = "drop"
  )

# Annual source/sink area time series (Model A classified at P≥0.5)
annual_source_area <- cp %>%
  group_by(Year, cell) %>%
  summarise(
    area_mha          = first(area_mha),
    frac_source_months = mean(source_prob_A >= 0.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    source_area = area_mha * frac_source_months,
    sink_area   = area_mha * (1 - frac_source_months)
  ) %>%
  group_by(Year) %>%
  summarise(
    `Weak-source` = sum(source_area, na.rm = TRUE),
    `Weak-sink`   = sum(sink_area,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(`Weak-source`, `Weak-sink`), names_to = "class", values_to = "area_mha")

# Cell-level annual fluxes 2025 (all three approaches)
cell_flux_2025 <- cp %>%
  filter(Year == 2025) %>%
  group_by(cell, x, y) %>%
  summarise(
    annual_continuous_gC_m2_yr  = sum(flux_continuous,  na.rm = TRUE),
    annual_balanced_gC_m2_yr    = sum(flux_balanced,     na.rm = TRUE),
    annual_constrained_gC_m2_yr = sum(flux_constrained,  na.rm = TRUE),
    .groups = "drop"
  )

flux_map_abs_max <- max(abs(range(
  c(cell_flux_2025$annual_continuous_gC_m2_yr,
    cell_flux_2025$annual_balanced_gC_m2_yr,
    cell_flux_2025$annual_constrained_gC_m2_yr),
  na.rm = TRUE)))

# Seasonal area-weighted flux by EcoType (all three approaches)
flux_seasonal <- cp %>%
  group_by(month, EcoType_plot) %>%
  summarise(
    Continuous  = weighted.mean(flux_continuous,  area_mha, na.rm = TRUE),
    Balanced    = weighted.mean(flux_balanced,     area_mha, na.rm = TRUE),
    Constrained = weighted.mean(flux_constrained,  area_mha, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(Continuous, Balanced, Constrained),
               names_to = "approach", values_to = "flux") %>%
  mutate(approach = factor(approach,
    levels = c("Continuous", "Balanced", "Constrained")))

# Annual budget by EcoType — continuous approach
budget_by_ecotype_year <- cp %>%
  mutate(tg_cell = gC_m2_yr_to_tg_ch4(flux_continuous, area_mha)) %>%
  group_by(Year, EcoType_plot) %>%
  summarise(annual_tg = sum(tg_cell, na.rm = TRUE), .groups = "drop")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1  Source-probability model  (8.0 × 7.0 in)
# A: P(source) density by EcoType (Model A)
# B: Annual mean P(source) map, 2025 (Model A)
# C: Annual exchange-class map, 2025 (Model A, P≥0.5) + inset area time series
# ─────────────────────────────────────────────────────────────────────────────

p1a <- cp %>%
  group_by(cell, EcoType_plot) %>%
  summarise(mean_prob = mean(source_prob_A, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = mean_prob, fill = EcoType_plot, color = EcoType_plot)) +
  geom_density(alpha = 0.28, linewidth = 0.7) +
  geom_vline(xintercept = 0.5, color = "grey30", linetype = "dotted", linewidth = 0.8) +
  scale_fill_manual(values  = ecotype_colors) +
  scale_color_manual(values = ecotype_colors) +
  labs(title = "A. P(source) by Ecosystem Type",
       x = "Mean P(source)", y = "Density", fill = NULL, color = NULL) +
  fig_theme +
  theme(legend.position   = c(0.18, 0.75),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size   = unit(0.35, "cm"),
        legend.text       = element_text(size = 6.5))

p1b <- annual_cell_A_2025 %>%
  ggplot(aes(x = x, y = y, fill = mean_source_prob_A)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  coord_equal(expand = FALSE) +
  prob_scale +
  labs(title = "B. Annual Mean P(source), 2025",
       x = "Longitude", y = "Latitude") +
  map_theme +
  theme(legend.position   = c(0.08, 0.28),
        legend.direction  = "vertical",
        legend.key.width  = unit(0.4, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.text       = element_text(size = 6.5),
        legend.title      = element_text(size = 7),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA))

# Build exchange-class raster for 2025
class_raster <- template
names(class_raster) <- "class"
values(class_raster) <- NA_real_
class_raster[annual_cell_A_2025$cell] <- if_else(
  annual_cell_A_2025$annual_exchange_class == "Weak-source", 2, 1)
class_map_df <- as.data.frame(class_raster, xy = TRUE, na.rm = TRUE) %>%
  mutate(annual_exchange_class = factor(class,
    levels = c(1, 2), labels = c("Weak-sink", "Weak-source")))

p1c_base <- class_map_df %>%
  ggplot(aes(x = x, y = y, fill = annual_exchange_class)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  coord_equal(expand = FALSE) +
  scale_fill_manual(values   = c("Weak-sink" = sink_col, "Weak-source" = source_col),
                    na.value = "transparent") +
  labs(title = "C. Annual Exchange Class, 2025  (P(source) ≥ 0.5)",
       x = "Longitude", y = "Latitude", fill = NULL) +
  map_theme +
  theme(legend.position  = "top",
        legend.direction = "horizontal",
        legend.key.size  = unit(0.4, "cm"),
        legend.text      = element_text(size = 7))

p1d <- annual_source_area %>%
  ggplot(aes(x = Year, y = area_mha / 100, color = class)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_color_manual(values = c("Weak-sink" = sink_col, "Weak-source" = source_col)) +
  labs(title = "D. Source/Sink Area, 2000–2025",
       x = "Year", y = "Area (100 Mha)", color = NULL) +
  fig_theme

p1c_inset <- ggdraw() +
  draw_plot(p1c_base) +
  draw_plot(p1d + theme(plot.title       = element_blank(),
                        axis.title       = element_text(size = 7),
                        axis.text        = element_text(size = 6),
                        legend.position  = "none",
                        plot.background  = element_rect(fill = "transparent", color = NA),
                        panel.background = element_rect(fill = "transparent")),
            x = 0.16, y = 0.08, width = 0.211, height = 0.237) +
  draw_label("D.", x = 0.16, y = 0.317, hjust = 0, vjust = 0, size = 8, fontface = "bold")

fig1 <- plot_grid(
  plot_grid(p1a, p1b, ncol = 2, rel_widths = c(0.7, 1.3)),
  p1c_inset,
  ncol = 1, rel_heights = c(1, 1.1)
)

ggsave(file.path(figure_dir, "Fig1_source_probability.png"),
  fig1, width = 8.0, height = 7.0, units = "in", dpi = 300, bg = "white")
message("Fig 1 written.")

# ─────────────────────────────────────────────────────────────────────────────
# Fig S1  RF Model Performance  (9.0 × 9.5 in)
# A: Probability calibration — P(source) model (OOB → isotonic)
# B: Stage 1 classification skill metrics
# C: Stage 1 variable importance (top 8 predictors)
# ─────────────────────────────────────────────────────────────────────────────

pS1a <- probability_calibration_skill %>%
  ggplot(aes(x = mean_isotonic_cal_prob, y = observed_source_fraction)) +
  geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(size = n), color = "#009688", alpha = 0.85) +
  geom_line(color = "#009688", linewidth = 0.8) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_size_continuous(range = c(1.5, 5)) +
  labs(title = "A. Probability Calibration  P(source)",
       x = "Isotonic-calibrated P(source)", y = "Observed source fraction", size = "n months") +
  fig_theme

# Classification skill — long format
skill_long <- comparison_class %>%
  pivot_longer(
    cols      = any_of(c("auc", "brier_skill_score", "tjur_r2",
                          "accuracy", "sensitivity", "specificity")),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    metric = factor(metric,
      levels = c("auc", "brier_skill_score", "tjur_r2",
                 "accuracy", "sensitivity", "specificity"),
      labels = c("AUC", "Brier Skill Score", "Tjur R²",
                 "Accuracy", "Sensitivity", "Specificity")),
    group = if_else(metric %in% c("AUC", "Brier Skill Score", "Tjur R²"),
                    "Discrimination", "Classification")
  )

pS1b <- skill_long %>%
  ggplot(aes(x = value, y = metric)) +
  geom_col(width = 0.55, alpha = 0.85, fill = "#009688") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_text(aes(label = sprintf("%.3f", value)), hjust = -0.15, size = 3.2) +
  scale_x_continuous(limits = c(0, 1.15), expand = c(0, 0),
                     breaks = seq(0, 1, by = 0.25)) +
  facet_grid(group ~ ., scales = "free_y", space = "free") +
  labs(title = "B. Stage 1 Classification Skill",
       x = "Value", y = NULL) +
  fig_theme + theme(legend.position = "none")

# Variable importance — top 8 predictors
pS1c <- rf_class_importance %>%
  arrange(desc(importance)) %>%
  head(8) %>%
  mutate(predictor = fct_reorder(predictor, importance)) %>%
  ggplot(aes(x = importance, y = predictor)) +
  geom_col(fill = "#009688", alpha = 0.85) +
  labs(title = "C. Stage 1 Variable Importance",
       x = "Permutation importance", y = NULL) +
  fig_theme + theme(legend.position = "none")

figS1 <- plot_grid(
  pS1a,
  plot_grid(pS1b, pS1c, ncol = 2),
  ncol = 1, rel_heights = c(1, 1.1)
)

ggsave(file.path(figure_dir, "FigS1_model_performance.png"),
  figS1, width = 9.0, height = 9.5, units = "in", dpi = 300, bg = "white")
message("Fig S1 written.")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2  Magnitude models  (8.0 × 8.0 in)
# A: Obs vs fitted — Weak-sink
# B: Obs vs fitted — Weak-source
# C: Seasonal flux cycle by EcoType (three approaches)
# ─────────────────────────────────────────────────────────────────────────────

make_mag_plot <- function(model_label, panel_tag) {
  n_obs <- nrow(filter(magnitude_fitted_values, magnitude_model == model_label))
  magnitude_fitted_values %>%
    filter(magnitude_model == model_label, !is.na(EcoType)) %>%
    ggplot(aes(x = monthly_flux_gC_m2_month, y = fitted_flux_gC_m2_month, color = EcoType)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "grey35",
                linetype = "dashed", linewidth = 0.5) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = ecotype_colors) +
    labs(title = paste0(panel_tag, ". ", model_label, " model  (n = ", n_obs, ")"),
         x = "Observed (g C m⁻² mo⁻¹)",
         y = "Fitted (g C m⁻² mo⁻¹)",
         color = NULL) +
    fig_theme
}

p2a <- make_mag_plot("Weak-sink",   "A")
p2b <- make_mag_plot("Weak-source", "B")

p2c <- flux_seasonal %>%
  ggplot(aes(x = month, y = flux, color = EcoType_plot, linetype = approach)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 1.3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = ecotype_colors) +
  scale_linetype_manual(
    values = c("Continuous" = "solid", "Balanced" = "dashed", "Constrained" = "dotted"),
    labels = c("Continuous"  = "Continuous",
               "Balanced"    = "Dichotomous",
               "Constrained" = "All-Sink")
  ) +
  facet_wrap(~EcoType_plot, scales = "free_y", ncol = 2) +
  labs(title = "C. Seasonal Flux Cycle by Ecosystem Type",
       x = "Month", y = "Flux (g C m⁻² mo⁻¹)",
       color = NULL, linetype = "Approach") +
  fig_theme + guides(color = "none")

shared_ecotype_legend <- get_legend(
  p2a + theme(legend.position = "top",
              legend.key.size = unit(0.4, "cm"),
              legend.text     = element_text(size = 8))
)

fig2 <- plot_grid(
  shared_ecotype_legend,
  plot_grid(p2a + theme(legend.position = "none"),
            p2b + theme(legend.position = "none"),
            ncol = 2),
  p2c,
  ncol = 1, rel_heights = c(0.07, 1, 1.4)
)

ggsave(file.path(figure_dir, "Fig2_magnitude_models.png"),
  fig2, width = 8.0, height = 8.0, units = "in", dpi = 300, bg = "white")
message("Fig 2 written.")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 3  Budget  (7.5 × 9.5 in)
# A: Mean budget bar chart — all approaches
# B: Annual time series — three RF approaches
# C: GMB threshold sensitivity (P(source))
# ─────────────────────────────────────────────────────────────────────────────

p3a <- comparison_budget %>%
  mutate(approach_wrap = str_wrap(approach, 28)) %>%
  ggplot(aes(x = reorder(approach_wrap, mean_tg_ch4_yr),
             y = mean_tg_ch4_yr, fill = approach)) +
  geom_col(width = 0.65, alpha = 0.85, show.legend = FALSE) +
  geom_errorbar(aes(ymin = min_tg_ch4_yr, ymax = max_tg_ch4_yr),
                width = 0.25, linewidth = 0.6) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = gmb_col,
             linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  annotate("text", x = -Inf, y = gmb_soil_sink_tg_ch4_yr - 1.5,
           label = "GMB (−35)", hjust = 0, size = 3.0, color = gmb_col) +
  scale_fill_manual(values = approach_colors) +
  labs(title = "A. Mean Annual Budget (2000–2025)",
       x = NULL, y = "Net exchange (Tg CH₄ yr⁻¹)") +
  coord_flip() +
  fig_theme

p3b <- annual_budget_long %>%
  ggplot(aes(x = Year, y = annual_tg_ch4_yr, color = approach)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -34,
           fill = gmb_col, alpha = 0.10) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = gmb_col,
             linetype = "dashed", linewidth = 0.7) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_color_manual(values = approach_colors) +
  annotate("text", x = 2001, y = gmb_soil_sink_tg_ch4_yr + 1.8,
           label = "GMB reference (−35 Tg yr⁻¹)",
           hjust = 0, size = 3.0, color = gmb_col) +
  labs(title = "B. Annual Net Exchange, 2000–2025",
       x = "Year", y = "Net exchange (Tg CH₄ yr⁻¹)", color = NULL) +
  fig_theme + theme(legend.position = "none")

p3c <- gmb_sensitivity %>%
  ggplot(aes(x = threshold, y = mean_tg)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -34,
           fill = gmb_col, alpha = 0.10) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = gmb_col,
             linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_line(color = "black", linewidth = 1.1) +
  geom_point(color = "black", size = 1.5) +
  {
    budget_floor_val <- min(gmb_sensitivity$mean_tg, na.rm = TRUE)
    floor_threshold  <- gmb_sensitivity$threshold[which.min(gmb_sensitivity$mean_tg)]
    list(
      geom_vline(xintercept = floor_threshold, color = "black",
                 linetype = "dotted", linewidth = 0.9),
      annotate("label", x = floor_threshold, y = budget_floor_val,
               label = sprintf("Budget floor\n%.1f Tg (P = %.2f)", budget_floor_val, floor_threshold),
               hjust = -0.1, vjust = 0, size = 3.2, color = "black",
               fill = "white", label.size = 0)
    )
  } +
  scale_x_continuous(breaks = seq(0.3, 1.0, by = 0.1)) +
  labs(title = "C. Threshold Sensitivity  P(source)",
       x = "P(source) threshold",
       y = "Mean annual budget (Tg CH₄ yr⁻¹)") +
  fig_theme

fig3 <- plot_grid(p3a, p3b, p3c, ncol = 1, rel_heights = c(1.1, 1, 1))

ggsave(file.path(figure_dir, "Fig3_budget.png"),
  fig3, width = 7.5, height = 9.5, units = "in", dpi = 300, bg = "white")
message("Fig 3 written.")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4  Spatial flux maps, 2025  (9.0 × 9.5 in)
# A: Continuous  (P(source))
# B: Dichotomous  (P(source) ≥ 0.5)
# C: GMB-Dichotomous  (P(source) ≥ P*)
# ─────────────────────────────────────────────────────────────────────────────

p4a <- make_flux_map(
  cell_flux_2025, "annual_continuous_gC_m2_yr",
  "A. Annual Net Flux, 2025  —  Continuous  (P(source))",
  flux_map_abs_max
)

p4b <- make_flux_map(
  cell_flux_2025, "annual_balanced_gC_m2_yr",
  "B. Annual Net Flux, 2025  —  Dichotomous  (P(source) ≥ 0.5)",
  flux_map_abs_max
)

p4c <- make_flux_map(
  cell_flux_2025, "annual_constrained_gC_m2_yr",
  "C. Annual Net Flux, 2025  —  All-Sink  (all upland cells → sink flux)",
  flux_map_abs_max
)

fig4 <- plot_grid(p4a, p4b, p4c, ncol = 1)

ggsave(file.path(figure_dir, "Fig4_spatial_maps.png"),
  fig4, width = 9.0, height = 9.5, units = "in", dpi = 300, bg = "white")
message("Fig 4 written.")

# ─────────────────────────────────────────────────────────────────────────────
# Standalone budget time series (for external use)
# ─────────────────────────────────────────────────────────────────────────────

ggsave(file.path(figure_dir, "annual_budget_time_series_2000_2025.png"),
  p3a + labs(title = NULL),
  width = 7.0, height = 3.5, units = "in", dpi = 300, bg = "white")

message("All RF figures written to: ", figure_dir)
