# 14_Global_SpatialUpscalingFiguresRF.R
#
# Standalone figure script for the RF spatial-projection OUTPUT (the global
# upscaling grid produced by 13_Global_SpatialUpscalingRF.R). Loads
# pre-computed outputs from that script and regenerates the upscaling
# figures without re-running the projection.
#
# Model-performance/diagnostic figures (calibration, Stage 1 classification
# skill, variable importance, observed-vs-OOB-fitted magnitude) live in a
# separate script, 12b_Model_FIGURES_RF.R, since those depend only on
# 12_SourceProp_MagnitudeModels.R and have nothing to do with the spatial
# grid — this script no longer produces them as their OWN standalone
# figures. It DOES still reassemble the original, pre-split
# Fig2_magnitude_models.png (3 panels: A/B obs-vs-OOB-fitted Weak-sink/
# Weak-source from 12_SourceProp_MagnitudeModels.R, C seasonal flux cycle
# from this script's own upscaling output) for anything downstream that
# still expects that exact filename — this is the only figure in this
# script that reads from 12's output rather than 13's.
#
# Produces:
#   Fig1_source_probability.png      — P(source) density, maps, source/sink area
#   Fig2_seasonal_flux_cycle.png     — Seasonal flux cycle by Ecosystem Type (three approaches)
#   Fig2_magnitude_models.png        — Original 3-panel figure (A/B from 12, C from this script)
#   Fig3_budget.png                  — Annual time series, bar chart, GMB sensitivity
#   Fig4_spatial_maps.png            — Per-approach flux maps for 2025
#   annual_budget_time_series_2000_2025.png
#   annual_dichotomous_exchange_class_area.csv         — Weak-source/
#     Weak-sink area (Mha) by year, 2000-2025, dichotomous P*=0.5 rule
#   annual_dichotomous_exchange_class_area_summary.csv — mean/SD/CV%/min/max
#     of the above, for interannual-stability claims
#
# Prerequisites (written by 13_Global_SpatialUpscalingRF.R):
#   <output_dir>/monthly_cell_predictions_2000_2025.rds
#   <output_dir>/era5_template.tif
#   <output_dir>/spatial_projection_parameters.csv
#   <output_dir>/annual_budget_three_approaches.csv
#   <output_dir>/annual_budget_long.csv
#   <output_dir>/budget_summary_three_approaches.csv
#   <output_dir>/comparison_budget_all_approaches.csv
#   <output_dir>/gmb_threshold_sensitivity.csv
#
# Optional prerequisite (written by 12_SourceProp_MagnitudeModels.R; if
# absent, Fig2_magnitude_models.png is skipped with a message but everything
# else in this script still runs):
#   <output_dir>/magnitude_model_fitted_values.csv

library(tidyverse)
library(terra)
library(cowplot)

# ─────────────────────────────────────────────────────────────────────────────
# Paths  (must match 13_Global_SpatialUpscalingRF.R)
# ─────────────────────────────────────────────────────────────────────────────

# When sourced from 13_Global_SpatialUpscalingRF.R, output_dir and
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

# A SpatRaster's underlying data lives in an external C++ pointer that does
# NOT survive R session save/restore (e.g. RStudio's "restore .RData on
# startup", or re-running this script in a session where the raster was
# created earlier and then the workspace was reloaded). inherits(template,
# "SpatRaster") stays TRUE even when that pointer has gone stale/NULL, so a
# class-only check silently reuses a broken object and any later raster
# operation fails with an Rcpp error such as: NULL value passed as symbol
# address. Actually exercising the pointer via terra::ncell(), inside try(),
# catches this and forces a fresh read from disk when needed.
template_is_valid <- exists("template") &&
  inherits(template, "SpatRaster") &&
  !inherits(try(terra::ncell(template), silent = TRUE), "try-error")
if (!template_is_valid) {
  template <- rast(file.path(output_dir, "era5_template.tif"))
}

# Scalar parameters: always re-read from CSV for consistency
params              <- read.csv(file.path(output_dir, "spatial_projection_parameters.csv"))
P_star              <- params$P_star
binary_threshold    <- params$binary_threshold
gmb_soil_sink_tg_ch4_yr <- params$gmb_target_tg_ch4_yr

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
if (!exists("gmb_sensitivity"))
  gmb_sensitivity <- read.csv(file.path(output_dir, "gmb_threshold_sensitivity.csv"))

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
    strip.background  = element_rect(fill = "black"),
    strip.text        = element_text(size = 11, color = "white"),
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

# Dichotomous annual exchange-class area, ALL years (generalizes
# annual_cell_A_2025's 2025-only classification to every year: a cell's
# ANNUAL MEAN P(source) >= 0.5 assigns the WHOLE cell to Weak-source for
# that year, Weak-sink otherwise -- the same "dichotomous, P* = 0.5" rule
# the panel-D map applies for 2025 alone). This is a coarser, different
# split than annual_source_area above (which area-weights by the FRACTION
# of source-classified months instead of assigning a cell wholesale) --
# written out here, and not just plotted, so claims about the interannual
# stability of the dichotomous source/sink split (e.g. coefficient of
# variation across 2000-2025) can be checked exactly from a CSV rather than
# estimated by eye from the figure.
annual_cell_A_dichotomous <- cp %>%
  group_by(Year, cell) %>%
  summarise(
    area_mha            = first(area_mha),
    mean_source_prob_A   = mean(source_prob_A, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(annual_exchange_class = if_else(mean_source_prob_A >= 0.5, "Weak-source", "Weak-sink"))

annual_dichotomous_area <- annual_cell_A_dichotomous %>%
  group_by(Year, annual_exchange_class) %>%
  summarise(area_mha = sum(area_mha, na.rm = TRUE), .groups = "drop")

annual_dichotomous_area_summary <- annual_dichotomous_area %>%
  group_by(annual_exchange_class) %>%
  summarise(
    mean_area_mha = mean(area_mha, na.rm = TRUE),
    sd_area_mha   = sd(area_mha, na.rm = TRUE),
    cv_pct        = 100 * sd_area_mha / mean_area_mha,
    min_area_mha  = min(area_mha, na.rm = TRUE),
    max_area_mha  = max(area_mha, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(annual_dichotomous_area, file.path(output_dir, "annual_dichotomous_exchange_class_area.csv"), row.names = FALSE)
write.csv(annual_dichotomous_area_summary, file.path(output_dir, "annual_dichotomous_exchange_class_area_summary.csv"), row.names = FALSE)

message(sprintf(
  "Dichotomous exchange-class area (2000-2025): Weak-source mean = %.0f Mha (CV = %.2f%%), Weak-sink mean = %.0f Mha (CV = %.2f%%)",
  annual_dichotomous_area_summary$mean_area_mha[annual_dichotomous_area_summary$annual_exchange_class == "Weak-source"],
  annual_dichotomous_area_summary$cv_pct[annual_dichotomous_area_summary$annual_exchange_class == "Weak-source"],
  annual_dichotomous_area_summary$mean_area_mha[annual_dichotomous_area_summary$annual_exchange_class == "Weak-sink"],
  annual_dichotomous_area_summary$cv_pct[annual_dichotomous_area_summary$annual_exchange_class == "Weak-sink"]
))

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
# Figure 1  Source-probability model  (7.5 × 9.5 in — fits a manuscript page)
# A: Annual mean P(source) map, 2025 (Model A), with a narrow zonal-mean
#    marginal panel on the right sharing A's latitude axis — own row, top
# B: P(source) density by EcoType (Model A) — middle row, left
# C: Latitudinal gradient of P(source) by EcoType (Model A) — middle row, right
# D: Annual exchange-class map, 2025 (Model A, P≥0.5) + inset area time series
#    (E) — own row, bottom
# E: Source/Sink area time series, 2000–2025 — embedded inset within D
# ─────────────────────────────────────────────────────────────────────────────

p1a <- cp %>%
  group_by(cell, EcoType_plot) %>%
  summarise(mean_prob = mean(source_prob_A, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = mean_prob, fill = EcoType_plot, color = EcoType_plot)) +
  geom_density(alpha = 0.28, linewidth = 0.7) +
  geom_vline(xintercept = 0.5, color = "grey30", linetype = "dotted", linewidth = 0.8) +
  scale_fill_manual(values  = ecotype_colors) +
  scale_color_manual(values = ecotype_colors) +
  labs(title = "B.",
       x = "Mean P(source)", y = "Density", fill = NULL, color = NULL) +
  fig_theme +
  theme(legend.position   = c(0.18, 0.75),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size   = unit(0.35, "cm"),
        legend.text       = element_text(size = 6.5))

# Zonal (latitude-only) mean P(source), pooled across all EcoTypes.
lat_bin_width <- 5
zonal_mean_data <- annual_cell_A_2025 %>%
  mutate(lat_bin = floor(y / lat_bin_width) * lat_bin_width + lat_bin_width / 2) %>%
  group_by(lat_bin) %>%
  summarise(
    n_cells = n(),
    mean_p  = mean(mean_source_prob_A, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 5) %>%
  mutate(mean_p = pmin(pmax(mean_p, 0), 1))

# The zonal-mean strip is drawn INSIDE p1b's own panel (extending the
# longitude axis to the right to make room for it) rather than as a
# separate ggplot object glued on with cowplot. coord_equal() gives the
# map a fixed 1:1 aspect ratio, so its rendered panel doesn't necessarily
# fill the full height of whatever cell cowplot allocates to it — that
# mismatch is what broke plot_grid(align="h", axis="tb") alignment even
# after matching title rows. Drawing the strip in the same coordinate
# system as the map guarantees its vertical (latitude) axis is identical
# to the map's, since both share the same un-expanded lat_range.
lon_range <- range(annual_cell_A_2025$x)
lat_range <- range(annual_cell_A_2025$y)
margin_gap   <- diff(lon_range) * 0.03
margin_width <- diff(lon_range) * 0.20
margin_x0 <- lon_range[2] + margin_gap
margin_x1 <- margin_x0 + margin_width
prob_to_x <- function(p) margin_x0 + p * margin_width

zonal_mean_data <- zonal_mean_data %>%
  mutate(x_line = prob_to_x(mean_p))

p1b <- annual_cell_A_2025 %>%
  ggplot(aes(x = x, y = y, fill = mean_source_prob_A)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  annotate("rect", xmin = margin_x0, xmax = margin_x1,
           ymin = lat_range[1], ymax = lat_range[2],
           fill = "white", color = "grey50", linewidth = 0.3) +
  geom_segment(data = data.frame(xv = prob_to_x(0.5)),
               aes(x = xv, xend = xv, y = lat_range[1], yend = lat_range[2]),
               inherit.aes = FALSE, color = "grey40",
               linetype = "dotted", linewidth = 0.5) +
  geom_path(data = zonal_mean_data, aes(x = x_line, y = lat_bin),
            inherit.aes = FALSE, color = "grey15", linewidth = 0.7) +
  annotate("text", x = c(prob_to_x(0), prob_to_x(1)),
           y = lat_range[1] + 0.03 * diff(lat_range),
           label = c("0", "1"), size = 2, color = "grey30") +
  annotate("text", x = margin_x1 - 0.015 * diff(lon_range),
           y = mean(lat_range), label = "Zonal mean P(source)",
           angle = 90, size = 2.3, fontface = "italic", color = "grey20") +
  coord_equal(xlim = c(lon_range[1], margin_x1), ylim = lat_range, expand = FALSE) +
  scale_x_continuous(breaks = pretty(lon_range, n = 5)) +
  prob_scale +
  labs(title = "A.",
       x = "Longitude", y = "Latitude") +
  map_theme +
  theme(legend.position   = c(0.045, 0.28),
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
  labs(title = "D.",
       x = "Longitude", y = "Latitude", fill = NULL) +
  map_theme +
  # Legend embedded inside the panel (in the open ocean south of Africa)
  # instead of an external "bottom" legend. An external legend added a row
  # of height to D that A/B/C/E don't have, which is what threw off the
  # row-to-row size matching with map A — an inside legend keeps D's full
  # plot the same height as a bare map.
  theme(legend.position   = c(0.56, 0.05),
        legend.direction  = "horizontal",
        legend.key.size   = unit(0.4, "cm"),
        legend.text       = element_text(size = 7),
        legend.background = element_rect(fill = "transparent", color = NA))

p1d <- annual_source_area %>%
  ggplot(aes(x = Year, y = area_mha / 100, color = class)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_color_manual(values = c("Weak-sink" = sink_col, "Weak-source" = source_col)) +
  labs(title = "E.",
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
            x = 0.10, y = 0.14, width = 0.211, height = 0.237) +
  draw_label("E.", x = 0.10, y = 0.377, hjust = 0, vjust = 0, size = 8, fontface = "bold")

# C: Latitudinal gradient of P(source) by EcoType — supports the manuscript
# claim that tropical forest/humid grassland cells carry the highest
# P(source), while boreal forest, temperate shrubland, and arid regions are
# net sinks/near-neutral. Built from the same per-cell mean P(source) as
# panel B, binned into 5° latitude bands per EcoType with mean ± SE per
# band. Bands with too few cells (< 5) are dropped as unstable. coord_flip()
# keeps latitude vertical, matching the map orientation in panel A. No
# legend on this panel — the EcoType colors are already keyed in panel B.
# (lat_bin_width is defined earlier, alongside panel A's zonal-mean margin.)
cell_mean_prob <- cp %>%
  group_by(cell, y, EcoType_plot) %>%
  summarise(mean_prob = mean(source_prob_A, na.rm = TRUE), .groups = "drop")

lat_gradient_data <- cell_mean_prob %>%
  mutate(lat_bin = floor(y / lat_bin_width) * lat_bin_width + lat_bin_width / 2) %>%
  group_by(lat_bin, EcoType_plot) %>%
  summarise(
    n_cells = n(),
    mean_p  = mean(mean_prob, na.rm = TRUE),
    se_p    = sd(mean_prob, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 5)

p1e <- lat_gradient_data %>%
  ggplot(aes(x = lat_bin, y = mean_p, color = EcoType_plot, fill = EcoType_plot)) +
  geom_hline(yintercept = 0.5, color = "grey30", linetype = "dotted", linewidth = 0.8) +
  geom_ribbon(aes(ymin = pmax(mean_p - 2 * se_p, 0), ymax = pmin(mean_p + 2 * se_p, 1)),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  coord_flip(ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(-90, 90, by = 30)) +
  scale_color_manual(values = ecotype_colors) +
  scale_fill_manual(values = ecotype_colors) +
  labs(title = "C.",
       x = "Latitude (°)", y = "Mean P(source)", color = NULL, fill = NULL) +
  fig_theme + theme(legend.position = "none")

# Three rows: A (map, with zonal-mean strip drawn inside its own panel)
# alone on top, B+C side by side in the middle, D (map + embedded inset)
# alone on the bottom. A and D are both full-width world maps built from
# the same template/cell footprint, so they share the same fixed
# coord_equal() aspect ratio — giving their rows equal rel_heights makes
# the two maps render at the same size and keeps their longitude
# gridlines lined up down the page. B and C are shrunk to free up that
# space.
fig1 <- plot_grid(
  p1b,
  plot_grid(p1a, p1e, ncol = 2, rel_widths = c(1, 1)),
  p1c_inset,
  ncol = 1, rel_heights = c(1.1, 0.75, 1.1)
)

# 7.5 x 9.5 in fits within a standard manuscript page (single full-width
# figure, letter/A4 with margins) at 300 dpi.
ggsave(file.path(figure_dir, "FIGURE4_source_probability.png"),
  fig1, width = 7.5, height = 9.5, units = "in", dpi = 300, bg = "white")
message("Fig 1 written.")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2  Seasonal flux cycle by Ecosystem Type  (7.5 × 6.0 in)
# (three approaches; model-performance figures — calibration, Stage 1 skill,
# variable importance, obs-vs-fitted magnitude — live in
# 12b_Model_FIGURES_RF.R instead)
# ─────────────────────────────────────────────────────────────────────────────

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
  labs(title = "",
       x = "Month", y = "Flux (g C m⁻² mo⁻¹)",
       color = NULL, linetype = "Approach") +
  fig_theme + guides(color = "none")

ggsave(file.path(figure_dir, "Fig2_seasonal_flux_cycle.png"),
  p2c, width = 7.5, height = 6.0, units = "in", dpi = 300, bg = "white")
message("Fig 2 (seasonal flux cycle) written.")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 (original, pre-split layout)  Fig2_magnitude_models.png  (8.0 × 8.0 in)
# A: Obs vs OOB-fitted — Weak-sink        (from 12_SourceProp_MagnitudeModels.R)
# B: Obs vs OOB-fitted — Weak-source      (from 12_SourceProp_MagnitudeModels.R)
# C: Seasonal flux cycle by Ecosystem Type (from this script's own upscaling output)
#
# Kept for anything downstream that still expects this exact filename/layout.
# The split-out versions (Fig2_magnitude_model_fit.png in 12b, just A/B; and
# Fig2_seasonal_flux_cycle.png above, just C) remain the primary outputs —
# this reassembles both into the original combined figure without
# duplicating any fitting logic, by reading 12's already-written CSV.
# ─────────────────────────────────────────────────────────────────────────────

magnitude_fitted_values_file <- file.path(output_dir, "magnitude_model_fitted_values.csv")
if (file.exists(magnitude_fitted_values_file)) {

  magnitude_fitted_values <- read.csv(magnitude_fitted_values_file) %>%
    mutate(magnitude_model = recode(magnitude_model,
                                    "Weak sink"   = "Weak-sink",
                                    "Weak source" = "Weak-source"))

  # Duplicated from 12_SourceProp_MagnitudeModels.R / 12b_Model_FIGURES_RF.R
  # (kept in sync manually, same pattern used between 12/12b/13/14/19).
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
      labs(title = paste0(panel_tag),
           x = "Observed (g C m⁻² mo⁻¹)",
           y = "Fitted (g C m⁻² mo⁻¹, OOB)",
           color = NULL) +
      fig_theme
  }

  p2a <- make_mag_plot("Weak-sink",   "A")
  p2b <- make_mag_plot("Weak-source", "B")

  shared_ecotype_legend <- get_legend(
    p2a + theme(legend.position = "top",
                legend.key.size = unit(0.4, "cm"),
                legend.text     = element_text(size = 8))
  )

  p2c_labeled <- p2c + labs(title = paste0("C. ", p2c$labels$title))

  fig2_combined <- plot_grid(
    shared_ecotype_legend,
    plot_grid(p2a + theme(legend.position = "none"),
              p2b + theme(legend.position = "none"),
              ncol = 2),
    p2c_labeled,
    ncol = 1, rel_heights = c(0.07, 1, 1.4)
  )

  ggsave(file.path(figure_dir, "FIGURE5_magnitude_models.png"),
    fig2_combined, width = 8.0, height = 8.0, units = "in", dpi = 300, bg = "white")
  message("Fig 2 (original combined magnitude_models layout) written.")

} else {
  message("Skipping Fig2_magnitude_models.png: ", magnitude_fitted_values_file,
          " not found — run 12_SourceProp_MagnitudeModels.R first.")
}

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
  labs(title = "A.",
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
  labs(title = "B. ",
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
  labs(title = "C.",
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
  "A. ",
  flux_map_abs_max
)

p4b <- make_flux_map(
  cell_flux_2025, "annual_balanced_gC_m2_yr",
  "B.",
  flux_map_abs_max
)

p4c <- make_flux_map(
  cell_flux_2025, "annual_constrained_gC_m2_yr",
  "C. ",
  flux_map_abs_max
)

fig4 <- plot_grid(p4a, p4b, p4c, ncol = 1)

ggsave(file.path(figure_dir, "FIGURE6_spatial_maps.png"),
  fig4, width = 9.0, height = 9.5, units = "in", dpi = 300, bg = "white")
message("Fig 4 written.")

# ─────────────────────────────────────────────────────────────────────────────
# Standalone budget time series (for external use)
# ─────────────────────────────────────────────────────────────────────────────

ggsave(file.path(figure_dir, "annual_budget_time_series_2000_2025.png"),
  p3a + labs(title = NULL),
  width = 7.0, height = 3.5, units = "in", dpi = 300, bg = "white")

message("All RF figures written to: ", figure_dir)

