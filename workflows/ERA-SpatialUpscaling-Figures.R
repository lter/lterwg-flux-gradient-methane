# ERA-SpatialUpscaling-Figures.R
#
# Standalone figure script for the ERA5-Land spatial upscaling pipeline.
# Loads pre-computed outputs from ERA-SpatialUpscaling-Monthly.R and
# regenerates all four publication figures without re-fitting any models
# or re-projecting any raster layers.
#
# Prerequisites (written by the main script on first run):
#   <output_dir>/monthly_cell_predictions_2000_2025.rds
#   <output_dir>/era5_template.tif
#   <output_dir>/model_parameters.csv
#   <output_dir>/annual_expected_flux_2000_2025.csv
#   <output_dir>/magnitude_model_fitted_values.csv
#   <output_dir>/probability_calibration_skill.csv
#   <output_dir>/annual_source_sink_map_cells_2025.csv
#   <output_dir>/annual_source_sink_area_2000_2025.csv

library(tidyverse)
library(terra)
library(cowplot)

# ─────────────────────────────────────────────────────────────────────────────
# Paths  (must match ERA-SpatialUpscaling-Monthly.R)
# ─────────────────────────────────────────────────────────────────────────────

# When sourced from ERA-SpatialUpscaling-Monthly.R, output_dir and figure_dir
# are already defined in the calling environment — use them directly and skip
# the directory detection below.
if (!exists("output_dir")) {
  spatial_dir <- Sys.getenv(
    "MONTHLY_UPSCALING_DIR",
    unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly"
  )

  # Which trial to use — set UPSCALING_TRIAL_SUBDIR to a specific sub-directory
  # name (e.g. "ERA5_land_MODIS_WAD2M_lmer_2000_2025"), or leave unset to
  # auto-select the most recently modified ERA5 trial directory.
  trial_subdir <- Sys.getenv("UPSCALING_TRIAL_SUBDIR", unset = "")

  if (nchar(trial_subdir) == 0) {
    trial_dirs <- list.dirs(spatial_dir, recursive = FALSE, full.names = TRUE)
    trial_dirs <- trial_dirs[grepl("ERA5", basename(trial_dirs))]
    if (length(trial_dirs) == 0)
      stop("No ERA5 trial directories found under ", spatial_dir,
           ". Run ERA-SpatialUpscaling-Monthly.R first.")
    trial_subdir <- basename(trial_dirs[which.max(file.info(trial_dirs)$mtime)])
    message("Auto-selected trial: ", trial_subdir)
  }

  output_dir <- file.path(spatial_dir, trial_subdir)
}

if (!exists("figure_dir")) {
  figure_dir <- file.path(output_dir, "figures")
}
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Load saved outputs  (skipped when already in memory from the main script)
# ─────────────────────────────────────────────────────────────────────────────

if (!exists("monthly_cell_predictions")) {
  message("Loading saved outputs from ", output_dir)
  monthly_cell_predictions <- readRDS(
    file.path(output_dir, "monthly_cell_predictions_2000_2025.rds"))
} else {
  message("Using monthly_cell_predictions already in memory.")
}

if (!exists("template") || !inherits(template, "SpatRaster")) {
  template <- rast(file.path(output_dir, "era5_template.tif"))
}

# Scalar params always re-read from the saved CSV (cheap, ensures consistency)
params <- read.csv(file.path(output_dir, "model_parameters.csv"))
source_probability_threshold            <- params$source_probability_threshold
source_probability_threshold_percentile <- params$source_probability_threshold_percentile
gmb_soil_sink_tg_ch4_yr                 <- params$gmb_soil_sink_tg_ch4_yr

# Tabular outputs: read from disk when not already in memory
if (!exists("annual_expected_flux_2000_2025"))
  annual_expected_flux_2000_2025 <- read.csv(
    file.path(output_dir, "annual_expected_flux_2000_2025.csv"))
if (!exists("magnitude_fitted_values"))
  magnitude_fitted_values <- read.csv(
    file.path(output_dir, "magnitude_model_fitted_values.csv"))
if (!exists("probability_calibration_skill"))
  probability_calibration_skill <- read.csv(
    file.path(output_dir, "probability_calibration_skill.csv"))
if (!exists("annual_cell_class_2025"))
  annual_cell_class_2025 <- read.csv(
    file.path(output_dir, "annual_source_sink_map_cells_2025.csv"))
if (!exists("annual_source_sink_area"))
  annual_source_sink_area <- read.csv(
    file.path(output_dir, "annual_source_sink_area_2000_2025.csv"))
if (!exists("class_probability_model_skill"))
  class_probability_model_skill <- read.csv(
    file.path(output_dir, "class_probability_model_skill.csv"))

# ─────────────────────────────────────────────────────────────────────────────
# Constants and helpers
# ─────────────────────────────────────────────────────────────────────────────

source_probability_thresholds <- seq(0.50, 0.95, by = 0.05)
high_threshold                <- 0.95

# Convert g C m⁻² month⁻¹ × area (Mha) to Tg CH4 yr⁻¹
# g C → g CH4: × (16/12); g CH4 → Tg CH4: × 1e6 km² × 1e8 cm²/km² × 1e-15 Tg/g
gC_m2_yr_to_tg_ch4 <- function(flux_gC_m2_yr, area_mha) {
  flux_gC_m2_yr * area_mha * 0.01333333333333333  # g C m⁻² × Mha × (16/12) × 1e10 m²/Mha × 1e-12 Tg/g
}

model_colors <- c(
  "Binary: hard-threshold + lmer conditional magnitude" = "#1B9E77",
  "Continuous: P(source) × source mag + (1-P) × sink mag" = "#7570B3"
)
binary_color     <- model_colors["Binary: hard-threshold + lmer conditional magnitude"]
continuous_color <- model_colors["Continuous: P(source) × source mag + (1-P) × sink mag"]

ecotype_colors <- c(
  "Forest"    = "#1B9E77",
  "Grassland" = "#D95F02",
  "Shrubland" = "#7570B3",
  "Arid"      = "#C49A2C"
)

# ─────────────────────────────────────────────────────────────────────────────
# Shared ggplot theme
# ─────────────────────────────────────────────────────────────────────────────

fig_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor    = element_blank(),
    plot.title          = element_text(face = "bold", size = 12),
    plot.tag            = element_text(face = "bold", size = 12),
    legend.position     = "bottom",
    legend.title        = element_text(size = 11),
    legend.text         = element_text(size = 10),
    strip.background    = element_rect(fill = "grey92"),
    strip.text          = element_text(size = 11),
    axis.title          = element_text(size = 11),
    axis.text           = element_text(size = 10)
  )

map_theme <- fig_theme +
  theme(panel.grid = element_blank(), legend.key.width = unit(1.2, "cm"))

sink_col   <- "#2166AC"
source_col <- "#B2182B"

flux_scale <- function(lim, name = "g C m⁻² yr⁻¹") {
  scale_fill_gradient2(
    low = sink_col, mid = "#D4C08A", high = source_col,
    midpoint = 0, limits = c(-lim, lim),
    breaks = pretty(c(-lim, lim), n = 5),
    name = name, na.value = "transparent"
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

mcp <- monthly_cell_predictions %>%
  mutate(EcoType_plot = if_else(is_arid == 1L, "Arid", as.character(EcoType)))

# Annual classification raster
annual_class_raster <- template
names(annual_class_raster) <- "annual_exchange_class"
values(annual_class_raster) <- NA_real_
annual_class_raster[annual_cell_class_2025$cell] <- if_else(
  annual_cell_class_2025$annual_exchange_class == "Weak source", 2, 1)
annual_class_map_df <- as.data.frame(annual_class_raster, xy = TRUE, na.rm = TRUE) %>%
  mutate(annual_exchange_class = factor(annual_exchange_class,
    levels = c(1, 2), labels = c("Weak sink", "Weak source")))

# Cell-level annual fluxes 2025
cell_flux_2025 <- monthly_cell_predictions %>%
  filter(Year == 2025) %>%
  mutate(flux_cont = source_probability * predicted_source_flux_gC_m2_month +
           (1 - source_probability) * predicted_sink_flux_gC_m2_month) %>%
  group_by(cell, x, y) %>%
  summarise(
    annual_binary_gC_m2_yr     = sum(expected_flux_gC_m2_month, na.rm = TRUE),
    annual_continuous_gC_m2_yr = sum(flux_cont,                 na.rm = TRUE),
    .groups = "drop"
  )
flux_map_abs_max <- max(abs(range(
  c(cell_flux_2025$annual_binary_gC_m2_yr, cell_flux_2025$annual_continuous_gC_m2_yr),
  na.rm = TRUE)))

# Seasonal P(source) 2025
prob_seasonal_2025 <- monthly_cell_predictions %>%
  filter(Year == 2025) %>%
  mutate(season = case_when(
    month %in% c(12, 1, 2) ~ "DJF",
    month %in% c(6, 7, 8)  ~ "JJA",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(season)) %>%
  group_by(cell, x, y, season) %>%
  summarise(mean_source_probability = mean(source_probability, na.rm = TRUE),
            .groups = "drop")

# Threshold sensitivity
threshold_sensitivity <- map_dfr(source_probability_thresholds, function(thresh) {
  monthly_cell_predictions %>%
    mutate(
      flux_t = if_else(source_probability >= thresh,
                       predicted_source_flux_gC_m2_month,
                       predicted_sink_flux_gC_m2_month),
      tg_t   = gC_m2_yr_to_tg_ch4(flux_t, area_mha)
    ) %>%
    group_by(Year) %>%
    summarise(annual_tg = sum(tg_t, na.rm = TRUE), .groups = "drop") %>%
    mutate(threshold = thresh)
})
thresh_summary <- threshold_sensitivity %>%
  group_by(threshold) %>%
  summarise(mean_tg = mean(annual_tg), min_tg = min(annual_tg), max_tg = max(annual_tg),
            .groups = "drop")

# Budget by EcoType
budget_by_ecotype_year <- mcp %>%
  mutate(tg_cell = gC_m2_yr_to_tg_ch4(expected_flux_gC_m2_month, area_mha)) %>%
  group_by(Year, EcoType_plot) %>%
  summarise(annual_tg = sum(tg_cell, na.rm = TRUE), .groups = "drop")

# Seasonal flux by EcoType
flux_seasonal <- mcp %>%
  mutate(
    flux_bin  = expected_flux_gC_m2_month,
    flux_cont = source_probability * predicted_source_flux_gC_m2_month +
      (1 - source_probability) * predicted_sink_flux_gC_m2_month
  ) %>%
  group_by(month, EcoType_plot) %>%
  summarise(
    Binary     = weighted.mean(flux_bin,  area_mha, na.rm = TRUE),
    Continuous = weighted.mean(flux_cont, area_mha, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(Binary, Continuous), names_to = "approach", values_to = "flux")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1  Source-probability model  (6.5 × 8.0 in)
# A: P(source) distribution by EcoType
# B: Annual mean P(source) map, 2025
# C: Annual source/sink area time series, 2000-2025
# ─────────────────────────────────────────────────────────────────────────────

p1a <- mcp %>%
  group_by(cell, EcoType_plot) %>%
  summarise(mean_prob = mean(source_probability, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = mean_prob, fill = EcoType_plot, color = EcoType_plot)) +
  geom_density(alpha = 0.28, linewidth = 0.7) +
  geom_vline(xintercept = source_probability_threshold, color = "grey30",
             linetype = "dotted", linewidth = 0.8) +
  scale_fill_manual(values  = ecotype_colors) +
  scale_color_manual(values = ecotype_colors) +
  labs(title = "A. P(source) by Ecosystem Type",
       x = "Mean P(source)", y = "Density",
       fill = NULL, color = NULL) +
  fig_theme +
  theme(legend.position      = c(0.18, 0.75),
        legend.background    = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size      = unit(0.35, "cm"),
        legend.text          = element_text(size = 6.5))

p1b <- annual_cell_class_2025 %>%
  ggplot(aes(x = x, y = y, fill = mean_source_probability)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  coord_equal(expand = FALSE) +
  prob_scale +
  labs(title = "B. Continuous P(source), 2025", x = "Longitude", y = "Latitude") +
  map_theme +
  theme(legend.position   = c(0.08, 0.28),
        legend.direction  = "vertical",
        legend.key.width  = unit(0.4, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.text       = element_text(size = 6.5),
        legend.title      = element_text(size = 7),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA))

p1c <- annual_class_map_df %>%
  ggplot(aes(x = x, y = y, fill = annual_exchange_class)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  coord_equal(expand = FALSE) +
  scale_fill_manual(values = c("Weak sink" = sink_col, "Weak source" = source_col),
                    labels = c("Weak sink" = "Weak-sink", "Weak source" = "Weak-source"),
                    na.value = "transparent") +
  labs(title = "C. Binary P(source), 2025", x = "Longitude", y = "Latitude",
       fill = NULL) +
  map_theme +
  theme(plot.margin      = margin(0, 0, 0, 0),
        legend.position  = "top",
        legend.direction = "horizontal",
        legend.key.size  = unit(0.4, "cm"),
        legend.text      = element_text(size = 7))

p1d <- annual_source_sink_area %>%
  ggplot(aes(x = Year, y = mean_monthly_area_mha / 100, color = selected_exchange_class)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_color_manual(values = c("Weak sink" = sink_col, "Weak source" = source_col),
                     labels = c("Weak sink" = "Weak-sink", "Weak source" = "Weak-source")) +
  labs(title = "D. Annual Source/Sink Area, 2000-2025",
       x = "Year", y = "Area (100 Mha)", color = NULL) +
  fig_theme

# Embed D (time series) as inset in bottom-left of C (exchange class map)
p1c_inset <- ggdraw() +
  draw_plot(p1c) +
  draw_plot(p1d + theme(plot.title = element_blank(),
                        axis.title = element_text(size = 7),
                        axis.text  = element_text(size = 6),
                        legend.position = "none",
                        plot.background = element_rect(fill = "transparent", color = NA),
                        panel.background = element_rect(fill = "transparent")),
            x = 0.13, y = 0.08, width = 0.234, height = 0.263) +
  draw_label("D.", x = 0.13, y = 0.345, hjust = 0, vjust = 0, size = 8, fontface = "bold")

fig1 <- plot_grid(
  plot_grid(p1a, p1b, ncol = 2, rel_widths = c(0.7, 1.3)),
  p1c_inset,
  ncol = 1, rel_heights = c(1, 1.1)
)

ggsave(file.path(figure_dir, "Fig1_source_probability.png"),
  fig1, width = 8.0, height = 7.0, units = "in", dpi = 300, bg = "white")

message("Fig 1 written.")

# ─────────────────────────────────────────────────────────────────────────────
# Supplement Fig S1  Probability model performance  (6.5 × 7.0 in, 2×2)
# A: Calibration curve
# B: Classification performance metrics
# C: Annual mean P(source) map, 2025
# D: P(source) distribution by EcoType
# ─────────────────────────────────────────────────────────────────────────────

ps <- class_probability_model_skill

pS1a <- probability_calibration_skill %>%
  ggplot(aes(x = mean_calibrated_source_probability, y = observed_source_fraction)) +
  geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = "dashed",
              linewidth = 0.5) +
  geom_point(aes(size = n), color = source_col, alpha = 0.85) +
  geom_line(color = source_col, linewidth = 0.8) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_size_continuous(range = c(1.5, 5)) +
  labs(title = "A. Probability Calibration",
       x = "Calibrated P(source)", y = "Observed fraction", size = "n months") +
  fig_theme

pS1b <- tibble(
  metric = c("AUC", "Tjur R²", "Dev. explained", "Accuracy",
             "Sensitivity", "Specificity"),
  value  = c(ps$auc, ps$tjur_r2, ps$deviance_explained,
             ps$accuracy, ps$sensitivity_source, ps$specificity_sink),
  group  = c("Discrimination", "Discrimination", "Discrimination",
             "Classification", "Classification", "Classification")
) %>%
  mutate(metric = factor(metric,
    levels = c("AUC", "Tjur R²", "Dev. explained",
               "Accuracy", "Sensitivity", "Specificity"))) %>%
  ggplot(aes(x = value, y = metric, fill = group)) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f", value)), hjust = -0.1, size = 3.2) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_x_continuous(limits = c(0, 1.12), expand = c(0, 0)) +
  scale_fill_manual(values = c("Discrimination" = "#7570B3",
                                "Classification" = "#1B9E77")) +
  labs(title = sprintf("B. Model Performance  (threshold = %.2f)", ps$threshold),
       x = "Value", y = NULL, fill = NULL) +
  fig_theme + theme(legend.position = "bottom")

pS1c <- annual_cell_class_2025 %>%
  ggplot(aes(x = x, y = y, fill = mean_source_probability)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  coord_equal(expand = FALSE) +
  prob_scale +
  labs(title = "C. Annual Mean P(source), 2025", x = "Longitude", y = "Latitude") +
  map_theme

pS1d <- mcp %>%
  group_by(cell, EcoType_plot) %>%
  summarise(mean_prob = mean(source_probability, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = mean_prob, fill = EcoType_plot, color = EcoType_plot)) +
  geom_density(alpha = 0.28, linewidth = 0.7) +
  geom_vline(xintercept = source_probability_threshold, color = "grey30",
             linetype = "dotted", linewidth = 0.8) +
  annotate("text", x = min(source_probability_threshold + 0.06, 0.90), y = Inf,
           label = sprintf("threshold = %.2f", source_probability_threshold),
           hjust = 0, vjust = 1.4, size = 3.0, color = "grey30") +
  scale_fill_manual(values  = ecotype_colors) +
  scale_color_manual(values = ecotype_colors) +
  labs(title = "D. P(source) Distribution by Ecosystem Type",
       x = "Mean P(source)", y = "Density", fill = NULL, color = NULL) +
  fig_theme

figS1 <- plot_grid(
  plot_grid(pS1a, pS1b, ncol = 2),
  plot_grid(pS1c, pS1d, ncol = 2),
  ncol = 1, rel_heights = c(1, 1)
)

ggsave(file.path(figure_dir, "FigS1_probability_model_performance.png"),
  figS1, width = 6.5, height = 7.0, units = "in", dpi = 300, bg = "white")

message("Fig S1 (supplement) written.")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2  Magnitude models  (6.5 × 8.0 in)
# A: Sink model obs vs fitted | B: Source model obs vs fitted
# C: Seasonal flux cycle by EcoType (2×2 facet, binary vs continuous)
# ─────────────────────────────────────────────────────────────────────────────

make_mag_plot <- function(model_label, panel_tag) {
  n_obs <- sum(magnitude_fitted_values$magnitude_model == model_label)
  magnitude_fitted_values %>%
    filter(magnitude_model == model_label, !is.na(EcoType)) %>%
    ggplot(aes(x = monthly_flux_gC_m2_month, y = fitted_flux_gC_m2_month,
               color = EcoType)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "grey35", linetype = "dashed",
                linewidth = 0.5) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = ecotype_colors) +
    labs(title = paste0(panel_tag, ". ", model_label, " model  (n = ", n_obs, ")"),
         x = "Observed (g C m⁻² mo⁻¹)", y = "Fitted (g C m⁻² mo⁻¹)",
         color = NULL) +
    fig_theme
}

p2a <- make_mag_plot("Weak sink",   "A")
p2b <- make_mag_plot("Weak source", "B")

p2c <- flux_seasonal %>%
  ggplot(aes(x = month, y = flux, color = EcoType_plot, linetype = approach)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.85) + geom_point(size = 1.3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = ecotype_colors) +
  scale_linetype_manual(values = c("Binary" = "solid", "Continuous" = "dashed")) +
  facet_wrap(~EcoType_plot, scales = "free_y", ncol = 2) +
  labs(title = "C. Seasonal Flux Cycle by Ecosystem Type",
       x = "Month", y = "Flux (g C m⁻² mo⁻¹)",
       color = NULL, linetype = "Approach") +
  fig_theme +
  guides(color = "none")

# Shared EcoType legend extracted from p2a
shared_ecotype_legend <- get_legend(
  p2a + theme(legend.position = "top",
              legend.key.size = unit(0.4, "cm"),
              legend.text = element_text(size = 8))
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
# Figure 3  Budget  (6.5 × 8.5 in)
# A: Annual net exchange, both approaches
# B: Binary budget stacked by EcoType
# C: Threshold sensitivity (data-driven + 0.95 annotated)
# ─────────────────────────────────────────────────────────────────────────────

p3a <- annual_expected_flux_2000_2025 %>%
  ggplot(aes(x = Year, y = annual_net_exchange_tg_ch4_yr, color = magnitude_model)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -35,
           fill = sink_col, alpha = 0.10) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = sink_col,
             linetype = "dashed", linewidth = 0.7) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  annotate("text", x = 2025, y = gmb_soil_sink_tg_ch4_yr - 1.5,
           label = "GMB sink (–35)", hjust = 1, size = 3.2, color = sink_col) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_color_manual(values = model_colors) +
  labs(title = "A. Annual Net Exchange, 2000–2025",
       x = "Year", y = "Net exchange (Tg CH₄ yr⁻¹)", color = NULL) +
  fig_theme

p3b <- budget_by_ecotype_year %>%
  ggplot(aes(x = Year, y = annual_tg, fill = EcoType_plot)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_area(alpha = 0.82, position = "stack") +
  geom_line(aes(color = EcoType_plot), position = "stack",
            linewidth = 0.35, show.legend = FALSE) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_fill_manual(values  = ecotype_colors) +
  scale_color_manual(values = ecotype_colors) +
  labs(title = "B. Binary Budget by Ecosystem Type",
       x = "Year", y = "Net exchange (Tg CH₄ yr⁻¹)", fill = NULL) +
  fig_theme

p3c <- threshold_sensitivity %>%
  ggplot(aes(x = threshold, y = annual_tg, group = Year)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -36, ymax = -35,
           fill = sink_col, alpha = 0.10) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = gmb_soil_sink_tg_ch4_yr, color = sink_col,
             linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = source_probability_threshold, color = "grey40",
             linetype = "dotted", linewidth = 0.9) +
  geom_vline(xintercept = high_threshold, color = "grey40",
             linetype = "longdash", linewidth = 0.7) +
  geom_line(alpha = 0.15, linewidth = 0.4, color = as.character(binary_color)) +
  geom_ribbon(data = thresh_summary,
              aes(y = mean_tg, ymin = min_tg, ymax = max_tg, group = 1),
              fill = as.character(binary_color), alpha = 0.17, color = NA) +
  geom_line(data = thresh_summary, aes(y = mean_tg, group = 1),
            color = as.character(binary_color), linewidth = 1.1) +
  geom_point(data = thresh_summary, aes(y = mean_tg, group = 1),
             color = as.character(binary_color), size = 2) +
  annotate("text", x = source_probability_threshold + 0.005, y = Inf,
           label = sprintf("data-driven\n%.2f", source_probability_threshold),
           hjust = 0, vjust = 1.2, size = 3.1, color = "grey30") +
  annotate("text", x = high_threshold + 0.005, y = Inf,
           label = "high\n0.95", hjust = 0, vjust = 1.2, size = 3.1, color = "grey30") +
  scale_x_continuous(breaks = source_probability_thresholds) +
  labs(title = "C. Binary Approach: Threshold Sensitivity",
       x = "P(source) threshold", y = "Net exchange (Tg CH₄ yr⁻¹)") +
  fig_theme

fig3 <- plot_grid(p3a, p3b, p3c, ncol = 1, rel_heights = c(1, 1, 1))

ggsave(file.path(figure_dir, "Fig3_budget.png"),
  fig3, width = 6.5, height = 8.5, units = "in", dpi = 300, bg = "white")

message("Fig 3 written.")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4  Spatial flux maps  (6.5 × 9.0 in)
# A: Annual exchange class, 2025  (binary)
# B: Annual net flux, 2025  (binary)
# C: Annual net flux, 2025  (continuous)
# ─────────────────────────────────────────────────────────────────────────────

p4a <- annual_class_map_df %>%
  ggplot(aes(x = x, y = y, fill = annual_exchange_class)) +
  geom_tile(width = res(template)[1], height = res(template)[2]) +
  coord_equal(expand = FALSE) +
  scale_fill_manual(values = c("Weak sink" = sink_col, "Weak source" = source_col),
                    na.value = "transparent") +
  labs(title = "A. Annual Exchange Class, 2025  (binary)",
       x = "Longitude", y = "Latitude", fill = NULL) +
  map_theme

p4b <- make_flux_map(cell_flux_2025, "annual_binary_gC_m2_yr",
                     "A. Annual Net Flux, 2025  (binary)", flux_map_abs_max)

p4c <- make_flux_map(cell_flux_2025, "annual_continuous_gC_m2_yr",
                     "B. Annual Net Flux, 2025  (continuous)", flux_map_abs_max)

fig4 <- plot_grid(p4b, p4c, ncol = 1, rel_heights = c(1, 1))

ggsave(file.path(figure_dir, "Fig4_spatial_maps.png"),
  fig4, width = 6.5, height = 6.5, units = "in", dpi = 300, bg = "white")

message("Fig 4 written.")

# ─────────────────────────────────────────────────────────────────────────────
# Standalone outputs also used elsewhere
# ─────────────────────────────────────────────────────────────────────────────

ggsave(file.path(figure_dir, "annual_expected_flux_time_series_2000_2025.png"),
  p3a + labs(title = NULL), width = 6.5, height = 3.5, units = "in", dpi = 300)

message("All figures written to ", figure_dir)
