# Model-performance figures — Stage 1 (source-probability) and Stage 2
# (flux-magnitude) Random Forest models.
#
# Standalone figure script covering ONLY model fit/skill diagnostics:
# probability calibration, Stage 1 classification skill, Stage 1 variable
# importance, and Stage 2 observed-vs-OOB-fitted magnitude. Loads outputs
# already written by 12_SourceProp_MagnitudeModels.R — no dependency on the
# spatial projection grid (13_Global_SpatialUpscalingRF.R / cell_preds).
#
# For upscaling-OUTPUT figures instead (P(source) maps, exchange-class maps,
# seasonal flux cycle, annual budget, spatial flux maps), see
# 14_Global_SpatialUpscalingFiguresRF.R.
#
# Produces:
#   FigS1_model_performance.png  — Calibration, Stage 1 skill, variable importance
#   Fig2_magnitude_model_fit.png — Obs vs OOB-fitted, Weak-sink / Weak-source
#
# Prerequisites (written by 12_SourceProp_MagnitudeModels.R):
#   <output_dir>/probability_calibration_skill.csv
#   <output_dir>/comparison_class_skill.csv
#   <output_dir>/rf_class_variable_importance.csv
#   <output_dir>/magnitude_model_fitted_values.csv

library(tidyverse)
library(cowplot)

# ─────────────────────────────────────────────────────────────────────────────
# Paths  (must match 12_SourceProp_MagnitudeModels.R's output_dir)
# ─────────────────────────────────────────────────────────────────────────────

# When sourced from 12_SourceProp_MagnitudeModels.R, output_dir and
# figure_dir are already defined — use them directly.
if (!exists("output_dir")) {
  rf_dir     <- Sys.getenv("MONTHLY_RF_DIR",
    unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF")
  output_dir <- file.path(rf_dir, "OUTPUT")
}
if (!exists("figure_dir")) {
  figure_dir <- file.path(dirname(output_dir), "FIGURES")
}
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  file.path(output_dir, "probability_calibration_skill.csv"),
  file.path(output_dir, "comparison_class_skill.csv"),
  file.path(output_dir, "rf_class_variable_importance.csv"),
  file.path(output_dir, "magnitude_model_fitted_values.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0)
  stop("Missing required inputs (run 12_SourceProp_MagnitudeModels.R first): ",
       paste(missing_files, collapse = ", "))

# ─────────────────────────────────────────────────────────────────────────────
# Load saved outputs  (skipped when already in memory from 12)
# ─────────────────────────────────────────────────────────────────────────────

probability_calibration_skill <- read.csv(file.path(output_dir, "probability_calibration_skill.csv"))
comparison_class <- read.csv(file.path(output_dir, "comparison_class_skill.csv"))
rf_class_importance <- read.csv(file.path(output_dir, "rf_class_variable_importance.csv"))
magnitude_fitted_values <- read.csv(file.path(output_dir, "magnitude_model_fitted_values.csv")) %>%
  mutate(magnitude_model = recode(magnitude_model,
                                  "Weak sink"   = "Weak-sink",
                                  "Weak source" = "Weak-source"))

# ─────────────────────────────────────────────────────────────────────────────
# Shared aesthetics
# ─────────────────────────────────────────────────────────────────────────────

ecotype_colors <- c(
  "Forest"    = "#1B7837",
  "Grassland" = "#D9B86C",
  "Shrubland" = "#C2A5CF",
  "Arid"      = "#D95F02"
)

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
  labs(title = "A.",
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
  labs(title = "B.",
       x = "Value", y = NULL) +
  fig_theme + theme(legend.position = "none")

# Variable importance — top 8 predictors
pS1c <- rf_class_importance %>%
  arrange(desc(importance)) %>%
  head(8) %>%
  mutate(predictor = fct_reorder(predictor, importance)) %>%
  ggplot(aes(x = importance, y = predictor)) +
  geom_col(fill = "#009688", alpha = 0.85) +
  labs(title = "C. ",
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
# Figure 2  Magnitude model fit  (8.0 × 4.5 in)
# A: Obs vs OOB-fitted — Weak-sink
# B: Obs vs OOB-fitted — Weak-source
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
    labs(#title = paste0(panel_tag, ". ", model_label, " model  (n = ", n_obs, ")"),
      title = paste0(panel_tag),
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

fig2 <- plot_grid(
  shared_ecotype_legend,
  plot_grid(p2a + theme(legend.position = "none"),
            p2b + theme(legend.position = "none"),
            ncol = 2),
  ncol = 1, rel_heights = c(0.1, 1)
)

ggsave(file.path(figure_dir, "Fig2_magnitude_model_fit.png"),
  fig2, width = 8.0, height = 4.5, units = "in", dpi = 300, bg = "white")
message("Fig 2 (magnitude model fit) written.")

message("Model-performance figures written to: ", figure_dir)

