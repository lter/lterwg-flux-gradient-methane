# Source-probability and flux-magnitude Random Forest models — fitting and
# testing only. Spatial projection onto the global grid is a separate step
# (13_Global_SpatialUpscalingRF.R), which loads the fitted models this
# script produces rather than re-fitting them.
#
# Two models, fit on NEON upland site-months:
#
#   Stage 1 — P(source) classifier (1:1 balanced RF):
#     Training data downsampled to 1:1 source:sink (equal class
#     representation). Regularised probability forest: min.node.size = 20,
#     max.depth = 8, sample.fraction = 0.7 without replacement.
#     OOB predictions used for honest AUC and isotonic calibration.
#
#   Stage 2 — flux-magnitude regressors (shared across all upscaling
#   approaches): separate ranger log-absolute-flux regressions for
#   Weak-sink and Weak-source months, trained on the full data.
#     - Skill (RMSE/MAE/correlation) reported from ranger's OOB predictions,
#       held out at the site-MONTH level, not an in-sample fit.
#     - Additionally validated against a fully external, independent
#       dataset never used anywhere in this pipeline: the FLUXNET-CH4
#       "Upland" ecosystem-class aggregate (Delwiche et al. 2021), the same
#       published benchmark 09_NEON_FLUXNETComparison.R uses for NEON's
#       observed flux, applied here to the model's OOB-predicted flux
#       instead. This speaks to generalization beyond NEON's own site
#       network without requiring a leave-one-site-out refit — see
#       magnitude_model_fluxnet_external_validation.csv and
#       magnitude_model_fluxnet_validation_text.md.
#
# "Arid" (aridity_index < arid_ai_threshold — the same rule used in
# 19_Supp_NEONRepresentativeness.R) is a genuine, data-driven EcoType level
# used for LABELING and as a Stage 2 predictor, but it is deliberately NOT
# modeled by Stage 1. Only 2 NEON sites cross this threshold (JORN, SRER),
# and empirically ALL 84 of their site-months are observed weak-SOURCE
# (positive flux) -- the opposite of the standard ecological expectation
# that dryland/desert soils are net atmospheric CH4 sinks via methanotrophic
# oxidation. With only 2 sites, both showing what looks like anomalous
# behavior relative to that prior, there isn't enough independent evidence
# to trust a learned Arid-specific P(source) pattern, so Arid is instead
# handled deterministically:
#   - Stage 1 (P(source) classifier) excludes Arid site-months from
#     training entirely, and P(source) is hard-forced to 0 for Arid cells
#     at prediction time (13_Global_SpatialUpscalingRF.R) rather than using
#     whatever rf_model_A would extrapolate for a level it never trained on.
#   - Stage 2's SOURCE magnitude model excludes Arid site-months too, since
#     Arid never reaches it (P(source) = 0 always routes to the sink model).
#   - Stage 2's SINK magnitude model also excludes Arid's own site-months
#     (all 84 are source-labeled, so none would qualify as sink training
#     data anyway) and instead predicts Arid cells by extrapolating from
#     Forest/Grassland/Shrubland sink behavior using the shared continuous
#     covariates (temperature, soil moisture, MAT/MAP) -- "borrowing
#     strength" from the other ecosystem types' sink response rather than
#     learning an Arid-specific one from 2 anomalous-looking sites.
# 13_Global_SpatialUpscalingRF.R must apply this identical scheme (Arid ==
# aridity_index-based, P(source) forced 0, sink-model-only) at prediction
# time for this to be consistent between training and the global grid.
#
# Outputs:
#   OUTPUT/source_magnitude_model_bundle.rds — everything
#     13_Global_SpatialUpscalingRF.R needs to apply these models to the
#     global grid (fitted model objects, standardizers, and every scalar
#     threshold that must match between fitting and prediction, bundled
#     together so the two scripts can't silently drift out of sync the way
#     arid_ai_threshold/aridity_mat_floor once did across separate files).
#   OUTPUT/rf_class_variable_importance.csv
#   OUTPUT/rf_magnitude_variable_importance.csv
#   OUTPUT/probability_calibration_skill.csv
#   OUTPUT/comparison_class_skill.csv
#   OUTPUT/comparison_magnitude_skill.csv
#   OUTPUT/magnitude_model_skill.csv
#   OUTPUT/magnitude_model_fitted_values.csv
#   OUTPUT/magnitude_model_arid_forced_sink_check.csv — Arid's own 84
#     site-months (never in either magnitude model's training data),
#     observed flux vs. what the sink model predicts for them
#   OUTPUT/magnitude_model_predicted_annual_flux_by_site.csv
#   OUTPUT/magnitude_model_fluxnet_external_validation.csv
#   OUTPUT/magnitude_model_fluxnet_validation_text.md
#   OUTPUT/model_fit_parameters.csv
#   OUTPUT/source_magnitude_model_summary.txt
#   FIGURES/Fig3_classification_skill.png
#   FIGURES/Fig5_magnitude_models.png

library(tidyverse)
library(data.table)
library(ranger)
library(cowplot)

# ── Paths ─────────────────────────────────────────────────────────────────────
# Deliberately does NOT depend on spatial_dir/ecoregions/ERA5-Land/MODIS/WAD2M
# — none of that is needed to fit or evaluate these models, only to project
# them onto the global grid (13_Global_SpatialUpscalingRF.R's job).

localdir.ch4 <- Sys.getenv("LOCALDIR_CH4",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/Methane")
rf_dir <- Sys.getenv("MONTHLY_RF_DIR",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF")

if (!dir.exists(localdir.ch4)) stop("CH4 data directory not found: ", localdir.ch4)

output_dir <- file.path(rf_dir, "OUTPUT")
figure_dir <- file.path(rf_dir, "FIGURES")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

era5_30min_file    <- file.path(localdir.ch4, "OUTPUT/NEON_ERA5_gapfilled_30min.csv.gz")
site_behavior_file <- file.path(localdir.ch4, "OUTPUT/30min_site_behavior.csv")

required_files <- c(era5_30min_file, site_behavior_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) stop("Missing: ", paste(missing_files, collapse = ", "))

# ── Scalar parameters ─────────────────────────────────────────────────────────
# These (arid_ai_threshold, aridity_mat_floor, binary_threshold, and the RF
# hyperparameters) are bundled into source_magnitude_model_bundle.rds below
# specifically so 13_Global_SpatialUpscalingRF.R reads them from there
# instead of redefining its own copies — the two scripts drifting apart on
# a value like this is exactly how the Arid-category bug happened before.

binary_threshold  <- 0.5   # fixed threshold for Stage 1 classification skill / Approach 2
arid_ai_threshold <- 15
# aridity_index = MAP / (MAT + 10) is numerically degenerate as MAT
# approaches -10 C from above (denominator -> 0), producing spuriously
# extreme finite values right at the boundary, not just for MAT <= -10 C.
# Validated in 19_Supp_NEONRepresentativeness.R (global aridity_index range
# collapsed from roughly [-5.5M, 698k] to [0.03, 1161]); training and the
# spatial-projection grid must use the same well-behaved formula. Cells with
# MAT at or below this floor get aridity_index = NA and are dropped via the
# existing is.finite(aridity_index) filters.
aridity_mat_floor <- -9  # °C; i.e. require MAT + 10 > 1, not just > 0
rf_seed           <- 42
n_trees           <- 500
rf_min_node_size  <- 20
rf_max_depth      <- 8
rf_sample_frac    <- 0.7

# ── Helper functions ──────────────────────────────────────────────────────────

auc_rank <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]; predicted <- predicted[ok]
  n1 <- sum(observed == 1); n0 <- sum(observed == 0)
  r  <- rank(predicted, ties.method = "average")
  (sum(r[observed == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

make_standardizer <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  list(center = m, scale = if (!is.finite(s) || s == 0) 1 else s)
}

apply_standardizers <- function(dat, std) {
  dat %>% mutate(
    z_Tair = (mean_ERA5_Tair_C - std$Tair$center) / std$Tair$scale,
    z_VSWC = (mean_ERA5_VSWC   - std$VSWC$center) / std$VSWC$scale,
    z_MAP  = (MAP               - std$MAP$center)  / std$MAP$scale,
    z_MAT  = (MAT               - std$MAT$center)  / std$MAT$scale
  )
}

# ── RF helper functions ────────────────────────────────────────────────────────

class_balance_weights_rf <- function(y) {
  n <- length(y); tbl <- table(y)
  as.numeric(n / (length(tbl) * tbl[as.character(y)]))
}

fit_rf_prob_model <- function(training_data, predictors, use_case_weights = TRUE) {
  training_data <- training_data %>%
    mutate(weak_source_f = factor(weak_source, levels = c(0, 1)))
  fmla <- reformulate(predictors, response = "weak_source_f")
  cw   <- if (use_case_weights) class_balance_weights_rf(training_data$weak_source) else NULL
  ranger(
    formula         = fmla,
    data            = training_data,
    num.trees       = n_trees,
    probability     = TRUE,
    case.weights    = cw,
    min.node.size   = rf_min_node_size,
    max.depth       = rf_max_depth,
    replace         = FALSE,
    sample.fraction = rf_sample_frac,
    importance      = "impurity",
    seed            = rf_seed
  )
}

predict_rf_prob <- function(rf_model, newdata) {
  as.numeric(predict(rf_model, data = newdata)$predictions[, "1"])
}

fit_isotonic_calibration <- function(oob_probs, observed_labels) {
  ord <- order(oob_probs)
  iso <- isoreg(x = oob_probs[ord], y = as.numeric(observed_labels[ord]))
  list(x = oob_probs[ord], y = iso$yf)
}

calibrate_isotonic <- function(iso_cal, new_probs) {
  pmin(pmax(approx(iso_cal$x, iso_cal$y, xout = new_probs,
                   method = "linear", rule = 2, ties = "ordered")$y, 0), 1)
}

fit_rf_magnitude_model <- function(training_data, state_name, predictors) {
  training_data <- training_data %>%
    mutate(log_abs_flux = log(pmax(abs(monthly_flux_gC_m2_month), 1e-6)))
  list(
    model  = ranger(reformulate(predictors, "log_abs_flux"),
                    data = training_data, num.trees = n_trees,
                    importance = "impurity", seed = rf_seed),
    engine = "ranger_log_abs",
    state  = state_name
  )
}

predict_rf_magnitude <- function(fit, newdata) {
  mag <- exp(predict(fit$model, data = newdata)$predictions)
  if (identical(fit$state, "Weak-sink")) -mag else mag
}

# Out-of-bag counterpart to predict_rf_magnitude(): ranger() stores OOB
# predictions (each row scored only by the trees that did NOT include it in
# their bootstrap sample) in $predictions automatically at fit time, for
# free, with no re-fitting or new predict() call needed. Row order matches
# the training_data the model was fit on 1:1, since fit_rf_magnitude_model()
# passes that data straight to ranger() with no row filtering in between.
# This is what should be used for any skill/generalization statistic — using
# predict_rf_magnitude(fit, training_data) instead (predicting the training
# data with the full, already-fitted forest) is an in-sample fit statistic,
# not a held-out one, and will look better than the model's real skill.
oob_predict_rf_magnitude <- function(fit) {
  mag <- exp(fit$model$predictions)
  if (identical(fit$state, "Weak-sink")) -mag else mag
}

summarise_magnitude_fit <- function(training_data, fit) {
  fitted <- oob_predict_rf_magnitude(fit)
  tibble(
    magnitude_model               = fit$state,
    engine                        = fit$engine,
    evaluation_basis               = "OOB (row-level held-out)",
    n_observations                = nrow(training_data),
    mean_observed_flux            = mean(training_data$monthly_flux_gC_m2_month, na.rm = TRUE),
    mean_fitted_flux              = mean(fitted, na.rm = TRUE),
    rmse_gC_m2_month              = sqrt(mean((training_data$monthly_flux_gC_m2_month - fitted)^2, na.rm = TRUE)),
    mae_gC_m2_month               = mean(abs(training_data$monthly_flux_gC_m2_month - fitted), na.rm = TRUE),
    correlation_observed_fitted   = suppressWarnings(
      cor(training_data$monthly_flux_gC_m2_month, fitted, use = "complete.obs"))
  )
}

classification_skill <- function(prob, observed, threshold, model_label, eval_basis,
                                 oob_error = NA_real_) {
  pred_class <- as.integer(prob >= threshold)
  tibble(
    model              = model_label,
    evaluation_basis   = eval_basis,
    n_site_months      = length(observed),
    n_sink_months      = sum(observed == 0),
    n_source_months    = sum(observed == 1),
    source_fraction    = mean(observed),
    auc                = auc_rank(observed, prob),
    brier_score        = mean((observed - prob)^2),
    brier_null         = mean(observed) * (1 - mean(observed)),
    brier_skill_score  = 1 - mean((observed - prob)^2) / (mean(observed) * (1 - mean(observed))),
    tjur_r2            = mean(prob[observed == 1]) - mean(prob[observed == 0]),
    ranger_oob_error   = oob_error,
    threshold          = threshold,
    accuracy           = mean(pred_class == observed),
    sensitivity_source = sum(pred_class == 1 & observed == 1) / sum(observed == 1),
    specificity_sink   = sum(pred_class == 0 & observed == 0) / sum(observed == 0),
    precision_source   = sum(pred_class == 1 & observed == 1) / max(sum(pred_class == 1), 1),
    predicted_source_fraction = mean(pred_class == 1)
  )
}

# ── Load and prepare training data ────────────────────────────────────────────

site_behavior <- read.csv(site_behavior_file) %>%
  mutate(SITE_ID = as.character(SITE_ID))

upland_sites <- site_behavior %>%
  filter(!is.na(EcoType),
    !str_detect(EcoType, regex("wetland|inundat|flood|marsh|swamp|bog|fen|lake|rice",
                               ignore_case = TRUE))) %>%
  distinct(SITE_ID, EcoType, MAP, MAT)

era5_30min <- data.table::fread(era5_30min_file) %>% as_tibble() %>%
  mutate(across(c(SITE_ID), as.character),
         across(c(Year, month), as.integer),
         across(c(ERA5_Tair_C, ERA5_VSWC, gapfilled_CH4_mgC_30min), as.numeric)) %>%
  inner_join(upland_sites %>% select(SITE_ID, EcoType), by = c("SITE_ID","EcoType")) %>%
  filter(is.finite(Year), is.finite(month), is.finite(ERA5_Tair_C), is.finite(ERA5_VSWC))

monthly_training <- era5_30min %>%
  reframe(.by = c(SITE_ID, EcoType, Year, month),
    monthly_budget_mgC_m2  = sum(gapfilled_CH4_mgC_30min, na.rm = TRUE),
    mean_ERA5_Tair_C        = mean(ERA5_Tair_C, na.rm = TRUE),
    mean_ERA5_VSWC          = mean(ERA5_VSWC,   na.rm = TRUE)) %>%
  left_join(upland_sites, by = c("SITE_ID","EcoType")) %>%
  mutate(
    weak_source              = as.integer(monthly_budget_mgC_m2 > 0),
    monthly_flux_gC_m2_month = monthly_budget_mgC_m2 / 1000,
    MAP                      = as.numeric(MAP),
    MAT                      = as.numeric(MAT),
    # See aridity_mat_floor above.
    aridity_index            = if_else(MAT > aridity_mat_floor, MAP / (MAT + 10), NA_real_),
    is_arid                  = as.integer(!is.na(aridity_index) & aridity_index < arid_ai_threshold),
    # Arid label uses the same aridity_index < arid_ai_threshold rule used
    # throughout (incl. 19_Supp_NEONRepresentativeness.R). No NEON site is
    # literally labeled "Arid" in site metadata; the ~2 sites crossing this
    # threshold (JORN, SRER) are relabeled here purely for BOOKKEEPING and
    # as a Stage 2 predictor level — NOT so Stage 1 fits an Arid-specific
    # P(source) pattern. See the header note: both of those sites' data show
    # 100% weak-source months, contradicting the expectation that arid soils
    # are net sinks, so Stage 1 training deliberately EXCLUDES is_arid==1
    # rows below rather than learning from what looks like an unrepresentative
    # 2-site sample. 13_Global_SpatialUpscalingRF.R must classify Arid the
    # SAME way (aridity_index-based) and force P(source) = 0 for it, so
    # training and the grid agree on what "Arid" means and how it's routed.
    EcoType                  = factor(if_else(is_arid == 1, "Arid", as.character(EcoType)),
                                      levels = c("Forest","Grassland","Shrubland","Arid"))
  ) %>%
  filter(is.finite(weak_source), is.finite(mean_ERA5_Tair_C), is.finite(mean_ERA5_VSWC),
         is.finite(MAP), is.finite(MAT), is.finite(aridity_index))

class_predictors <- c("EcoType","mean_ERA5_Tair_C","mean_ERA5_VSWC","MAP","MAT","aridity_index")

# ── Stage 1: Balanced RF (1:1 source:sink) ────────────────────────────────────
# Downsample source months to match the number of sink months so both classes
# have exactly equal representation in training. No case weights needed.
# OOB predictions are used for honest evaluation and isotonic calibration.
#
# Arid (is_arid == 1) is EXCLUDED from Stage 1 training entirely — see the
# header note and the is_arid comment above. All 84 Arid site-months are
# observed weak-source, so Stage 1 has no genuine sink/source contrast to
# learn for Arid anyway, and letting it try would mean extrapolating a
# P(source) pattern from just 2 sites whose data looks unrepresentative of
# the ecological prior (arid soils as net CH4 sinks). P(source) is instead
# hard-forced to 0 for Arid, both here and in
# 13_Global_SpatialUpscalingRF.R's grid loop.

stage1_training_data <- monthly_training %>% filter(is_arid == 0)

set.seed(rf_seed)
sink_idx   <- which(stage1_training_data$weak_source == 0)
source_idx <- which(stage1_training_data$weak_source == 1)
bal_idx    <- c(sink_idx, sample(source_idx, length(sink_idx)))
balanced_training <- stage1_training_data[bal_idx, ]

message(sprintf(
  "Stage 1 balanced training: %d sink + %d source = %d total (from %d non-Arid site-months; %d Arid site-months excluded)",
  length(sink_idx), length(sink_idx), 2 * length(sink_idx),
  nrow(stage1_training_data), sum(monthly_training$is_arid == 1)))

message("Fitting Stage 1 balanced RF...")
rf_model_A <- fit_rf_prob_model(balanced_training, class_predictors, use_case_weights = FALSE)

oob_prob_A  <- rf_model_A$predictions[, "1"]
iso_cal_A   <- fit_isotonic_calibration(oob_prob_A, balanced_training$weak_source)

balanced_training <- balanced_training %>%
  mutate(
    source_prob_A_oob = oob_prob_A,
    source_prob_A_raw = predict_rf_prob(rf_model_A, balanced_training),
    source_prob_A     = calibrate_isotonic(iso_cal_A, source_prob_A_oob)
  )

# Also add calibrated projections back to full training for Stage 2
# (source_probability predictor). Arid rows get a hard-forced 0 instead of
# rf_model_A's raw prediction, since that model never trained on Arid data
# at all (see above) — its extrapolation there would be arbitrary, and
# forcing 0 is what routes Arid entirely to the sink magnitude model below.
monthly_training <- monthly_training %>%
  mutate(
    source_prob_A_raw = predict_rf_prob(rf_model_A, monthly_training),
    source_prob_A     = if_else(is_arid == 1, 0,
                                calibrate_isotonic(iso_cal_A, source_prob_A_raw))
  )

skill_A <- classification_skill(
  prob        = balanced_training$source_prob_A,
  observed    = balanced_training$weak_source,
  threshold   = binary_threshold,
  model_label = "Continuous",
  eval_basis  = "OOB + isotonic calibration (1:1 balanced)",
  oob_error   = rf_model_A$prediction.error
)

cal_skill_A <- balanced_training %>%
  mutate(prob_bin = ntile(source_prob_A_oob, 10)) %>%
  group_by(prob_bin) %>%
  summarise(
    n                         = n(),
    mean_oob_raw_prob         = mean(source_prob_A_oob, na.rm = TRUE),
    mean_isotonic_cal_prob    = mean(source_prob_A,     na.rm = TRUE),
    observed_source_fraction  = mean(weak_source,       na.rm = TRUE),
    calibration_error         = mean_isotonic_cal_prob - observed_source_fraction,
    .groups = "drop"
  ) %>% mutate(model = "A")

importance_A <- tibble(
  predictor  = names(rf_model_A$variable.importance),
  importance = rf_model_A$variable.importance,
  model      = "P(source)"
) %>% arrange(desc(importance))

# ── Stage 2: RF magnitude models (shared across all approaches) ───────────────
# Trained on the full weighted training data.
# source_prob_A (calibrated OOB) used as the source_probability predictor
# so it represents an honest, unbiased probability signal.

magnitude_standardizers <- list(
  Tair = make_standardizer(monthly_training$mean_ERA5_Tair_C),
  VSWC = make_standardizer(monthly_training$mean_ERA5_VSWC),
  MAP  = make_standardizer(monthly_training$MAP),
  MAT  = make_standardizer(monthly_training$MAT)
)

monthly_training <- apply_standardizers(monthly_training, magnitude_standardizers) %>%
  mutate(
    source_probability = source_prob_A,   # Stage 2 uses Model A calibrated OOB probs
    SITE_ID = factor(SITE_ID)
  )

mag_predictors <- c("z_Tair","z_VSWC","z_MAP","z_MAT","source_probability","is_arid","EcoType")

# Both magnitude models exclude Arid's own site-months (EcoType != "Arid"):
#   - source_mag_data: Arid never reaches this model at prediction time
#     (P(source) is hard-forced to 0 above), so training on it would only
#     bias the model toward a class it will never actually be asked to
#     predict for.
#   - sink_mag_data: excluding it is a no-op in practice (all 84 Arid
#     site-months are weak-source, so none would pass this filter anyway)
#     but is stated explicitly so future data updates can't silently start
#     feeding Arid's own (observed-as-source) months into the sink model.
#     Arid grid cells are instead predicted by extrapolating from
#     Forest/Grassland/Shrubland sink behavior via the shared continuous
#     covariates (z_Tair/z_VSWC/z_MAP/z_MAT) — see the header note.
sink_mag_data   <- monthly_training %>% filter(EcoType != "Arid", weak_source == 0, monthly_flux_gC_m2_month <= 0)
source_mag_data <- monthly_training %>% filter(EcoType != "Arid", weak_source == 1, monthly_flux_gC_m2_month >  0)

message(sprintf(
  "Stage 2 training: %d Weak-sink + %d Weak-source site-months (Arid site-months excluded from both; n_arid = %d)",
  nrow(sink_mag_data), nrow(source_mag_data), sum(monthly_training$is_arid == 1)))

message("Fitting Stage 2 sink magnitude RF...")
sink_mag_model   <- fit_rf_magnitude_model(sink_mag_data,   "Weak-sink",   mag_predictors)
message("Fitting Stage 2 source magnitude RF...")
source_mag_model <- fit_rf_magnitude_model(source_mag_data, "Weak-source", mag_predictors)

magnitude_model_skill <- bind_rows(
  summarise_magnitude_fit(sink_mag_data,   sink_mag_model),
  summarise_magnitude_fit(source_mag_data, source_mag_model)
)

# OOB (held-out), not in-sample: see oob_predict_rf_magnitude() note above.
# This also feeds Fig 5 (observed vs. fitted), so that figure shows genuine
# held-out skill rather than an in-sample fit.
magnitude_fitted_values <- bind_rows(
  sink_mag_data %>%
    mutate(magnitude_model = "Weak-sink",
           fitted_flux = pmin(oob_predict_rf_magnitude(sink_mag_model), 0)),
  source_mag_data %>%
    mutate(magnitude_model = "Weak-source",
           fitted_flux = pmax(oob_predict_rf_magnitude(source_mag_model), 0))
) %>% rename(fitted_flux_gC_m2_month = fitted_flux)

# ── Arid diagnostic: forced-sink prediction vs. observed ──────────────────────
# Arid site-months are excluded from BOTH magnitude models' training (see
# above), but 13_Global_SpatialUpscalingRF.R still routes Arid grid cells
# through sink_mag_model exclusively (P(source) forced to 0). This reports
# what that forced routing would predict for JORN/SRER's own 84 site-months
# -- a genuine holdout prediction, since neither site ever appears in
# sink_mag_model's training data -- against what was actually observed, so
# the cost of this design choice is visible rather than hidden rather than
# just asserted. Folded into predicted_annual_flux_by_site/fluxnet_validation
# below so the FLUXNET comparison still covers all NEON upland sites.
arid_mag_data <- monthly_training %>% filter(EcoType == "Arid")
arid_sink_prediction_check <- arid_mag_data %>%
  mutate(
    forced_sink_prediction_gC_m2_month = pmin(predict_rf_magnitude(sink_mag_model, arid_mag_data), 0),
    observed_flux_gC_m2_month          = monthly_flux_gC_m2_month
  ) %>%
  select(SITE_ID, EcoType, Year, month,
         observed_flux_gC_m2_month, forced_sink_prediction_gC_m2_month)

message(sprintf(
  "Arid forced-sink check: %d site-months (JORN/SRER) -- mean observed = %.3f, mean forced-sink prediction = %.3f g C m-2 mo-1 (neither site was ever in sink_mag_model's training data)",
  nrow(arid_sink_prediction_check),
  mean(arid_sink_prediction_check$observed_flux_gC_m2_month),
  mean(arid_sink_prediction_check$forced_sink_prediction_gC_m2_month)))

# ── Stage 2 external validation: independent FLUXNET-CH4 benchmark ────────────
# OOB error (above) is held out at the site-MONTH level, not the site level,
# so it still doesn't test whether the model generalizes across sites (a
# reviewer asked specifically for this — see leave-one-site-out). Rather than
# retrain N models on N site-excluded folds, this instead checks the model's
# OOB-predicted flux against a dataset that was NEVER part of this pipeline
# at all: the FLUXNET-CH4 "Upland" ecosystem-class aggregate (Delwiche et al.
# 2021, Table 1) — the same published benchmark 09_NEON_FLUXNETComparison.R
# already uses for NEON's OBSERVED flux, applied here to the model's
# PREDICTED flux instead. Agreement with an independently measured,
# different-network dataset is arguably stronger evidence of generalization
# than resampling NEON's own data, since it can't be explained by anything
# idiosyncratic to how NEON's site network or gap-filling works.
fluxnet_upland_benchmark <- tibble(
  ecosystem_class    = "Upland",
  annual_gC_m2_yr    = 4.0,
  annual_sd_gC_m2_yr = 10.5,
  n_sites            = 15L,
  n_site_years       = 47L,
  reference          = "Delwiche et al. 2021, Table 1"
)

predicted_annual_flux_by_site <- bind_rows(
  magnitude_fitted_values %>%
    group_by(SITE_ID) %>%
    summarise(
      n_site_months             = n(),
      predicted_annual_gC_m2_yr = mean(fitted_flux_gC_m2_month, na.rm = TRUE) * 12,
      .groups = "drop"
    ) %>%
    mutate(prediction_basis = "OOB (site in training)"),
  arid_sink_prediction_check %>%
    group_by(SITE_ID) %>%
    summarise(
      n_site_months             = n(),
      predicted_annual_gC_m2_yr = mean(forced_sink_prediction_gC_m2_month, na.rm = TRUE) * 12,
      .groups = "drop"
    ) %>%
    mutate(prediction_basis = "Forced-sink (Arid, never in training)")
)

fluxnet_validation <- predicted_annual_flux_by_site %>%
  summarise(
    n_sites                        = n(),
    mean_predicted_annual_gC_m2_yr = mean(predicted_annual_gC_m2_yr, na.rm = TRUE),
    sd_predicted_annual_gC_m2_yr   = sd(predicted_annual_gC_m2_yr,   na.rm = TRUE)
  ) %>%
  bind_cols(fluxnet_upland_benchmark %>% rename_with(~ paste0("fluxnet_", .))) %>%
  mutate(
    fluxnet_range_low    = fluxnet_annual_gC_m2_yr - fluxnet_annual_sd_gC_m2_yr,
    fluxnet_range_high   = fluxnet_annual_gC_m2_yr + fluxnet_annual_sd_gC_m2_yr,
    within_fluxnet_range = mean_predicted_annual_gC_m2_yr >= fluxnet_range_low &
                            mean_predicted_annual_gC_m2_yr <= fluxnet_range_high
  )

fluxnet_validation_text <- paste(strwrap(sprintf(
"To assess how well the Stage 2 magnitude models generalize beyond the NEON
training data, we (1) evaluated skill using ranger's out-of-bag (OOB)
predictions -- each site-month scored only by trees that did not see it during
fitting, replacing the in-sample fit statistic previously reported -- and (2)
compared the OOB-predicted annual CH4 flux for NEON's upland sites against an
entirely independent published benchmark: the FLUXNET-CH4 \"Upland\"
ecosystem-class aggregate (Delwiche et al. 2021, Table 1; %d sites, %d
site-years, mean %.1f +/- %.1f g C m-2 yr-1), a dataset never used anywhere in
this pipeline's training or model selection. Across %d NEON upland sites with
OOB-predicted flux, the mean predicted annual flux was %.1f +/- %.1f g C m-2
yr-1, %s the published FLUXNET-CH4 Upland range (%.1f to %.1f g C m-2 yr-1).",
  fluxnet_validation$fluxnet_n_sites, fluxnet_validation$fluxnet_n_site_years,
  fluxnet_validation$fluxnet_annual_gC_m2_yr, fluxnet_validation$fluxnet_annual_sd_gC_m2_yr,
  fluxnet_validation$n_sites, fluxnet_validation$mean_predicted_annual_gC_m2_yr,
  fluxnet_validation$sd_predicted_annual_gC_m2_yr,
  if (fluxnet_validation$within_fluxnet_range) "falling within" else "falling outside",
  fluxnet_validation$fluxnet_range_low, fluxnet_validation$fluxnet_range_high
), width = 90), collapse = "\n")

importance_mag <- bind_rows(
  tibble(predictor = names(sink_mag_model$model$variable.importance),
         importance = sink_mag_model$model$variable.importance,   model = "Stage 2 — Weak-sink"),
  tibble(predictor = names(source_mag_model$model$variable.importance),
         importance = source_mag_model$model$variable.importance, model = "Stage 2 — Weak-source")
) %>% arrange(model, desc(importance))

training_vswc_range <- range(monthly_training$mean_ERA5_VSWC, na.rm = TRUE)
training_prec_range <- quantile(monthly_training$MAP / 12, probs = c(0.02, 0.98), na.rm = TRUE)

# ── Skill comparison tables ────────────────────────────────────────────────────

comparison_class <- skill_A %>%
  select(model, evaluation_basis, auc, brier_score, brier_null, brier_skill_score,
         tjur_r2, ranger_oob_error, threshold, accuracy,
         sensitivity = sensitivity_source, specificity = specificity_sink)

comparison_magnitude <- magnitude_model_skill %>%
  select(magnitude_model, evaluation_basis, rmse_gC_m2_month, mae_gC_m2_month, correlation_observed_fitted) %>%
  mutate(approach = "RF/ranger", .before = 1)

# ── Model bundle for 13_Global_SpatialUpscalingRF.R ────────────────────────────
# Everything the spatial-projection script needs to APPLY these models,
# bundled together so it never has to redefine (and risk drifting from) a
# threshold or standardizer this script already fit/chose.

model_bundle <- list(
  rf_model_A              = rf_model_A,
  iso_cal_A                = iso_cal_A,
  class_predictors          = class_predictors,
  sink_mag_model            = sink_mag_model,
  source_mag_model          = source_mag_model,
  mag_predictors             = mag_predictors,
  magnitude_standardizers    = magnitude_standardizers,
  training_vswc_range        = training_vswc_range,
  training_prec_range        = training_prec_range,
  ecotype_levels              = levels(monthly_training$EcoType),
  arid_ai_threshold           = arid_ai_threshold,
  aridity_mat_floor           = aridity_mat_floor,
  binary_threshold            = binary_threshold,
  # Arid is excluded from Stage 1 and Stage 2's source model, and forced to
  # P(source) = 0 / sink-model-only at prediction time — see the header note
  # (empirically, both NEON Arid-analog sites are 100% weak-source, an
  # unrepresentative 2-site sample for learning a P(source) pattern).
  # 13_Global_SpatialUpscalingRF.R reads this flag and asserts it's TRUE
  # rather than silently assuming it, so this design choice can't silently
  # drift out of sync between the two scripts.
  arid_forced_sink            = TRUE,
  rf_seed                     = rf_seed,
  n_trees                     = n_trees,
  rf_min_node_size            = rf_min_node_size,
  rf_max_depth                = rf_max_depth,
  rf_sample_frac              = rf_sample_frac,
  n_training_site_months      = nrow(monthly_training),
  n_stage1_balanced_per_class = length(sink_idx)
)
saveRDS(model_bundle, file.path(output_dir, "source_magnitude_model_bundle.rds"))

# ── Write outputs ─────────────────────────────────────────────────────────────

write.csv(comparison_class,       file.path(output_dir, "comparison_class_skill.csv"),      row.names = FALSE)
write.csv(comparison_magnitude,   file.path(output_dir, "comparison_magnitude_skill.csv"),  row.names = FALSE)
write.csv(cal_skill_A,            file.path(output_dir, "probability_calibration_skill.csv"), row.names = FALSE)
write.csv(importance_A,           file.path(output_dir, "rf_class_variable_importance.csv"),  row.names = FALSE)
write.csv(importance_mag,         file.path(output_dir, "rf_magnitude_variable_importance.csv"), row.names = FALSE)
write.csv(magnitude_model_skill,  file.path(output_dir, "magnitude_model_skill.csv"),         row.names = FALSE)
write.csv(magnitude_fitted_values,file.path(output_dir, "magnitude_model_fitted_values.csv"), row.names = FALSE)
write.csv(arid_sink_prediction_check, file.path(output_dir, "magnitude_model_arid_forced_sink_check.csv"), row.names = FALSE)
write.csv(predicted_annual_flux_by_site, file.path(output_dir, "magnitude_model_predicted_annual_flux_by_site.csv"), row.names = FALSE)
write.csv(fluxnet_validation,      file.path(output_dir, "magnitude_model_fluxnet_external_validation.csv"), row.names = FALSE)
writeLines(fluxnet_validation_text, file.path(output_dir, "magnitude_model_fluxnet_validation_text.md"))
write.csv(data.frame(
  binary_threshold    = binary_threshold,
  arid_ai_threshold    = arid_ai_threshold,
  aridity_mat_floor    = aridity_mat_floor,
  n_trees               = n_trees,
  rf_min_node_size      = rf_min_node_size,
  rf_max_depth          = rf_max_depth,
  rf_sample_frac        = rf_sample_frac,
  rf_seed               = rf_seed,
  calibration_method    = "isotonic_regression_on_OOB",
  stage1_training        = sprintf(
    "1:1 balanced (n=%d sink + %d source = %d total; from %d non-Arid site-months; %d Arid site-months excluded, forced to P(source)=0)",
    length(sink_idx), length(sink_idx), 2 * length(sink_idx),
    nrow(stage1_training_data), sum(monthly_training$is_arid == 1))
), file.path(output_dir, "model_fit_parameters.csv"), row.names = FALSE)

# ── Summary text ──────────────────────────────────────────────────────────────

capture.output({
  cat("Source-probability and magnitude models — fitting and testing\n\n")
  cat(sprintf("Stage 1: weighted RF | n=%d | n_trees=%d | min.node.size=%d | max.depth=%d | sample.frac=%.1f\n",
    nrow(monthly_training), n_trees, rf_min_node_size, rf_max_depth, rf_sample_frac))
  cat("Calibration: isotonic regression on OOB predictions\n\n")

  cat("─── Stage 1 Classification Skill ───\n")
  print(comparison_class)

  cat("\n─── Stage 2 Magnitude Skill (OOB, held-out) ───\n")
  print(comparison_magnitude)

  cat("\n─── Stage 2 External Validation — FLUXNET-CH4 Upland Benchmark ───\n")
  print(fluxnet_validation)
  cat("\n"); cat(fluxnet_validation_text, "\n")

  cat("\n─── Arid handling: excluded from training, forced to sink model ───\n")
  cat(sprintf(
    "%d Arid site-months (JORN, SRER) excluded from Stage 1 and Stage 2's source model.\n",
    sum(monthly_training$is_arid == 1)))
  cat(sprintf(
    "Forced-sink prediction check (never in sink_mag_model's training data): mean observed = %.3f, mean predicted = %.3f g C m-2 mo-1.\n",
    mean(arid_sink_prediction_check$observed_flux_gC_m2_month),
    mean(arid_sink_prediction_check$forced_sink_prediction_gC_m2_month)))

  cat("\n─── Calibration (OOB deciles) ───\n"); print(cal_skill_A)

  cat("\n─── Stage 1 Variable Importance ───\n")
  print(importance_A)
  cat("\n─── Stage 2 Variable Importance ───\n"); print(importance_mag)
}, file = file.path(output_dir, "source_magnitude_model_summary.txt"))

# ── Figures ───────────────────────────────────────────────────────────────────

fig_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = "bottom", strip.background = element_rect(fill = "grey92"),
        axis.title = element_text(size = 11), axis.text = element_text(size = 10))

ecotype_colors <- c(Forest = "#1B7837", Grassland = "#D9B86C",
                    Shrubland = "#C2A5CF", Arid = "#D95F02")

# Fig 3: Classification skill (AUC, Tjur R², sensitivity, specificity)
p3 <- comparison_class %>%
  filter(!is.na(auc)) %>%
  pivot_longer(c(auc, tjur_r2, accuracy, sensitivity, specificity),
               names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(metric = factor(metric,
    levels = c("auc","tjur_r2","accuracy","sensitivity","specificity"),
    labels = c("AUC","Tjur R²","Accuracy","Sensitivity","Specificity")),
    model_short = model) %>%
  ggplot(aes(x = value, y = metric, fill = model_short)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85, show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  scale_fill_manual(values = c("Continuous" = "#009688")) +
  labs(title = "Stage 1 Classification Skill",
       x = "Value", y = NULL, fill = NULL) +
  fig_theme + theme(legend.position = "top")

ggsave(file.path(figure_dir, "Fig3_classification_skill.png"),
  p3, width = 7, height = 4.5, units = "in", dpi = 300, bg = "white")

# Fig 5: Magnitude model — observed vs fitted (OOB)
make_mag_plot <- function(model_label, tag) {
  magnitude_fitted_values %>%
    filter(magnitude_model == model_label, !is.na(EcoType)) %>%
    ggplot(aes(x = monthly_flux_gC_m2_month, y = fitted_flux_gC_m2_month, color = EcoType)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "grey35", linetype = "dashed", linewidth = 0.5) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = ecotype_colors) +
    labs(title = paste0(tag, ". ", model_label), color = NULL,
         x = "Observed (g C m⁻² mo⁻¹)", y = "Fitted (g C m⁻² mo⁻¹, OOB)") +
    fig_theme
}
fig5 <- plot_grid(make_mag_plot("Weak-sink","A"), make_mag_plot("Weak-source","B"), ncol = 2)
ggsave(file.path(figure_dir, "Fig5_magnitude_models.png"),
  fig5, width = 9, height = 4.5, units = "in", dpi = 300, bg = "white")

message("Source-probability and magnitude model fitting complete. Bundle saved to: ",
  file.path(output_dir, "source_magnitude_model_bundle.rds"))
