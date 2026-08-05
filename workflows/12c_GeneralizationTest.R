# Leave-one-site-out (LOSO) cross-validation — Stage 2 flux-magnitude RF
# models.
#
# 12_SourceProp_MagnitudeModels.R already reports two generalization
# statistics for Stage 2: ranger's OOB skill (each site-MONTH scored only by
# trees that didn't see it) and an external comparison against the
# published FLUXNET-CH4 "Upland" benchmark. Neither directly answers a
# reviewer's specific question: how well does the model do on a NEON SITE it
# has never seen at all? OOB held-out rows still come from sites the model
# trained on; the FLUXNET comparison is an aggregate, not a site-by-site
# test. This script closes that gap with an actual leave-one-site-out test,
# without touching the production models 13_Global_SpatialUpscalingRF.R
# uses.
#
# Design — nested, to avoid leakage:
#   For each NEON upland site in turn, EVERYTHING is refit on the other
#   sites only:
#     1. Stage 1 (balanced P(source) classifier) is refit on the training
#        sites, with its own isotonic calibration fit on ITS OWN OOB
#        predictions (never touching the held-out site).
#     2. That fold's calibrated Stage 1 model generates the
#        source_probability feature for both the training sites and the
#        held-out site — so Stage 2's "generalization" score is never
#        contaminated by a Stage 1 model that already learned something
#        about the held-out site's typical behavior.
#     3. Stage 2 standardizers (z-scoring) are fit on the training sites
#        only, then applied to the held-out site.
#     4. Stage 2 sink/source magnitude RFs are refit on the training sites
#        and used to predict the held-out site's months.
#   This means N+1 Stage-1 fits and up to 2N Stage-2 fits (N = number of
#   NEON upland sites) — slower than 12_SourceProp_MagnitudeModels.R's
#   single fit, but each fold is a small RF (500 trees, a few hundred rows)
#   so the whole loop typically finishes in a few minutes.
#
# This script is purely diagnostic. The fold-specific models are discarded
# after scoring; source_magnitude_model_bundle.rds (written by
# 12_SourceProp_MagnitudeModels.R, used by 13_Global_SpatialUpscalingRF.R)
# is untouched and still reflects models fit on ALL sites — refitting it on
# a held-out-site basis would defeat the purpose of using all available
# NEON data operationally.
#
# Arid handling mirrors 12_SourceProp_MagnitudeModels.R's production design
# exactly, per-fold: only 2 NEON sites (JORN, SRER) cross the aridity_index
# < arid_ai_threshold "Arid" definition, and both are empirically 100%
# weak-source — an unrepresentative sample for a learned P(source) pattern.
# So in every fold, Stage 1 and the source-magnitude model are refit
# EXCLUDING Arid site-months entirely, and source_probability is forced to 0
# for Arid rows rather than taken from the fold's Stage 1 model. When the
# HELD-OUT site itself is Arid (JORN or SRER), its test months are scored
# against the fold's sink-magnitude model directly (labeled "Weak-sink
# (forced, Arid)"), regardless of their own observed sign — matching how
# 13_Global_SpatialUpscalingRF.R actually treats Arid cells, rather than
# testing a probabilistic routing the production pipeline never uses. This
# is still a genuine holdout: that fold's sink model never saw ANY Arid
# site's months (not even the one still in "training", since Arid is
# excluded from all Stage 2 training regardless of fold).
#
# Outputs:
#   OUTPUT/magnitude_model_loso_predictions.csv       — every held-out
#     site-month, observed vs LOSO-predicted flux
#   OUTPUT/magnitude_model_loso_skill_overall.csv      — RMSE/MAE/bias/r by
#     magnitude_model (Weak-sink / Weak-source)
#   OUTPUT/magnitude_model_loso_skill_by_ecotype.csv
#   OUTPUT/magnitude_model_loso_skill_by_site.csv       — per-site skill,
#     useful for spotting sites the models systematically fail on
#   OUTPUT/magnitude_model_generalization_comparison.csv — LOSO vs OOB
#     side by side (OOB pulled from 12_SourceProp_MagnitudeModels.R's
#     output if present)
#   OUTPUT/magnitude_model_generalization_summary.txt  — narrative,
#     numbers regenerated from the tables above on every run
#   FIGURES/Fig_LOSO_magnitude_models.png               — observed vs
#     LOSO-predicted, Weak-sink / Weak-source

library(tidyverse)
library(data.table)
library(ranger)
library(cowplot)

# ── Paths ─────────────────────────────────────────────────────────────────────
# Same raw inputs as 12_SourceProp_MagnitudeModels.R — this script rebuilds
# monthly_training itself rather than loading the bundle, since a
# leave-one-site-out test needs to refit everything per fold anyway.

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
# Static site attributes (SITE_ID, EcoType, MAP, MAT) from 05_NEON_FluxAnalysis.R's
# actively-regenerated summary, not the orphaned OUTPUT/30min_site_behavior.csv.
site_attributes_file <- file.path(localdir.ch4, "OUTPUT/NEON_scale_annual_budget_summary.csv")

required_files <- c(era5_30min_file, site_attributes_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) stop("Missing: ", paste(missing_files, collapse = ", "))

# ── Scalar parameters ─────────────────────────────────────────────────────────
# Deliberately duplicated from 12_SourceProp_MagnitudeModels.R rather than
# read from source_magnitude_model_bundle.rds — that bundle holds models fit
# on ALL sites, and loading its scalars is harmless, but loading its MODELS
# would defeat this script's entire purpose. Keep these in sync manually if
# 12's values ever change (same pattern used between 12/13/19).

binary_threshold  <- 0.5
arid_ai_threshold <- 15
aridity_mat_floor <- -9  # see 12_SourceProp_MagnitudeModels.R for rationale
rf_seed           <- 42
n_trees           <- 500
rf_min_node_size  <- 20
rf_max_depth      <- 8
rf_sample_frac    <- 0.7

# Minimum training rows required to fit a fold-specific Stage 2 magnitude
# model at all. Below this, a 500-tree RF fit is numerically possible but
# the resulting skill numbers are not meaningful (single-digit training
# rows) — such folds are skipped and reported rather than silently
# included.
min_fold_training_rows <- 10

# ── Helper functions ──────────────────────────────────────────────────────────
# Duplicated from 12_SourceProp_MagnitudeModels.R (kept in sync manually,
# same pattern used between 12/13/19).

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
    importance      = "none",
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
                    importance = "none", seed = rf_seed),
    engine = "ranger_log_abs",
    state  = state_name
  )
}

predict_rf_magnitude <- function(fit, newdata) {
  mag <- exp(predict(fit$model, data = newdata)$predictions)
  if (identical(fit$state, "Weak-sink")) -mag else mag
}

# ── Load and prepare training data ────────────────────────────────────────────
# Identical construction to 12_SourceProp_MagnitudeModels.R.

site_attributes <- read.csv(site_attributes_file) %>%
  mutate(SITE_ID = as.character(SITE_ID))

upland_sites <- site_attributes %>%
  filter(!is.na(EcoType),
    !str_detect(EcoType, regex("wetland|inundat|flood|marsh|swamp|bog|fen|lake|rice|crop|agri",
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
    aridity_index            = if_else(MAT > aridity_mat_floor, MAP / (MAT + 10), NA_real_),
    is_arid                  = as.integer(!is.na(aridity_index) & aridity_index < arid_ai_threshold),
    # Arid label is derived from site-level MAT/MAP climate normals, not
    # from the flux target itself, so computing it once up front (rather
    # than per LOSO fold) introduces no leakage.
    EcoType                  = factor(if_else(is_arid == 1, "Arid", as.character(EcoType)),
                                      levels = c("Forest","Grassland","Shrubland","Arid"))
  ) %>%
  filter(is.finite(weak_source), is.finite(mean_ERA5_Tair_C), is.finite(mean_ERA5_VSWC),
         is.finite(MAP), is.finite(MAT), is.finite(aridity_index))

class_predictors <- c("EcoType","mean_ERA5_Tair_C","mean_ERA5_VSWC","MAP","MAT","aridity_index")
mag_predictors   <- c("z_Tair","z_VSWC","z_MAP","z_MAT","source_probability","is_arid","EcoType")

sites         <- sort(unique(monthly_training$SITE_ID))
n_arid_sites  <- monthly_training %>% filter(is_arid == 1) %>% distinct(SITE_ID) %>% nrow()

message(sprintf(
  "Leave-one-site-out CV: %d sites, %d site-months total (%d Arid sites: %s)",
  length(sites), nrow(monthly_training), n_arid_sites,
  paste(monthly_training %>% filter(is_arid == 1) %>% distinct(SITE_ID) %>% pull(SITE_ID),
        collapse = ", ")))

# ── One LOSO fold: refit Stage 1 + Stage 2 excluding one site, predict it ─────

run_loso_fold <- function(train_data, test_data, fold_site) {
  # is_arid is a static per-SITE flag (derived from site-level MAT/MAP
  # normals, not from month-to-month flux), so every row of a held-out
  # site's test_data shares the same value.
  is_arid_fold <- isTRUE(test_data$is_arid[1] == 1)

  # Stage 1 refit: Arid site-months are excluded from training entirely,
  # mirroring 12_SourceProp_MagnitudeModels.R's production design (see the
  # header note — both NEON Arid-analog sites are 100% weak-source, an
  # unrepresentative sample for a learned P(source) pattern). This holds
  # regardless of which site is held out.
  stage1_train_f <- train_data %>% filter(is_arid == 0)

  set.seed(rf_seed)
  sink_idx_f   <- which(stage1_train_f$weak_source == 0)
  source_idx_f <- which(stage1_train_f$weak_source == 1)
  if (length(sink_idx_f) < 2 || length(source_idx_f) < 2) {
    message(sprintf("  %s: too few non-Arid sink/source months in training fold (%d sink, %d source) to balance Stage 1 — fold skipped.",
      fold_site, length(sink_idx_f), length(source_idx_f)))
    return(NULL)
  }
  bal_idx_f        <- c(sink_idx_f, sample(source_idx_f, length(sink_idx_f)))
  balanced_train_f <- stage1_train_f[bal_idx_f, ]

  rf_model_A_f <- fit_rf_prob_model(balanced_train_f, class_predictors, use_case_weights = FALSE)
  iso_cal_A_f  <- fit_isotonic_calibration(rf_model_A_f$predictions[, "1"], balanced_train_f$weak_source)

  # source_probability: forced to 0 for Arid rows (never modeled by Stage
  # 1), calibrated Stage-1 prediction otherwise — mirrors production.
  train_data <- train_data %>%
    mutate(source_probability = if_else(is_arid == 1, 0,
      calibrate_isotonic(iso_cal_A_f, predict_rf_prob(rf_model_A_f, train_data))))
  test_data <- test_data %>%
    mutate(source_probability = if_else(is_arid == 1, 0,
      calibrate_isotonic(iso_cal_A_f, predict_rf_prob(rf_model_A_f, test_data))))

  std_f <- list(
    Tair = make_standardizer(train_data$mean_ERA5_Tair_C),
    VSWC = make_standardizer(train_data$mean_ERA5_VSWC),
    MAP  = make_standardizer(train_data$MAP),
    MAT  = make_standardizer(train_data$MAT)
  )
  train_data <- apply_standardizers(train_data, std_f)
  test_data  <- apply_standardizers(test_data,  std_f)

  # Stage 2 training: both magnitude models exclude Arid's own site-months,
  # mirroring production (source model never reaches Arid; sink model
  # extrapolates to Arid from Forest/Grassland/Shrubland instead of
  # learning from Arid's own — entirely source-labeled — months).
  sink_train_f   <- train_data %>% filter(EcoType != "Arid", weak_source == 0, monthly_flux_gC_m2_month <= 0)
  source_train_f <- train_data %>% filter(EcoType != "Arid", weak_source == 1, monthly_flux_gC_m2_month >  0)

  sink_model_f <- if (nrow(sink_train_f) >= min_fold_training_rows) {
    fit_rf_magnitude_model(sink_train_f, "Weak-sink", mag_predictors)
  } else {
    message(sprintf("  %s: only %d non-Arid Weak-sink training months (< %d) — sink model not fit for this fold.",
      fold_site, nrow(sink_train_f), min_fold_training_rows))
    NULL
  }

  if (is_arid_fold) {
    # Held-out site IS Arid (JORN or SRER): production always routes Arid
    # cells through the sink model regardless of the site's own observed
    # sign, so score it the same way — a genuine holdout, since this fold's
    # sink model never saw ANY Arid site's months (Arid is excluded from
    # Stage 2 training in every fold, not just this one).
    if (is.null(sink_model_f)) {
      message(sprintf("  %s: Arid fold has no sink model to force predictions through — fold skipped.", fold_site))
      return(NULL)
    }
    return(test_data %>%
      mutate(magnitude_model = "Weak-sink (forced, Arid)",
             predicted_flux_gC_m2_month = pmin(predict_rf_magnitude(sink_model_f, test_data), 0)) %>%
      select(SITE_ID, EcoType, Year, month, is_arid, weak_source, magnitude_model,
             monthly_flux_gC_m2_month, predicted_flux_gC_m2_month))
  }

  # Held-out site is NOT Arid: route each test month through whichever
  # model matches its own observed sign, as before.
  source_model_f <- if (nrow(source_train_f) >= min_fold_training_rows) {
    fit_rf_magnitude_model(source_train_f, "Weak-source", mag_predictors)
  } else {
    message(sprintf("  %s: only %d non-Arid Weak-source training months (< %d) — source model not fit for this fold.",
      fold_site, nrow(source_train_f), min_fold_training_rows))
    NULL
  }

  sink_test_f   <- test_data %>% filter(weak_source == 0, monthly_flux_gC_m2_month <= 0)
  source_test_f <- test_data %>% filter(weak_source == 1, monthly_flux_gC_m2_month >  0)

  preds <- list()
  if (nrow(sink_test_f) > 0 && !is.null(sink_model_f)) {
    preds$sink <- sink_test_f %>%
      mutate(magnitude_model = "Weak-sink",
             predicted_flux_gC_m2_month = pmin(predict_rf_magnitude(sink_model_f, sink_test_f), 0))
  } else if (nrow(sink_test_f) > 0) {
    message(sprintf("  %s: %d Weak-sink test months skipped (no fold-specific sink model).", fold_site, nrow(sink_test_f)))
  }
  if (nrow(source_test_f) > 0 && !is.null(source_model_f)) {
    preds$source <- source_test_f %>%
      mutate(magnitude_model = "Weak-source",
             predicted_flux_gC_m2_month = pmax(predict_rf_magnitude(source_model_f, source_test_f), 0))
  } else if (nrow(source_test_f) > 0) {
    message(sprintf("  %s: %d Weak-source test months skipped (no fold-specific source model).", fold_site, nrow(source_test_f)))
  }

  if (length(preds) == 0) return(NULL)
  bind_rows(preds) %>%
    select(SITE_ID, EcoType, Year, month, is_arid, weak_source, magnitude_model,
           monthly_flux_gC_m2_month, predicted_flux_gC_m2_month)
}

# ── Run all folds ───────────────────────────────────────────────────────────

loso_predictions <- map(seq_along(sites), function(i) {
  s <- sites[i]
  message(sprintf("Fold %d/%d: holding out %s...", i, length(sites), s))
  run_loso_fold(
    train_data = monthly_training %>% filter(SITE_ID != s),
    test_data  = monthly_training %>% filter(SITE_ID == s),
    fold_site  = s
  )
}) %>% bind_rows()

message(sprintf("LOSO complete: %d held-out site-months scored across %d sites.",
  nrow(loso_predictions), n_distinct(loso_predictions$SITE_ID)))

# ── Skill summaries ───────────────────────────────────────────────────────────

# Lin's concordance correlation coefficient (CCC) and its precision/accuracy
# decomposition. CCC = r * C_b (accuracy factor <= 1); slope_pred_obs < 1 flags
# compression of predictions toward the mean, which r alone cannot detect.
ccc_components <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]; predicted <- predicted[ok]
  mo <- mean(observed); mp <- mean(predicted)
  vo <- mean((observed - mo)^2); vp <- mean((predicted - mp)^2)
  cov_op <- mean((observed - mo) * (predicted - mp))
  r   <- if (vo > 0 && vp > 0) cov_op / sqrt(vo * vp) else NA_real_
  ccc <- 2 * cov_op / (vo + vp + (mo - mp)^2)
  list(
    ccc            = ccc,
    accuracy_cb    = if (!is.na(r) && r != 0) ccc / r else NA_real_,
    slope_pred_obs = if (vo > 0) cov_op / vo else NA_real_
  )
}

summarise_loso <- function(df, ...) {
  df %>%
    group_by(...) %>%
    summarise(
      n_observations                 = n(),
      n_sites                        = n_distinct(SITE_ID),
      rmse_gC_m2_month                = sqrt(mean((monthly_flux_gC_m2_month - predicted_flux_gC_m2_month)^2, na.rm = TRUE)),
      mae_gC_m2_month                 = mean(abs(monthly_flux_gC_m2_month - predicted_flux_gC_m2_month), na.rm = TRUE),
      bias_gC_m2_month                = mean(predicted_flux_gC_m2_month - monthly_flux_gC_m2_month, na.rm = TRUE),
      correlation_observed_predicted  = suppressWarnings(
        cor(monthly_flux_gC_m2_month, predicted_flux_gC_m2_month, use = "complete.obs")),
      ccc_observed_predicted          = ccc_components(monthly_flux_gC_m2_month, predicted_flux_gC_m2_month)$ccc,
      ccc_accuracy_cb                 = ccc_components(monthly_flux_gC_m2_month, predicted_flux_gC_m2_month)$accuracy_cb,
      slope_predicted_obs             = ccc_components(monthly_flux_gC_m2_month, predicted_flux_gC_m2_month)$slope_pred_obs,
      .groups = "drop"
    ) %>%
    mutate(evaluation_basis = "Leave-one-site-out (unseen site)", .after = 1)
}

loso_skill_overall    <- summarise_loso(loso_predictions, magnitude_model)
loso_skill_by_ecotype <- summarise_loso(loso_predictions, magnitude_model, EcoType)
loso_skill_by_site    <- summarise_loso(loso_predictions, magnitude_model, SITE_ID, EcoType)

# ── Compare against 12_SourceProp_MagnitudeModels.R's OOB skill, if present ───

oob_skill_file <- file.path(output_dir, "magnitude_model_skill.csv")
oob_skill <- if (file.exists(oob_skill_file)) read.csv(oob_skill_file) else NULL
if (is.null(oob_skill))
  message("Note: ", oob_skill_file, " not found — run 12_SourceProp_MagnitudeModels.R first ",
          "for a side-by-side OOB-vs-LOSO comparison. Proceeding with LOSO-only outputs.")

generalization_comparison <- loso_skill_overall %>%
  select(magnitude_model, evaluation_basis, n_observations, n_sites,
         rmse_gC_m2_month, mae_gC_m2_month,
         correlation_observed_predicted)

if (!is.null(oob_skill)) {
  generalization_comparison <- bind_rows(
    generalization_comparison,
    oob_skill %>%
      select(magnitude_model, evaluation_basis, n_observations,
             rmse_gC_m2_month, mae_gC_m2_month,
             correlation_observed_predicted = correlation_observed_fitted)
  )
}

# ── Narrative summary (regenerated from the tables above on every run) ────────

narrative_lines <- c(sprintf(
  "Leave-one-site-out (LOSO) cross-validation refits Stage 1 and Stage 2 on all NEON upland sites except one, predicts that held-out site's months, and repeats for every site (%d folds, %d sites total). %d of those sites cross the aridity_index < %d Arid threshold; mirroring 12_SourceProp_MagnitudeModels.R's production design, Arid site-months are excluded from Stage 1 and the source-magnitude model in EVERY fold (not just when an Arid site is held out), and when an Arid site itself is held out, its test months are scored against the fold's sink model directly (labeled \"Weak-sink (forced, Arid)\" in the predictions table), matching how 13_Global_SpatialUpscalingRF.R actually treats Arid cells.",
  length(sites), length(sites), n_arid_sites, arid_ai_threshold
))

for (ml in c("Weak-sink", "Weak-source")) {
  loso_row <- loso_skill_overall %>% filter(magnitude_model == ml)
  if (nrow(loso_row) == 0) next
  oob_row <- if (!is.null(oob_skill)) oob_skill %>% filter(magnitude_model == ml) else NULL
  if (!is.null(oob_row) && nrow(oob_row) > 0) {
    narrative_lines <- c(narrative_lines, sprintf(
      "%s: LOSO RMSE = %.2f g C m-2 mo-1 (r = %.2f, n = %d months across %d unseen sites), versus %.2f g C m-2 mo-1 for the same-site OOB estimate reported by 12_SourceProp_MagnitudeModels.R (%.1fx higher) -- the gap between these two numbers is the honest cost of extrapolating to a site the model never saw at all, as opposed to a held-out month from a site it did.",
      ml, loso_row$rmse_gC_m2_month, loso_row$correlation_observed_predicted,
      loso_row$n_observations, loso_row$n_sites,
      oob_row$rmse_gC_m2_month, loso_row$rmse_gC_m2_month / oob_row$rmse_gC_m2_month
    ))
  } else {
    narrative_lines <- c(narrative_lines, sprintf(
      "%s: LOSO RMSE = %.2f g C m-2 mo-1 (r = %.2f, n = %d months across %d unseen sites). OOB comparison unavailable -- run 12_SourceProp_MagnitudeModels.R first for the same-site benchmark.",
      ml, loso_row$rmse_gC_m2_month, loso_row$correlation_observed_predicted,
      loso_row$n_observations, loso_row$n_sites
    ))
  }
}

arid_loso       <- loso_skill_by_ecotype %>% filter(EcoType == "Arid")
arid_loso_sites <- loso_predictions %>% filter(EcoType == "Arid") %>% distinct(SITE_ID) %>% nrow()
if (nrow(arid_loso) > 0) {
  narrative_lines <- c(narrative_lines, sprintf(
    "Arid specifically (\"Weak-sink (forced, Arid)\" in the predictions table): LOSO RMSE = %.2f g C m-2 mo-1 across %d Arid site(s) held out in turn, scored against a sink model that never saw any Arid site's months (in this fold or any other) -- interpret with caution given only %d NEON sites cross the Arid threshold at all, so this reflects how well non-Arid sink behavior extrapolates to Arid conditions, not a within-class generalization estimate.",
    mean(arid_loso$rmse_gC_m2_month, na.rm = TRUE), arid_loso_sites, n_arid_sites
  ))
}

loso_narrative_text <- paste(
  vapply(narrative_lines, function(x) paste(strwrap(x, width = 90), collapse = "\n"), character(1)),
  collapse = "\n\n")

# ── Write outputs ─────────────────────────────────────────────────────────────

write.csv(loso_predictions,           file.path(output_dir, "magnitude_model_loso_predictions.csv"), row.names = FALSE)
write.csv(loso_skill_overall,         file.path(output_dir, "magnitude_model_loso_skill_overall.csv"), row.names = FALSE)
write.csv(loso_skill_by_ecotype,      file.path(output_dir, "magnitude_model_loso_skill_by_ecotype.csv"), row.names = FALSE)
write.csv(loso_skill_by_site,         file.path(output_dir, "magnitude_model_loso_skill_by_site.csv"), row.names = FALSE)
write.csv(generalization_comparison,  file.path(output_dir, "magnitude_model_generalization_comparison.csv"), row.names = FALSE)
writeLines(loso_narrative_text,       file.path(output_dir, "magnitude_model_generalization_summary.txt"))

# ── Figure: observed vs LOSO-predicted, Weak-sink / Weak-source ───────────────

fig_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = "bottom", strip.background = element_rect(fill = "grey92"),
        axis.title = element_text(size = 11), axis.text = element_text(size = 10))

ecotype_colors <- c(Forest = "#1B7837", Grassland = "#D9B86C",
                    Shrubland = "#C2A5CF", Arid = "#D95F02")

# Only "Weak-sink"/"Weak-source" (non-Arid) predictions are plotted here.
# Arid's forced-sink LOSO check ("Weak-sink (forced, Arid)") is a different
# kind of comparison — non-Arid sink behavior extrapolated to Arid
# conditions, not within-class generalization — and is reported separately
# in loso_skill_by_ecotype.csv and the narrative summary instead.
make_loso_plot <- function(model_label, tag) {
  df    <- loso_predictions %>% filter(magnitude_model == model_label, !is.na(EcoType))
  skill <- loso_skill_overall %>% filter(magnitude_model == model_label)
  df %>%
    ggplot(aes(x = monthly_flux_gC_m2_month, y = predicted_flux_gC_m2_month, color = EcoType)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "grey35", linetype = "dashed", linewidth = 0.5) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = ecotype_colors) +
    labs(title = sprintf("%s. %s model  (n = %d, %d unseen sites, RMSE = %.2f)",
                         tag, model_label, nrow(df), skill$n_sites[1], skill$rmse_gC_m2_month[1]),
         color = NULL,
         x = "Observed (g C m⁻² mo⁻¹)", y = "LOSO-predicted (g C m⁻² mo⁻¹, unseen site)") +
    fig_theme
}

fig_loso <- plot_grid(
  make_loso_plot("Weak-sink",   "A"),
  make_loso_plot("Weak-source", "B"),
  ncol = 2
)

ggsave(file.path(figure_dir, "Fig_LOSO_magnitude_models.png"),
  fig_loso, width = 9, height = 4.5, units = "in", dpi = 300, bg = "white")

message("Leave-one-site-out generalization test complete. Outputs written to: ", output_dir)
message(loso_narrative_text)
