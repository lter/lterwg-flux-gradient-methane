# Independent retrospective audit; does not replace production flux products.
library(data.table)
library(ranger)
library(lubridate)
library(jsonlite)
library(ggplot2)
repo <- normalizePath(Sys.getenv("GF_AUDIT_REPO", unset = "."))
source(file.path(repo, "functions/gapfill_block_audit.R"))
root <- Sys.getenv("LOCALDIR_CH4", unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE")
out <- Sys.getenv("GF_AUDIT_OUT_DIR", unset = file.path(root, "Gapfill_Audit/blocked_20261021"))
cache <- "/private/tmp/methane_gapfill_audit_cache"
local_logs <- file.path(repo, "run_logs/20261021_gapfill_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "checkpoints"), showWarnings = FALSE)
dir.create(local_logs, recursive = TRUE, showWarnings = FALSE)
trees <- as.integer(Sys.getenv("GF_AUDIT_TREES", unset = "500"))
repeats <- as.integer(Sys.getenv("GF_AUDIT_REPEATS", unset = "3"))
scenarios <- strsplit(Sys.getenv("GF_AUDIT_SCENARIOS",
  unset = "random_rows,1d,7d,30d,90d,site_year"), ",", fixed = TRUE)[[1]]
modes <- c("full_data_climatology", "fold_climatology", "fold_loo_climatology", "no_climatology")
stopifnot(trees >= 50, repeats >= 1,
          all(scenarios %in% c("random_rows", "1d", "7d", "30d", "90d", "site_year")))
signature <- list(trees = trees, repeats = repeats, scenarios = scenarios, modes = modes,
  seed = 42L, design_seed = 20261021L, mtry_pop = 4L, mtry_site = 5L,
  ranger = as.character(packageVersion("ranger")),
  code = unname(tools::md5sum(c(file.path(repo, "functions/gapfill_block_audit.R"),
                             file.path(repo, "workflows/check/NEON.ERA5.BlockedGapfillAudit.R")))))
input_names <- c("NEON_ERA5_gapfilled_30min.csv.gz", "30min_ch4_model_data.csv")
for (name in input_names) if (!file.exists(file.path(cache, name)))
  stop("Stage input locally with run_logs/20261019_noEC/stage_input.py: ", name)
signature$inputs <- lapply(input_names, function(n)
  fromJSON(file.path(cache, paste0(n, ".json"))))
write_json(signature, file.path(out, "audit_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
status <- list(status = "running", started = format(Sys.time(), tz = "America/New_York"),
               output_dir = out, completed_jobs = 0L)
save_status <- function() {
  write_json(status, file.path(local_logs, "audit_status.json"), pretty = TRUE, auto_unbox = TRUE)
  write_json(status, file.path(out, "audit_status.json"), pretty = TRUE, auto_unbox = TRUE)
}
save_status()
options(error = function() {
  status$status <- "failed"; status$error <- geterrmessage(); save_status()
  quit(save = "no", status = 1L)
})
atomic_rds <- function(x, file) {
  tmp <- tempfile(fileext = ".rds"); on.exit(unlink(tmp))
  saveRDS(x, tmp)
  if (!file.copy(tmp, file, overwrite = TRUE)) stop("Cannot save checkpoint: ", file)
}
weather <- c("ERA5_Tair_C", "ERA5_Tsoil_C", "ERA5_VSWC", "ERA5_VSWC_deep",
  "ERA5_precip_mm", "ERA5_SW_Wm2", "ERA5_VPD_kPa", "tair_roll7", "precip_roll7", "Tair_VSWC")
predictors <- c(weather, "climatological_CH4_mgC_30min", "sin_hour", "cos_hour",
                "doy", "month", "season", "EcoType")
prepared <- file.path(cache, "prepared_observations.rds")
prep_signature <- list(inputs = signature$inputs, version = 1L)
saved <- if (file.exists(prepared)) readRDS(prepared) else NULL
if (!is.null(saved) && identical(saved$signature, prep_signature)) {
  d <- saved$d; raw <- saved$raw
} else {
  message("Reading observed targets and weather-only features from the frozen production output.")
  cols <- c("SITE_ID", "Year", "time.rounded", "month", "doy", "hour_num", "season", "EcoType",
            weather, "CH4_mgC_30min", "observed_flux", "climatological_CH4_mgC_30min")
  d <- fread(file.path(cache, input_names[1]), select = cols)
  d <- d[observed_flux == TRUE & is.finite(CH4_mgC_30min)]
  d[, `:=`(SITE_ID = as.character(SITE_ID), season = as.character(season),
            EcoType = as.character(EcoType), time.rounded = as.POSIXct(time.rounded, tz = "UTC"))]
  raw <- fread(file.path(cache, input_names[2]),
               select = c("SITE_ID", "time.rounded", "CH4_mgC_30min", "season"))
  raw[, time.rounded := ymd_hms(time.rounded, tz = "UTC", truncated = 3, quiet = TRUE)]
  raw <- raw[is.finite(CH4_mgC_30min) & !is.na(time.rounded)]
  raw[, `:=`(SITE_ID = as.character(SITE_ID), season = as.character(season),
               month = month(time.rounded), hour_num = hour(time.rounded) + minute(time.rounded)/60)]
  stopifnot(!anyDuplicated(raw[, .(SITE_ID, time.rounded)]),
            !anyDuplicated(d[, .(SITE_ID, time.rounded)]))
  joined <- raw[d, on = .(SITE_ID, time.rounded), .(raw_y = x.CH4_mgC_30min, donor_season = x.season)]
  stopifnot(all(is.finite(joined$raw_y)), max(abs(joined$raw_y - d$CH4_mgC_30min)) < 1e-10)
  d[, clim_donor_season := joined$donor_season]
  d[, `:=`(sin_hour = sin(2*pi*hour_num/24), cos_hour = cos(2*pi*hour_num/24))]
  setorder(d, SITE_ID, time.rounded); d[, row_id := .I]
  saveRDS(list(signature = prep_signature, d = d, raw = raw), prepared)
}
message("Observed evaluation rows: ", nrow(d), "; raw climatology donors: ", nrow(raw))
full_clim <- audit_fit_climatology(raw)
full_values <- audit_apply_climatology(d, full_clim)
available <- is.finite(d$climatological_CH4_mgC_30min)
stopifnot(max(abs(full_values$climatological_CH4_mgC_30min[available] -
                   d$climatological_CH4_mgC_30min[available])) < 1e-8)
message("Reconstructed full-data climatology agrees with the production predictor.")
design <- audit_make_design(d, scenarios, repeats)
blocks <- design$blocks
if (!nrow(blocks)) stop("No eligible test blocks.")
block_export <- copy(blocks); block_export[, n_test_obs := lengths(row_ids)]; block_export[, row_ids := NULL]
fwrite(block_export, file.path(out, "test_block_design.csv"))
fwrite(design$eligibility, file.path(out, "block_eligibility.csv"))
gaps <- d[, .(gap_days = pmax(as.numeric(diff(time.rounded), units = "days") - 1/48, 0)), by = SITE_ID]
gaps <- gaps[gap_days > 1e-8]
gap_summary <- gaps[, .(n_internal_gaps = .N, median_days = median(gap_days),
  p90_days = quantile(gap_days, .9), p99_days = quantile(gap_days, .99), max_days = max(gap_days)), by = SITE_ID]
fwrite(gap_summary, file.path(out, "observed_internal_gap_lengths.csv"))
eco_levels <- sort(unique(na.omit(d$EcoType)))
seasons <- c("Winter", "Spring", "Summer", "Autumn")
predictions <- list(); diagnostics <- list(); jobs <- 0L
for (scenario in scenarios) for (rep in seq_len(repeats)) {
  scenario_value <- scenario
  repetition_value <- rep
  selected <- blocks[blocks$scenario == scenario_value &
                       blocks$repetition == repetition_value]
  if (!nrow(selected)) next
  test_ids <- sort(unique(unlist(selected$row_ids)))
  test <- copy(d[row_id %in% test_ids]); train <- copy(d[!row_id %in% test_ids])
  stopifnot(!length(intersect(train$row_id, test$row_id)))
  # Exclude every raw flux donor in each test interval, not only rows that
  # happened to survive model eligibility/weather filtering.
  raw_keep <- rep(TRUE, nrow(raw))
  if (scenario == "random_rows") {
    withheld <- raw[test, on = .(SITE_ID, time.rounded), which = TRUE]
    raw_keep[withheld[!is.na(withheld)]] <- FALSE
  } else for (i in seq_len(nrow(selected))) {
    b <- selected[i]
    raw_keep <- raw_keep & !(raw$SITE_ID == b$SITE_ID & raw$time.rounded >= b$start & raw$time.rounded < b$end)
  }
  train_raw <- raw[raw_keep]
  stopifnot(!nrow(train_raw[test, on = .(SITE_ID, time.rounded), nomatch = 0L]))
  fold_clim <- audit_fit_climatology(train_raw)
  for (mode in modes) {
    jobs <- jobs + 1L; job <- paste(scenario, rep, mode, sep = "_")
    checkpoint <- file.path(out, "checkpoints", paste0(job, ".rds"))
    status$current <- job; status$completed_jobs <- jobs - 1L; save_status()
    prior <- if (file.exists(checkpoint)) readRDS(checkpoint) else NULL
    if (!is.null(prior) && identical(prior$signature, signature)) {
      message("Reusing ", job); predictions[[job]] <- prior$prediction; diagnostics[[job]] <- prior$diagnostic
      next
    }
    message(format(Sys.time()), " FIT ", job, " (", nrow(train), " train; ", nrow(test), " held out)")
    tr <- copy(train); te <- copy(test)
    has_clim <- mode != "no_climatology"
    if (has_clim) {
      lookup <- if (mode == "full_data_climatology") full_clim else fold_clim
      cv_train <- audit_apply_climatology(tr, lookup, leave_one_out = mode == "fold_loo_climatology")
      cv_test <- audit_apply_climatology(te, lookup)
      tr[, climatological_CH4_mgC_30min := cv_train$climatological_CH4_mgC_30min]
      te[, climatological_CH4_mgC_30min := cv_test$climatological_CH4_mgC_30min]
    } else {
      cv_test <- data.table(climatology_level = rep("not_used", nrow(te)), climatology_donors = 0L)
    }
    base_ok <- function(x) is.finite(x$ERA5_Tair_C) & is.finite(x$ERA5_VSWC) &
      is.finite(x$hour_num) & is.finite(x$doy) & !is.na(x$season) & !is.na(x$EcoType)
    tr_ok <- base_ok(tr); te_ok <- base_ok(te)
    if (has_clim) {
      tr_ok <- tr_ok & is.finite(tr$climatological_CH4_mgC_30min)
      te_ok <- te_ok & is.finite(te$climatological_CH4_mgC_30min)
    }
    tr <- tr[tr_ok]
    for (name in c("season", "EcoType")) {
      lev <- if (name == "season") seasons else eco_levels
      set(tr, j = name, value = factor(tr[[name]], levels = lev))
      set(te, j = name, value = factor(te[[name]], levels = lev))
    }
    site_levels <- sort(unique(tr$SITE_ID))
    tr[, SITE_ID := factor(SITE_ID, levels = site_levels)]
    te[, SITE_ID := factor(SITE_ID, levels = site_levels)]
    p <- if (has_clim) predictors else setdiff(predictors, "climatological_CH4_mgC_30min")
    stopifnot(nrow(tr) >= 500)
    pop <- ranger(reformulate(p, "CH4_mgC_30min"), data = tr,
      num.trees = trees, mtry = 4L, min.node.size = 5L, sample.fraction = .7,
      num.threads = 4L, seed = 42L, write.forest = TRUE)
    site <- ranger(reformulate(c(p, "SITE_ID"), "CH4_mgC_30min"), data = tr,
      num.trees = trees, mtry = 5L, min.node.size = 5L, sample.fraction = .7,
      num.threads = 4L, seed = 42L, write.forest = TRUE)
    counts <- table(tr$SITE_ID)
    n_site_train <- as.integer(counts[as.character(test$SITE_ID)]); n_site_train[is.na(n_site_train)] <- 0L
    route_site <- te_ok & n_site_train >= 50L & !is.na(te$SITE_ID)
    pred <- rep(NA_real_, nrow(te))
    if (any(te_ok)) pred[te_ok] <- predict(pop, data = te[te_ok, ..p])$predictions
    pp <- c(p, "SITE_ID")
    if (any(route_site)) pred[route_site] <- predict(site, data = te[route_site, ..pp])$predictions
    result <- test[, .(row_id, SITE_ID, time.rounded, Year, CH4_mgC_30min)]
    result[, `:=`(scenario = scenario, repetition = rep, mode = mode, predicted = pred,
      prediction_route = ifelse(!te_ok, "unsupported_features", ifelse(route_site, "site_rf", "population_rf")),
      n_site_training = n_site_train, climatology_level = cv_test$climatology_level,
      climatology_donors = cv_test$climatology_donors)]
    diagnostic <- data.table(job, scenario, repetition = rep, mode,
      n_train_labels = nrow(tr), n_train_climatology_donors = if (mode == "full_data_climatology") nrow(raw) else nrow(train_raw),
      n_test = nrow(test), n_scored = sum(is.finite(pred)), n_site_route = sum(route_site))
    atomic_rds(list(signature = signature, prediction = result, diagnostic = diagnostic), checkpoint)
    predictions[[job]] <- result; diagnostics[[job]] <- diagnostic
    rm(pop, site); gc(verbose = FALSE)
  }
}
all <- rbindlist(predictions); diagnostic <- rbindlist(diagnostics)
fwrite(all, file.path(out, "held_out_predictions.csv.gz"))
fwrite(diagnostic, file.path(out, "fold_diagnostics.csv"))
metrics <- all[, audit_score(CH4_mgC_30min, predicted), by = .(scenario, mode)]
by_site <- all[, audit_score(CH4_mgC_30min, predicted), by = .(scenario, mode, SITE_ID)]
by_block <- all[, audit_score(CH4_mgC_30min, predicted), by = .(scenario, mode, repetition, SITE_ID)]
by_block[, block_observed_sign_agreement := sign(observed_sum_gC_m2) == sign(predicted_sum_gC_m2)]
fwrite(metrics, file.path(out, "metrics_by_gap.csv")); fwrite(by_site, file.path(out, "metrics_by_site.csv"))
fwrite(by_block, file.path(out, "metrics_by_block.csv"))
set.seed(20261021)
site_errors <- all[is.finite(predicted), .(n = .N, sse = sum((predicted-CH4_mgC_30min)^2),
  error_sum = sum(predicted-CH4_mgC_30min)), by = .(scenario, mode, SITE_ID)]
ci <- site_errors[, {
  draws <- replicate(500, { ii <- sample.int(.N, .N, replace = TRUE)
    c(rmse = sqrt(sum(sse[ii])/sum(n[ii])), bias = sum(error_sum[ii])/sum(n[ii])) })
  list(n_sites = .N, rmse_low = quantile(draws[1,], .025), rmse_high = quantile(draws[1,], .975),
       bias_low = quantile(draws[2,], .025), bias_high = quantile(draws[2,], .975))
}, by = .(scenario, mode)]
fwrite(ci, file.path(out, "site_cluster_bootstrap_intervals.csv"))
metrics[, scenario := factor(scenario, levels = scenarios)]
plot <- ggplot(metrics, aes(scenario, rmse, color = mode, group = mode)) +
  geom_point() + geom_line() + theme_bw(base_size = 11) +
  labs(x = "Artificial gap", y = "Held-out RMSE (mg C m⁻² per 30 min)", color = NULL,
       caption = "Only hidden observed intervals can be scored. This is not validation of whole-gap or annual budgets.") +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))
ggsave(file.path(out, "blocked_gapfill_rmse.png"), plot, width = 10, height = 5, dpi = 200)
report <- c("# Blocked gapfilling audit", "",
  "Retrospective reconstruction, not forecasting. Models can use observations outside the hidden interval in other years.",
  "Full-data climatology is a deliberately leaky control. Fold climatology excludes every hidden-block flux donor.",
  "Fold LOO additionally removes each training row's own response from its climatological predictor. No-climatology is a weather/site-only sensitivity.",
  "Weather-only cached features are held fixed. All forests are refit without test labels; mtry 4/5, 500 trees (unless a pilot override), seed 42, no audit-driven tuning.",
  "Minimum eight observations per eligible block; eligibility tables disclose excluded sparse sites. Different gap lengths can score different cohorts.",
  "Block sums cover hidden observed timestamps only, not the unobserved remainder of a block. They cannot establish annual budget truth.",
  "Confidence intervals resample sites. Repetitions are not independent new measurements. These tests do not validate FG instrumentation or storage corrections.", "",
  capture.output(print(metrics)), "", "See metrics_by_site.csv, metrics_by_block.csv, fold_diagnostics.csv and site_cluster_bootstrap_intervals.csv.")
writeLines(report, file.path(out, "README.md"))
status$status <- "complete"; status$completed_jobs <- jobs
status$finished <- format(Sys.time(), tz = "America/New_York"); save_status()
message("AUDIT COMPLETE: ", out)
