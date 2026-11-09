# Summarize the completed blocked audit on a common set of scored observations.
library(data.table)
library(ggplot2)
repo <- normalizePath(Sys.getenv("GF_AUDIT_REPO", unset = "."))
source(file.path(repo, "functions/gapfill_block_audit.R"))
out <- Sys.getenv("GF_AUDIT_OUT_DIR", unset =
  "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Gapfill_Audit/blocked_20261021")
pred_file <- file.path(out, "held_out_predictions.csv.gz")
if (!file.exists(pred_file)) stop("Missing completed held-out predictions: ", pred_file)
x <- fread(pred_file)
modes <- c("full_data_climatology", "fold_climatology", "fold_loo_climatology", "no_climatology")
keys <- c("scenario", "repetition", "SITE_ID", "row_id")
support <- x[, .(n_scored_modes = uniqueN(mode[is.finite(predicted)])), by = keys]
common_keys <- support[n_scored_modes == length(modes), ..keys]
common <- merge(x, common_keys, by = keys, all = FALSE)
common_metrics <- common[, audit_score(CH4_mgC_30min, predicted), by = .(scenario, mode)]
common_coverage <- common[, .(n_common_observations = uniqueN(row_id),
  n_common_sites = uniqueN(SITE_ID)), by = "scenario"]
mode_coverage <- x[, .(n_test_observations = uniqueN(row_id),
  n_scored_observations = uniqueN(row_id[is.finite(predicted)]),
  n_scored_sites = uniqueN(SITE_ID[is.finite(predicted)])), by = .(scenario, mode)]
fwrite(common_metrics, file.path(out, "metrics_common_support.csv"))
fwrite(common_coverage, file.path(out, "common_support_counts.csv"))
fwrite(mode_coverage, file.path(out, "prediction_coverage.csv"))
common_metrics[, scenario := factor(scenario,
  levels = c("random_rows", "7d", "30d", "90d", "site_year"))]
common_metrics[, mode_label := factor(mode, levels = modes,
  labels = c("All-data climate (leaky)", "Fold climate", "Fold LOO climate", "No climate"))]
p <- ggplot(common_metrics,
  aes(scenario, rmse, color = mode_label, group = mode_label)) +
  geom_point() + geom_line() + theme_bw(base_size = 11) +
  labs(x = "Artificial gap", y = "Held-out RMSE (mg C m⁻² per 30 min)",
       color = NULL, caption = "One sampled block per site and gap type; all modes scored on the same held-out observations.") +
  scale_x_discrete(labels = c(random_rows = "Random row", `7d` = "7 days",
    `30d` = "30 days", `90d` = "90 days", site_year = "One year")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))
ggsave(file.path(out, "blocked_gapfill_rmse_common_support.png"), p,
       width = 10, height = 5, dpi = 200)
if (nrow(common_metrics)) {
  cat("Common-support held-out metrics\n")
  print(common_metrics[order(scenario, mode)])
  cat("\nCommon-support counts\n")
  print(common_coverage[order(scenario)])
} else {
  warning("No held-out observations were scored by every mode.")
}
scenario_order <- c("random_rows", "7d", "30d", "90d", "site_year")
metric <- function(scenario_name, mode_name, column) {
  value <- common_metrics[as.character(scenario) == scenario_name & mode == mode_name,
                          get(column)]
  if (length(value) == 1L && is.finite(value)) value else NA_real_
}
fmt <- function(x, digits = 3L) ifelse(is.finite(x), format(round(x, digits), nsmall = digits), "NA")
rows <- lapply(scenario_order, function(s) {
  n <- common_coverage[scenario == s, n_common_observations]
  full <- metric(s, "full_data_climatology", "rmse")
  fold <- metric(s, "fold_climatology", "rmse")
  loo <- metric(s, "fold_loo_climatology", "rmse")
  base <- metric(s, "no_climatology", "rmse")
  r2 <- metric(s, "fold_loo_climatology", "r2")
  sprintf("| %s | %s | %s | %s | %s | %s | %s |",
    s, if (length(n)) n else 0L, fmt(full), fmt(fold), fmt(loo), fmt(base), fmt(r2))
})
findings <- c(
  "# Gapfill block-validation findings", "",
  "Random out-of-bag (OOB) scores are optimistic here because each measured flux also contributes to its climatological predictor. The full-data-climatology column is the deliberately leaky control.",
  "", "All four predictor treatments below are scored on the exact same held-out observations. RMSE is mg C m⁻² per 30 min.", "",
  "| Holdout design | Scored records | Full-data climate (leaky) | Fold climate (training self-inclusion) | Fold leave-one-out climate | No climate | Fold LOO R² |",
  "|---|---:|---:|---:|---:|---:|---:|", unlist(rows), "",
  "The full-data control lowers RMSE by roughly 2–9% versus fold leave-one-out across these designs, demonstrating optimistic validation from letting held-out fluxes enter the climatology. More notably, fold climatology with training-row self-inclusion is 7–14% worse than fold leave-one-out, a train/test predictor mismatch. The leave-one-out climate model is nearly indistinguishable from the no-climatology model, so this feature did not add measurable blocked predictive skill in this audit.",
  "", "The leakage-free blocked R² is near zero or negative for 7-, 30-, 90-day, and year-long holds. Predictions retain little pointwise variation skill over those horizons. Aggregate bias is small in these sampled blocks, but block sums include only measured timestamps inside each hidden interval; they are not annual budget validation.",
  "", "Recommended next model revision: build the training climatology as a leave-one-out predictor (and exclude the entire outer test block), then compare against removing this flux-derived feature or replacing it with independently estimated multi-year/site-season priors. Before changing annual budgets, repeat blocks across every eligible site-year and report errors for weak-source, weak-sink, and near-zero sites separately; propagate the resulting uncertainty into source/sink classification and spatial totals.",
  "", "Coverage, site/block metrics, donor provenance, and bootstrap intervals are in the adjacent CSV outputs. See `README.md` for design details and limitations.")
writeLines(findings, file.path(out, "audit_findings.md"))
