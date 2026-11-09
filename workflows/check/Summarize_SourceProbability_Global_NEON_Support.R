library(data.table)
rfdir <- Sys.getenv("MONTHLY_RF_DIR", unset="/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF")
out <- file.path(rfdir,"OUTPUT")
year_tag <- "2000-2010-2020-2025"
tag <- paste0("years_",year_tag,"_months_1-2-3-4-5-6-7-8-9-10-11-12")
univariate <- fread(file.path(out,paste0("source_probability_global_neon_univariate_support_",tag,".csv")))
joint <- fread(file.path(out,paste0("source_probability_global_neon_joint_support_",tag,".csv")))
weighted <- fread(file.path(out,"source_probability_supported_weighted_loso_summary.csv"))
fwrite(univariate,file.path(out,paste0("source_probability_global_neon_univariate_support_",year_tag,".csv")))
fwrite(joint,file.path(out,paste0("source_probability_global_neon_joint_support_",year_tag,".csv")))
lines <- c("# Global-versus-NEON predictor support and weighted validation", "",
  sprintf("The classifier's six Stage-1 inputs were compared within ecosystem: EcoType plus monthly ERA5-Land temperature and soil water, static MAP/MAT, and derived aridity index. The global sample covers all months of four representative projection years (%s), with 10,000 cell-month records per ecosystem when available. Global area summaries use cell area weights.", year_tag),
  "", "NEON tests transfer to held-out sites; its source prevalence is not used as the global prior.",
  "", "## Joint support", "", capture.output(print(joint)), "",
  "The multivariate flag uses Euclidean nearest-neighbor distance across five numeric inputs, scaled by within-ecosystem NEON IQR. The threshold is the 95th percentile of NEON leave-one-site-out distances. It is a transparent screening rule, not a calibrated probability of extrapolation.",
  "", "## Univariate support", "", capture.output(print(univariate)), "",
  "## Supported-only, site-held-out weighting", "", capture.output(print(weighted)), "",
  "The domain-density ratios are estimated fold by fold without the held-out site's predictors or outcomes, within the supported portion only. Global-domain mass uses area weights; NEON-domain mass gives each training site equal total weight. Ratios are capped at 10 before normalization. Results are averages of held-out-site metrics, with no shrubland reweighting.",
  "", "These weighted validation metrics describe only the supported global covariate subset. They cannot validate the excluded global area or establish global source prevalence. Effective sample sizes and capped-weight fractions show material instability, especially for grassland; treat weighted calibration as a sensitivity diagnostic.")
writeLines(lines,file.path(out,paste0("source_probability_global_neon_support_",year_tag,".md")))
print(joint); print(weighted)
