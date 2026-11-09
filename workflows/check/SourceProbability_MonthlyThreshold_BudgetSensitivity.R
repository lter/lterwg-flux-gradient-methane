#!/usr/bin/env Rscript
# Apply monthly probability thresholds to the saved global cell-month predictions.
# Includes uniform 0.5, 0.8, and 0.85 scenarios plus a forest-only nested
# threshold with 0.5 fallbacks elsewhere. Uniform 0.85 is an explicit transfer
# assumption: nested NEON validation supports it for forest only. This diagnostic
# does not alter the production bundle or the continuous probability-weighted flux.
# Rerun after NEON gapfilled fluxes and site summaries are recalculated, then
# rerun step 13 before using refreshed fluxes for the global budget comparison.
suppressPackageStartupMessages(library(data.table))
rfdir <- Sys.getenv("MONTHLY_RF_DIR",unset="/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF")
out <- file.path(rfdir,"OUTPUT")
x <- readRDS(file.path(out,"monthly_cell_predictions_2000_2025.rds"))
setDT(x)
required <- c("Year","EcoType","source_prob_A","predicted_sink_flux","predicted_source_flux","area_mha")
stopifnot(all(required %in% names(x)))
gmb <- fread(file.path(out,"gmb_threshold_sensitivity.csv"))
calc <- function(thresholds,label) {
  z <- copy(x)
  z[,threshold:=unname(thresholds[as.character(EcoType)])]
  z[!is.finite(threshold),threshold:=.5]
  z[,flux:=fifelse(source_prob_A>=threshold,predicted_source_flux,predicted_sink_flux)]
  z[,tg:=flux*area_mha*0.0133333333333333]
  annual <- z[,.(annual_tg=sum(tg,na.rm=TRUE)),by=.(Year,EcoType)]
  annual[,scenario:=label]
  annual[]
}
all_ecos <- unique(as.character(x$EcoType))
uniform05 <- setNames(rep(.5,length(all_ecos)),all_ecos)
uniform08 <- setNames(rep(.8,length(all_ecos)),all_ecos)
uniform085 <- setNames(rep(.85,length(all_ecos)),all_ecos)
monthly <- setNames(rep(.5,length(all_ecos)),all_ecos)
monthly["Forest"] <- .85
annual <- rbindlist(list(calc(uniform05,"uniform_0.50"),
  calc(uniform08,"uniform_0.80"),
  calc(uniform085,"uniform_0.85"),
  calc(monthly,"forest_nested_0.85_other_ecosystems_fallback_0.50")))
fwrite(annual,file.path(out,"source_probability_monthly_threshold_annual_budgets_by_year_ecotype.csv"))
summary <- annual[,.(mean_annual_tg=mean(annual_tg),sd_annual_tg=sd(annual_tg),
  n_years=uniqueN(Year)),by=.(scenario,EcoType)]
annual_total <- annual[,.(annual_tg=sum(annual_tg)),by=.(Year,scenario)]
overall <- annual_total[,.(mean_annual_tg=mean(annual_tg),sd_annual_tg=sd(annual_tg),
  n_years=uniqueN(Year)),by=scenario]
fwrite(summary,file.path(out,"source_probability_monthly_threshold_annual_budgets_summary.csv"))
fwrite(overall,file.path(out,"source_probability_monthly_threshold_annual_budgets_overall.csv"))
print(summary); print(overall); print(gmb[threshold %in% c(.5,.8,.85),])
