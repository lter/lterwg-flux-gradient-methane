# 02b_VAL_TotalFlux.R
#
# Combines ensemble gradient flux (FG_mean) with single-point storage flux to
# produce a total CH4 flux (FG_total) comparable to EC_mean, for the three
# validation towers (SE-Sto, SE-Svb, US-Uaf). Adapted from 04_NEON_TotalFlux.R
# in the main workflow: same full_join-by-time.rounded-and-add pattern, but
# the storage input is VAL_StorageFlux_SinglePoint.R's single-point
# approximation (in lterwg-flux-gradient) rather than the NEON-style
# multi-level column integration -- see that script for why (the multi-level
# method never clears its 3-level minimum at SE-Sto/SE-Svb).
#
# Runs as step 02b: right after 02_VAL_SensorHeightPairs.R (which produces
# SITEval_DATA_FILTERED_RSHP_EnSEMBLE.Rdata) and before 05_VAL_FluxAnalysis.R
# / 07_VAL_Supplement.R, both of which load this step's output and use
# FG_total in place of FG_mean wherever they compare against EC_mean (EC
# inherently includes storage; FG does not until this step has run).
# (03_VAL_DielAnalysis.R / 04_VAL_Figures.R were removed from the workflow --
# see 03_VAL_DielAnalysis.R for why.)
#
# Fallback rule (per project decision): when no single-point storage
# estimate is available for a given half-hour (concentration missing, gap
# >1 hr at that site), FG_total falls back to FG_mean alone rather than
# propagating NA. storage_added flags which case applies, so downstream
# analyses can subset to storage-corrected rows only if a stricter
# comparison is wanted.
#
# Units (per project decision): SITE_storage_flux_1pt.csv reports
# storage_flux_1pt in nmol CH4 m-2 s-1; FG_mean/EC_mean are treated as
# umol CH4 m-2 s-1 per 05_VAL_FluxAnalysis.R, so storage is divided by 1000
# before being added. Caveat: after this conversion storage is a small
# correction relative to FG_mean (median contribution well under 1% at all
# three sites) -- worth double-checking this units assumption against the
# raw NEON/AmeriFlux source before relying on it for the paper, since it
# makes the storage correction close to negligible.
#
# SE-Sto note: dry_air_concentration for SE-Sto was previously ~100x too
# small due to two raw-unit bugs (see apply_site_unit_corrections() in
# flow.validation.storage.R, lterwg-flux-gradient repo); this is now fixed
# upstream, so FG_total for SE-Sto reflects the corrected storage flux as
# long as Storage_Flux/SE-Sto_storage_flux.csv has been regenerated since
# the fix.
#
# Output: SITEval_DATA_FILTERED_RSHP_EnSEMBLE_TotalFlux.Rdata (list
# SITEval_DATA_FILTERED_RSHPc_H_total, by site) saved alongside
# SITEval_DATA_FILTERED_RSHP_EnSEMBLE.Rdata in Validation_Sites/.

library(tidyverse)

load(fs::path(localdir.ch4, "SITEval_DATA_FILTERED_RSHP_EnSEMBLE.Rdata"))
load(fs::path(localdir.ch4, "Storage_Flux/SinglePoint/Validation_storage_singlepoint.RData"))

site.list <- names(SITEval_DATA_FILTERED_RSHPc_H)

SITEval_DATA_FILTERED_RSHPc_H_total <- list()

for (site in site.list) {

  flux.site <- SITEval_DATA_FILTERED_RSHPc_H[[site]]

  flux.storage.site <- Validation_storage_singlepoint[[site]] %>%
    mutate(storage_flux_filled_umol = storage_flux_1pt / 1000) %>%
    select(time.rounded, storage_flux_filled_umol)

  merged <- flux.site %>%
    full_join(flux.storage.site, by = "time.rounded") %>%
    mutate(
      storage_added = !is.na(storage_flux_filled_umol),
      FG_total = case_when(
        !is.na(FG_mean) & storage_added ~ FG_mean + storage_flux_filled_umol,
        !is.na(FG_mean)                 ~ FG_mean,
        TRUE                            ~ NA_real_
      )
    )

  SITEval_DATA_FILTERED_RSHPc_H_total[[site]] <- merged

  pct_storage <- round(100 * mean(merged$storage_added[!is.na(merged$FG_mean)]), 1)
  print(paste0("Done with ", site, " (", pct_storage, "% of FG rows got a storage term)"))
}

save(SITEval_DATA_FILTERED_RSHPc_H_total,
     file = fs::path(localdir.ch4, "SITEval_DATA_FILTERED_RSHP_EnSEMBLE_TotalFlux.Rdata"))
