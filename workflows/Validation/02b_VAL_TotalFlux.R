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
# Units (corrected Oct 2026): SITE_storage_flux_1pt.csv reports storage_flux_1pt in nmol CH4 m-2 s-1 at
# all three sites (SE-Sto's ppm profile is converted to ppb upstream in flow.validation.storage.R), so storage
# is divided by 1000. FG_mean and EC_mean are converted to umol CH4 m-2 s-1 per site with val_ch4_to_umol
# below (SE-Svb and US-Uaf were in nmol); before this, nmol FG was added to umol storage at those two sites,
# which is why storage looked negligible there.
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

# Validation-site CH4 flux units (checked Oct 2026 against the data): the harmonized SE-Svb and US-Uaf
# profiles are in ppb, so their FG_mean is nmol CH4 m-2 s-1, and their EC fluxes (ch4_flux_nmolm2s, CH4.flux)
# are nmol too; SE-Sto's profiles and EC (Fch4) are in ppm / umol. Everything is converted to umol CH4 m-2 s-1.
if (!exists("val_ch4_to_umol")) val_ch4_to_umol <- c("SE-Sto" = 1, "SE-Svb" = 1e-3, "US-Uaf" = 1e-3)

site.list <- names(SITEval_DATA_FILTERED_RSHPc_H)

SITEval_DATA_FILTERED_RSHPc_H_total <- list()

for (site in site.list) {

  if (!site %in% names(val_ch4_to_umol)) stop("No CH4 unit factor for validation site ", site)
  flux.site <- SITEval_DATA_FILTERED_RSHPc_H[[site]] %>%
    mutate(FG_mean = FG_mean * val_ch4_to_umol[[site]],
           EC_mean = EC_mean * val_ch4_to_umol[[site]],
           flux_units = "umol m-2 s-1")

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
