# 04_NEON_TotalFlux.R
# Total CH4 flux = ensemble gradient flux + storage flux.
# Storage flux (lterwg-flux-gradient: flow.neon.storage.R -> flow.neon.storage.compile.R), one file per NEON period:
#   NEON_Storage_Flux/NEON_storage.RData                 Aug 2021 - Jun 2024
#   NEON_Storage_Flux_202407_202608/NEON_storage.RData   Jul 2024 - Aug 2026 (Oct 2026 run)
# Both are stacked per site (CH4 only); where a half-hour appears in both, the newer period is kept.

storage.files <- c(file.path(localdir, "NEON_Storage_Flux", "NEON_storage.RData"),
                   file.path(localdir, "NEON_Storage_Flux_202407_202608", "NEON_storage.RData"))

storage.ch4 <- list()
for (k in seq_along(storage.files)) {
  e <- new.env(); load(storage.files[k], envir = e)
  for (s in names(e$NEON_storage)) {
    storage.ch4[[s]] <- rbind(storage.ch4[[s]],
      e$NEON_storage[[s]] %>% filter(gas == "CH4") %>%
        transmute(time.rounded = time_halfhour_storage_end, storage_flux_filled, period = k))
  }
  rm(e); invisible(gc())
}
storage.ch4 <- lapply(storage.ch4, function(x) x %>% arrange(time.rounded, desc(period)) %>%
                        distinct(time.rounded, .keep_all = TRUE) %>% select(-period))

load(fs::path(localdir.ch4, "SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata"))
site.list <- names(SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE)

SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage <- list()
for (site in site.list) {
  flux.site <- SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]]
  flux.storage.site <- storage.ch4[[site]]
  if (is.null(flux.storage.site)) flux.storage.site <- data.frame(time.rounded = as.POSIXct(character(0), tz = "UTC"), storage_flux_filled = numeric(0))
  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage[[site]] <- flux.site %>%
    full_join(flux.storage.site, by = "time.rounded") %>%
    mutate(flux_total = case_when(!is.na(FG_ENSEMBLE_RSHP) ~ FG_ENSEMBLE_RSHP + storage_flux_filled, TRUE ~ NA))
  print(paste0("Done with ", site))
}

save(SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage,
     file = fs::path(localdir.ch4, "SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata"))
