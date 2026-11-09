# 04_NEON_TotalFlux.R
# Total CH4 flux = ensemble gradient flux + storage flux.
# Storage flux: v3b release (lterwg-flux-gradient, recomputed Nov 2026), one CSV per site covering Aug 2021 - Aug 2026:
#   NEON_GradientFlux_EDI_2021_2026_v3b/<SITE>/<SITE>_STORAGE_FLUX_30min.csv   (timestamps in UTC)
# CH4 rows only; storage_flux_filled (native sign, gap-filled) is joined on the interval-end label
# time_halfhour_storage_end, as before. (Previously: NEON_Storage_Flux/ and NEON_Storage_Flux_202407_202608/ RData.)
library(data.table)
# storage.dir may also be several folders in time order (e.g. the two per-period EDI folders on Bouchet,
# NEON_GradientFlux_EDI_202108_202406_v3b and _202407_202608_v3b); a half-hour present in both keeps the later folder.
if (!exists("storage.dir")) storage.dir <- file.path(localdir, "NEON_GradientFlux_EDI_2021_2026_v3b")

read_storage_ch4 <- function(site) {
  f <- file.path(storage.dir, site, paste0(site, "_STORAGE_FLUX_30min.csv"))
  f <- f[file.exists(f)]
  if (!length(f)) { message("  no storage file for ", site); return(NULL) }
  d <- rbindlist(lapply(seq_along(f), function(k)
    fread(f[k], select = c("gas", "time_halfhour_storage_end", "storage_flux_filled"), showProgress = FALSE)[gas == "CH4"][, period := k]))
  d <- d[, .(time.rounded = as.POSIXct(as.character(time_halfhour_storage_end), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
             storage_flux_filled = as.numeric(storage_flux_filled), period)]
  setorder(d, time.rounded, -period)
  d <- unique(d[!is.na(time.rounded)], by = "time.rounded")[, period := NULL]
  as.data.frame(d)
}

load(fs::path(localdir.ch4, "SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata"))
site.list <- names(SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE)

SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage <- list()
for (site in site.list) {
  flux.site <- SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]]
  flux.storage.site <- read_storage_ch4(site)
  if (is.null(flux.storage.site)) flux.storage.site <- data.frame(time.rounded = as.POSIXct(character(0), tz = "UTC"), storage_flux_filled = numeric(0))
  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage[[site]] <- flux.site %>%
    full_join(flux.storage.site, by = "time.rounded") %>%
    mutate(flux_total = case_when(!is.na(FG_ENSEMBLE_RSHP) ~ FG_ENSEMBLE_RSHP + storage_flux_filled, TRUE ~ NA))
  print(paste0("Done with ", site))
}

save(SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage,
     file = fs::path(localdir.ch4, "SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata"))
