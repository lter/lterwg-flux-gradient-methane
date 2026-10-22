# Source conventions in the current validation formatters.
# SE-Sto: ppm profiles / umol fluxes; SE-Svb and US-Uaf: ppb / nmol.
validation_flux_to_umol <- function(df, site) {
  native_unit <- c("SE-Sto"="umol", "SE-Svb"="nmol", "US-Uaf"="nmol")
  if (length(site) != 1L || !site %in% names(native_unit)) stop("Unknown validation site: ", site)
  if ("CH4_flux_unit" %in% names(df)) {
    if (anyNA(df$CH4_flux_unit) || any(df$CH4_flux_unit != "umol m-2 s-1")) {
      stop("Unexpected CH4_flux_unit for ", site)
    }
    return(df) # Already normalized: never divide again.
  }
  if ("FG_total" %in% names(df)) {
    stop("Unmarked total-flux input for ", site, "; rerun 02b_VAL_TotalFlux.R from the native ensemble.")
  }
  columns <- intersect(c("FG_mean", "EC_mean"), names(df))
  if (!all(c("FG_mean", "EC_mean") %in% columns)) stop("Missing gradient or EC flux for ", site)
  divisor <- if (native_unit[[site]] == "nmol") 1000 else 1
  df[columns] <- lapply(df[columns], function(x) x / divisor)
  df$CH4_flux_unit <- rep("umol m-2 s-1", nrow(df))
  df
}
