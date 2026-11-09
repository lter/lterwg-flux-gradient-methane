# Extend the monthly NEON soil-moisture covariates to cover the current fluxes.
# Same aggregation as depreciaded/flow.soilmoisture.NEON.R: 30-minute SWS,
# verticalPosition 501, all soil plots, with no additional QC selection.
# Newer months may be provisional; record each source URL, release and checksum.
library(data.table)
library(neonUtilities)

root <- Sys.getenv("LOCALDIR_CH4", "/Volumes/MaloneLab/Research/FluxGradient/METHANE")
repo <- Sys.getenv("DIRREPO_CH4", getwd())
cache <- Sys.getenv("SOIL_CACHE_DIR", file.path(dirname(repo), "NEON_SoilMoisture_cache"))
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "OUTPUT"), recursive = TRUE, showWarnings = FALSE)
for (f in unique(c(path.expand("~/.Renviron"), file.path(repo, ".Renviron"))))
  if (file.exists(f)) readRenviron(f)
token <- Sys.getenv("NEON_TOKEN")
if (!nzchar(token)) stop("NEON_TOKEN is required. Configure it in the git-ignored .Renviron file.")

soil_file <- file.path(root, "Soildata_YearMon.Rdata")
load(soil_file)
old <- as.data.table(Site_SoilData)
flux <- new.env()
load(file.path(root, "SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata"), flux)
needed <- rbindlist(lapply(names(flux$SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage), function(s) {
  d <- flux$SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage[[s]]
  months <- unique(format(d$time.rounded[is.finite(d$flux_total)], "%Y-%m"))
  if (!length(months)) return(data.table(Site = character(), YearMon = character()))
  data.table(Site = s, YearMon = months)
}))
rm(flux); invisible(gc())
start_month <- Sys.getenv("SOIL_START_MONTH")
end_month <- Sys.getenv("SOIL_END_MONTH")
if (nzchar(start_month) && nzchar(end_month)) {
  metadata_file <- file.path(Sys.getenv("LOCALDIR_FLUXGRADIENT", dirname(root)),
    "Ameriflux_NEON field-sites.csv")
  sites <- unique(read.csv(metadata_file)$Site_Id.NEON)
  months <- format(seq(as.Date(paste0(start_month, "-01")),
    as.Date(paste0(end_month, "-01")), by = "month"), "%Y-%m")
  needed <- unique(rbindlist(list(needed, CJ(Site = sites, YearMon = months))))
}
missing <- needed[!old[is.finite(VSWCMean)], on = .(Site, YearMon)]
message(nrow(missing), " site-months need soil-moisture data")
new_rows <- list(); manifest <- list(); availability <- list()
api <- get("getAPI", asNamespace("neonUtilities"))
product_response <- api("https://data.neonscience.org/api/v0/products/DP1.00094.001", token)
if (is.null(product_response) || product_response$status_code != 200)
  stop("Cannot read NEON soil-moisture availability")
product <- jsonlite::fromJSON(httr::content(product_response, "text", encoding = "UTF-8"))
available_sites <- product$data$siteCodes

download_cached <- function(url) {
  name <- basename(sub("[?].*$", "", url))
  f <- file.path(cache, name)
  if (!file.exists(f)) {
    tmp <- paste0(f, ".part")
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE, method = "libcurl")
    if (!file.rename(tmp, f)) stop("Cannot finalize ", f)
  }
  f
}

for (s in unique(missing$Site)) {
  months <- sort(missing[Site == s, YearMon])
  si <- match(s, available_sites$siteCode)
  published <- if (is.na(si)) character() else available_sites$availableMonths[[si]]
  months <- intersect(months, published)
  if (!length(months)) { message(s, ": no published measurements for the missing months"); next }
  message("Soil moisture: ", s, " ", min(months), " to ", max(months))
  qfile <- file.path(cache, paste0(s, "_", min(months), "_", max(months), "_query.rds"))
  if (file.exists(qfile)) {
    query <- readRDS(qfile)
  } else {
    query <- neonUtilities::queryFiles(
      dpID = "DP1.00094.001", site = s, startdate = min(months), enddate = max(months),
      package = "basic", timeIndex = "30", tabl = "SWS_30_minute", metadata = FALSE,
      include.provisional = TRUE, token = token)
    if (!is.null(query)) saveRDS(query, qfile)
  }
  if (is.null(query)) {
    stop("NEON file query failed for ", s, ". Check NEON_TOKEN and API availability; original soil data unchanged.")
  }
  files <- unique(as.data.table(query$filesall), by = "urlbase")
  files <- files[grepl("[.]501[.]", urlbase) & grepl("SWS_30_minute", urlbase)]
  files[, file_month := regmatches(urlbase, regexpr("[0-9]{4}-[0-9]{2}", urlbase))]
  files <- files[file_month %in% months]
  message("  ", nrow(files), " sensor files")
  downloads <- parallel::mclapply(files$urlbase, download_cached,
    mc.cores = if (.Platform$OS.type == "unix") 4L else 1L)
  if (any(vapply(downloads, inherits, logical(1), "try-error")))
    stop("One or more soil-moisture downloads failed for ", s)
  site_data <- list()
  for (i in seq_len(nrow(files))) {
    url <- files$urlbase[i]
    name <- basename(sub("[?].*$", "", url))
    parts <- strsplit(name, ".", fixed = TRUE)[[1]]
    month <- regmatches(name, regexpr("[0-9]{4}-[0-9]{2}", name))
    if (!month %in% months) next
    f <- downloads[[i]]
    expected <- files$md5[i]
    if (!is.na(expected) && nzchar(expected) &&
        !identical(tolower(unname(tools::md5sum(f))), tolower(expected)))
      stop("Checksum mismatch: ", name)
    d <- if (grepl("[.]parquet$", f)) as.data.table(arrow::read_parquet(f)) else fread(f)
    stopifnot(all(c("startDateTime", "VSWCMean", "VSWCMinimum", "VSWCMaximum", "VSWCVariance") %in% names(d)))
    d[, `:=`(Site = s, domain = parts[2],
      YearMon = format(as.POSIXct(startDateTime, tz = "UTC"), "%Y-%m"))]
    site_data[[length(site_data) + 1L]] <- d[, .(Site, domain, YearMon, VSWCMean,
      VSWCMinimum, VSWCMaximum, VSWCVariance)]
    manifest[[length(manifest) + 1L]] <- data.table(Site = s, YearMon = month,
      url = url, md5 = expected, release = files$release[i])
  }
  if (length(site_data)) {
    d <- rbindlist(site_data)
    agg <- d[, .(VSWCMean = mean(VSWCMean, na.rm = TRUE),
      VSWCMin = min(VSWCMinimum, na.rm = TRUE), VSWCMax = max(VSWCMaximum, na.rm = TRUE),
      VSWCVar = mean(VSWCVariance, na.rm = TRUE)), by = .(YearMon, Site, domain)]
    for (cc in c("VSWCMean", "VSWCMin", "VSWCMax", "VSWCVar"))
      set(agg, which(!is.finite(agg[[cc]])), cc, NA_real_)
    new_rows[[s]] <- agg
    availability[[s]] <- data.table(Site = s, YearMon = months)[agg[, .(Site, YearMon, VSWCMean)],
      on = .(Site, YearMon)][, status := fifelse(is.finite(VSWCMean), "available", "missing_or_nonfinite")]
  } else availability[[s]] <- data.table(Site = s, YearMon = months, status = "no_501_files")
}

manifest_file <- file.path(root, "OUTPUT/NEON_soil_moisture_download_manifest.csv")
new_manifest <- rbindlist(manifest, fill = TRUE)
if (nrow(new_manifest)) {
  prior_manifest <- if (file.exists(manifest_file) && file.info(manifest_file)$size > 0)
    fread(manifest_file) else data.table()
  full_manifest <- unique(rbindlist(list(prior_manifest, new_manifest), fill = TRUE), by = "url")
  fwrite(full_manifest, manifest_file)
}
added <- rbindlist(new_rows, fill = TRUE)
if (nrow(added)) {
  # Preserve all existing values except the missing/nonfinite months being filled.
  combined <- rbindlist(list(added, old), use.names = TRUE, fill = TRUE)
  combined <- unique(combined, by = c("Site", "YearMon"))
  setorder(combined, Site, YearMon)
  stopifnot(!anyDuplicated(combined[, .(Site, YearMon)]))
  archive <- file.path(root, "archive_20261019_soil_moisture")
  dir.create(archive, showWarnings = FALSE)
  backup <- file.path(archive, basename(soil_file))
  if (!file.exists(backup) && !file.copy(soil_file, backup)) stop("Cannot archive original soil data")
  Site_SoilData <- as.data.frame(combined)
  # Existing seasonal aggregates cover a different period and remain unchanged.
  # Only Site_SoilData (the monthly table) is consumed by the active workflow.
  tmp <- tempfile(fileext = ".Rdata")
  save(Site_SoilData, Site_SoilData_Season, file = tmp)
  if (!file.copy(tmp, soil_file, overwrite = TRUE)) stop("Cannot save updated monthly soil data")
  unlink(tmp)
}
coverage <- merge(needed, as.data.table(Site_SoilData)[, .(Site, YearMon, VSWCMean)],
  by = c("Site", "YearMon"), all.x = TRUE)
coverage[, soil_available := is.finite(VSWCMean)]
fwrite(coverage, file.path(root, "OUTPUT/NEON_soil_moisture_flux_coverage.csv"))
message("Added/replaced ", nrow(added), " monthly rows; ", sum(!coverage$soil_available),
  " flux site-months still have no finite soil-moisture measurement. Fluxes are retained.")
