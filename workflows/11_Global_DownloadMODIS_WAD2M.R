# Download/process dynamic land-cover and inundation inputs for spatial CH4 upscaling.
#
# Outputs used by ERA-SpatialUpscaling-Monthly.R:
#   DATA/modis_mcd12c1_processed/MODIS_MCD12C1_ecotype_<year>.tif
#   DATA/wad2m/WAD2M_wetlands_2000-2020_05deg_Ver2.0.nc
#
# Notes:
# - MODIS MCD12C1.061 HDF files are large (~1.2 GB each) and LP DAAC data
#   access requires NASA Earthdata authorization for full downloads.
# - WAD2M is public on Zenodo, but that server can return temporary 504 errors;
#   this script retries the download.

library(tidyverse)
library(jsonlite)
library(terra)

spatial_dir <- Sys.getenv(
  "MONTHLY_UPSCALING_DIR",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly"
)

data_dir <- file.path(spatial_dir, "DATA")
era5_land_dir <- file.path(data_dir, "era5_land_monthly")
modis_hdf_dir <- file.path(data_dir, "modis_mcd12c1_hdf")
modis_processed_dir <- file.path(data_dir, "modis_mcd12c1_processed")
wad2m_dir <- file.path(data_dir, "wad2m")
dir.create(modis_hdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(modis_processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(wad2m_dir, recursive = TRUE, showWarnings = FALSE)

years <- 2001:2024

get_cmr_mcd12c1 <- function(year) {
  cmr_url <- paste0(
    "https://cmr.earthdata.nasa.gov/search/granules.json?",
    "short_name=MCD12C1&version=061&page_size=1&temporal=",
    year, "-01-01T00:00:00Z,", year, "-12-31T23:59:59Z"
  )
  cmr_file <- tempfile(fileext = ".json")
  download.file(cmr_url, cmr_file, quiet = TRUE, mode = "wb")
  cmr <- fromJSON(cmr_file, simplifyVector = FALSE)
  if (length(cmr$feed$entry) == 0) {
    return(tibble(Year = year, title = NA_character_, url = NA_character_))
  }
  entry <- cmr$feed$entry[[1]]
  data_url <- entry$links %>%
    keep(~ !is.null(.x$href) && str_detect(.x$href, "\\.hdf(?:\\?|$)")) %>%
    pluck(1, "href", .default = NA_character_)
  tibble(Year = year, title = entry$title, url = data_url)
}

manifest_file <- file.path(modis_hdf_dir, "MCD12C1_CMR_manifest.csv")
manifest_is_current <- file.exists(manifest_file) && {
  existing_manifest <- read.csv(manifest_file)
  all(years %in% existing_manifest$Year) &&
    all(!is.na(existing_manifest$url[match(years, existing_manifest$Year)]))
}
if (manifest_is_current) {
  message("Using existing MODIS MCD12C1 CMR manifest.")
  modis_manifest <- existing_manifest
} else {
  message("Building or refreshing MODIS MCD12C1 CMR manifest.")
  modis_manifest <- map_dfr(years, get_cmr_mcd12c1)
  write.csv(modis_manifest, manifest_file, row.names = FALSE)
}

is_hdf4_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 4) return(FALSE)
  signature <- readBin(path, what = "raw", n = 4)
  identical(signature, as.raw(c(0x0e, 0x03, 0x13, 0x01)))
}

download_modis_hdf <- identical(tolower(Sys.getenv("DOWNLOAD_MODIS_HDF", unset = "true")), "true")
if (download_modis_hdf) {
  netrc_file <- path.expand(Sys.getenv("EARTHDATA_NETRC", unset = "~/.netrc"))
  if (!file.exists(netrc_file) ||
      !any(str_detect(readLines(netrc_file, warn = FALSE), "machine\\s+urs\\.earthdata\\.nasa\\.gov"))) {
    stop(
      "MODIS downloads require NASA Earthdata credentials. Create ", netrc_file,
      " containing:\n",
      "machine urs.earthdata.nasa.gov login YOUR_USERNAME password YOUR_PASSWORD\n",
      "then restrict its permissions with: chmod 600 ", netrc_file,
      "\nYou may instead set EARTHDATA_NETRC to another credential-file path."
    )
  }
  cookie_file <- file.path(modis_hdf_dir, ".urs_cookies")

  for (i in seq_len(nrow(modis_manifest))) {
    row <- modis_manifest[i, ]
    if (is.na(row$url)) next
    out_file <- file.path(modis_hdf_dir, paste0(row$title, ".hdf"))
    if (is_hdf4_file(out_file)) {
      message("Skipping existing MODIS HDF: ", basename(out_file))
      next
    }
    part_file <- paste0(out_file, ".part")
    if (file.exists(part_file)) unlink(part_file)
    message("Downloading MODIS HDF for ", row$Year, " ...")
    status <- system2(
      "curl",
      c(
        "--fail", "--show-error", "--location",
        "--netrc-file", shQuote(netrc_file),
        "--cookie", shQuote(cookie_file), "--cookie-jar", shQuote(cookie_file),
        "--output", shQuote(part_file), shQuote(row$url)
      )
    )
    if (!identical(status, 0L) || !is_hdf4_file(part_file)) {
      if (file.exists(part_file)) unlink(part_file)
      stop(
        "MODIS download failed for ", row$Year,
        ". Confirm your Earthdata username/password and authorize LP DAAC data access ",
        "in your Earthdata account. curl exit status: ", status
      )
    }
    if (!file.rename(part_file, out_file)) stop("Could not finalize downloaded file: ", out_file)
  }
}

read_modis_igbp <- function(hdf_file) {
  datasets <- sds(hdf_file)
  dataset_descriptions <- vapply(seq_along(datasets), function(i) {
    paste(names(datasets[[i]]), collapse = " ")
  }, character(1))
  selected <- which(str_detect(dataset_descriptions, "^Majority_Land_Cover_Type_1$"))[1]
  if (is.na(selected)) {
    stop("Could not find IGBP majority land-cover layer in ", hdf_file)
  }
  datasets[[selected]]
}

classify_modis_to_ecotype <- function(igbp) {
  # IGBP classes: 0 water; 1-5 forest; 6-7 shrubland; 8-10 savanna/grassland;
  # 11 wetlands; 12 cropland; 13 urban; 14 cropland mosaic; 15 snow/ice; 16 barren.
  #
  # This is the ONLY place non-upland land cover gets excluded from the
  # global upscaling grid (13_Global_SpatialUpscalingRF.R's ecotype_lookup
  # inner_join drops any code not in its lookup table, so mapping a class to
  # NA here is what actually removes it from the grid downstream):
  #   Class  0 (Water)              -> NA (excluded)
  #   Class 11 (Permanent Wetlands) -> NA (excluded)
  #   Class 13 (Urban and Built-up) -> NA (excluded)
  #   Class 15 (Snow/Ice)           -> NA (excluded)
  #   Class 12 (Croplands) and 14 (Cropland/Natural Vegetation Mosaic) ->
  #     code 1 (Cropland); NOT excluded here, but 13_Global_SpatialUpscalingRF.R's
  #     ecotype_lookup filters EcoType != "Cropland" right after the join, so
  #     these are excluded one step downstream instead.
  # Note this only screens each 0.05deg pixel's MAJORITY land cover -- a
  # pixel that's mostly upland but partly flooded (floodplain forest,
  # scattered wetlands) stays in the grid and is instead area-down-weighted
  # by WAD2M's continuous inundation fraction in 13_Global_SpatialUpscalingRF.R,
  # not excluded here or there (see that script's WAD2M loading comment for
  # why a hard cutoff there would wrongly drop legitimate upland cells too).
  #
  # Both barren / desert classes are assigned code 5 (Arid) so that hyper-arid
  # cells receive sink-only treatment instead of inheriting Shrubland source
  # propensity:
  #   Class  7 (Open Shrublands) — covers the Sahara, Arabian Peninsula, and
  #     Central Asian deserts (the geographically dominant desert class).
  #   Class 16 (Barren) — pure rock/sand/tundra with essentially no vegetation.
  # Class 6 (Closed Shrublands) stays as Shrubland (code 4) — those cells have
  # meaningful canopy cover and NEON shrubland observations are applicable.
  subst(
    igbp,
    from = 0:16,
    to = c(NA, 2, 2, 2, 2, 2, 4, 5, 3, 3, 3, NA, 1, NA, 1, NA, 5)
  )
}

era5_template_file <- file.path(era5_land_dir, "era5_land_monthly_2001.nc")
if (file.exists(era5_template_file)) {
  era5_template <- rast(era5_template_file)[[1]]
} else {
  era5_template <- NULL
}

hdf_files <- list.files(modis_hdf_dir, pattern = "^MCD12C1\\..*\\.hdf$", full.names = TRUE)
for (hdf_file in hdf_files) {
  if (!is_hdf4_file(hdf_file)) next
  year <- as.integer(str_match(basename(hdf_file), "A([0-9]{4})001")[, 2])
  out_file <- file.path(modis_processed_dir, sprintf("MODIS_MCD12C1_ecotype_%s.tif", year))
  if (file.exists(out_file)) next
  message("Processing MODIS land cover for ", year)
  ecotype <- classify_modis_to_ecotype(read_modis_igbp(hdf_file))
  names(ecotype) <- "ecotype_code"
  if (!is.null(era5_template)) {
    crs(ecotype) <- crs(era5_template)
    ecotype <- resample(ecotype, era5_template, method = "near")
  }
  writeRaster(ecotype, out_file, overwrite = TRUE, datatype = "INT2S", NAflag = -9999)
}

wad2m_zip <- file.path(data_dir, "WAD2M_wetlands_2000-2020_05deg_Ver2.0.nc.zip")
wad2m_nc <- file.path(wad2m_dir, "WAD2M_wetlands_2000-2020_05deg_Ver2.0.nc")
wad2m_url <- "https://zenodo.org/records/5553187/files/WAD2M_wetlands_2000-2020_05deg_Ver2.0.nc.zip?download=1"

if (!file.exists(wad2m_nc)) {
  if (!file.exists(wad2m_zip) || file.info(wad2m_zip)$size < 1e6) {
    message("Downloading WAD2M monthly inundation/wetland fraction.")
    for (attempt in 1:5) {
      status <- try(
        download.file(wad2m_url, wad2m_zip, mode = "wb", quiet = FALSE, timeout = 180),
        silent = TRUE
      )
      if (!inherits(status, "try-error") && file.exists(wad2m_zip) && file.info(wad2m_zip)$size > 1e6) break
      Sys.sleep(10 * attempt)
    }
  }
  if (file.exists(wad2m_zip) && file.info(wad2m_zip)$size > 1e6) {
    utils::unzip(wad2m_zip, exdir = wad2m_dir)
  } else {
    warning("WAD2M download did not complete. The spatial workflow will use zero inundation fraction until this file is available.")
  }
}

message("Done. MODIS processed files: ", length(list.files(modis_processed_dir, pattern = "\\.tif$")))
message("WAD2M NetCDF present: ", file.exists(wad2m_nc))
