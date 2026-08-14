# NEON site metadata table with MODIS land cover and WAD2M inundation.
#
# Builds a per-site attribute table for the NEON CH4 flux-gradient sites that
# combines:
#   1. NEON/AmeriFlux metadata  -> IGBP class + site description (name + IGBP
#      vegetation description), read from the same "Ameriflux_NEON field-sites.csv"
#      used by 05_NEON_FluxAnalysis.R.
#   2. MODIS MCD12C1 land cover  -> the RAW IGBP majority class at each site
#      (read from the same HDFs as 11_Global_DownloadMODIS_WAD2M.R) AND the
#      project's derived 5-class ecotype (from the processed ecotype rasters).
#   3. WAD2M inundation           -> the range (min / max / mean) of monthly
#      inundated fraction (Fw) at each site's 0.5-deg cell, 2000-2020.
#   4. A flag indicating whether the site is included in the upland CH4
#      source/sink models (12_SourceProp_MagnitudeModels.R / 13_Global_
#      SpatialUpscalingRF.R). Those models are trained/applied on upland
#      ecotypes only; Wetland and Cropland sites are excluded (see
#      classify_modis_to_ecotype() and ecotype_lookup in scripts 11 and 13).
#
# Output: OUTPUT/NEON_site_metadata_landcover_inundation.csv
#
# Column naming follows the repo convention (snake_case; unit after last "_").

library(tidyverse)
library(terra)

# ---- Paths (env-overridable, matching the rest of the workflow) --------------

localdir <- Sys.getenv(
  "LOCALDIR_FLUXGRADIENT",
  unset = "/Volumes/MaloneLab/Research/FluxGradient"
)
localdir.ch4 <- Sys.getenv(
  "LOCALDIR_CH4",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE"
)
spatial_dir <- Sys.getenv(
  "MONTHLY_UPSCALING_DIR",
  unset = file.path(localdir.ch4, "Upscaling_Monthly")
)

metadata_file       <- file.path(localdir, "Ameriflux_NEON field-sites.csv")
modis_hdf_dir       <- file.path(spatial_dir, "DATA", "modis_mcd12c1_hdf")
modis_processed_dir <- file.path(spatial_dir, "DATA", "modis_mcd12c1_processed")
wad2m_nc            <- file.path(spatial_dir, "DATA", "wad2m",
                                 "WAD2M_wetlands_2000-2020_05deg_Ver2.0.nc")
site_filter_dir     <- file.path(localdir.ch4, "NEON_GradientFlux_Data_Filter")

output_dir <- file.path(localdir.ch4, "OUTPUT")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(output_dir, "NEON_site_metadata_landcover_inundation.csv")

if (!file.exists(metadata_file)) {
  stop("NEON/AmeriFlux metadata file not found: ", metadata_file)
}

# ---- Lookup tables -----------------------------------------------------------

# MODIS MCD12C1 IGBP majority land-cover code (0-16) -> abbreviation + name.
igbp_code_lookup <- tibble::tribble(
  ~modis_igbp_code, ~modis_igbp_class, ~modis_igbp_description,
  0L,  "WAT", "Water Bodies",
  1L,  "ENF", "Evergreen Needleleaf Forests",
  2L,  "EBF", "Evergreen Broadleaf Forests",
  3L,  "DNF", "Deciduous Needleleaf Forests",
  4L,  "DBF", "Deciduous Broadleaf Forests",
  5L,  "MF",  "Mixed Forests",
  6L,  "CSH", "Closed Shrublands",
  7L,  "OSH", "Open Shrublands",
  8L,  "WSA", "Woody Savannas",
  9L,  "SAV", "Savannas",
  10L, "GRA", "Grasslands",
  11L, "WET", "Permanent Wetlands",
  12L, "CRO", "Croplands",
  13L, "URB", "Urban and Built-up Lands",
  14L, "CVM", "Cropland/Natural Vegetation Mosaics",
  15L, "SNO", "Permanent Snow and Ice",
  16L, "BSV", "Barren or Sparsely Vegetated"
)

# IGBP abbreviation -> full description (used to fill the site "description"
# from the metadata file's IGBP abbreviation when it lacks a description column).
igbp_abbrev_lookup <- igbp_code_lookup %>%
  transmute(igbp_class = modis_igbp_class, igbp_description = modis_igbp_description) %>%
  distinct(igbp_class, .keep_all = TRUE)

# Project's 5-class ecotype (processed MODIS rasters; script 11 classify_modis_to_ecotype).
ecotype_code_lookup <- tibble::tribble(
  ~modis_ecotype_code, ~modis_ecotype,
  1L, "Cropland",
  2L, "Forest",
  3L, "Grassland",
  4L, "Shrubland",
  5L, "Arid"
)

# ---- Helpers -----------------------------------------------------------------

# Find the first column whose name matches any of the patterns (case-insensitive),
# tolerant to header-spelling variants in the metadata CSV.
pick_col <- function(df, patterns) {
  nm <- names(df)
  for (p in patterns) {
    hit <- which(grepl(p, nm, ignore.case = TRUE))
    if (length(hit) > 0) return(nm[hit[1]])
  }
  NA_character_
}

mode_int <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_integer_)
  as.integer(names(sort(table(x), decreasing = TRUE))[1])
}

# Read the MODIS IGBP majority land-cover layer from an MCD12C1 HDF
# (identical selection to read_modis_igbp() in 11_Global_DownloadMODIS_WAD2M.R).
read_modis_igbp_layer <- function(hdf_file) {
  datasets <- terra::sds(hdf_file)
  desc <- vapply(seq_along(datasets),
                 function(i) paste(names(datasets[[i]]), collapse = " "),
                 character(1))
  sel <- which(stringr::str_detect(desc, "^Majority_Land_Cover_Type_1$"))[1]
  if (is.na(sel)) stop("Majority_Land_Cover_Type_1 not found in ", hdf_file)
  datasets[[sel]]
}

# ---- 1. NEON/AmeriFlux metadata ----------------------------------------------

meta_raw <- read.csv(metadata_file, check.names = TRUE, stringsAsFactors = FALSE)

col_site  <- pick_col(meta_raw, c("Site_Id\\.NEON", "NEON", "^site.?id"))
col_name  <- pick_col(meta_raw, c("^Name$", "Site.?Name", "field_site_name"))
col_igbp  <- pick_col(meta_raw, c("Vegetation.*Abbreviation.*IGBP", "IGBP.*Abbrev", "^IGBP$"))
col_igbpd <- pick_col(meta_raw, c("Vegetation.*Description.*IGBP", "IGBP.*Descr"))
col_lat   <- pick_col(meta_raw, c("Latitude", "^lat$"))
col_lon   <- pick_col(meta_raw, c("Longitude", "^lon$", "^long$"))
col_kop   <- pick_col(meta_raw, c("Climate.*Koeppen", "Koeppen", "Koppen"))

if (any(is.na(c(col_site, col_igbp, col_lat, col_lon)))) {
  stop("Could not locate required metadata columns. Found: ",
       "site=", col_site, " igbp=", col_igbp, " lat=", col_lat, " lon=", col_lon)
}

site_meta <- tibble(
  site_id          = as.character(meta_raw[[col_site]]),
  site_name        = if (!is.na(col_name))  as.character(meta_raw[[col_name]]) else NA_character_,
  igbp_class       = toupper(trimws(as.character(meta_raw[[col_igbp]]))),
  igbp_description = if (!is.na(col_igbpd)) as.character(meta_raw[[col_igbpd]]) else NA_character_,
  koppen_class     = if (!is.na(col_kop))   as.character(meta_raw[[col_kop]]) else NA_character_,
  latitude         = suppressWarnings(as.numeric(meta_raw[[col_lat]])),
  longitude        = suppressWarnings(as.numeric(meta_raw[[col_lon]]))
) %>%
  filter(!is.na(site_id), site_id != "") %>%
  distinct(site_id, .keep_all = TRUE) %>%
  # Fill IGBP description from abbreviation when the file lacks its own.
  left_join(igbp_abbrev_lookup, by = "igbp_class", suffix = c("", "_lu")) %>%
  mutate(
    igbp_description = coalesce(igbp_description, igbp_description_lu),
    site_description = case_when(
      !is.na(site_name) & !is.na(igbp_description) ~ paste0(site_name, " - ", igbp_description),
      !is.na(site_name)                            ~ site_name,
      TRUE                                         ~ igbp_description
    ),
    # Project ecotype from the IGBP abbreviation (same mapping as 05_NEON_FluxAnalysis.R).
    ecotype = case_when(
      igbp_class %in% c("ENF", "DBF", "MF", "EBF", "SAV") ~ "Forest",
      igbp_class == "WET"                                  ~ "Wetland",
      igbp_class == "GRA"                                  ~ "Grassland",
      igbp_class %in% c("CVM", "CRO")                      ~ "Cropland",
      igbp_class == "OSH"                                  ~ "Shrubland",
      TRUE                                                 ~ NA_character_
    )
  ) %>%
  select(-igbp_description_lu)

# Restrict to the NEON gradient-flux study sites (QC'd site folders); fall back
# to all NEON sites in the metadata if that directory is unavailable.
study_sites <- if (dir.exists(site_filter_dir)) {
  list.dirs(site_filter_dir, full.names = FALSE, recursive = FALSE) %>%
    stringr::str_subset("^[A-Z]{4}$")
} else {
  character(0)
}
if (length(study_sites) > 0) {
  site_meta <- filter(site_meta, site_id %in% study_sites)
}

# Point extraction requires valid coordinates; drop (and report) any site
# missing lat/lon so the extraction points stay row-aligned with the table.
missing_coords <- site_meta %>% filter(!is.finite(latitude) | !is.finite(longitude))
if (nrow(missing_coords) > 0) {
  warning("Dropping ", nrow(missing_coords), " site(s) without coordinates: ",
          paste(missing_coords$site_id, collapse = ", "))
}
site_meta <- site_meta %>%
  filter(is.finite(latitude), is.finite(longitude)) %>%
  arrange(site_id)

message("Sites in table: ", nrow(site_meta))

pts <- terra::vect(
  as.data.frame(site_meta[, c("longitude", "latitude")]),
  geom = c("longitude", "latitude"),
  crs  = "EPSG:4326"
)

# ---- 2a. RAW MODIS IGBP majority class at each site --------------------------

hdf_files <- list.files(modis_hdf_dir, pattern = "^MCD12C1\\..*\\.hdf$", full.names = TRUE)
if (length(hdf_files) == 0) {
  warning("No MODIS MCD12C1 HDFs found in ", modis_hdf_dir,
          " -- raw IGBP class will be NA. Run 11_Global_DownloadMODIS_WAD2M.R first.")
  modis_igbp_site <- tibble(site_id = site_meta$site_id, modis_igbp_code = NA_integer_)
} else {
  message("Extracting raw MODIS IGBP class from ", length(hdf_files), " yearly HDF(s)...")
  igbp_by_year <- purrr::map(hdf_files, function(f) {
    r <- read_modis_igbp_layer(f)
    as.integer(terra::extract(r, pts, ID = FALSE)[[1]])
  })
  igbp_mat <- do.call(cbind, igbp_by_year)         # sites x years
  modis_igbp_site <- tibble(
    site_id         = site_meta$site_id,
    modis_igbp_code = apply(igbp_mat, 1, mode_int)  # majority class across years
  )
}
modis_igbp_site <- modis_igbp_site %>%
  left_join(igbp_code_lookup, by = "modis_igbp_code")

# ---- 2b. Processed 5-class ecotype at each site ------------------------------

tif_files <- list.files(modis_processed_dir,
                        pattern = "^MODIS_MCD12C1_ecotype_[0-9]{4}\\.tif$",
                        full.names = TRUE)
if (length(tif_files) == 0) {
  warning("No processed MODIS ecotype rasters in ", modis_processed_dir,
          " -- modis_ecotype will be NA.")
  modis_ecotype_site <- tibble(site_id = site_meta$site_id, modis_ecotype_code = NA_integer_)
} else {
  message("Extracting processed MODIS ecotype from ", length(tif_files), " raster(s)...")
  eco_by_year <- purrr::map(tif_files, function(f) {
    as.integer(terra::extract(terra::rast(f), pts, ID = FALSE)[[1]])
  })
  eco_mat <- do.call(cbind, eco_by_year)
  modis_ecotype_site <- tibble(
    site_id            = site_meta$site_id,
    modis_ecotype_code = apply(eco_mat, 1, mode_int)
  )
}
modis_ecotype_site <- modis_ecotype_site %>%
  left_join(ecotype_code_lookup, by = "modis_ecotype_code")

# ---- 3. WAD2M inundation-fraction range at each site -------------------------

if (file.exists(wad2m_nc)) {
  message("Extracting WAD2M inundation fraction (Fw) time series per site...")
  wad2m <- terra::rast(wad2m_nc)                    # 252 monthly layers, 2000-2020
  fw <- terra::extract(wad2m, pts, ID = FALSE)      # sites x 252
  fw <- as.matrix(fw)
  inund_site <- tibble(
    site_id                  = site_meta$site_id,
    inundation_fraction_min  = apply(fw, 1, function(v) suppressWarnings(min(v, na.rm = TRUE))),
    inundation_fraction_max  = apply(fw, 1, function(v) suppressWarnings(max(v, na.rm = TRUE))),
    inundation_fraction_mean = apply(fw, 1, function(v) mean(v, na.rm = TRUE))
  ) %>%
    mutate(across(starts_with("inundation_fraction_"), ~ ifelse(is.finite(.x), .x, NA_real_)),
           inundation_fraction_range = inundation_fraction_max - inundation_fraction_min)
} else {
  warning("WAD2M NetCDF not found: ", wad2m_nc, " -- inundation columns will be NA.")
  inund_site <- tibble(
    site_id = site_meta$site_id,
    inundation_fraction_min = NA_real_, inundation_fraction_max = NA_real_,
    inundation_fraction_mean = NA_real_, inundation_fraction_range = NA_real_
  )
}

# ---- 4. Assemble table + upland-model inclusion flag -------------------------

# The upland source/sink models (12/13) are trained and applied on upland
# ecotypes only. classify_modis_to_ecotype() (script 11) maps IGBP Water,
# Permanent Wetlands, Urban, and Snow/Ice to NA, and ecotype_lookup (script 13)
# drops Cropland -- so Forest, Grassland, and Shrubland are the included upland
# ecotypes. Arid is a climate-derived subset of Shrubland (only JORN & SRER
# cross the aridity threshold) that is routed through the sink-magnitude model.
upland_ecotypes <- c("Forest", "Grassland", "Shrubland")

site_table <- site_meta %>%
  left_join(modis_igbp_site,    by = "site_id") %>%
  left_join(modis_ecotype_site, by = "site_id") %>%
  left_join(inund_site,         by = "site_id") %>%
  mutate(
    in_upland_source_sink_model = ecotype %in% upland_ecotypes,
    upland_model_note = case_when(
      ecotype == "Wetland"  ~ "Excluded: wetland ecotype (not an upland site)",
      ecotype == "Cropland" ~ "Excluded: cropland ecotype (dropped by ecotype_lookup)",
      site_id %in% c("JORN", "SRER") ~ "Included (Arid subset; routed to sink-magnitude model)",
      ecotype %in% upland_ecotypes   ~ "Included (upland ecotype)",
      TRUE ~ "Excluded: ecotype undetermined"
    )
  ) %>%
  select(
    site_id, site_name, site_description,
    igbp_class, igbp_description, ecotype, koppen_class,
    latitude, longitude,
    modis_igbp_code, modis_igbp_class, modis_igbp_description,
    modis_ecotype_code, modis_ecotype,
    inundation_fraction_min, inundation_fraction_max,
    inundation_fraction_mean, inundation_fraction_range,
    in_upland_source_sink_model, upland_model_note
  ) %>%
  arrange(site_id)

readr::write_csv(site_table, output_file)

message("Wrote ", nrow(site_table), " sites to ", output_file)
message("Included in upland source/sink models: ",
        sum(site_table$in_upland_source_sink_model, na.rm = TRUE), " sites")
print(site_table, n = nrow(site_table), width = Inf)
