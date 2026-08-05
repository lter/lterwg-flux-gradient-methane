# Compare ERA5-gapfilled NEON gradient CH4 fluxes with published reference rates.
#
# Published ecosystem-class values are annual CH4 fluxes from Delwiche et al. (2021),
# Table 1, converted from g C m-2 yr-1 to mg C m-2 d-1.
# NEON values are ERA5-gapfilled annual gradient/tower CH4 budgets from
# OUTPUT/NEON_ERA5_gapfilled_annual_budget_by_year.csv, converted to
# mg C m-2 d-1 for comparison with FLUXNET, process-model, and chamber rates.

library(tidyverse)
library(ggplot2)
library(scales)
library(readxl)
library(patchwork)

localdir.ch4 <- Sys.getenv(
  "LOCALDIR_CH4",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/Methane"
)
localdir <- Sys.getenv(
  "LOCALDIR_FLUXGRADIENT",
  unset = "/Volumes/MaloneLab/Research/FluxGradient"
)
localdir.smud <- Sys.getenv(
  "LOCALDIR_SMUD",
  unset = "/Volumes/MaloneLab/Research/SMUD"
)
if (dir.exists(localdir.ch4)) {
  setwd(localdir.ch4)
}

message("Working directory: ", getwd())
if (!dir.exists("OUTPUT"))  dir.create("OUTPUT",  recursive = TRUE)
if (!dir.exists("FIGURES")) dir.create("FIGURES", recursive = TRUE)
message("OUTPUT writable: ",  file.access("OUTPUT",  mode = 2) == 0)
message("FIGURES writable: ", file.access("FIGURES", mode = 2) == 0)

site_behavior_file <- "OUTPUT/30min_site_behavior.csv"
era5_annual_budget_file <- "OUTPUT/NEON_ERA5_gapfilled_annual_budget_by_year.csv"
era5_mean_annual_budget_file <- "OUTPUT/NEON_ERA5_gapfilled_mean_annual_budget.csv"

if (!file.exists(site_behavior_file)) {
  stop("Missing OUTPUT/30min_site_behavior.csv. Run flow.30min.analysis.R first.")
}

if (!file.exists(era5_annual_budget_file)) {
  stop(
    "Missing ", era5_annual_budget_file, ". Run NEON.ERA5.HalfHourlyGapfill.R first ",
    "to create ERA5-gapfilled annual gradient flux budgets."
  )
}

if (!file.exists(era5_mean_annual_budget_file)) {
  stop(
    "Missing ", era5_mean_annual_budget_file, ". Run NEON.ERA5.HalfHourlyGapfill.R first."
  )
}

behavior_levels <- c("Weak-sink", "Fluctuating", "Weak-source")
behavior_labels <- c(
  "Weak-sink"   = "Weak sink",
  "Fluctuating" = "Fluctuating",
  "Weak-source" = "Weak source"
)
behavior_colors <- c(
  "Weak-sink"   = "#2166AC",
  "Fluctuating" = "grey35",
  "Weak-source" = "#B2182B"
)

base_plot_size <- 11
panel_title_size <- 14
axis_title_size <- 12
axis_text_size <- 10
legend_text_size <- 10
legend_title_size <- 11

# Annual mean fluxes and SD from FLUXNET-CH4 ecosystem classes.
# Units in the source are g C m-2 yr-1.
fluxnet_reference <- tribble(
  ~comparison_group, ~ecosystem_class, ~annual_gC_m2_yr, ~annual_sd_gC_m2_yr, ~n_sites, ~n_site_years, ~reference_id,
  "Published FLUXNET-CH4", "Urban",      46.5,  5.6,  1,  1, "Delwiche et al. 2021, Sect. 3.1.3",
  "Published FLUXNET-CH4", "Marsh",      40.8, 20.7, 10, 42, "Delwiche et al. 2021, Table 1",
  "Published FLUXNET-CH4", "Lake",       28.2, 33.4,  2,  4, "Delwiche et al. 2021, Table 1",
  "Published FLUXNET-CH4", "Swamp",      26.4, 19.9,  6, 15, "Delwiche et al. 2021, Table 1",
  "Published FLUXNET-CH4", "Fen",        20.5, 16.0,  8, 40, "Delwiche et al. 2021, Table 1",
  "Published FLUXNET-CH4", "Rice",       14.4,  8.8,  7, 20, "Delwiche et al. 2021, Table 1",
  "Published FLUXNET-CH4", "Mangrove",   11.1,  0.5,  1,  3, "Delwiche et al. 2021, Table 1",
  "Published FLUXNET-CH4", "Bog",        10.5,  6.4,  7, 32, "Delwiche et al. 2021, Table 1",
  "Published FLUXNET-CH4", "Drained",     6.3,  7.1,  7, 20, "Delwiche et al. 2021, Table 1",
  "Published FLUXNET-CH4", "Upland",      4.0, 10.5, 15, 47, "Delwiche et al. 2021, Table 1",
  "Published FLUXNET-CH4", "Wet tundra",  3.8,  1.8, 11, 39, "Delwiche et al. 2021, Table 1",
  "Published FLUXNET-CH4", "Salt marsh",  2.9,  4.7,  5, 10, "Delwiche et al. 2021, Table 1"
) %>%
  mutate(
    daily_mgC_m2_day = annual_gC_m2_yr * 1000 / 365,
    daily_sd_mgC_m2_day = annual_sd_gC_m2_yr * 1000 / 365,
    daily_low_mgC_m2_day = daily_mgC_m2_day - daily_sd_mgC_m2_day,
    daily_high_mgC_m2_day = daily_mgC_m2_day + daily_sd_mgC_m2_day,
    source_type = "Published ecosystem class"
  )

# Process-based soil CH4 uptake rates from MeMo v1.0.
# Source units are mg CH4 m-2 yr-1. Values are converted to mg C m-2 d-1,
# with negative sign to match the flux convention where uptake is negative.
process_model_uptake <- tribble(
  ~model, ~ecosystem_class, ~uptake_mgCH4_m2_yr, ~uptake_sd_mgCH4_m2_yr, ~reference_id,
  "MeMo v1.0", "Tropical deciduous forest", 602,  63, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Open shrubland",            518, 134, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Temperate broadleaf forest",512,  82, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Savanna",                   500, 132, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Dense shrubland",           481,  90, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Grassland/steppe",          392, 110, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Temperate needleleaf forest",347, 90, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Tropical evergreen forest", 332,  45, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Temperate deciduous forest",321, 70, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Boreal deciduous forest",   282, 117, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Boreal evergreen forest",   269,  94, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Mixed forest",              182,  82, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Tundra",                    176, 143, "Murguia-Flores et al. 2018, Table 9",
  "MeMo v1.0", "Polar desert/rock/ice",     105,  48, "Murguia-Flores et al. 2018, Table 9"
) %>%
  mutate(
    comparison_group = "Process-based soil CH4 model",
    daily_mgC_m2_day = -uptake_mgCH4_m2_yr * 12 / 16 / 365,
    daily_sd_mgC_m2_day = uptake_sd_mgCH4_m2_yr * 12 / 16 / 365,
    daily_low_mgC_m2_day = daily_mgC_m2_day - daily_sd_mgC_m2_day,
    daily_high_mgC_m2_day = daily_mgC_m2_day + daily_sd_mgC_m2_day,
    annual_gC_m2_yr = daily_mgC_m2_day * 365 / 1000,
    annual_sd_gC_m2_yr = daily_sd_mgC_m2_day * 365 / 1000,
    n_sites = NA_integer_,
    n_site_years = NA_integer_,
    source_type = "Process-based model"
  ) %>%
  arrange(daily_mgC_m2_day) %>%
  mutate(plot_rank = rev(seq_len(dplyr::n())))

# Published soil-chamber CH4 flux rates. Values are converted to mg C m-2 d-1
# with positive values indicating net CH4 emission and negative values
# indicating net atmospheric CH4 uptake by soil.
soil_chamber_reference <- tribble(
  ~comparison_group, ~ecosystem_class, ~daily_mgC_m2_day, ~daily_low_mgC_m2_day, ~daily_high_mgC_m2_day, ~n_sites, ~n_site_years, ~reference_id, ~conversion_note,
  "Published soil chamber", "Wetland soil chambers", 75.0, 75.0, 75.0, NA_integer_, NA_integer_,
  "Whalen 2005 review: wetland CH4 emissions commonly about 100 mg CH4 m-2 d-1",
  "100 mg CH4 m-2 d-1 * 12/16 = 75 mg C m-2 d-1",
  "Published soil chamber", "Temperate/subarctic forest soil chambers", -1.50, -2.25, -0.75, 2L, NA_integer_,
  "Adamsen and King 1993: forest soil CH4 uptake generally 1-3 mg CH4 m-2 d-1",
  "-1 to -3 mg CH4 m-2 d-1 * 12/16 = -0.75 to -2.25 mg C m-2 d-1",
  "Published soil chamber", "Rural forest soil chambers", -1.26, -1.26, -1.26, NA_integer_, NA_integer_,
  "Groffman and Pouyat 2009: rural forest CH4 uptake 1.68 mg CH4 m-2 d-1",
  "-1.68 mg CH4 m-2 d-1 * 12/16 = -1.26 mg C m-2 d-1",
  "Published soil chamber", "Desert soil chambers", -0.435, -0.69, -0.18, NA_integer_, NA_integer_,
  "Striegl et al. 1992: central 50% of desert soil CH4 uptake rates 0.24-0.92 mg CH4 m-2 d-1",
  "-0.24 to -0.92 mg CH4 m-2 d-1 * 12/16 = -0.18 to -0.69 mg C m-2 d-1",
  "Published soil chamber", "Grassland soil chambers", -0.548, -0.822, -0.274, NA_integer_, NA_integer_,
  "Mosier et al. 1991: aerobic grassland/soil CH4-C uptake about 1-3 kg CH4-C ha-1 yr-1",
  "-1 to -3 kg CH4-C ha-1 yr-1 = -0.274 to -0.822 mg C m-2 d-1",
  "Published soil chamber", "Cropland soil chambers", -0.1125, -0.3225, 0.1425, NA_integer_, NA_integer_,
  "Danish farmland chamber study: CH4 flux ranged -0.43 to 0.19 mg CH4 m-2 d-1, mean -0.15 mg CH4 m-2 d-1",
  "mg CH4 m-2 d-1 * 12/16 = mg C m-2 d-1"
) %>%
  mutate(
    annual_gC_m2_yr = daily_mgC_m2_day * 365 / 1000,
    annual_sd_gC_m2_yr = NA_real_,
    source_type = "Soil chamber literature"
  )

nmol_ch4_s_to_mg_c_day <- 12.011e-6 * 86400
umol_ch4_s_to_mg_c_day <- 12.011e-3 * 86400
min_non_fluxnet_obs_per_year <- 100

summarise_non_fluxnet_tower_years <- function() {
  validation_dir <- file.path(localdir, "Validation_Sites")
  meta_file <- file.path(localdir, "metadata_validation.csv")

  site_metadata <- if (file.exists(meta_file)) {
    read.csv(meta_file, check.names = FALSE) %>%
      transmute(
        site_id = SITE_ID,
        latitude = suppressWarnings(as.numeric(LATITUDE)),
        longitude = suppressWarnings(as.numeric(LONGITUDE)),
        igbp = IGBP
      )
  } else {
    tibble(site_id = character(), latitude = numeric(), longitude = numeric(), igbp = character())
  }

  tower_metadata <- tribble(
    ~site_id, ~ecosystem_class, ~location, ~reference_id, ~full_reference, ~flux_units_original,
    "SE-Sto", "Subarctic mire", "Stordalen Mire, Abisko, Sweden", "Lakomiec et al. 2021", "Lakomiec, P., Holst, J., Friborg, T., Crill, P., Rakos, N., Kljun, N., Olsson, P.-O., Eklundh, L., Persson, A., and Rinne, J. (2021). Field-scale CH4 emission at a subarctic mire with heterogeneous permafrost thaw status. Biogeosciences, 18, 5811-5830. https://doi.org/10.5194/bg-18-5811-2021", "umol CH4 m-2 s-1",
    "SE-Svb", "Managed boreal forest", "Svartberget, Krycklan Catchment, Sweden", "Chi et al. 2020", "Chi, J., Nilsson, M. B., Laudon, H., Lindroth, A., Wallerman, J., Fransson, J. E. S., Kljun, N., Lundmark, T., Ottosson Lofvenius, M., and Peichl, M. (2020). The Net Landscape Carbon Balance-Integrating terrestrial and aquatic carbon fluxes in a managed boreal forest landscape in Sweden. Global Change Biology, 26, 2353-2367. https://doi.org/10.1111/gcb.14983", "nmol CH4 m-2 s-1",
    "US-Uaf", "Poorly drained black spruce forest over permafrost", "University of Alaska Fairbanks, Alaska, USA", "Iwata et al. 2015", "Iwata, H., Harazono, Y., Ueyama, M., Sakabe, A., Nagano, H., Kosugi, Y., Takahashi, K., and Kim, Y. (2015). Methane exchange in a poorly-drained black spruce forest over permafrost observed using the eddy covariance technique. Agricultural and Forest Meteorology, 214-215, 157-168. https://doi.org/10.1016/j.agrformet.2015.08.252", "nmol CH4 m-2 s-1"
  ) %>%
    left_join(site_metadata, by = "site_id")

  tower_years <- list()

  se_sto_file <- file.path(validation_dir, "SE-Sto", "SE-Sto_gas_fluxes_30min.csv")
  if (file.exists(se_sto_file)) {
    tower_years[["SE-Sto"]] <- read.csv(se_sto_file, check.names = TRUE) %>%
      transmute(
        site_id = "SE-Sto",
        timestamp = as.POSIXct(
          paste(substr(as.character(date), 1, 10), as.character(time)),
          format = "%Y-%m-%d %H:%M:%S",
          tz = "UTC"
        ),
        flux_original = suppressWarnings(as.numeric(Fch4_1_1_1)),
        flux_daily_mgC_m2_day = flux_original * umol_ch4_s_to_mg_c_day
      )
  }

  se_svb_file <- file.path(validation_dir, "SE-Svb", "CH4_SE_SVB_FLUX+PROFILE_2019.csv")
  if (file.exists(se_svb_file)) {
    tower_years[["SE-Svb"]] <- read.csv(se_svb_file, check.names = TRUE) %>%
      transmute(
        site_id = "SE-Svb",
        timestamp = as.POSIXct(timestamp, format = "%d-%b-%Y %H:%M:%S", tz = "UTC"),
        flux_original = suppressWarnings(as.numeric(ch4_flux_nmolm2s_85m)),
        flux_daily_mgC_m2_day = flux_original * nmol_ch4_s_to_mg_c_day
      )
  }

  us_uaf_file <- file.path(validation_dir, "US-Uaf", "US-Uaf CH4_concentration.csv")
  if (file.exists(us_uaf_file)) {
    tower_years[["US-Uaf"]] <- read.csv(us_uaf_file, check.names = TRUE) %>%
      transmute(
        site_id = "US-Uaf",
        timestamp = as.POSIXct(Datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
        flux_original = suppressWarnings(as.numeric(CH4.flux)),
        flux_daily_mgC_m2_day = flux_original * nmol_ch4_s_to_mg_c_day
      )
  }

  if (length(tower_years) == 0) {
    return(tibble())
  }

  bind_rows(tower_years) %>%
    filter(is.finite(flux_daily_mgC_m2_day), !is.na(timestamp)) %>%
    mutate(Year = as.integer(format(timestamp, "%Y"))) %>%
    reframe(
      .by = c(site_id, Year),
      n_observations = dplyr::n(),
      daily_mgC_m2_day = mean(flux_daily_mgC_m2_day, na.rm = TRUE),
      daily_median_mgC_m2_day = median(flux_daily_mgC_m2_day, na.rm = TRUE)
    ) %>%
    filter(n_observations >= min_non_fluxnet_obs_per_year) %>%
    left_join(tower_metadata, by = "site_id")
}

non_fluxnet_tower_years <- summarise_non_fluxnet_tower_years()

non_fluxnet_tower_reference <- non_fluxnet_tower_years %>%
  reframe(
    .by = c(
      site_id, ecosystem_class, location, latitude, longitude, igbp,
      reference_id, full_reference, flux_units_original
    ),
    comparison_group = "Towers outside FLUXNET-CH4",
    CH4_behavior = "Weak-source",
    annual_gC_m2_yr = median(daily_mgC_m2_day, na.rm = TRUE) * 365 / 1000,
    annual_sd_gC_m2_yr = sd(daily_mgC_m2_day, na.rm = TRUE) * 365 / 1000,
    daily_low_mgC_m2_day = quantile(daily_mgC_m2_day, 0.25, na.rm = TRUE),
    daily_high_mgC_m2_day = quantile(daily_mgC_m2_day, 0.75, na.rm = TRUE),
    daily_mgC_m2_day = median(daily_mgC_m2_day, na.rm = TRUE),
    n_sites = 1L,
    n_site_years = dplyr::n(),
    n_observations = sum(n_observations, na.rm = TRUE),
    measurement_years = paste(range(Year, na.rm = TRUE), collapse = "-"),
    source_type = "Tower literature outside FLUXNET-CH4"
  )

additional_draft_tower_reference <- tribble(
  ~comparison_group, ~site_id, ~ecosystem_class, ~location, ~latitude, ~longitude, ~igbp, ~measurement_years, ~daily_mgC_m2_day, ~daily_low_mgC_m2_day, ~daily_high_mgC_m2_day, ~n_sites, ~n_site_years, ~reference_id, ~full_reference, ~flux_units_original, ~conversion_note, ~CH4_behavior,
  "Towers outside FLUXNET-CH4", "US-Blo-literature", "Ponderosa pine plantation", "Blodgett Forest, Sierra Nevada, California, USA", 38.90, -120.63, "ENF", "2007", -2.5 * 12 / 16, -2.5 * 12 / 16, -2.5 * 12 / 16, 1L, NA_integer_, "Smeets et al. 2009", "Smeets, C. J. P. P., Holzinger, R., Vigano, I., Goldstein, A. H., and Rockmann, T. (2009). Eddy covariance methane measurements at a Ponderosa pine plantation in California. Atmospheric Chemistry and Physics, 9, 8365-8375. https://doi.org/10.5194/acp-9-8365-2009", "mg CH4 m-2 d-1", "Daily mean whole-period downward CH4 flux of 2.5 mg CH4 m-2 d-1 converted as -2.5 * 12/16 = mg C m-2 d-1.", "Weak-sink",
  "Towers outside FLUXNET-CH4", "SE-Nor-literature", "Boreal mixed pine-spruce forest", "Norunda research station, central Sweden", 60 + 5 / 60, 17 + 29 / 60, "MF", "2010-2011", mean(c(1.48, 4.57)) * umol_ch4_s_to_mg_c_day / 3600, 1.48 * umol_ch4_s_to_mg_c_day / 3600, 4.57 * umol_ch4_s_to_mg_c_day / 3600, 1L, NA_integer_, "Sundqvist et al. 2015", "Sundqvist, E., Molder, M., Crill, P., Kljun, N., and Lindroth, A. (2015). Methane exchange in a boreal forest estimated by gradient method. Tellus B: Chemical and Physical Meteorology, 67, 26688. https://doi.org/10.3402/tellusb.v67.26688", "umol CH4 m-2 h-1", "Mean tower-method emissions of 1.48-4.57 umol CH4 m-2 h-1 converted to mg C m-2 d-1 with umol CH4 * 12.011e-3 * 24.", "Weak-source",
  "Towers outside FLUXNET-CH4", "CA-Hal-literature", "Temperate mixed forest", "Haliburton Forest and Wildlife Reserve, central Ontario, Canada", 45 + 17 / 60 + 11 / 3600, -(78 + 32 / 60 + 19 / 3600), "MF", "2011", -2.7 * nmol_ch4_s_to_mg_c_day, -(2.7 + 0.13) * nmol_ch4_s_to_mg_c_day, -(2.7 - 0.13) * nmol_ch4_s_to_mg_c_day, 1L, NA_integer_, "Wang et al. 2013", "Wang, J. M., Murphy, J. G., Geddes, J. A., Winsborough, C. L., Basiliko, N., and Thomas, S. C. (2013). Methane fluxes measured by eddy covariance and static chamber techniques at a temperate forest in central Ontario, Canada. Biogeosciences, 10, 4371-4382. https://doi.org/10.5194/bg-10-4371-2013", "nmol CH4 m-2 s-1", "Average EC uptake flux of -2.7 +/- 0.13 nmol CH4 m-2 s-1 converted with nmol CH4 * 12.011e-6 * 86400.", "Weak-sink"
) %>%
  mutate(
    annual_gC_m2_yr = daily_mgC_m2_day * 365 / 1000,
    annual_sd_gC_m2_yr = NA_real_,
    n_observations = NA_integer_,
    source_type = "Tower literature outside FLUXNET-CH4"
  )

additional_draft_chamber_reference <- tribble(
  ~comparison_group, ~site_id, ~ecosystem_class, ~location, ~latitude, ~longitude, ~igbp, ~measurement_years, ~daily_mgC_m2_day, ~daily_low_mgC_m2_day, ~daily_high_mgC_m2_day, ~n_sites, ~n_site_years, ~reference_id, ~full_reference, ~flux_units_original, ~conversion_note,
  "Published soil chamber", "SE-Nor-chambers", "Boreal forest soil chambers", "Norunda research station, central Sweden", 60 + 5 / 60, 17 + 29 / 60, "MF", "2010", -10 * umol_ch4_s_to_mg_c_day / 3600, -10 * umol_ch4_s_to_mg_c_day / 3600, -10 * umol_ch4_s_to_mg_c_day / 3600, 1L, NA_integer_, "Sundqvist et al. 2015", "Sundqvist, E., Molder, M., Crill, P., Kljun, N., and Lindroth, A. (2015). Methane exchange in a boreal forest estimated by gradient method. Tellus B: Chemical and Physical Meteorology, 67, 26688. https://doi.org/10.3402/tellusb.v67.26688", "umol CH4 m-2 h-1", "Soil chamber uptake of around -10 umol CH4 m-2 h-1 converted to mg C m-2 d-1 with umol CH4 * 12.011e-3 * 24."
) %>%
  mutate(
    annual_gC_m2_yr = daily_mgC_m2_day * 365 / 1000,
    annual_sd_gC_m2_yr = NA_real_,
    n_observations = NA_integer_,
    CH4_behavior = NA_character_,
    source_type = "Soil chamber literature"
  )

# === SMUD Global CH4 Compilation ===
# Source: Studies_and_Fluxes_cleaned.xlsx (Malone Lab, /Volumes/MaloneLab/Research/SMUD)
# Annual column is already in mg C m-2 d-1 (daily mean rates).
# Groffman 2009 and Striegl 1992 are excluded because numeric values from those papers
# are already hard-coded in soil_chamber_reference above, preventing double-counting.
smud_duplicate_authors <- c("Groffman", "Striegl")

classify_smud_behavior <- function(flux_daily) {
  case_when(
    flux_daily < -0.05 ~ "Weak-sink",
    flux_daily >  0.05 ~ "Weak-source",
    TRUE               ~ "Fluctuating"
  )
}

# Returns a list-column: one character vector of applicable ecosystem_group
# categories per input row. Most SMUD ecosystem labels map to exactly one
# category; Tundra is duplicated into both Shrubland and Grassland/Savanna
# (same policy as assign_reference_ecosystem_groups() below, and the same
# choice confirmed for the literature reference data), so each individual
# SMUD tundra measurement contributes to both panels' medians rather than
# being forced into a single "best guess" bucket.
map_smud_ecosystem <- function(ecosystem) {
  purrr::map(ecosystem, function(e) {
    cat <- dplyr::case_when(
      e %in% c("Forest", "Rainforest", "Woodland") ~ "Forest",
      e %in% c("Grassland", "Savanna")             ~ "Grassland/Savanna",
      e %in% c("Shrub", "Shrubland")                ~ "Shrubland",
      e == "Agriculture"                             ~ "Cropland",
      e == "Tundra"                                  ~ "MULTI_TUNDRA",
      e %in% c("Desert", "Bare")                     ~ "Arid",
      e == "Wetland"                                  ~ "Wetland/Lake",
      e == "Urban"                                    ~ "Urban",
      TRUE                                             ~ "Other"
    )
    if (identical(cat, "MULTI_TUNDRA")) c("Shrubland", "Grassland/Savanna") else cat
  })
}

smud_file <- file.path(localdir.smud, "Studies_and_Fluxes_cleaned.xlsx")

smud_raw <- if (file.exists(smud_file)) {
  read_xlsx(smud_file, sheet = "Sheet1") %>%
    filter(
      !is.na(Annual),
      is.finite(Annual),
      !str_detect(coalesce(Paper_author, ""), paste(smud_duplicate_authors, collapse = "|"))
    ) %>%
    mutate(
      ecosystem_group = map_smud_ecosystem(Ecosystem),
      CH4_behavior    = classify_smud_behavior(Annual),
      is_ec           = str_detect(coalesce(Flux_method, ""), "Eddy covariance"),
      study_id        = paste0(Paper_author, " ", Paper_year)
    ) %>%
    tidyr::unnest_longer(ecosystem_group) %>%
    filter(ecosystem_group != "Other")
} else {
  warning("SMUD file not found: ", smud_file, ". SMUD data will not be included.")
  tibble()
}

smud_ec_reference <- if (nrow(smud_raw) > 0) {
  smud_raw %>%
    filter(is_ec) %>%
    group_by(ecosystem_group) %>%
    summarise(
      comparison_group      = "SMUD global compilation (EC)",
      ecosystem_class       = paste0("SMUD ", first(ecosystem_group), " EC"),
      n_sites               = n_distinct(study_id),
      n_site_years          = n(),
      daily_mgC_m2_day      = median(Annual, na.rm = TRUE),
      daily_low_mgC_m2_day  = quantile(Annual, 0.25, na.rm = TRUE),
      daily_high_mgC_m2_day = quantile(Annual, 0.75, na.rm = TRUE),
      annual_sd_gC_m2_yr    = sd(Annual, na.rm = TRUE) * 365 / 1000,
      latitude              = median(Latitude, na.rm = TRUE),
      longitude             = median(Longitude, na.rm = TRUE),
      measurement_years     = paste(range(Paper_year, na.rm = TRUE), collapse = "-"),
      .groups               = "drop"
    ) %>%
    mutate(
      CH4_behavior      = classify_smud_behavior(daily_mgC_m2_day),
      annual_gC_m2_yr   = daily_mgC_m2_day * 365 / 1000,
      site_id           = NA_character_,
      location          = NA_character_,
      igbp              = NA_character_,
      n_observations    = NA_integer_,
      flux_units_original = "mg C m-2 d-1 (pre-converted in SMUD compilation)",
      full_reference    = "SMUD global CH4 compilation: Studies_and_Fluxes_cleaned.xlsx, Malone Lab",
      reference_id      = "SMUD global compilation (Studies_and_Fluxes_cleaned.xlsx)",
      source_type       = "SMUD tower (EC)"
    )
} else {
  tibble()
}

smud_chamber_reference <- if (nrow(smud_raw) > 0) {
  smud_raw %>%
    filter(!is_ec) %>%
    group_by(ecosystem_group) %>%
    summarise(
      comparison_group      = "SMUD global compilation (chambers)",
      ecosystem_class       = paste0("SMUD ", first(ecosystem_group), " chambers"),
      n_sites               = n_distinct(study_id),
      n_site_years          = n(),
      daily_mgC_m2_day      = median(Annual, na.rm = TRUE),
      daily_low_mgC_m2_day  = quantile(Annual, 0.25, na.rm = TRUE),
      daily_high_mgC_m2_day = quantile(Annual, 0.75, na.rm = TRUE),
      annual_sd_gC_m2_yr    = sd(Annual, na.rm = TRUE) * 365 / 1000,
      latitude              = median(Latitude, na.rm = TRUE),
      longitude             = median(Longitude, na.rm = TRUE),
      measurement_years     = paste(range(Paper_year, na.rm = TRUE), collapse = "-"),
      .groups               = "drop"
    ) %>%
    mutate(
      CH4_behavior      = classify_smud_behavior(daily_mgC_m2_day),
      annual_gC_m2_yr   = daily_mgC_m2_day * 365 / 1000,
      site_id           = NA_character_,
      location          = NA_character_,
      igbp              = NA_character_,
      n_observations    = NA_integer_,
      flux_units_original = "mg C m-2 d-1 (pre-converted in SMUD compilation)",
      full_reference    = "SMUD global CH4 compilation: Studies_and_Fluxes_cleaned.xlsx, Malone Lab",
      reference_id      = "SMUD global compilation (Studies_and_Fluxes_cleaned.xlsx)",
      source_type       = "SMUD soil chamber"
    )
} else {
  tibble()
}

if (nrow(smud_raw) > 0) {
  message(sprintf(
    "SMUD loaded: %d EC rows in %d ecosystem groups; %d chamber rows in %d ecosystem groups.",
    nrow(smud_raw %>% filter(is_ec)),    nrow(smud_ec_reference),
    nrow(smud_raw %>% filter(!is_ec)),   nrow(smud_chamber_reference)
  ))
}

draft_candidate_sources <- tribble(
  ~source_type, ~reference_id, ~ecosystem_class, ~location, ~measurement_type, ~status, ~doi_or_url, ~notes,
  "Tower/EC/gradient candidate", "Shoemaker et al. 2014", "Temperate forest", "Harvard Forest, Massachusetts, USA", "Tower/ecosystem-scale CH4", "Identified in draft; numeric flux still needs full-text/table extraction", "https://doi.org/10.1002/2013GL058691", "Publisher PDF was blocked by browser challenge during this pass.",
  "Tower/EC/gradient candidate", "Hill and Vargas 2022", "Temperate tidal salt marsh", "Delaware, USA", "Plot and ecosystem CH4", "Identified in draft; numeric flux still needs full-text/table extraction", "https://doi.org/10.1029/2022JG006943", "Useful wetland/salt-marsh benchmark if not already represented by FLUXNET-CH4.",
  "Tower/EC/gradient candidate", "Werner et al. 2003", "Mixed temperate/boreal lowland and wetland forest", "Northern Wisconsin, USA", "Tall-tower CH4", "Identified in draft; numeric flux still needs full-text/table extraction", "https://doi.org/10.1046/j.1365-2486.2003.00670.x", "Regional/tall-tower footprint rather than site-scale EC.",
  "Tower/EC/gradient candidate", "Schrier-Uijl et al. 2010", "Peat grassland", "Netherlands", "EC and chambers", "Identified in draft; numeric flux still needs full-text/table extraction", "https://doi.org/10.1016/j.agrformet.2010.05.005", "Potential bridge between chamber and EC estimates.",
  "Tower/EC/gradient candidate", "Zhang et al. 2012", "Permafrost ecosystem", "Qinghai-Tibetan Plateau, China", "Chambers scaled to EC", "Identified in draft; numeric flux still needs full-text/table extraction", "https://doi.org/10.1111/j.1365-2486.2011.02587.x", "Potential tundra/permafrost benchmark.",
  "Tower/EC/gradient candidate", "Yu et al. 2013", "Alpine wetland", "Tibetan Plateau, China", "EC and manual/automated chambers", "Identified in draft; numeric flux still needs full-text/table extraction", "https://doi.org/10.1016/j.envpol.2013.06.018", "Potential wetland benchmark.",
  "Tower/EC/gradient candidate", "Zhao et al. 2019", "Small ponds", "Boreal/temperate ponds", "Flux-gradient and EC methods", "Identified in draft; numeric flux still needs full-text/table extraction", "https://doi.org/10.1016/j.agrformet.2019.05.032", "Aquatic source; may belong in a pond/lake group rather than tower forest comparison.",
  "Chamber/soil candidate", "Keller et al. 1983", "Forest soils", NA_character_, "Soil chambers", "Identified in draft; numeric flux still needs full-text/table extraction", NA_character_, "Older chamber benchmark.",
  "Chamber/soil candidate", "Keller et al. 1990", "Tropical/agricultural development soils", "Central Panama", "Soil chambers", "Identified in draft; numeric flux still needs full-text/table extraction", NA_character_, "Land-use contrast benchmark.",
  "Chamber/soil candidate", "Steudler et al. 1989", "Temperate forest soils", NA_character_, "Soil chambers", "Identified in draft; numeric flux still needs full-text/table extraction", NA_character_, "N-fertilization/forest soil uptake benchmark.",
  "Chamber/soil candidate", "Yavitt et al. 1990", "Temperate forest soils", NA_character_, "Soil chambers", "Identified in draft; numeric flux still needs full-text/table extraction", NA_character_, "Site-level forest soil uptake/source contrast.",
  "Chamber/soil candidate", "Yavitt et al. 1995", "Northern hardwood forest", NA_character_, "Soil chambers", "Identified in draft; numeric flux still needs full-text/table extraction", NA_character_, "Forest ecosystem methane dynamics.",
  "Chamber/soil candidate", "Crill 1991", "Temperate woodland soil", NA_character_, "Soil chambers", "Identified in draft; numeric flux still needs full-text/table extraction", NA_character_, "Soil CH4 uptake benchmark.",
  "Chamber/soil candidate", "Whalen et al. 1991", "Taiga soils", NA_character_, "Soil chambers", "Identified in draft; numeric flux still needs full-text/table extraction", NA_character_, "Boreal/taiga soil uptake benchmark.",
  "Chamber/soil candidate", "Seiler et al. 1984", "Tropical soils and termite nests", "Tropical regions", "Soil chambers", "Identified in draft; numeric flux still needs full-text/table extraction", NA_character_, "May be too broad for site-level comparison.",
  "Chamber/soil candidate", "Scharffe et al. 1990", "Tropical soils/savanna", "Guayana Shield, Venezuela", "Soil chambers", "Identified in draft; numeric flux still needs full-text/table extraction", NA_character_, "Tropical benchmark.",
  "Chamber/soil candidate", "Teh et al. 2005", "Humid tropical forest soils", NA_character_, "Soil chambers", "Identified in draft; numeric flux still needs full-text/table extraction", NA_character_, "Tropical soil process benchmark.",
  "Chamber/soil candidate", "von Fischer and Hedin 2002", "Forest soils", NA_character_, "Isotope pool dilution/chambers", "Identified in draft; numeric flux still needs full-text/table extraction", "https://doi.org/10.1029/2001GB001448", "Production and consumption rates, not only net flux.",
  "Chamber/soil candidate", "Angle et al. 2017", "Wetland oxygenated soils", NA_character_, "Soil incubations/chambers", "Identified in draft; numeric flux still needs full-text/table extraction", "https://doi.org/10.1038/ncomms15617", "Process-oriented wetland source.",
  "Chamber/soil synthesis candidate", "Treat et al. 2014", "Permafrost peatlands", "Pan-Arctic", "Synthesis", "Identified in draft; numeric flux still needs synthesis extraction", NA_character_, "Synthesis may provide broad wetland/permafrost benchmark.",
  "Chamber/soil synthesis candidate", "Treat et al. 2018", "Permafrost-region ecosystems", "Pan-Arctic", "Synthesis", "Identified in draft; numeric flux still needs synthesis extraction", NA_character_, "Nongrowing-season benchmark.",
  "Chamber/soil synthesis candidate", "Turetsky et al. 2014", "Wetlands", "Global", "Synthesis", "Identified in draft; numeric flux still needs synthesis extraction", "https://doi.org/10.1111/gcb.12580", "Large wetland synthesis; likely best for broad wetland comparison rather than site count.",
  "Chamber/soil synthesis candidate", "Smith et al. 2000", "Northern European soils", "Northern Europe", "Comparison/review", "Identified in draft; numeric flux still needs synthesis extraction", NA_character_, "Useful regional soil uptake/emission benchmark."
)

draft_candidate_extracted_info <- tribble(
  ~reference_id, ~extraction_status, ~plot_ready, ~ecosystem_class, ~location, ~measurement_type, ~reported_metric, ~reported_units, ~daily_mgC_m2_day, ~daily_low_mgC_m2_day, ~daily_high_mgC_m2_day, ~conversion_note, ~full_reference, ~notes,
  "Shoemaker et al. 2014", "Blocked/full text needed", FALSE, "Temperate evergreen forest", "Howland Forest, Maine, USA", "Ecosystem-scale tower flux", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Shoemaker, J. K., Keenan, T. F., Hollinger, D. Y., and Richardson, A. D. (2014). Forest ecosystem changes from annual methane source to sink depending on late summer water balance. Geophysical Research Letters, 41, 673-679. https://doi.org/10.1002/2013GL058691", "Accessible abstract reports neutral-to-source behavior in 2011 and small sink behavior in 2012, but not numeric rates.",
  "Hill and Vargas 2022", "Relative comparison extracted; absolute flux still needs full text/table", FALSE, "Temperate tidal salt marsh", "Delaware, USA", "EC and chamber comparison", "Chamber-upscaled estimates underestimated CH4 emissions by 69% relative to eddy covariance", "percent difference", NA_real_, NA_real_, NA_real_, NA_character_, "Hill, A. C., and Vargas, R. (2022). Methane and carbon dioxide fluxes in a temperate tidal salt marsh: Comparisons between plot and ecosystem measurements. Journal of Geophysical Research: Biogeosciences, 127. https://doi.org/10.1029/2022JG006943", "Useful scale-comparison source, but accessible abstract does not provide absolute CH4 flux rate.",
  "Werner et al. 2003", "Blocked/full text needed", FALSE, "Mixed temperate/boreal lowland and wetland forest", "Northern Wisconsin, USA", "Tall-tower CH4 exchange", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Werner, C., Davis, K., Bakwin, P., Yi, C., Hurst, D., and Lock, L. (2003). Regional-scale measurements of CH4 exchange from a tall tower over a mixed temperate/boreal lowland and wetland forest. Global Change Biology, 9, 1251-1261. https://doi.org/10.1046/j.1365-2486.2003.00670.x", "Likely relevant for regional tower benchmark; numeric rate needs full text/table.",
  "Schrier-Uijl et al. 2010", "Blocked/full text needed", FALSE, "Heterogeneous grass ecosystem on peat", "Netherlands", "EC and chamber comparison", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Schrier-Uijl, A. P., Kroon, P. S., Hensen, A., Leffelaar, P. A., Berendse, F., and Veenendaal, E. M. (2010). Comparison of chamber and eddy covariance-based CO2 and CH4 emission estimates in a heterogeneous grass ecosystem on peat. Agricultural and Forest Meteorology, 150, 825-831. https://doi.org/10.1016/j.agrformet.2010.05.005", "Relevant scale-comparison source; numeric rate needs full text/table.",
  "Zhang et al. 2012", "Blocked/full text needed", FALSE, "Permafrost ecosystem", "Siberian permafrost region", "Chamber upscaling to EC/model", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Zhang, Y., Sachs, T., Li, C., and Boike, J. (2012). Upscaling methane fluxes from closed chambers to eddy covariance based on a permafrost biogeochemistry integrated model. Global Change Biology, 18, 1428-1440. https://doi.org/10.1111/j.1365-2486.2011.02587.x", "Relevant chamber-to-EC scaling source; numeric comparison value needs full text/table.",
  "Yu et al. 2013", "Blocked/full text needed", FALSE, "Alpine wetland", "Tibetan Plateau, China", "EC and chamber comparison", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Yu, L., Wang, H., Wang, G., Song, W., Huang, Y., Li, S.-G., Liang, N., Tang, Y., and He, J.-S. (2013). A comparison of methane emission measurements using eddy covariance and manual and automated chamber-based techniques in Tibetan Plateau alpine wetland. Environmental Pollution, 181, 81-90. https://doi.org/10.1016/j.envpol.2013.06.018", "Relevant wetland scale-comparison source; numeric rate needs full text/table.",
  "Zhao et al. 2019", "Blocked/full text needed", FALSE, "Small ponds", NA_character_, "Flux-gradient and EC method evaluation", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Zhao, J., Zhang, M., Xiao, W., Wang, W., Zhang, Z., Yu, Z., Xiao, Q., Cao, Z., Xu, J., Zhang, X., Liu, S., and Lee, X. (2019). An evaluation of the flux-gradient and the eddy covariance method to measure CH4, CO2, and H2O fluxes from small ponds. Agricultural and Forest Meteorology, 275, 255-264. https://doi.org/10.1016/j.agrformet.2019.05.032", "Methodologically relevant to gradient fluxes, but aquatic/pond source and numeric rate needs full text/table.",
  "Keller et al. 1983", "Blocked/full text needed", FALSE, "Forest soils", NA_character_, "Soil chambers", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Keller, M., Goreau, T. J., Wofsy, S. C., Kaplan, W. A., and McElroy, M. B. (1983). Production of nitrous oxide and consumption of methane by forest soils. Geophysical Research Letters, 10, 1156-1159. https://doi.org/10.1029/GL010i012p01156", "Found as historical chamber/source reference; numeric rate needs full text/table.",
  "Keller et al. 1990", "Blocked/full text needed", FALSE, "Central Panama soils", "Central Panama", "Soil chambers", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Keller, M., Mitre, M. E., and Stallard, R. F. (1990). Consumption of atmospheric methane in soils of central Panama: Effects of agricultural development. Global Biogeochemical Cycles, 4, 21-27. https://doi.org/10.1029/GB004i001p00021", "Found as historical chamber/source reference; numeric rate needs full text/table.",
  "Steudler et al. 1989", "Blocked/full text needed", FALSE, "Temperate forest soils", NA_character_, "Soil chambers", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Steudler, P. A., Bowden, R. D., Melillo, J. M., and Aber, J. D. (1989). Influence of nitrogen fertilization on methane uptake in temperate forest soils. Nature, 341, 314-316. https://doi.org/10.1038/341314a0", "Found as historical chamber/source reference; numeric rate needs full text/table.",
  "Yavitt et al. 1990", "Blocked/full text needed", FALSE, "Temperate forest soils", NA_character_, "Soil chambers", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Yavitt, J. B., Downey, D. M., Lang, G. E., and Sexstone, A. J. (1990). Methane consumption in two temperate forest soils. Biogeochemistry, 9, 39-52. https://doi.org/10.1007/BF00002716", "Found as historical chamber/source reference; numeric rate needs full text/table.",
  "Yavitt et al. 1995", "Blocked/full text needed", FALSE, "Northern hardwood ecosystem", NA_character_, "Soil chambers", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Yavitt, J. B., Fahey, T. J., and Simmons, J. A. (1995). Methane and carbon dioxide dynamics in a northern hardwood ecosystem. Soil Science Society of America Journal, 59, 796-804. https://doi.org/10.2136/sssaj1995.03615995005900030023x", "Found as historical chamber/source reference; numeric rate needs full text/table.",
  "Crill 1991", "Blocked/full text needed", FALSE, "Temperate woodland soil", NA_character_, "Soil chambers", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Crill, P. M. (1991). Seasonal patterns of methane uptake and carbon dioxide release by a temperate woodland soil. Global Biogeochemical Cycles, 5, 319-334. https://doi.org/10.1029/91GB02466", "Found as historical chamber/source reference; numeric rate needs full text/table.",
  "Whalen et al. 1991", "Blocked/full text needed", FALSE, "Taiga", NA_character_, "Soil chambers", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Whalen, S. C., Reeburgh, W. S., and Kizer, K. S. (1991). Methane consumption and emission by taiga. Global Biogeochemical Cycles, 5, 261-273. https://doi.org/10.1029/91GB01303", "Found as historical chamber/source reference; numeric rate needs full text/table.",
  "Seiler et al. 1984", "Blocked/full text needed", FALSE, "Tropical soils and termite nests", NA_character_, "Field chambers", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Seiler, W., Conrad, R., and Scharffe, D. (1984). Field studies of methane emission from termite nests into the atmosphere and measurements of methane uptake by tropical soils. Journal of Atmospheric Chemistry, 1, 171-186. https://doi.org/10.1007/BF00053839", "Found as historical chamber/source reference; numeric rate needs full text/table.",
  "Scharffe et al. 1990", "Blocked/full text needed", FALSE, "Northern Guayana Shield soils", "Venezuela", "Soil flux chambers", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Scharffe, D., Hao, W. M., Donoso, L., Crutzen, P. J., and Sanhueza, E. (1990). Soil fluxes and atmospheric concentration of CO and CH4 in the northern part of the Guayana Shield, Venezuela. Journal of Geophysical Research, 95, 22475. https://doi.org/10.1029/JD095iD13p22475", "Found as historical chamber/source reference; numeric rate needs full text/table.",
  "Teh et al. 2005", "Blocked/full text needed", FALSE, "Humid tropical forest soils", NA_character_, "Soil incubations/process rates", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Teh, Y. A., Silver, W. L., and Conrad, M. E. (2005). Oxygen effects on methane production and oxidation in humid tropical forest soils. Global Change Biology, 11, 1283-1297. https://doi.org/10.1111/j.1365-2486.2005.00983.x", "Mechanistic production/oxidation source; numeric net flux needs full text/table.",
  "von Fischer and Hedin 2002", "Gross rates extracted; not directly net-flux comparable", FALSE, "Diverse soils", "17 field sites", "Field-based isotope pool dilution", "Gross production 0.04-930; gross consumption 0.1-9.2; mean production in dry oxic soils 0.15", "mg CH4-C m-2 d-1", NA_real_, NA_real_, NA_real_, "Already in CH4-C daily units; gross production/consumption rates are not plotted as net fluxes.", "von Fischer, J. C., and Hedin, L. O. (2002). Separating methane production and consumption with a field-based isotope pool dilution technique. Global Biogeochemical Cycles, 16, 8-1-8-13. https://doi.org/10.1029/2001GB001448", "Important evidence for co-occurring production and consumption, but not a net ecosystem/chamber flux point.",
  "Angle et al. 2017", "Mechanistic contribution extracted; not an absolute flux", FALSE, "Freshwater wetland", "Old Woman Creek, Ohio, USA", "Soil geochemistry/metagenomics/process attribution", "Up to 80% of methane fluxes attributed to methanogenesis in oxygenated soils", "percent contribution", NA_real_, NA_real_, NA_real_, NA_character_, "Angle, J. C., Morin, T. H., Solden, L. M., Narrowe, A. B., Smith, G. J., Borton, M. A., et al. (2017). Methanogenesis in oxygenated soils is a substantial fraction of wetland methane emissions. Nature Communications, 8, 1567. https://doi.org/10.1038/s41467-017-01753-4", "Supports mechanism behind scale mismatch; does not provide a standalone net flux rate for the comparison figure.",
  "Treat et al. 2014", "Process-production source; not directly net-flux comparable", FALSE, "Alaskan permafrost peats", "Alaska, USA", "Laboratory peat production", NA_character_, NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Treat, C. C., Wollheim, W. M., Varner, R. K., Grandy, A. S., Talbot, J., and Frolking, S. (2014). Temperature and peat type control CO2 and CH4 production in Alaskan permafrost peats. Global Change Biology, 20, 2674-2686. https://doi.org/10.1111/gcb.12572", "Useful for mechanisms; requires full text to extract production rates and is not a net field flux.",
  "Treat et al. 2018", "Annual synthesis ranges extracted", FALSE, "Northern wetlands and uplands", "Temperate, boreal, and tundra regions", "Synthesis of nongrowing-season and annual CH4 fluxes", "Annual wetland emissions ranged 0.9-78; upland median annual flux 0.0 +/- 0.2", "g CH4 m-2 yr-1", NA_real_, 0.9 * 12 / 16 * 1000 / 365, 78 * 12 / 16 * 1000 / 365, "Annual g CH4 m-2 yr-1 converted to mg C m-2 d-1 as g CH4 * 12/16 * 1000/365.", "Treat, C. C., Bloom, A. A., and Marushchak, M. E. (2018). Nongrowing season methane emissions-a significant component of annual emissions across northern ecosystems. Global Change Biology, 24, 3331-3343. https://doi.org/10.1111/gcb.14137", "Useful broad synthesis benchmark, but not plotted because the accessible value is a multi-ecosystem annual range rather than a site/source-class daily estimate.",
  "Turetsky et al. 2014", "Controls and dataset size extracted; numeric summary flux still needs full text/table", FALSE, "Northern, temperate, and subtropical wetlands", "71 wetland sites", "Synthesis of instantaneous wetland CH4 fluxes", "Approximately 19,000 instantaneous measurements from 71 wetland sites", "count", NA_real_, NA_real_, NA_real_, NA_character_, "Turetsky, M. R., Kotowska, A., Bubier, J., Dise, N. B., Crill, P., Hornibrook, E. R. C., et al. (2014). A synthesis of methane emissions from 71 northern, temperate, and subtropical wetlands. Global Change Biology, 20, 2183-2197. https://doi.org/10.1111/gcb.12580", "Accessible abstract did not report a numeric flux summary; full table/database extraction needed.",
  "Smith et al. 2000", "Source identified; numeric rate still needs full text/table", FALSE, "Northern European soils", "Northern Europe", "Soil methane oxidation synthesis/comparison", "Range and statistical distribution of oxidation rates reported, but numeric summary not exposed in accessible abstract", NA_character_, NA_real_, NA_real_, NA_real_, NA_character_, "Smith, K. A., Dobbie, K. E., Ball, B. C., Bakken, L. R., Sitaula, B. K., Hansen, S., et al. (2000). Oxidation of atmospheric methane in Northern European soils, comparison with other ecosystems, and uncertainties in the global terrestrial sink. Global Change Biology, 6, 791-803. https://doi.org/10.1046/j.1365-2486.2000.00356.x", "The earlier related Ball et al. 1997 paper reports CH4 oxidation rates from 0 to 2.5 mg m-2 d-1, but that value was not assigned to Smith et al. 2000 without full-text confirmation."
)

# Returns a list-column: one character vector of applicable ecosystem_group
# categories per input ecosystem_class string. Most reference classes map to
# exactly one category; a few composite/ambiguous classes are duplicated
# across every category their site composition actually spans, so a single
# published/measured flux estimate contributes to each relevant panel rather
# than being forced into one "best guess" bucket:
#  - "Upland" (Delwiche et al. 2021 FLUXNET-CH4 aggregate, n=15 sites: 6
#    needleleaf + 1 mixed forest, 2 alpine meadow + 1 grassland + 1 tundra,
#    3 cropland, 1 urban) -> Forest, Grassland/Savanna, Cropland, Urban.
#  - "Drained" (Delwiche et al. 2021 aggregate, n=7: former wetlands now used
#    as grassland n=3 or cropland n=3) -> Grassland/Savanna, Cropland.
#  - Any "tundra" class -> Shrubland, Grassland/Savanna.
assign_reference_ecosystem_groups <- function(ecosystem_class) {
  purrr::map(ecosystem_class, function(ec) {
    if (is.na(ec)) return(NA_character_)
    if (str_detect(ec, regex("^upland$", ignore_case = TRUE))) {
      return(c("Forest", "Grassland/Savanna", "Cropland", "Urban"))
    }
    if (str_detect(ec, regex("^drained$", ignore_case = TRUE))) {
      return(c("Grassland/Savanna", "Cropland"))
    }
    if (str_detect(ec, regex("tundra", ignore_case = TRUE))) {
      return(c("Shrubland", "Grassland/Savanna"))
    }
    if (str_detect(ec, regex("marsh|fen|swamp|bog|mire|wetland|mangrove|lake", ignore_case = TRUE))) {
      return("Wetland/Lake")
    }
    if (str_detect(ec, regex("rice|cropland", ignore_case = TRUE))) return("Cropland")
    if (str_detect(ec, regex("forest|black spruce", ignore_case = TRUE))) return("Forest")
    if (str_detect(ec, regex("desert|rock|ice", ignore_case = TRUE))) return("Arid")
    if (str_detect(ec, regex("shrubland", ignore_case = TRUE))) return("Shrubland")
    if (str_detect(ec, regex("grassland|steppe|savanna", ignore_case = TRUE))) return("Grassland/Savanna")
    if (str_detect(ec, regex("urban", ignore_case = TRUE))) return("Urban")
    "Other"
  })
}

# Condensed to align with the categories used in the RF upscaling model
# (12_SourceProp_MagnitudeModels.R / 13_Global_SpatialUpscalingRF.R: Forest,
# Grassland, Shrubland, plus "Arid" as an aridity_index-based split of
# Shrubland — see assign_neon_ecosystem() below). Wetland and Lake are
# merged (both are inundated/aquatic CH4 sources); Urban, Wetland/Lake, and
# Cropland are kept as additional comparison-only categories not used in
# the upscaling itself.
ecosystem_levels <- c(
  "Urban", "Wetland/Lake", "Cropland", "Forest", "Shrubland", "Grassland/Savanna", "Arid"
)

# "All ecosystems" is a summary facet appended after the 7 real categories --
# with facet_wrap(ncol = 3) that's 8 panels in a 3x3 grid, using the one
# empty grid slot left over. It is NOT used for per-record ecosystem
# assignment (assign_reference_ecosystem_groups(), assign_neon_ecosystem(),
# upland_eco_groups all ignore it) -- it's added only when building
# bar_source_data for comparison_figure, as a deduplicated pooled copy of
# every other panel's data. Keeping it as the last level of ecosystem_levels
# means it's automatically included wherever comparison_figure factors
# ecosystem_group by that vector, without affecting any other dataset.
# "All Upland" is a second summary facet, pooling just the non-inundated
# categories (everything except Wetland/Lake) -- same idea as the density
# plot's "upland" pooling below, but broken out by data source instead of by
# tower vs. chamber. With this and "All ecosystems" added, the 7 real
# categories + 2 summaries = 9 panels, filling the 3x3 grid exactly.
ecosystem_levels <- c(ecosystem_levels, "All ecosystems", "All Upland")

# Upland (non-inundated) ecosystem groups — excludes Wetland/Lake. Used both
# for the "All Upland" summary facet in comparison_figure and for the
# tower-vs-chamber density plot below.
upland_eco_groups <- c("Urban", "Cropland", "Forest", "Shrubland",
                       "Grassland/Savanna", "Arid")

# Checks whether ANY of a reference class's applicable categories (see
# assign_reference_ecosystem_groups()) falls in upland_eco_groups -- used for
# pooled "is this upland at all" filters that don't need a per-category
# facet split (the density plot, and the "All Upland" summary facet below).
is_any_upland <- function(ecosystem_class) {
  purrr::map_lgl(
    assign_reference_ecosystem_groups(ecosystem_class),
    function(cats) any(cats %in% upland_eco_groups, na.rm = TRUE)
  )
}

# Same aridity_index definition used throughout the RF pipeline
# (12_SourceProp_MagnitudeModels.R, 13_Global_SpatialUpscalingRF.R,
# 19_Supp_NEONRepresentativeness.R): aridity_index = MAP / (MAT + 10), with
# a floor on MAT to avoid the formula blowing up as MAT+10 -> 0. No NEON
# site is literally labeled "Arid" in site metadata (see the note by
# assign_neon_ecosystem() below) -- this is what identifies the ~2 sites
# (JORN, SRER; both raw EcoType "Shrubland") that should be shown under the
# Arid comparison panel instead.
arid_ai_threshold <- 15
aridity_mat_floor <- -9  # deg C; see 12_SourceProp_MagnitudeModels.R for rationale

site_behavior <- read.csv(site_behavior_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    CH4_behavior = factor(CH4_behavior, levels = behavior_levels),
    CH4_gradient_behavior = factor(CH4_gradient_behavior, levels = behavior_levels),
    MAP = as.numeric(MAP),
    MAT = as.numeric(MAT),
    aridity_index = if_else(MAT > aridity_mat_floor, MAP / (MAT + 10), NA_real_),
    is_arid = as.integer(!is.na(aridity_index) & aridity_index < arid_ai_threshold)
  )

era5_annual_budget <- read.csv(era5_annual_budget_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    Year = as.integer(Year),
    annual_budget_gC_m2_yr = as.numeric(annual_budget_gC_m2_yr),
    daily_mgC_m2_day = annual_budget_gC_m2_yr * 1000 / 365,
    model_only_daily_mgC_m2_day = model_only_annual_budget_gC_m2_yr * 1000 / 365
  ) %>%
  left_join(
    site_behavior %>% dplyr::select(SITE_ID, CH4_behavior, CH4_gradient_behavior, EcoType, is_arid),
    by = "SITE_ID"
  ) %>%
  filter(
    is.finite(annual_budget_gC_m2_yr),
    is.finite(daily_mgC_m2_day)
  )

era5_mean_annual_budget <- read.csv(era5_mean_annual_budget_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    era5_annual_behavior = factor(era5_annual_behavior, levels = behavior_levels)
  )

# NEON's own EcoType field only ever takes: Cropland, Forest, Grassland,
# Shrubland, Wetland (confirmed directly from OUTPUT/30min_site_behavior.csv)
# -- no NEON site is literally "Tundra," "Desert," "Lake," or "Urban," so
# this stays a simple scalar mapping (no multi-category duplication needed
# here, unlike assign_reference_ecosystem_groups() for the literature/model
# data). Desert/Bare and Wetland/Lake branches are kept for robustness in
# case EcoType ever gains those values.
#
# is_arid (see site_behavior above) takes priority over the raw text label:
# JORN and SRER are both raw EcoType "Shrubland" but cross the
# aridity_index < arid_ai_threshold definition used throughout the RF
# pipeline, so they're shown under the Arid panel here instead of
# Shrubland -- the same reclassification 12_SourceProp_MagnitudeModels.R
# applies for training. Empirically both sites are 100% observed
# weak-source (see 12's header note), which is the whole reason this
# comparison is worth showing: it lets the Arid panel visibly disagree with
# the "arid soils are net CH4 sinks" prior instead of that disagreement
# being invisible inside the Shrubland panel's pooled statistics.
assign_neon_ecosystem <- function(ecotype, is_arid = 0L) {
  case_when(
    is_arid == 1 ~ "Arid",
    ecotype == "Grassland" ~ "Grassland/Savanna",
    ecotype == "Wetland" ~ "Wetland/Lake",
    ecotype %in% c("Desert", "Bare") ~ "Arid",
    ecotype %in% ecosystem_levels ~ ecotype,
    TRUE ~ "Other"
  )
}

coverage_summary <- era5_annual_budget %>%
  summarise(
    coverage_median = median(observed_coverage, na.rm = TRUE),
    coverage_mean = mean(observed_coverage, na.rm = TRUE),
    coverage_max = max(observed_coverage, na.rm = TRUE)
  )

coverage_note <- paste0(
  "Current ERA5 annual budgets have low observed half-hourly coverage ",
  "(median ", percent(coverage_summary$coverage_median, accuracy = 0.1),
  ", mean ", percent(coverage_summary$coverage_mean, accuracy = 0.1),
  ", max ", percent(coverage_summary$coverage_max, accuracy = 0.1),
  "), so NEON annual values should be interpreted as model-gapfilled estimates."
)

era5_site_summary <- era5_annual_budget %>%
  left_join(
    era5_mean_annual_budget %>%
      dplyr::select(SITE_ID, era5_annual_behavior, mean_observed_coverage),
    by = "SITE_ID"
  ) %>%
  mutate(
    era5_annual_behavior = factor(era5_annual_behavior, levels = behavior_levels),
    gradient_behavior_for_plot = coalesce(era5_annual_behavior, CH4_gradient_behavior, CH4_behavior)
  ) %>%
  group_by(SITE_ID, gradient_behavior_for_plot, CH4_gradient_behavior, CH4_behavior, EcoType, is_arid) %>%
  summarise(
    n_years = n(),
    mean_observed_coverage = mean(observed_coverage, na.rm = TRUE),
    median_n_observed_30min_per_year = median(n_observed, na.rm = TRUE),
    min_n_observed_30min_per_year = min(n_observed, na.rm = TRUE),
    daily_mgC_m2_day = median(daily_mgC_m2_day, na.rm = TRUE),
    daily_q25_mgC_m2_day = quantile(daily_mgC_m2_day, 0.25, na.rm = TRUE),
    daily_q75_mgC_m2_day = quantile(daily_mgC_m2_day, 0.75, na.rm = TRUE),
    annual_gC_m2_yr = median(annual_budget_gC_m2_yr, na.rm = TRUE),
    annual_q25_gC_m2_yr = quantile(annual_budget_gC_m2_yr, 0.25, na.rm = TRUE),
    annual_q75_gC_m2_yr = quantile(annual_budget_gC_m2_yr, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    gradient_behavior_for_plot = factor(gradient_behavior_for_plot, levels = behavior_levels),
    ecosystem_group = assign_neon_ecosystem(EcoType, is_arid),
    era5_gradient_class = recode(
      as.character(gradient_behavior_for_plot),
      "Weak-sink" = "NEON ERA5 gapfilled sink",
      "Fluctuating" = "NEON ERA5 gapfilled fluctuating",
      "Weak-source" = "NEON ERA5 gapfilled source"
    )
  )

era5_class_summary <- era5_site_summary %>%
  group_by(era5_gradient_class, gradient_behavior_for_plot) %>%
  summarise(
    comparison_group = "NEON ERA5 gapfilled gradient flux",
    ecosystem_class = first(era5_gradient_class),
    n_sites = n_distinct(SITE_ID),
    n_site_years = sum(n_years, na.rm = TRUE),
    daily_median_mgC_m2_day = median(daily_mgC_m2_day, na.rm = TRUE),
    daily_low_mgC_m2_day = quantile(daily_mgC_m2_day, 0.25, na.rm = TRUE),
    daily_high_mgC_m2_day = quantile(daily_mgC_m2_day, 0.75, na.rm = TRUE),
    daily_sd_mgC_m2_day = sd(daily_mgC_m2_day, na.rm = TRUE),
    reference_id = "This analysis: ERA5-gapfilled annual gradient flux budgets from OUTPUT/NEON_ERA5_gapfilled_annual_budget_by_year.csv, converted to daily rates",
    source_type = "NEON ERA5 gapfilled gradient",
    .groups = "drop"
  ) %>%
  mutate(
    CH4_behavior = gradient_behavior_for_plot,
    daily_mgC_m2_day = daily_median_mgC_m2_day,
    annual_gC_m2_yr = daily_mgC_m2_day * 365 / 1000,
    annual_sd_gC_m2_yr = daily_sd_mgC_m2_day * 365 / 1000
  )

comparison_columns <- c(
  "comparison_group", "ecosystem_class", "CH4_behavior", "annual_gC_m2_yr",
  "annual_sd_gC_m2_yr", "daily_mgC_m2_day", "daily_low_mgC_m2_day",
  "daily_high_mgC_m2_day", "n_sites", "n_site_years", "reference_id",
  "source_type", "site_id", "location", "latitude", "longitude", "igbp",
  "measurement_years", "n_observations", "flux_units_original", "full_reference"
)

comparison_table <- bind_rows(
  fluxnet_reference %>%
    mutate(CH4_behavior = NA_character_) %>%
    dplyr::select(any_of(comparison_columns)),
  soil_chamber_reference %>%
    mutate(CH4_behavior = NA_character_) %>%
    dplyr::select(any_of(comparison_columns)),
  additional_draft_chamber_reference %>%
    dplyr::select(any_of(comparison_columns)),
  smud_chamber_reference %>%
    dplyr::select(any_of(comparison_columns)),
  process_model_uptake %>%
    mutate(CH4_behavior = NA_character_) %>%
    dplyr::select(any_of(comparison_columns)),
  non_fluxnet_tower_reference %>%
    dplyr::select(any_of(comparison_columns)),
  additional_draft_tower_reference %>%
    dplyr::select(any_of(comparison_columns)),
  smud_ec_reference %>%
    dplyr::select(any_of(comparison_columns)),
  era5_class_summary %>%
    dplyr::select(any_of(comparison_columns))
) %>%
  arrange(desc(daily_mgC_m2_day)) %>%
  mutate(
    source_type = factor(
      source_type,
      levels = c(
        "Published ecosystem class",
        "Soil chamber literature",
        "SMUD soil chamber",
        "Process-based model",
        "Tower literature outside FLUXNET-CH4",
        "SMUD tower (EC)",
        "NEON ERA5 gapfilled gradient"
      )
    ),
    # List-column: most rows get exactly one ecosystem_group, but composite
    # reference classes (Upland, Drained, Tundra) get several -- see
    # assign_reference_ecosystem_groups(). unnest_longer() below duplicates
    # those rows' flux estimate into every applicable panel.
    ecosystem_group = purrr::map2(
      as.character(source_type), ecosystem_class,
      function(st, ec) {
        if (identical(st, "NEON ERA5 gapfilled gradient")) return("Across NEON sites")
        assign_reference_ecosystem_groups(ec)[[1]]
      }
    )
  ) %>%
  tidyr::unnest_longer(ecosystem_group)

write.csv(comparison_table, "OUTPUT/CH4_flux_FLUXNET_NEON_comparison_values.csv", row.names = FALSE)
write.csv(era5_annual_budget, "OUTPUT/NEON_ERA5_annual_gradient_flux_rates_for_FLUXNET_comparison.csv", row.names = FALSE)
write.csv(era5_site_summary, "OUTPUT/NEON_ERA5_site_median_daily_gradient_flux_classes.csv", row.names = FALSE)
write.csv(process_model_uptake, "OUTPUT/process_model_upland_CH4_uptake_values.csv", row.names = FALSE)
write.csv(soil_chamber_reference, "OUTPUT/soil_chamber_CH4_flux_reference_values.csv", row.names = FALSE)
write.csv(additional_draft_chamber_reference, "OUTPUT/draft_reference_soil_chamber_CH4_flux_values.csv", row.names = FALSE)
write.csv(non_fluxnet_tower_reference, "OUTPUT/non_FLUXNET_CH4_tower_flux_reference_values.csv", row.names = FALSE)
write.csv(non_fluxnet_tower_years, "OUTPUT/non_FLUXNET_CH4_tower_flux_reference_years.csv", row.names = FALSE)
write.csv(additional_draft_tower_reference, "OUTPUT/draft_reference_non_FLUXNET_tower_CH4_flux_values.csv", row.names = FALSE)
write.csv(draft_candidate_sources, "OUTPUT/draft_reference_CH4_candidate_sources_for_extraction.csv", row.names = FALSE)
write.csv(draft_candidate_extracted_info, "OUTPUT/draft_reference_CH4_extracted_candidate_info.csv", row.names = FALSE)
if (nrow(smud_raw) > 0) {
  write.csv(smud_raw, "OUTPUT/SMUD_CH4_flux_all_filtered_rows.csv", row.names = FALSE)
  write.csv(smud_ec_reference, "OUTPUT/SMUD_CH4_EC_tower_reference_by_ecosystem.csv", row.names = FALSE)
  write.csv(smud_chamber_reference, "OUTPUT/SMUD_CH4_chamber_reference_by_ecosystem.csv", row.names = FALSE)
}

# ── Summary table: medians by data source and state class ─────────────────────

neon_for_table <- era5_site_summary %>%
  transmute(
    source_label     = "NEON ERA5",
    daily_mgC_m2_day = daily_mgC_m2_day,
    behavior = case_when(
      daily_mgC_m2_day < 0 ~ "Weak-sink",
      daily_mgC_m2_day > 0 ~ "Weak-source",
      TRUE                  ~ "Fluctuating"
    )
  )

ref_for_table <- comparison_table %>%
  filter(source_type != "NEON ERA5 gapfilled gradient") %>%
  mutate(
    source_label = case_when(
      source_type == "Published ecosystem class"            ~ "FLUXNET-CH4",
      source_type %in% c("Soil chamber literature",
                         "SMUD soil chamber")              ~ "Chambers",
      source_type == "Process-based model"                 ~ "Process model",
      source_type %in% c("Tower literature outside FLUXNET-CH4",
                         "SMUD tower (EC)")                ~ "Upland towers",
      TRUE                                                 ~ as.character(source_type)
    ),
    behavior = case_when(
      daily_mgC_m2_day < 0 ~ "Weak-sink",
      daily_mgC_m2_day > 0 ~ "Weak-source",
      TRUE                  ~ "Fluctuating"
    )
  ) %>%
  dplyr::select(source_label, daily_mgC_m2_day, behavior)

source_behavior_summary <- bind_rows(ref_for_table, neon_for_table) %>%
  filter(behavior %in% c("Weak-sink", "Weak-source"), is.finite(daily_mgC_m2_day)) %>%
  group_by(source_label, behavior) %>%
  summarise(
    n                 = n(),
    median_mgC_m2_day = median(daily_mgC_m2_day, na.rm = TRUE),
    q25_mgC_m2_day    = quantile(daily_mgC_m2_day, 0.25, na.rm = TRUE),
    q75_mgC_m2_day    = quantile(daily_mgC_m2_day, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    source_label = factor(source_label,
                          levels = c("FLUXNET-CH4", "Chambers", "Process model",
                                     "Upland towers", "NEON ERA5")),
    behavior     = factor(behavior, levels = c("Weak-sink", "Weak-source"))
  ) %>%
  arrange(source_label, behavior)

write.csv(source_behavior_summary,
          "OUTPUT/CH4_flux_medians_by_source_and_behavior.csv",
          row.names = FALSE)
message("Wrote OUTPUT/CH4_flux_medians_by_source_and_behavior.csv")

reference_lines <- c(
  "# References and Values Used",
  "",
  "## Published FLUXNET-CH4 ecosystem-class rates",
  "Delwiche, K. B., Knox, S. H., Malhotra, A., Fluet-Chouinard, E., McNicol, G., Feron, S., et al. (2021). FLUXNET-CH4: a global, multi-ecosystem dataset and analysis of methane seasonality from freshwater wetlands. Earth System Science Data, 13, 3607-3689. https://doi.org/10.5194/essd-13-3607-2021",
  "",
  "Values used from Delwiche et al. (2021) Table 1: salt marsh 2.9 +/- 4.7, wet tundra 3.8 +/- 1.8, upland 4.0 +/- 10.5, drained 6.3 +/- 7.1, bog 10.5 +/- 6.4, mangrove 11.1 +/- 0.5, rice 14.4 +/- 8.8, fen 20.5 +/- 16.0, swamp 26.4 +/- 19.9, lake 28.2 +/- 33.4, marsh 40.8 +/- 20.7 g C m-2 yr-1.",
  "Urban value used from Delwiche et al. (2021), Sect. 3.1.3: UK-LBT urban CH4 flux 46.5 +/- 5.6 g C m-2 yr-1.",
  "",
  "All published values were converted to mg C m-2 d-1 as: annual g C m-2 yr-1 * 1000 / 365.",
  "",
  "## Process-based upland soil CH4 uptake model rates",
  "Murguia-Flores, F., Arndt, S., Ganesan, A. L., Murray-Tortarolo, G., and Hornibrook, E. R. C. (2018). Soil Methanotrophy Model (MeMo v1.0): a process-based model to quantify global uptake of atmospheric methane by soil. Geoscientific Model Development, 11, 2009-2032. https://doi.org/10.5194/gmd-11-2009-2018",
  "",
  "Values used from Murguia-Flores et al. (2018) Table 9: tropical deciduous forest 602 +/- 63, open shrubland 518 +/- 134, temperate broadleaf forest 512 +/- 82, savanna 500 +/- 132, dense shrubland 481 +/- 90, grassland/steppe 392 +/- 110, temperate needleleaf forest 347 +/- 90, tropical evergreen forest 332 +/- 45, temperate deciduous forest 321 +/- 70, boreal deciduous forest 282 +/- 117, boreal evergreen forest 269 +/- 94, mixed forest 182 +/- 82, tundra 176 +/- 143, polar desert/rock/ice 105 +/- 48 mg CH4 m-2 yr-1.",
  "Process-model uptake values were converted to the NEON flux sign convention as: -1 * mg CH4 m-2 yr-1 * 12 / 16 / 365 = mg C m-2 d-1. Negative values indicate uptake.",
  "",
  "## Soil chamber CH4 flux reference rates",
  "Soil chamber values are included as literature benchmarks distinct from ecosystem-scale FLUXNET-CH4 and process-model estimates.",
  "Whalen (2005) reports that wetland CH4 emissions are commonly about 100 mg CH4 m-2 d-1; this was converted to 75 mg C m-2 d-1.",
  "Adamsen and King (1993) report temperate/subarctic forest soil CH4 uptake generally between 1 and 3 mg CH4 m-2 d-1; this was converted to -0.75 to -2.25 mg C m-2 d-1.",
  "Groffman and Pouyat (2009) report rural forest CH4 uptake of 1.68 mg CH4 m-2 d-1; this was converted to -1.26 mg C m-2 d-1.",
  "Striegl et al. (1992) report that the central 50% of desert soil CH4 uptake rates were 0.24 to 0.92 mg CH4 m-2 d-1; this was converted to -0.18 to -0.69 mg C m-2 d-1.",
  "Mosier et al. (1991) report aerobic soil CH4-C uptake of about 1 to 3 kg CH4-C ha-1 yr-1 across diverse ecosystems including grasslands; this was converted to -0.274 to -0.822 mg C m-2 d-1.",
  "A Danish farmland chamber study reported CH4 fluxes from -0.43 to 0.19 mg CH4 m-2 d-1, with mean -0.15 mg CH4 m-2 d-1; this was converted to -0.3225 to 0.1425 mg C m-2 d-1, with mean -0.1125 mg C m-2 d-1.",
  "Sundqvist et al. (2015) report Norunda automated soil-chamber uptake of around -10 umol CH4 m-2 h-1; this was converted to -2.88 mg C m-2 d-1.",
  "",
  "## Tower CH4 flux rates outside FLUXNET-CH4",
  "Additional tower benchmarks were extracted from the validation-tower source files associated with references included in the manuscript draft. Half-hourly tower CH4 fluxes were summarized to annual mean daily rates, years with fewer than 100 finite half-hourly observations were excluded, and site-level values are medians across retained years.",
  "Lakomiec et al. (2021), Stordalen Mire (SE-Sto), subarctic mire, Abisko, Sweden: source file SE-Sto_gas_fluxes_30min.csv, Fch4_1_1_1 in umol CH4 m-2 s-1, converted to mg C m-2 d-1.",
  "Chi et al. (2020), Svartberget/Krycklan (SE-Svb), managed boreal forest, Sweden: source file CH4_SE_SVB_FLUX+PROFILE_2019.csv, ch4_flux_nmolm2s_85m in nmol CH4 m-2 s-1, converted to mg C m-2 d-1.",
  "Iwata et al. (2015), US-Uaf, poorly drained black spruce forest over permafrost, Fairbanks, Alaska, USA: source file US-Uaf CH4_concentration.csv, CH4 flux in nmol CH4 m-2 s-1, converted to mg C m-2 d-1.",
  "Smeets et al. (2009), Blodgett Forest, ponderosa pine plantation, California, USA: whole-period daily mean downward CH4 flux of 2.5 mg CH4 m-2 d-1 converted to -1.875 mg C m-2 d-1.",
  "Sundqvist et al. (2015), Norunda research station, boreal mixed pine-spruce forest, Sweden: tower-method mean emissions of 1.48-4.57 umol CH4 m-2 h-1 converted to 0.43-1.32 mg C m-2 d-1, with midpoint 0.87 mg C m-2 d-1.",
  "Wang et al. (2013), Haliburton Forest and Wildlife Reserve, temperate forest, Ontario, Canada: EC uptake flux of -2.7 +/- 0.13 nmol CH4 m-2 s-1 converted to -2.80 mg C m-2 d-1.",
  "Tower coordinates, when available, came from metadata_validation.csv.",
  "",
  "## Candidate CH4 data sources from draft references",
  "Additional draft references with potential CH4 flux data were inventoried in OUTPUT/draft_reference_CH4_candidate_sources_for_extraction.csv. These rows are not plotted until a numeric CH4 flux rate, units, location, and extraction note are confirmed.",
  "Accessible information extracted from the candidate sources is recorded in OUTPUT/draft_reference_CH4_extracted_candidate_info.csv. This table separates plot-ready net fluxes from process-only values, relative comparisons, synthesis ranges, and sources that still require full-text/table extraction.",
  "Treat et al. (2018) provided broad synthesis ranges for annual wetland emissions and upland annual fluxes. These were not added to the plot because they summarize many ecosystem types and seasons rather than a single site/source class.",
  "von Fischer and Hedin (2002) and Angle et al. (2017) were retained as mechanistic evidence for co-occurring production/consumption and oxic methanogenesis, but not plotted because they report gross process rates or percent attribution rather than directly comparable net CH4 fluxes.",
  "",
  "## NEON ERA5-gapfilled gradient values",
  "NEON values are from this repository's ERA5-gapfilled annual gradient/tower flux output: OUTPUT/NEON_ERA5_gapfilled_annual_budget_by_year.csv and OUTPUT/NEON_ERA5_gapfilled_mean_annual_budget.csv.",
  "Annual ERA5 budgets retained observed half-hourly gradient/tower fluxes where available and filled missing half-hours with predictions from the ERA5-driven GAM in NEON.ERA5.HalfHourlyGapfill.R.",
  "Annual budgets in g C m-2 yr-1 were converted to daily rates as annual budget * 1000 / 365 = mg C m-2 d-1.",
  "NEON points in the figure are site-level medians across ERA5-gapfilled site-years. NEON class summaries are grouped by ERA5 annual behavior class when available, with gradient/total behavior classes from OUTPUT/30min_site_behavior.csv as fallbacks.",
  coverage_note,
  "",
  "## Output tables",
  "- OUTPUT/CH4_flux_FLUXNET_NEON_comparison_values.csv",
  "- OUTPUT/NEON_ERA5_annual_gradient_flux_rates_for_FLUXNET_comparison.csv",
  "- OUTPUT/NEON_ERA5_site_median_daily_gradient_flux_classes.csv",
  "- OUTPUT/process_model_upland_CH4_uptake_values.csv",
  "- OUTPUT/soil_chamber_CH4_flux_reference_values.csv",
  "- OUTPUT/draft_reference_soil_chamber_CH4_flux_values.csv",
  "- OUTPUT/non_FLUXNET_CH4_tower_flux_reference_values.csv",
  "- OUTPUT/non_FLUXNET_CH4_tower_flux_reference_years.csv",
  "- OUTPUT/draft_reference_non_FLUXNET_tower_CH4_flux_values.csv",
  "- OUTPUT/draft_reference_CH4_candidate_sources_for_extraction.csv",
  "- OUTPUT/draft_reference_CH4_extracted_candidate_info.csv"
)
writeLines(reference_lines, "OUTPUT/CH4_flux_FLUXNET_NEON_comparison_references.md")

source_label_levels <- c(
  "FLUXNET-CH4", "Chambers", "Process model", "Upland towers", "NEON ERA5"
)

# Row shown within each ecosystem facet panel. Ecosystem is now handled by
# facet_wrap(), so this only needs to separate sources vertically within a
# panel. Condensed from the earlier 8-row layout: Validation towers are
# folded into Upland towers (both are non-FLUXNET tower measurements), and
# NEON's three state-class rows (Weak-sink/Fluctuating/Weak-source) collapse
# into a single "NEON" row — sink vs. source is now shown by bar direction
# and color rather than by row.
source_row_levels <- c(
  "FLUXNET-CH4", "Chambers", "Process model", "Upland towers", "NEON"
)

# Fixed, shared weak/strong magnitude threshold (±10 mg C m-2 d-1) used for
# both sides, so the background shading is symmetric around zero. Replaces
# the earlier data-driven threshold (median |flux| on the sink side, ≈0.62).
# Note this fixed value still means most real NEON "weak-source" sites
# (source magnitudes are ~20x larger than sink magnitudes) may fall into the
# "strong-source" zone rather than "weak-source" — same tradeoff as before,
# now with a round, easy-to-communicate cutoff.
sink_threshold   <- 10
source_threshold <- 10
message(sprintf(
  "Weak/strong threshold (fixed, shared): %.3f mg C m-2 d-1",
  sink_threshold
))

# Individual literature/ecosystem-class estimates, one row per estimate. These
# feed into the condensed per-source bar chart below (bar_source_data).
reference_plot_data <- comparison_table %>%
  filter(source_type != "NEON ERA5 gapfilled gradient") %>%
  mutate(
    source_label = recode(
      as.character(source_type),
      "Published ecosystem class"            = "FLUXNET-CH4",
      "Soil chamber literature"              = "Chambers",
      "SMUD soil chamber"                    = "Chambers",
      "Process-based model"                  = "Process model",
      "Tower literature outside FLUXNET-CH4" = "Upland towers",
      "SMUD tower (EC)"                      = "Upland towers"
    ),
    source_label = factor(source_label, levels = source_label_levels),
    source_row = factor(as.character(source_label), levels = source_row_levels),
    ecosystem_group = factor(ecosystem_group, levels = ecosystem_levels)
  ) %>%
  filter(!is.na(source_row), !is.na(ecosystem_group), is.finite(daily_mgC_m2_day))

# SMUD aggregated by ecosystem × behavior class (sink / source shown separately)
smud_behavior_plot_data <- if (nrow(smud_raw) > 0) {
  smud_raw %>%
    mutate(
      smud_source    = if_else(is_ec, "SMUD EC", "SMUD chambers"),
      ecosystem_group = factor(ecosystem_group, levels = ecosystem_levels),
      CH4_behavior    = factor(CH4_behavior, levels = behavior_levels)
    ) %>%
    group_by(ecosystem_group, smud_source, CH4_behavior) %>%
    summarise(
      n_obs                 = n(),
      daily_mgC_m2_day      = median(Annual, na.rm = TRUE),
      daily_low_mgC_m2_day  = quantile(Annual, 0.25, na.rm = TRUE),
      daily_high_mgC_m2_day = quantile(Annual, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_base = length(ecosystem_levels) - as.integer(ecosystem_group) + 1,
      y_plot = y_base
    ) %>%
    filter(!is.na(y_base), is.finite(daily_mgC_m2_day))
} else {
  tibble()
}

# Validation towers (non-FLUXNET): one point per directly-measured site only.
# additional_draft_tower_reference (US-Blo, SE-Nor, CA-Hal) is already aggregated
# into the "Upland towers" symbol and should not appear here as well. These are
# folded into the "Upland towers" row below (condensed row set) rather than
# kept as a separate row.
validation_tower_points <- non_fluxnet_tower_reference %>%
  dplyr::select(site_id, ecosystem_class, daily_mgC_m2_day) %>%
  filter(is.finite(daily_mgC_m2_day)) %>%
  mutate(ecosystem_group = assign_reference_ecosystem_groups(ecosystem_class)) %>%
  tidyr::unnest_longer(ecosystem_group) %>%
  mutate(
    ecosystem_group = factor(ecosystem_group, levels = ecosystem_levels),
    source_row      = factor("Upland towers", levels = source_row_levels)
  ) %>%
  filter(!is.na(ecosystem_group))

neon_site_points <- era5_site_summary %>%
  mutate(
    ecosystem_group = factor(ecosystem_group, levels = ecosystem_levels),
    source_row      = factor("NEON", levels = source_row_levels)
  ) %>%
  filter(!is.na(source_row), is.finite(daily_mgC_m2_day))

axis_text_size <- 10

# ── One boxplot per data source per ecosystem panel ─────────────────────────
# Replaces the earlier point/errorbar-per-estimate design. Every source now
# contributes to a single long-format table (ecosystem_group, source_row,
# daily_mgC_m2_day); boxplots are then grouped by the SIGN of the flux within
# each (ecosystem_group, source_row) cell, so a row shows one box if all its
# values share a sign, or two boxes (one extending left, one right of zero) if
# both sink and source values are present — this is how NEON's own
# weak-sink/weak-source split is preserved without a dedicated row per class.
# Box color encodes sign only (blue = sink, red = source), matching the
# background shading, so no shape/pch legend is needed.
bar_source_data <- bind_rows(
  reference_plot_data %>%
    filter(source_label %in% c("FLUXNET-CH4", "Chambers", "Process model")) %>%
    transmute(ecosystem_group, source_row = as.character(source_label), daily_mgC_m2_day),
  reference_plot_data %>%
    filter(source_label == "Upland towers") %>%
    transmute(ecosystem_group, source_row = "Upland towers", daily_mgC_m2_day),
  validation_tower_points %>%
    transmute(ecosystem_group, source_row = "Upland towers", daily_mgC_m2_day),
  neon_site_points %>%
    transmute(ecosystem_group, source_row = "NEON", daily_mgC_m2_day)
) %>%
  filter(is.finite(daily_mgC_m2_day)) %>%
  mutate(
    source_row = factor(source_row, levels = source_row_levels),
    sign_class = case_when(
      daily_mgC_m2_day < 0 ~ "Sink",
      daily_mgC_m2_day > 0 ~ "Source",
      TRUE                 ~ "Fluctuating"
    )
  )

bar_source_data <- bar_source_data %>%
  mutate(sign_class = factor(sign_class, levels = c("Sink", "Fluctuating", "Source")))

# Shared builder for the two summary facets below ("All ecosystems", "All
# Upland"): pools the same four source rows, deduplicated on the ORIGINAL
# record identity (ecosystem_class / site_id), not the already-exploded
# ecosystem_group, so composite reference rows replicated across several
# ecosystem panels (FLUXNET's "Upland"/"Drained" aggregates, any "tundra"
# class -- see assign_reference_ecosystem_groups()) are only counted once
# rather than once per panel they appear in. When restrict_upland = TRUE,
# only rows whose exploded ecosystem_group falls in upland_eco_groups are
# pooled (mirrors the density plot's upland-only filtering).
build_summary_facet_data <- function(facet_label, restrict_upland = FALSE) {
  ref_data  <- reference_plot_data
  val_data  <- validation_tower_points
  neon_data <- neon_site_points
  if (restrict_upland) {
    ref_data  <- ref_data  %>% filter(as.character(ecosystem_group) %in% upland_eco_groups)
    val_data  <- val_data  %>% filter(as.character(ecosystem_group) %in% upland_eco_groups)
    neon_data <- neon_data %>% filter(as.character(ecosystem_group) %in% upland_eco_groups)
  }
  bind_rows(
    ref_data %>%
      filter(source_label %in% c("FLUXNET-CH4", "Chambers", "Process model")) %>%
      distinct(source_label, ecosystem_class, daily_mgC_m2_day, .keep_all = TRUE) %>%
      transmute(source_row = as.character(source_label), daily_mgC_m2_day),
    ref_data %>%
      filter(source_label == "Upland towers") %>%
      distinct(ecosystem_class, daily_mgC_m2_day, .keep_all = TRUE) %>%
      transmute(source_row = "Upland towers", daily_mgC_m2_day),
    val_data %>%
      distinct(site_id, daily_mgC_m2_day, .keep_all = TRUE) %>%
      transmute(source_row = "Upland towers", daily_mgC_m2_day),
    neon_data %>%
      transmute(source_row = "NEON", daily_mgC_m2_day)
  ) %>%
    filter(is.finite(daily_mgC_m2_day)) %>%
    mutate(
      ecosystem_group = factor(facet_label, levels = ecosystem_levels),
      source_row = factor(source_row, levels = source_row_levels),
      sign_class = case_when(
        daily_mgC_m2_day < 0 ~ "Sink",
        daily_mgC_m2_day > 0 ~ "Source",
        TRUE                 ~ "Fluctuating"
      ),
      sign_class = factor(sign_class, levels = c("Sink", "Fluctuating", "Source"))
    )
}

all_ecosystems_data <- build_summary_facet_data("All ecosystems")
all_upland_data     <- build_summary_facet_data("All Upland", restrict_upland = TRUE)

bar_source_data <- bind_rows(bar_source_data, all_ecosystems_data, all_upland_data)

sign_colors      <- c("Sink" = "#2166AC", "Fluctuating" = "grey35", "Source" = "#B2182B")
sign_colors_dark <- c("Sink" = "#0A3161", "Fluctuating" = "grey15", "Source" = "#7A1216")

comparison_figure <- ggplot() +
  # Four magnitude zones (repeated automatically in every facet panel). Zone
  # labels aren't drawn per-panel — with several ecosystem facets that was too much
  # repeated text — the boundaries are explained once in the plot subtitle
  # instead, and marked with dotted threshold lines below.
  annotate(  # strong-sink
    "rect",
    xmin = -Inf, xmax = -sink_threshold, ymin = -Inf, ymax = Inf,
    fill = "#2166AC", alpha = 0.30
  ) +
  annotate(  # weak-sink
    "rect",
    xmin = -sink_threshold, xmax = 0, ymin = -Inf, ymax = Inf,
    fill = "#BFDDF5", alpha = 0.52
  ) +
  annotate(  # weak-source
    "rect",
    xmin = 0, xmax = source_threshold, ymin = -Inf, ymax = Inf,
    fill = "#F5C0BD", alpha = 0.48
  ) +
  annotate(  # strong-source
    "rect",
    xmin = source_threshold, xmax = Inf, ymin = -Inf, ymax = Inf,
    fill = "#B2182B", alpha = 0.28
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.5) +
  geom_vline(
    xintercept = c(-sink_threshold, source_threshold),
    linetype = "dotted", color = "grey45", linewidth = 0.4
  ) +
  # One boxplot per (ecosystem, source, sign) — up to two boxes per row when a
  # source has both sink and source values in a given ecosystem panel. Boxes
  # sit at the same y position (source_row, no dodge) and separate naturally
  # along x because sink values are negative and source values are positive.
  geom_boxplot(
    data = bar_source_data,
    aes(x = daily_mgC_m2_day, y = source_row, fill = sign_class, color = sign_class),
    orientation = "y", width = 0.6, linewidth = 0.4,
    alpha = 0.85, outlier.size = 1.2, outlier.alpha = 0.6,
    position = position_identity()
  ) +
  # Validation tower sites overlaid as stars so their individual values are
  # visible within the "Upland towers" box (they're folded into that row's
  # boxplot above and would otherwise be indistinguishable from it).
  geom_point(
    data = validation_tower_points,
    aes(x = daily_mgC_m2_day, y = source_row, shape = "Validation tower site"),
    color = "black", fill = "black", size = 2.6, stroke = 0.8,
    position = position_jitter(height = 0.12, width = 0, seed = 20260525)
  ) +
  facet_wrap(~ ecosystem_group, ncol = 3) +
  scale_fill_manual(values = sign_colors, breaks = c("Sink", "Source"), name = "Flux direction") +
  scale_color_manual(values = sign_colors_dark, breaks = c("Sink", "Source"), guide = "none") +
  scale_shape_manual(values = c("Validation tower site" = 8), name = NULL) +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 0.01),
    breaks = c(-20, -5, -1, -0.1, 0, 0.1, 1, 10, 100),
    labels = c("-20", "-5", "-1", "-0.1", "0", "0.1", "1", "10", "100")
  ) +
  scale_y_discrete(limits = rev(source_row_levels), drop = FALSE) +
  labs(
  #  title = expression(bold("ERA5-gapfilled NEON gradient CH"[4]*" fluxes vs. published benchmarks")),
    x = expression("Daily CH"[4] * " flux (mg C m"^-2 * " d"^-1 * "; pseudo-log scale)"),
    y = NULL
  ) +
  theme_bw(base_size = base_plot_size) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    plot.title = element_text(face = "bold", size = panel_title_size),
    plot.margin = margin(t = 5, r = 15, b = 5, l = 5),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    axis.text.y = element_text(size = axis_text_size - 1),
    legend.title = element_text(size = legend_title_size),
    legend.text = element_text(size = legend_text_size),
    strip.background = element_rect(fill = "black", color = NA),
    strip.text = element_text(face = "bold", size = axis_title_size, color = "white"),
    panel.spacing = unit(1, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "FIGURES/NEON_FLUXNET_CH4_flux_comparison.png",
  plot = comparison_figure,
  width = 16,
  height = 15,
  units = "in",
  dpi = 300
)

message("Wrote FIGURES/NEON_FLUXNET_CH4_flux_comparison.png")

# Detailed caption, written alongside the figure rather than as an in-plot
# subtitle (the full explanation doesn't fit on the plot itself). Threshold
# value is pulled from sink_threshold so the text always matches the figure.
comparison_figure_caption <- paste0(
  "Figure. ERA5-gapfilled NEON gradient CH4 fluxes compared with published benchmarks. ",
  "Each panel is an ecosystem class; within a panel, rows are data sources: ",
  "FLUXNET-CH4 (Delwiche et al. 2021 published ecosystem-class values), Chambers ",
  "(soil chamber literature, including the SMUD compilation), Process model ",
  "(MeMo v1.0 process-based estimates), Upland towers (non-FLUXNET tower literature ",
  "and SMUD eddy-covariance towers; individual validation-tower site values are ",
  "additionally marked with black stars), and NEON (ERA5-gapfilled gradient/tower ",
  "CH4 budgets for individual NEON sites). Each row is drawn as a boxplot (median, ",
  "interquartile range, and whiskers to the data range) rather than a single point ",
  "estimate, so the spread of values contributing to each source/ecosystem cell is ",
  "visible. A row shows two boxes when a source has both sink (flux < 0) and source ",
  "(flux > 0) values within an ecosystem class -- most visibly for NEON, whose sites ",
  "split into weak-sink and weak-source behavior within several ecosystem classes. ",
  "Box and background shading color encodes flux direction only (blue = sink, ",
  "red = source), not data source. ",
  "Background shading divides each panel into four magnitude zones -- strong-sink, ",
  "weak-sink, weak-source, strong-source -- separated at a fixed, shared threshold ",
  "of ±", round(sink_threshold, 2), " mg C m-2 d-1. The same threshold magnitude is ",
  "used on both the sink and source sides so the shading is symmetric around zero; ",
  "because source fluxes are ",
  "typically an order of magnitude larger than sink fluxes (e.g., wetland CH4 ",
  "emission vs. upland CH4 uptake), most individual source-side estimates -- ",
  "including many NEON sites classified elsewhere in this analysis as ",
  "\"weak-source\" -- fall in the strong-source shaded zone here; the shading should ",
  "be read as a magnitude reference scale rather than a formal sink/source ",
  "classification. The x-axis uses a pseudo-log transform (linear near zero, ",
  "logarithmic at larger magnitudes) to accommodate the wide range of flux ",
  "magnitudes across ecosystem classes on a single scale."
)
writeLines(
  strwrap(comparison_figure_caption, width = 100),
  "FIGURES/NEON_FLUXNET_CH4_flux_comparison_caption.txt"
)
message("Wrote FIGURES/NEON_FLUXNET_CH4_flux_comparison_caption.txt")

# ── Density plot: NEON towers vs upland towers vs chambers ────────────────────

density_group_colors <- c(
  "Towers"   = "maroon",
  "Chambers" = "navy"
)

density_plot_data <- bind_rows(
  # Towers: NEON ERA5 site medians — upland only
  era5_site_summary %>%
    filter(ecosystem_group %in% upland_eco_groups) %>%
    transmute(flux = daily_mgC_m2_day, group = "Towers"),
  # Towers: site-year values from measured validation towers — upland only
  non_fluxnet_tower_years %>%
    filter(is_any_upland(ecosystem_class)) %>%
    transmute(flux = daily_mgC_m2_day, group = "Towers"),
  # Towers: literature summary values (all upland forests)
  additional_draft_tower_reference %>%
    transmute(flux = daily_mgC_m2_day, group = "Towers"),
  # Towers: individual SMUD EC studies — upland only
  smud_raw %>%
    filter(is_ec, ecosystem_group %in% upland_eco_groups) %>%
    transmute(flux = Annual, group = "Towers"),
  # Chambers: literature summary values — upland only
  soil_chamber_reference %>%
    filter(is_any_upland(ecosystem_class)) %>%
    transmute(flux = daily_mgC_m2_day, group = "Chambers"),
  additional_draft_chamber_reference %>%
    transmute(flux = daily_mgC_m2_day, group = "Chambers"),
  # Chambers: individual SMUD chamber studies — upland only
  smud_raw %>%
    filter(!is_ec, ecosystem_group %in% upland_eco_groups) %>%
    transmute(flux = Annual, group = "Chambers")
) %>%
  filter(is.finite(flux)) %>%
  mutate(group = factor(group, levels = names(density_group_colors)))

density_n_labels <- density_plot_data %>%
  count(group) %>%
  mutate(label = paste0(group, "\n(n = ", n, ")"))

# Center each zone label within the section actually visible on the plot.
# Zone midpoints must be computed in the SAME pseudo-log transformed space
# the x-axis is drawn in, not in raw data units — the x-axis is heavily
# log-compressed (sigma = 0.01), so e.g. the raw arithmetic midpoint of the
# weak-source zone (0 to 10) is 5, but 5 sits ~90% of the way through that
# zone visually (it's deep in the last decade, 1-10), which is why the
# earlier raw-mean version put labels hard against the zone edges instead of
# centered. Transforming, averaging, then inverting back to data units fixes
# this for all four zones. Strong-sink/strong-source still use the real data
# extent (density_flux_range) as their finite outer bound, since ±Inf has no
# transformed value to average.
density_flux_range <- range(density_plot_data$flux, na.rm = TRUE)
pslog <- pseudo_log_trans(sigma = 0.01)
zone_mid <- function(lo, hi) pslog$inverse(mean(pslog$transform(c(lo, hi))))
zone_label_x <- c(
  zone_mid(density_flux_range[1], -sink_threshold),   # strong-sink
  zone_mid(-sink_threshold, 0),                        # weak-sink
  zone_mid(0, source_threshold),                       # weak-source
  zone_mid(source_threshold, density_flux_range[2])    # strong-source
)

flux_density_figure <- ggplot(density_plot_data, aes(x = flux, fill = group, color = group)) +
  # Background shading matches comparison_figure's sink/source zone colors
  # and thresholds (sink_threshold == source_threshold — shared, balanced cutoff).
  annotate(  # strong-sink
    "rect", xmin = -Inf, xmax = -sink_threshold, ymin = -Inf, ymax = Inf,
    fill = "#2166AC", alpha = 0.30
  ) +
  annotate(  # weak-sink
    "rect", xmin = -sink_threshold, xmax = 0, ymin = -Inf, ymax = Inf,
    fill = "#BFDDF5", alpha = 0.52
  ) +
  annotate(  # weak-source
    "rect", xmin = 0, xmax = source_threshold, ymin = -Inf, ymax = Inf,
    fill = "#F5C0BD", alpha = 0.48
  ) +
  annotate(  # strong-source
    "rect", xmin = source_threshold, xmax = Inf, ymin = -Inf, ymax = Inf,
    fill = "#B2182B", alpha = 0.28
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.5) +
  geom_vline(
    xintercept = c(-sink_threshold, source_threshold),
    linetype = "dotted", color = "grey45", linewidth = 0.4
  ) +
  # Zone labels — feasible here (single panel) where it wasn't in the faceted
  # comparison_figure. Placed near the top of the plot, one per zone.
  annotate(
    "text",
    x = zone_label_x,
    y = Inf,
    label = c("Strong-sink", "Weak-sink", "Weak-source", "Strong-source"),
    vjust = 1.3, size = 3.2, fontface = "bold", color = "grey15"
  ) +
  geom_density(alpha = 0.35, linewidth = 1.0, trim = FALSE) +
  scale_x_continuous(
    trans  = pseudo_log_trans(sigma = 0.01),
    breaks = c(-20, -5, -1, -0.1, 0, 0.1, 1, 10, 100),
    labels = c("-20", "-5", "-1", "-0.1", "0", "0.1", "1", "10", "100")
  ) +
  scale_fill_manual(
    values = density_group_colors,
    labels = setNames(density_n_labels$label, density_n_labels$group),
    name   = "Data source"
  ) +
  scale_color_manual(
    values = density_group_colors,
    labels = setNames(density_n_labels$label, density_n_labels$group),
    name   = "Data source"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    legend.text      = element_text(size = 10),
    legend.title     = element_text(size = 11),
    panel.grid.minor = element_blank()
  ) +
  labs(
   # title = expression(bold("Upland CH"[4]*" flux distributions: towers vs chambers")),
    x     = expression("Daily CH"[4] * " flux (mg C m"^-2 * " d"^-1 * "; pseudo-log scale)"),
    y     = "Density"
  )

ggsave(
  "FIGURES/NEON_FLUXNET_CH4_flux_density.png",
  plot   = flux_density_figure,
  width  = 12,
  height = 7,
  units  = "in",
  dpi    = 300
)

message("Wrote FIGURES/NEON_FLUXNET_CH4_flux_density.png")

# ── Two-panel combined figure (A = comparison, B = density) ──────────────────

combined_figure <- (
  comparison_figure +
    labs(tag = "A.") +
    theme(plot.tag = element_text(size = 14, face = "bold"))
) /
  (
    flux_density_figure +
      labs(tag = "B.") +
      theme(plot.tag = element_text(size = 14, face = "bold"))
  ) +
  plot_layout(heights = c(2, 1)) &
  theme(
    axis.text        = element_text(size = 10),
    axis.title       = element_text(size = 12),
    legend.text      = element_text(size = 10),
    legend.title     = element_text(size = 11),
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 10),
    plot.caption     = element_text(size = 10),
    strip.text       = element_text(size = 11)
  )

ggsave(
  "FIGURES/FIGURE3_NEON_FLUXNET_CH4_flux_combined.png",
  plot   = combined_figure,
  width  = 10,
  height = 12,
  units  = "in",
  dpi    = 300
)

message("Wrote FIGURES/NEON_FLUXNET_CH4_flux_combined.png")

# Detailed caption for the two-panel combined figure. Built from the same
# sink_threshold value used in both panels so the text always matches what's
# drawn (both panels share the same threshold and zone colors).
combined_figure_caption <- paste0(
  "Figure. ERA5-gapfilled NEON gradient CH4 fluxes compared with published benchmarks, ",
  "by ecosystem class (A) and for upland ecosystems only (B). ",
  "(A) Each panel is an ecosystem class; within a panel, rows are data sources: ",
  "FLUXNET-CH4 (Delwiche et al. 2021 published ecosystem-class values), Chambers ",
  "(soil chamber literature, including the SMUD compilation), Process model ",
  "(MeMo v1.0 process-based estimates), Upland towers (non-FLUXNET tower literature ",
  "and SMUD eddy-covariance towers; individual validation-tower site values are ",
  "additionally marked with black stars), and NEON (ERA5-gapfilled gradient/tower ",
  "CH4 budgets for individual NEON sites). Each row is drawn as a boxplot (median, ",
  "interquartile range, and whiskers to the data range); a row shows two boxes when ",
  "a source has both sink (flux < 0) and source (flux > 0) values within an ",
  "ecosystem class -- most visibly for NEON, whose sites split into weak-sink and ",
  "weak-source behavior within several ecosystem classes. Box color encodes flux ",
  "direction only (blue = sink, red = source), not data source. ",
  "(B) Kernel density distributions of daily CH4 flux for upland (non-inundated) ",
  "ecosystem classes only, pooling all towers (maroon; NEON ERA5-gapfilled sites, ",
  "validation towers, and SMUD eddy-covariance towers) against all chambers (navy; ",
  "soil chamber literature and SMUD chamber studies), to compare the overall shape ",
  "and central tendency of tower- vs. chamber-based flux estimates independent of ",
  "ecosystem-class binning. ",
  "In both panels, background shading divides the flux axis into four magnitude ",
  "zones -- strong-sink, weak-sink, weak-source, strong-source -- separated at a ",
  "fixed, shared threshold of ±", round(sink_threshold, 2), " mg C m-2 d-1. The same ",
  "threshold magnitude is used on both the sink and source sides so the shading is ",
  "symmetric around zero; because ",
  "source fluxes are typically an order of magnitude larger than sink fluxes (e.g., ",
  "wetland CH4 emission vs. upland CH4 uptake), most individual source-side ",
  "estimates -- including many NEON sites classified elsewhere in this analysis as ",
  "\"weak-source\" -- fall in the strong-source shaded zone; the shading should be ",
  "read as a magnitude reference scale rather than a formal sink/source ",
  "classification. Both panels use a pseudo-log x-axis transform (linear near zero, ",
  "logarithmic at larger magnitudes) to accommodate the wide range of flux ",
  "magnitudes on a single scale."
)
writeLines(
  strwrap(combined_figure_caption, width = 100),
  "FIGURES/FIGURE3_NEON_FLUXNET_CH4_flux_combined_caption.txt"
)
message("Wrote FIGURES/NEON_FLUXNET_CH4_flux_combined_caption.txt")
