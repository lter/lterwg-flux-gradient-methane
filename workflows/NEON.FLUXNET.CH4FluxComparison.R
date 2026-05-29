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

localdir.ch4 <- Sys.getenv(
  "LOCALDIR_CH4",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/Methane"
)
if (dir.exists(localdir.ch4)) {
  setwd(localdir.ch4)
}

dir.create("OUTPUT", showWarnings = FALSE)
dir.create("FIGURES", showWarnings = FALSE)

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

behavior_levels <- c("Consistent sink", "Fluctuating", "Consistent source")
behavior_colors <- c(
  "Consistent sink" = "#2166AC",
  "Fluctuating" = "grey35",
  "Consistent source" = "#B2182B"
)

base_plot_size <- 18
panel_title_size <- 21
axis_title_size <- 19
axis_text_size <- 15
legend_text_size <- 17
legend_title_size <- 18

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

assign_reference_ecosystem <- function(ecosystem_class) {
  case_when(
    str_detect(ecosystem_class, regex("marsh|fen|swamp|bog|wetland|mangrove", ignore_case = TRUE)) ~ "Wetland",
    str_detect(ecosystem_class, regex("lake", ignore_case = TRUE)) ~ "Lake",
    str_detect(ecosystem_class, regex("rice|cropland", ignore_case = TRUE)) ~ "Cropland",
    str_detect(ecosystem_class, regex("^upland$|drained|desert|rock|ice", ignore_case = TRUE)) ~ "Upland/Desert",
    str_detect(ecosystem_class, regex("forest", ignore_case = TRUE)) ~ "Forest",
    str_detect(ecosystem_class, regex("shrubland", ignore_case = TRUE)) ~ "Shrubland",
    str_detect(ecosystem_class, regex("grassland|steppe|savanna", ignore_case = TRUE)) ~ "Grassland/Savanna",
    str_detect(ecosystem_class, regex("tundra", ignore_case = TRUE)) ~ "Tundra",
    str_detect(ecosystem_class, regex("urban", ignore_case = TRUE)) ~ "Urban",
    TRUE ~ "Other"
  )
}

ecosystem_levels <- c(
  "Urban", "Lake", "Wetland", "Cropland", "Forest", "Shrubland",
  "Grassland/Savanna", "Tundra", "Upland/Desert"
)

site_behavior <- read.csv(site_behavior_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    CH4_behavior = factor(CH4_behavior, levels = behavior_levels),
    CH4_gradient_behavior = factor(CH4_gradient_behavior, levels = behavior_levels)
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
    site_behavior %>% dplyr::select(SITE_ID, CH4_behavior, CH4_gradient_behavior, EcoType),
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

assign_neon_ecosystem <- function(ecotype) {
  case_when(
    ecotype == "Grassland" ~ "Grassland/Savanna",
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
  group_by(SITE_ID, gradient_behavior_for_plot, CH4_gradient_behavior, CH4_behavior, EcoType) %>%
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
    ecosystem_group = assign_neon_ecosystem(EcoType),
    era5_gradient_class = recode(
      as.character(gradient_behavior_for_plot),
      "Consistent sink" = "NEON ERA5 gapfilled sink",
      "Fluctuating" = "NEON ERA5 gapfilled fluctuating",
      "Consistent source" = "NEON ERA5 gapfilled source"
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

comparison_table <- bind_rows(
  fluxnet_reference %>%
    mutate(CH4_behavior = NA_character_) %>%
    dplyr::select(
      comparison_group, ecosystem_class, CH4_behavior, annual_gC_m2_yr,
      annual_sd_gC_m2_yr, daily_mgC_m2_day, daily_low_mgC_m2_day,
      daily_high_mgC_m2_day, n_sites, n_site_years, reference_id, source_type
    ),
  soil_chamber_reference %>%
    mutate(CH4_behavior = NA_character_) %>%
    dplyr::select(
      comparison_group, ecosystem_class, CH4_behavior, annual_gC_m2_yr,
      annual_sd_gC_m2_yr, daily_mgC_m2_day, daily_low_mgC_m2_day,
      daily_high_mgC_m2_day, n_sites, n_site_years, reference_id, source_type
    ),
  process_model_uptake %>%
    mutate(CH4_behavior = NA_character_) %>%
    dplyr::select(
      comparison_group, ecosystem_class, CH4_behavior, annual_gC_m2_yr,
      annual_sd_gC_m2_yr, daily_mgC_m2_day, daily_low_mgC_m2_day,
      daily_high_mgC_m2_day, n_sites, n_site_years, reference_id, source_type
    ),
  era5_class_summary %>%
    dplyr::select(
      comparison_group, ecosystem_class, CH4_behavior, annual_gC_m2_yr,
      annual_sd_gC_m2_yr, daily_mgC_m2_day, daily_low_mgC_m2_day,
      daily_high_mgC_m2_day, n_sites, n_site_years, reference_id, source_type
    )
) %>%
  arrange(desc(daily_mgC_m2_day)) %>%
  mutate(
    source_type = factor(
      source_type,
      levels = c(
        "Published ecosystem class",
        "Soil chamber literature",
        "Process-based model",
        "NEON ERA5 gapfilled gradient"
      )
    ),
    ecosystem_group = case_when(
      source_type == "NEON ERA5 gapfilled gradient" ~ "Across NEON sites",
      TRUE ~ assign_reference_ecosystem(ecosystem_class)
    )
  )

write.csv(comparison_table, "OUTPUT/CH4_flux_FLUXNET_NEON_comparison_values.csv", row.names = FALSE)
write.csv(era5_annual_budget, "OUTPUT/NEON_ERA5_annual_gradient_flux_rates_for_FLUXNET_comparison.csv", row.names = FALSE)
write.csv(era5_site_summary, "OUTPUT/NEON_ERA5_site_median_daily_gradient_flux_classes.csv", row.names = FALSE)
write.csv(process_model_uptake, "OUTPUT/process_model_upland_CH4_uptake_values.csv", row.names = FALSE)
write.csv(soil_chamber_reference, "OUTPUT/soil_chamber_CH4_flux_reference_values.csv", row.names = FALSE)

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
  "- OUTPUT/soil_chamber_CH4_flux_reference_values.csv"
)
writeLines(reference_lines, "OUTPUT/CH4_flux_FLUXNET_NEON_comparison_references.md")

source_offsets <- c(
  "FLUXNET-CH4" = 0,
  "Soil chambers" = 0,
  "Process model" = 0,
  "NEON ERA5" = 0
)

reference_plot_data <- comparison_table %>%
  filter(source_type != "NEON ERA5 gapfilled gradient") %>%
  mutate(
    source_label = recode(
      as.character(source_type),
      "Published ecosystem class" = "FLUXNET-CH4",
      "Soil chamber literature" = "Soil chambers",
      "Process-based model" = "Process model"
    ),
    source_label = factor(source_label, levels = c("FLUXNET-CH4", "Soil chambers", "Process model", "NEON ERA5")),
    ecosystem_group = factor(ecosystem_group, levels = ecosystem_levels),
    y_base = length(ecosystem_levels) - as.integer(ecosystem_group) + 1,
    y_plot = y_base + source_offsets[source_label]
  ) %>%
  filter(!is.na(y_plot), is.finite(daily_mgC_m2_day)) %>%
  reframe(
    .by = c(ecosystem_group, source_label, y_base, y_plot),
    n_estimates = dplyr::n(),
    daily_mgC_m2_day = median(daily_mgC_m2_day, na.rm = TRUE),
    daily_low_mgC_m2_day = min(daily_low_mgC_m2_day, na.rm = TRUE),
    daily_high_mgC_m2_day = max(daily_high_mgC_m2_day, na.rm = TRUE)
  )

neon_plot_data <- era5_site_summary %>%
  mutate(
    source_label = factor("NEON ERA5", levels = c("FLUXNET-CH4", "Soil chambers", "Process model", "NEON ERA5")),
    ecosystem_group = factor(ecosystem_group, levels = ecosystem_levels),
    y_base = length(ecosystem_levels) - as.integer(ecosystem_group) + 1,
    y_plot = y_base + source_offsets[source_label]
  ) %>%
  filter(!is.na(y_plot), is.finite(daily_mgC_m2_day)) %>%
  reframe(
    .by = c(ecosystem_group, source_label, y_base, y_plot),
    n_sites = dplyr::n(),
    prop_source_sites = mean(daily_mgC_m2_day > 0, na.rm = TRUE),
    daily_low_mgC_m2_day = quantile(daily_mgC_m2_day, 0.25, na.rm = TRUE),
    daily_high_mgC_m2_day = quantile(daily_mgC_m2_day, 0.75, na.rm = TRUE),
    daily_mgC_m2_day = median(daily_mgC_m2_day, na.rm = TRUE)
  ) %>%
  mutate(
    grouped_behavior = case_when(
      prop_source_sites >= 0.75 ~ "Consistent source",
      prop_source_sites <= 0.25 ~ "Consistent sink",
      TRUE ~ "Fluctuating"
    ),
    grouped_behavior = factor(grouped_behavior, levels = behavior_levels)
  )

neon_site_points <- era5_site_summary %>%
  mutate(
    ecosystem_group = factor(ecosystem_group, levels = ecosystem_levels),
    y_base = length(ecosystem_levels) - as.integer(ecosystem_group) + 1,
    y_plot = y_base,
    grouped_behavior = factor(gradient_behavior_for_plot, levels = behavior_levels)
  ) %>%
  filter(!is.na(y_plot), is.finite(daily_mgC_m2_day), !is.na(grouped_behavior))

ecosystem_axis <- tibble(
  ecosystem_group = factor(ecosystem_levels, levels = ecosystem_levels),
  y_base = length(ecosystem_levels) - seq_along(ecosystem_levels) + 1
)

source_shapes <- c(
  "FLUXNET-CH4" = 21,
  "Soil chambers" = 24,
  "Process model" = 22,
  "NEON ERA5" = 21
)
source_fills <- c(
  "FLUXNET-CH4" = NA,
  "Soil chambers" = "#C49A6C",
  "Process model" = "#CBC9E2",
  "NEON ERA5" = "#009E73"
)
source_outline <- c(
  "FLUXNET-CH4" = "grey15",
  "Soil chambers" = "#C49A6C",
  "Process model" = "#CBC9E2",
  "NEON ERA5" = "#009E73"
)
source_legend_fills <- source_fills
source_legend_fills["FLUXNET-CH4"] <- "white"
source_legend_fills["NEON ERA5"] <- "black"
source_legend_outline <- source_outline
source_legend_outline["NEON ERA5"] <- "black"

source_legend_data <- tibble(
  source_label = factor(names(source_shapes), levels = names(source_shapes)),
  daily_mgC_m2_day = 0,
  y_plot = 0
)

axis_text_size <-16

comparison_figure <- ggplot() +
  annotate(
    "rect",
    xmin = -20,
    xmax = 0,
    ymin = -Inf,
    ymax = Inf,
    fill = "#DCEEFF",
    alpha = 0.28
  ) +
  annotate(
    "text",
    x = -2.7,
    y = max(ecosystem_axis$y_base) + 0.48,
    label = "Upland Sink",
    color = "#2166AC",
    fontface = "bold",
    size = 6,
    hjust = 0
  ) +
  annotate(
    "rect",
    xmin = 0,
    xmax = 10,
    ymin = -Inf,
    ymax = Inf,
    fill = "#FDE0DD",
    alpha = 0.24
  ) +
  annotate(
    "text",
    x = 0.18,
    y = max(ecosystem_axis$y_base) + 0.48,
    label = "Upland Source",
    color = "#B2182B",
    fontface = "bold",
    size = 6,
    hjust = 0
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.5) +
  geom_hline(
    data = ecosystem_axis,
    aes(yintercept = y_base),
    color = "grey92",
    linewidth = 0.35
  ) +
  geom_errorbar(
    data = reference_plot_data %>% filter(source_label == "FLUXNET-CH4"),
    aes(
      x = daily_mgC_m2_day,
      y = y_plot,
      xmin = daily_low_mgC_m2_day,
      xmax = daily_high_mgC_m2_day
    ),
    orientation = "y",
    width = 0.07,
    linewidth = 0.85,
    color = source_outline["FLUXNET-CH4"],
    alpha = 0.78
  ) +
  geom_errorbar(
    data = reference_plot_data %>% filter(source_label == "Soil chambers"),
    aes(
      x = daily_mgC_m2_day,
      y = y_plot,
      xmin = daily_low_mgC_m2_day,
      xmax = daily_high_mgC_m2_day
    ),
    orientation = "y",
    width = 0.07,
    linewidth = 1.25,
    color = source_outline["Soil chambers"],
    alpha = 0.78
  ) +
  geom_errorbar(
    data = reference_plot_data %>% filter(source_label == "Process model"),
    aes(
      x = daily_mgC_m2_day,
      y = y_plot,
      xmin = daily_low_mgC_m2_day,
      xmax = daily_high_mgC_m2_day
    ),
    orientation = "y",
    width = 0.07,
    linewidth = 1.25,
    color = source_outline["Process model"],
    alpha = 0.78
  ) +
  geom_errorbar(
    data = neon_plot_data,
    aes(
      x = daily_mgC_m2_day,
      y = y_plot,
      xmin = daily_low_mgC_m2_day,
      xmax = daily_high_mgC_m2_day,
      color = grouped_behavior
    ),
    orientation = "y",
    width = 0.07,
    linewidth = 1.35,
    alpha = 0.85
  ) +
  geom_point(
    data = source_legend_data,
    aes(x = daily_mgC_m2_day, y = y_plot, shape = source_label),
    alpha = 0,
    size = 0,
    show.legend = TRUE
  ) +
  geom_point(
    data = reference_plot_data %>% filter(source_label == "FLUXNET-CH4"),
    aes(x = daily_mgC_m2_day, y = y_plot),
    shape = source_shapes["FLUXNET-CH4"],
    fill = source_fills["FLUXNET-CH4"],
    color = source_outline["FLUXNET-CH4"],
    size = 8.6,
    stroke = 1.0
  ) +
  geom_point(
    data = reference_plot_data %>% filter(source_label == "Soil chambers"),
    aes(x = daily_mgC_m2_day, y = y_plot),
    shape = source_shapes["Soil chambers"],
    fill = source_fills["Soil chambers"],
    color = source_outline["Soil chambers"],
    size = 8.6,
    stroke = 1.0
  ) +
  geom_point(
    data = reference_plot_data %>% filter(source_label == "Process model"),
    aes(x = daily_mgC_m2_day, y = y_plot),
    shape = source_shapes["Process model"],
    fill = source_fills["Process model"],
    color = source_outline["Process model"],
    size = 8.6,
    stroke = 1.0
  ) +
  geom_point(
    data = neon_site_points,
    aes(x = daily_mgC_m2_day, y = y_plot, color = grouped_behavior, fill = grouped_behavior),
    shape = source_shapes["NEON ERA5"],
    alpha = 0.3,
    size = 3.0,
    stroke = 0.8,
    position = position_jitter(height = 0.11, width = 0, seed = 20260525),
    show.legend = FALSE
  ) +
  geom_point(
    data = neon_plot_data,
    aes(x = daily_mgC_m2_day, y = y_plot, color = grouped_behavior, fill = grouped_behavior),
    shape = source_shapes["NEON ERA5"],
    alpha = 0.9,
    size = 8.6,
    stroke = 1.0
  ) +
  geom_point(
    data = tibble(
      daily_mgC_m2_day = 0,
      y_plot = 0,
      grouped_behavior = factor(behavior_levels, levels = behavior_levels)
    ),
    aes(x = daily_mgC_m2_day, y = y_plot, color = grouped_behavior),
    alpha = 0,
    size = 0,
    show.legend = TRUE
  ) +
  scale_shape_manual(values = source_shapes, breaks = names(source_shapes), name = "Data source") +
  scale_fill_manual(values = behavior_colors, breaks = behavior_levels, guide = "none", drop = FALSE, na.translate = FALSE) +
  scale_color_manual(values = behavior_colors, breaks = behavior_levels, name = "NEON class", drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 0.01),
    breaks = c(-20, -5, -1, -0.1, 0, 0.1, 1, 10, 100),
    labels = c("-20", "-5", "-1", "-0.1", "0", "0.1", "1", "10", "100")
  ) +
  scale_y_continuous(
    breaks = ecosystem_axis$y_base,
    labels = ecosystem_axis$ecosystem_group,
    expand = expansion(mult = c(0.05, 0.06))
  ) +
  labs(
    title = "ERA5-gapfilled NEON gradient CH4 fluxes compared with published benchmarks",
    subtitle = paste(
      "Rows align comparable ecosystem types; symbols distinguish FLUXNET-CH4, process models, soil chambers, and NEON ERA5.",
      "\nEach ecosystem/source pair is consolidated into one estimate; error bars show range or interquartile range. Negative values indicate CH4 uptake."
    ),
    x = expression("Daily CH"[4] * " flux (mg C m"^-2 * " d"^-1 * "; pseudo-log scale)"),
    y = NULL,
    caption = paste(
      "NEON values summarize site-level medians across ERA5-gapfilled annual gradient budgets converted to daily rates.",
      coverage_note
    )
  ) +
  guides(
    shape = guide_legend(
      override.aes = list(
        shape = unname(source_shapes[names(source_shapes)]),
        size = 9.6,
        alpha = 1,
        fill = unname(source_legend_fills[names(source_shapes)]),
        colour = unname(source_legend_outline[names(source_shapes)]),
        stroke = 1.2
      ),
      order = 1,
      nrow = 1
    ),
    color = guide_legend(
      override.aes = list(shape = 16, size = 9, alpha = 1, color = unname(behavior_colors[behavior_levels])),
      order = 2,
      nrow = 1
    ),
    fill = "none"
  ) +
  theme_bw(base_size = base_plot_size) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 25),
    plot.subtitle = element_text(size = 16),
    plot.caption = element_text(size = 14, color = "grey25"),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    legend.title = element_text(size = legend_title_size),
    legend.text = element_text(size = legend_text_size),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "FIGURES/NEON_FLUXNET_CH4_flux_comparison.png",
  plot = comparison_figure,
  width = 18,
  height = 10.5,
  units = "in",
  dpi = 300
)

message("Wrote FIGURES/NEON_FLUXNET_CH4_flux_comparison.png")
