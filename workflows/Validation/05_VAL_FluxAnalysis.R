# Validation tower CH4 flux products — 30-min, daily, and annual scales
#
# Adapted from flow.30min.analysis.R for SE-Sto, SE-Svb, and US-Uaf validation towers.
# Processes BOTH flux-gradient (FG_mean) and co-located eddy-covariance (EC_mean) fluxes
# through the same pipeline so their 30-min, daily, and annual estimates can be directly
# compared.
#
# Pipeline:
#   1. Load harmonised 30-min median ensemble data (SITEval_DATA_FILTERED_RSHPc_H)
#   2. Convert both FG and EC to mg C m⁻² 30min⁻¹
#   3. Standardised 30-min mean from balanced site × month × hour bins
#   4. Daily gap-fill lookup (site × month × hour median)
#   5. Annual budgets (FG and EC) by year
#   6. FG vs EC comparison at every scale
#   7. Outputs to OUTPUT/ (prefixed VAL_)
#
# Outputs written to localdir.val/OUTPUT/:
#   VAL_30min_ch4_model_data.csv             — all 30-min observations, both fluxes, covariates
#   VAL_30min_standardised_flux.csv          — balanced site-mean 30-min flux (FG and EC)
#   VAL_daily_flux_summary.csv               — daily mean FG and EC per site-year
#   VAL_annual_budgets_by_year.csv           — annual totals per site-year (FG and EC)
#   VAL_mean_annual_budget.csv               — mean ± SD annual budget (FG and EC)
#   VAL_FG_vs_EC_30min.csv                   — 30-min paired FG vs EC for scatterplots
#   VAL_FG_vs_EC_daily.csv                   — daily paired FG vs EC
#   VAL_FG_vs_EC_annual.csv                  — annual paired FG vs EC budget comparison
#   VAL_site_behavior.csv                    — source/sink classification (FG and EC)

library(tidyverse)
library(lubridate)

# ── Directories ───────────────────────────────────────────────────────────────
localdir.val <- Sys.getenv(
  "LOCALDIR_VAL",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/Validation_Sites"
)
DirRepo.ch4 <- Sys.getenv(
  "DIRREPO_CH4",
  unset = "/Users/sm3466/Library/CloudStorage/Dropbox-YSE/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-methane"
)

if (!dir.exists(localdir.val)) stop("Validation data directory not found: ", localdir.val)

setwd(localdir.val)
dir.create("OUTPUT",  showWarnings = FALSE, recursive = TRUE)
dir.create("FIGURES", showWarnings = FALSE, recursive = TRUE)

# ── Input: load or reconstruct SITEval_DATA_FILTERED_RSHPc_H ─────────────────
# Primary: SITEval_DATA_FILTERED_RSHP_EnSEMBLE.Rdata (output of flow.RSHP_VAL.R)
# Fallback: reconstruct from SITEval_DATA_FILTERED_CH4.Rdata + CCC_CH4.Rdata
#   using the same canopy filter (AA/AW only) + CCC > 0 + 30-min harmonisation
#   that flow.RSHP_VAL.R applies.

ensemble_file <- file.path(localdir.val, "SITEval_DATA_FILTERED_RSHP_EnSEMBLE.Rdata")

if (file.exists(ensemble_file)) {
  message("Loading ", ensemble_file)
  load(ensemble_file)   # provides SITEval_DATA_FILTERED_RSHPc_H
} else {
  message("SITEval_DATA_FILTERED_RSHP_EnSEMBLE.Rdata not found — reconstructing from raw files.")

  filtered_file <- file.path(localdir.val, "SITEval_DATA_FILTERED_CH4.Rdata")
  ccc_file      <- file.path(localdir.val, "CCC_CH4.Rdata")
  canopy_file   <- file.path(localdir.val, "Val_canopy.csv")

  missing <- c(filtered_file, ccc_file, canopy_file)[
    !file.exists(c(filtered_file, ccc_file, canopy_file))
  ]
  if (length(missing) > 0)
    stop("Cannot reconstruct SITEval_DATA_FILTERED_RSHPc_H — missing:\n",
         paste(missing, collapse = "\n"),
         "\nRun flow.filter.validation.R and flow.RSHP_VAL.R first.")

  load(filtered_file)   # → SITEval_DATA_FILTERED (list by site, 9-min, CH4 only)
  load(ccc_file)        # → CCC_VAL (CCC per site/approach/pair/season)

  canopy_info <- read.csv(canopy_file) %>%
    rename(Site = SITE_ID)

  # Apply canopy position filter (AA and AW only — same as flow.RSHP_VAL.R)
  SITEval_DATA_FILTEREDc <- purrr::imap(SITEval_DATA_FILTERED, function(df, site) {
    df %>%
      mutate(Site = site) %>%
      left_join(canopy_info, by = c("Site", "dLevelsAminusB")) %>%
      filter(Canopy_L1 %in% c("AA", "AW"))
  })

  # Join CCC and keep only pairs with CCC > 0
  SITEval_DATA_FILTEREDRSHPc <- purrr::imap(SITEval_DATA_FILTEREDc, function(df, site) {
    df %>%
      full_join(CCC_VAL, by = c("gas", "Approach", "dLevelsAminusB")) %>%
      filter(CCC > 0) %>%
      mutate(
        month  = format(timeEndA.local, "%m") %>% as.numeric(),
        Season = case_when(
          month %in% c(12, 1, 2) ~ "Winter",
          month %in% 3:5         ~ "Spring",
          month %in% 6:8         ~ "Summer",
          TRUE                   ~ "Autumn"
        ),
        hour = format(timeEndA.local, "%H")
      )
  })

  # Harmonise 9-min → 30-min medians (mirrors harmonize_val in flow.RSHP_VAL.R)
  harmonize_val_local <- function(tibble_list) {
    purrr::imap(tibble_list, function(df, site) {
      df %>%
        mutate(time.rounded = lubridate::round_date(timeEndA.local, unit = "30 minutes")) %>%
        reframe(
          .by = c(time.rounded, gas),
          FG_mean = median(FG_mean, na.rm = TRUE),
          EC_mean = mean(EC_mean,   na.rm = TRUE),
          Tair_C  = mean(Tair_C,    na.rm = TRUE),
          PAR     = mean(PAR,       na.rm = TRUE)
        ) %>%
        mutate(
          Month  = format(time.rounded, "%m") %>% as.numeric(),
          Season = case_when(
            Month %in% c(12, 1, 2) ~ "Winter",
            Month %in% 3:5         ~ "Spring",
            Month %in% 6:8         ~ "Summer",
            TRUE                   ~ "Autumn"
          )
        )
    })
  }

  SITEval_DATA_FILTERED_RSHPc_H <- harmonize_val_local(SITEval_DATA_FILTEREDRSHPc)
  message("Reconstruction complete. Sites: ",
          paste(names(SITEval_DATA_FILTERED_RSHPc_H), collapse = ", "))
}

# ── Constants ─────────────────────────────────────────────────────────────────
# FG_mean and EC_mean are in µmol CH4 m⁻² s⁻¹ (same units as NEON flux_total).
# Conversion to mg C m⁻² 30min⁻¹ follows NEON convention:
#   2 × 0.0000288872 × 1000  =  µmol/s × (1800 s/30min) × (12 µg C/µmol) / (1e6 µg/g) × (1e3 mg/g)
#
# If FG_mean is in different units at your sites, update flux_to_mgC_30min accordingly.
flux_to_mgC_30min <- 2 * 0.0000288872 * 1000   # mg C m⁻² 30min⁻¹ per (µmol m⁻² s⁻¹)
seconds_per_30min <- 30 * 60
ug_c_per_umol_c   <- 12.011
behavior_levels   <- c("Consistent sink", "Fluctuating", "Consistent source")
site_list         <- c("SE-Sto", "SE-Svb", "US-Uaf")

# Validation site metadata (coordinates from attr files)
site_metadata <- tibble(
  SITE_ID   = c("SE-Sto", "SE-Svb", "US-Uaf"),
  latitude  = c(64.18203, 64.25611, 64.8663),
  longitude = c(19.55654, 19.7745, -147.8555),
  EcoType   = c("Wetland", "Forest", "Wetland")
)

mg_c_30min_to_umol_c_s <- function(x) x * 1000 / ug_c_per_umol_c / seconds_per_30min

# ── Build 30-min data frame (FG + EC) ─────────────────────────────────────────
# Use CCC-filtered, harmonised 30-min ensemble (CCC > 0 per height pair).
# Filter to CH4 only; both FG and EC are kept as parallel columns.

ch4_30min <- purrr::imap_dfr(SITEval_DATA_FILTERED_RSHPc_H, function(df, site) {
  df %>%
    filter(gas == "CH4") %>%
    mutate(SITE_ID = site)
}) %>%
  mutate(
    time.rounded = as.POSIXct(time.rounded, tz = "UTC"),
    Date         = as.Date(time.rounded),
    Year         = as.integer(format(time.rounded, "%Y")),
    month        = as.integer(format(time.rounded, "%m")),
    doy          = as.integer(format(time.rounded, "%j")),
    hour_num     = as.numeric(format(time.rounded, "%H")) +
                   as.numeric(format(time.rounded, "%M")) / 60,
    hour_factor  = factor(as.integer(format(time.rounded, "%H")), levels = 0:23),
    sin_hour     = sin(2 * pi * hour_num / 24),
    cos_hour     = cos(2 * pi * hour_num / 24),
    season       = factor(coalesce(Season, "Summer"),
                          levels = c("Winter", "Spring", "Summer", "Autumn")),
    # Convert both fluxes to mg C m⁻² 30min⁻¹
    CH4_FG_mgC_30min = FG_mean * flux_to_mgC_30min,
    CH4_EC_mgC_30min = EC_mean * flux_to_mgC_30min,
    # log_PAR for GAM (matches NEON covariate)
    log_PAR      = log1p(pmax(PAR, 0)),
    VSWCMean     = NA_real_   # not available for validation sites; ERA5 VSWC used downstream
  ) %>%
  # Keep rows where at least FG is finite
  filter(is.finite(CH4_FG_mgC_30min))

# Quick sanity: show per-site counts
message("30-min observation counts per site:")
ch4_30min %>% count(SITE_ID) %>% print()

# Write the full 30-min model data (needed by Val.30min.Gapfill.R and Val.ERA5.HalfHourlyGapfill.R)
write.csv(ch4_30min, "OUTPUT/VAL_30min_ch4_model_data.csv", row.names = FALSE)

# ── 30-min standardised flux (balanced site × month × hour bin) ───────────────
# Mean FG and EC within each site × month × hour cell, computed only for cells
# that have ≥ 3 observations (matching NEON convention).

bin_min_obs <- 3L

site_month_hour_means <- ch4_30min %>%
  mutate(hour_int = as.integer(hour_num)) %>%
  group_by(SITE_ID, month, hour_int) %>%
  summarise(
    n_obs            = n(),
    FG_bin_mean      = mean(CH4_FG_mgC_30min, na.rm = TRUE),
    EC_bin_mean      = mean(CH4_EC_mgC_30min, na.rm = TRUE),
    Tair_bin_mean    = mean(Tair_C,           na.rm = TRUE),
    PAR_bin_mean     = mean(log_PAR,           na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_obs >= bin_min_obs)

# Standardised site-level 30-min flux: mean of all bin means (balanced)
site_30min_standardised <- site_month_hour_means %>%
  group_by(SITE_ID) %>%
  summarise(
    n_bins                      = n(),
    FG_30min_umolC_m2_s         = mg_c_30min_to_umol_c_s(mean(FG_bin_mean, na.rm = TRUE)),
    FG_30min_sd_umolC_m2_s      = mg_c_30min_to_umol_c_s(sd(FG_bin_mean,   na.rm = TRUE)),
    FG_30min_se_umolC_m2_s      = FG_30min_sd_umolC_m2_s / sqrt(n_bins),
    EC_30min_umolC_m2_s         = mg_c_30min_to_umol_c_s(mean(EC_bin_mean, na.rm = TRUE)),
    EC_30min_sd_umolC_m2_s      = mg_c_30min_to_umol_c_s(sd(EC_bin_mean,   na.rm = TRUE)),
    EC_30min_se_umolC_m2_s      = EC_30min_sd_umolC_m2_s / sqrt(n_bins),
    behavior_FG_30min           = if_else(mean(FG_bin_mean > 0) >= 0.75, "Consistent source",
                                   if_else(mean(FG_bin_mean > 0) <= 0.25, "Consistent sink",
                                           "Fluctuating")),
    behavior_EC_30min           = if_else(mean(EC_bin_mean > 0, na.rm = TRUE) >= 0.75, "Consistent source",
                                   if_else(mean(EC_bin_mean > 0, na.rm = TRUE) <= 0.25, "Consistent sink",
                                           "Fluctuating")),
    .groups = "drop"
  )

write.csv(site_30min_standardised, "OUTPUT/VAL_30min_standardised_flux.csv", row.names = FALSE)

# ── Daily gap-fill lookup ─────────────────────────────────────────────────────
# For each site-day: sum the 48 half-hourly values (using observed where available,
# filling missing slots from the site × month × hour lookup table).

# Lookup: site × month × hour medians
lookup_table <- ch4_30min %>%
  mutate(hour_int = as.integer(hour_num)) %>%
  group_by(SITE_ID, month, hour_int) %>%
  summarise(
    FG_lookup = median(CH4_FG_mgC_30min, na.rm = TRUE),
    EC_lookup = median(CH4_EC_mgC_30min, na.rm = TRUE),
    .groups = "drop"
  )

# Full grid of all half-hours in each observed day
all_days <- ch4_30min %>%
  distinct(SITE_ID, Date) %>%
  mutate(month = as.integer(format(Date, "%m")))

half_hour_grid <- all_days %>%
  mutate(
    hour_int = list(0:23),
    half_offset = list(c(0, 0.5))
  ) %>%
  unnest(hour_int) %>%
  unnest(half_offset) %>%
  mutate(
    hour_num_grid = hour_int + half_offset,
    time_rounded_grid = as.POSIXct(
      paste0(Date, " ", sprintf("%02d", hour_int), ":",
             if_else(half_offset == 0, "00", "30"), ":00"), tz = "UTC"
    )
  )

# Join observed half-hours, fill missing with lookup
daily_filled <- half_hour_grid %>%
  left_join(
    ch4_30min %>% dplyr::select(SITE_ID, time.rounded, CH4_FG_mgC_30min, CH4_EC_mgC_30min),
    by = c("SITE_ID", "time_rounded_grid" = "time.rounded")
  ) %>%
  left_join(lookup_table, by = c("SITE_ID", "month", "hour_int")) %>%
  mutate(
    FG_filled = coalesce(CH4_FG_mgC_30min, FG_lookup),
    EC_filled = coalesce(CH4_EC_mgC_30min, EC_lookup),
    FG_gapfilled = is.na(CH4_FG_mgC_30min),
    EC_gapfilled = is.na(CH4_EC_mgC_30min)
  ) %>%
  group_by(SITE_ID, Date) %>%
  summarise(
    n_slots          = n(),
    n_FG_observed    = sum(!FG_gapfilled, na.rm = TRUE),
    n_EC_observed    = sum(!EC_gapfilled, na.rm = TRUE),
    FG_coverage      = n_FG_observed / n_slots,
    EC_coverage      = n_EC_observed / n_slots,
    # Only compute daily sum for days with ≥ 25% coverage; rest NA
    FG_daily_mgC_m2  = if_else(FG_coverage >= 0.25, sum(FG_filled, na.rm = TRUE), NA_real_),
    EC_daily_mgC_m2  = if_else(EC_coverage >= 0.25, sum(EC_filled, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    Year  = as.integer(format(Date, "%Y")),
    month = as.integer(format(Date, "%m"))
  )

# Daily summary: site-level mean ± SD across all days
daily_summary <- daily_filled %>%
  group_by(SITE_ID) %>%
  summarise(
    n_days                  = sum(!is.na(FG_daily_mgC_m2)),
    FG_daily_mean_mgC_m2_day = mean(FG_daily_mgC_m2,  na.rm = TRUE),
    FG_daily_sd_mgC_m2_day   = sd(FG_daily_mgC_m2,    na.rm = TRUE),
    FG_daily_se_mgC_m2_day   = FG_daily_sd_mgC_m2_day / sqrt(n_days),
    EC_daily_mean_mgC_m2_day = mean(EC_daily_mgC_m2,  na.rm = TRUE),
    EC_daily_sd_mgC_m2_day   = sd(EC_daily_mgC_m2,    na.rm = TRUE),
    EC_daily_se_mgC_m2_day   = EC_daily_sd_mgC_m2_day / sqrt(n_days),
    behavior_FG_daily        = if_else(mean(FG_daily_mgC_m2 > 0, na.rm = TRUE) >= 0.75,
                                       "Consistent source",
                                 if_else(mean(FG_daily_mgC_m2 > 0, na.rm = TRUE) <= 0.25,
                                         "Consistent sink", "Fluctuating")),
    behavior_EC_daily        = if_else(mean(EC_daily_mgC_m2 > 0, na.rm = TRUE) >= 0.75,
                                       "Consistent source",
                                 if_else(mean(EC_daily_mgC_m2 > 0, na.rm = TRUE) <= 0.25,
                                         "Consistent sink", "Fluctuating")),
    .groups = "drop"
  )

write.csv(daily_filled,   "OUTPUT/VAL_daily_flux_all_days.csv",  row.names = FALSE)
write.csv(daily_summary,  "OUTPUT/VAL_daily_flux_summary.csv",   row.names = FALSE)

# ── Annual budgets by year (scaled from daily sums) ───────────────────────────
# Annual total = sum of filled daily values × (365 / n_days_with_data)
# This is the "scaled" approach used in the NEON pipeline.

annual_by_year <- daily_filled %>%
  filter(!is.na(FG_daily_mgC_m2) | !is.na(EC_daily_mgC_m2)) %>%
  group_by(SITE_ID, Year) %>%
  summarise(
    n_days_FG                    = sum(!is.na(FG_daily_mgC_m2)),
    n_days_EC                    = sum(!is.na(EC_daily_mgC_m2)),
    FG_annual_gC_m2_yr           = sum(FG_daily_mgC_m2, na.rm = TRUE) /
                                     max(n_days_FG, 1) * 365 / 1000,
    EC_annual_gC_m2_yr           = sum(EC_daily_mgC_m2, na.rm = TRUE) /
                                     max(n_days_EC, 1) * 365 / 1000,
    FG_annual_obs_coverage       = n_days_FG / 365,
    EC_annual_obs_coverage       = n_days_EC / 365,
    .groups = "drop"
  ) %>%
  # Only retain site-years with ≥ 90 days of data
  filter(n_days_FG >= 90 | n_days_EC >= 90) %>%
  mutate(
    sign_agree_FG_EC = sign(FG_annual_gC_m2_yr) == sign(EC_annual_gC_m2_yr)
  )

# Mean annual budget across available years
mean_annual_budget <- annual_by_year %>%
  group_by(SITE_ID) %>%
  summarise(
    n_years                       = n(),
    FG_mean_annual_gC_m2_yr       = mean(FG_annual_gC_m2_yr, na.rm = TRUE),
    FG_sd_annual_gC_m2_yr         = sd(FG_annual_gC_m2_yr,   na.rm = TRUE),
    FG_se_annual_gC_m2_yr         = FG_sd_annual_gC_m2_yr / sqrt(n_years),
    EC_mean_annual_gC_m2_yr       = mean(EC_annual_gC_m2_yr, na.rm = TRUE),
    EC_sd_annual_gC_m2_yr         = sd(EC_annual_gC_m2_yr,   na.rm = TRUE),
    EC_se_annual_gC_m2_yr         = EC_sd_annual_gC_m2_yr / sqrt(n_years),
    annual_difference_FG_minus_EC = FG_mean_annual_gC_m2_yr - EC_mean_annual_gC_m2_yr,
    sign_agree_annual             = sign(FG_mean_annual_gC_m2_yr) == sign(EC_mean_annual_gC_m2_yr),
    .groups = "drop"
  )

write.csv(annual_by_year,    "OUTPUT/VAL_annual_budgets_by_year.csv", row.names = FALSE)
write.csv(mean_annual_budget, "OUTPUT/VAL_mean_annual_budget.csv",    row.names = FALSE)

# ── Source/sink behavior classification ──────────────────────────────────────
# Using the 30-min observed fluxes (consistent with NEON.30min.Gapfill.R approach)

site_behavior <- ch4_30min %>%
  group_by(SITE_ID, month) %>%
  summarise(
    FG_monthly_source_frac = mean(CH4_FG_mgC_30min > 0, na.rm = TRUE),
    EC_monthly_source_frac = mean(CH4_EC_mgC_30min > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(SITE_ID) %>%
  summarise(
    n_months_FG_source  = sum(FG_monthly_source_frac >= 0.5, na.rm = TRUE),
    n_months_EC_source  = sum(EC_monthly_source_frac >= 0.5, na.rm = TRUE),
    FG_annual_behavior  = case_when(
      n_months_FG_source >= 9  ~ "Consistent source",
      n_months_FG_source <= 3  ~ "Consistent sink",
      TRUE ~ "Fluctuating"
    ),
    EC_annual_behavior  = case_when(
      n_months_EC_source >= 9  ~ "Consistent source",
      n_months_EC_source <= 3  ~ "Consistent sink",
      TRUE ~ "Fluctuating"
    ),
    behavior_agrees = FG_annual_behavior == EC_annual_behavior,
    .groups = "drop"
  ) %>%
  left_join(site_metadata, by = "SITE_ID") %>%
  left_join(mean_annual_budget %>%
              dplyr::select(SITE_ID, FG_mean_annual_gC_m2_yr, EC_mean_annual_gC_m2_yr,
                            annual_difference_FG_minus_EC, sign_agree_annual),
            by = "SITE_ID")

write.csv(site_behavior, "OUTPUT/VAL_site_behavior.csv", row.names = FALSE)

# ── FG vs EC comparison tables ────────────────────────────────────────────────

# 30-min: paired observations (both FG and EC finite)
fg_ec_30min <- ch4_30min %>%
  filter(is.finite(CH4_FG_mgC_30min), is.finite(CH4_EC_mgC_30min)) %>%
  dplyr::select(SITE_ID, time.rounded, Year, month, season, hour_num, doy,
                CH4_FG_mgC_30min, CH4_EC_mgC_30min, Tair_C, log_PAR) %>%
  mutate(
    FG_minus_EC = CH4_FG_mgC_30min - CH4_EC_mgC_30min,
    sign_agree  = sign(CH4_FG_mgC_30min) == sign(CH4_EC_mgC_30min)
  )

# 30-min agreement stats by site
fg_ec_30min_stats <- fg_ec_30min %>%
  group_by(SITE_ID) %>%
  summarise(
    n_pairs          = n(),
    r_pearson        = suppressWarnings(cor(CH4_FG_mgC_30min, CH4_EC_mgC_30min,
                                           use = "complete.obs")),
    rmse             = sqrt(mean(FG_minus_EC^2, na.rm = TRUE)),
    mae              = mean(abs(FG_minus_EC),    na.rm = TRUE),
    mean_bias        = mean(FG_minus_EC,         na.rm = TRUE),
    sign_accuracy    = mean(sign_agree,           na.rm = TRUE),
    .groups = "drop"
  )

# Daily: paired where both FG and EC daily values are available
fg_ec_daily <- daily_filled %>%
  filter(!is.na(FG_daily_mgC_m2), !is.na(EC_daily_mgC_m2)) %>%
  mutate(
    FG_minus_EC = FG_daily_mgC_m2 - EC_daily_mgC_m2,
    sign_agree  = sign(FG_daily_mgC_m2) == sign(EC_daily_mgC_m2)
  )

fg_ec_daily_stats <- fg_ec_daily %>%
  group_by(SITE_ID) %>%
  summarise(
    n_days       = n(),
    r_pearson    = suppressWarnings(cor(FG_daily_mgC_m2, EC_daily_mgC_m2, use = "complete.obs")),
    rmse         = sqrt(mean(FG_minus_EC^2, na.rm = TRUE)),
    mae          = mean(abs(FG_minus_EC),    na.rm = TRUE),
    mean_bias    = mean(FG_minus_EC,         na.rm = TRUE),
    sign_accuracy = mean(sign_agree,         na.rm = TRUE),
    .groups = "drop"
  )

# Annual
fg_ec_annual <- annual_by_year %>%
  filter(!is.na(FG_annual_gC_m2_yr), !is.na(EC_annual_gC_m2_yr)) %>%
  mutate(FG_minus_EC = FG_annual_gC_m2_yr - EC_annual_gC_m2_yr)

write.csv(fg_ec_30min,        "OUTPUT/VAL_FG_vs_EC_30min.csv",        row.names = FALSE)
write.csv(fg_ec_30min_stats,  "OUTPUT/VAL_FG_vs_EC_30min_stats.csv",  row.names = FALSE)
write.csv(fg_ec_daily,        "OUTPUT/VAL_FG_vs_EC_daily.csv",        row.names = FALSE)
write.csv(fg_ec_daily_stats,  "OUTPUT/VAL_FG_vs_EC_daily_stats.csv",  row.names = FALSE)
write.csv(fg_ec_annual,       "OUTPUT/VAL_FG_vs_EC_annual.csv",       row.names = FALSE)

# ── Summary report ────────────────────────────────────────────────────────────
writeLines(
  c(
    "# Validation Tower 30-min Analysis — FG vs EC",
    paste0("Generated: ", Sys.time()),
    "",
    "## 30-min FG vs EC agreement by site",
    capture.output(print(fg_ec_30min_stats)),
    "",
    "## Daily FG vs EC agreement by site",
    capture.output(print(fg_ec_daily_stats)),
    "",
    "## Annual budgets (scaled from daily, g C m-2 yr-1)",
    capture.output(print(mean_annual_budget %>%
                           dplyr::select(SITE_ID, FG_mean_annual_gC_m2_yr, EC_mean_annual_gC_m2_yr,
                                         annual_difference_FG_minus_EC, sign_agree_annual))),
    "",
    "## Site behavior classification (FG vs EC)",
    capture.output(print(site_behavior %>%
                           dplyr::select(SITE_ID, FG_annual_behavior, EC_annual_behavior,
                                         behavior_agrees)))
  ),
  "OUTPUT/VAL_30min_analysis_summary.txt"
)

message("flow.30min.analysis.VAL.R complete. Outputs in ", file.path(localdir.val, "OUTPUT"))
