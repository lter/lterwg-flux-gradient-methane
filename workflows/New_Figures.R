# make_figures.R ─────────────────────────────────────────────────────────────
# All 8 NEON CH4 figures using REAL data only.
#
# Data sources
# ─────────────────────────────────────────────────────────────────────────────
# LOCAL CSVs (in this repo — always available):
#   workflows/check/OUTPUT/NEON_total_flux_sign_summary_by_site.csv
#   workflows/check/OUTPUT/NEON_total_flux_gapfill_sign_bias_diagnostic.csv
#
# MALONE LAB DRIVE (mount at LOCALDIR_CH4 below before running):
#   SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata  → Figs 1, 2, 7
#   Soildata_YearMon.Rdata                                  → Fig 7
#   Ameriflux_NEON field-sites.csv                          → Fig 8 coords
#
# Figures produced
# ─────────────────────────────────────────────────────────────────────────────
#   1  Diel CH4 by season and behavior class
#   2  GAM population-level driver effects
#   3  Phase-space: prop_positive vs mean flux
#   4  Total vs gradient flux (storage contribution)
#   5  Annual budgets by site-year (model-gapfilled, with 95% CI)
#   6  Observed-only vs model-gapfilled annual budgets
#   7  Cliff's δ: environmental driver differences between behavior classes
#   8  Site map with Alaska / Hawai'i / Puerto Rico insets
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(sf)
  library(scales)
  library(mgcv)
  library(rnaturalearth)
  library(rnaturalearthdata)
})

# ══════════════════════════════════════════════════════════════════════════════
# 0. PATHS
# ══════════════════════════════════════════════════════════════════════════════
# Set this to where the MaloneLab drive is mounted on your machine
LOCALDIR_CH4 <- "/Volumes/MaloneLab/Research/FluxGradient/Methane"
LOCALDIR     <- "/Volumes/MaloneLab/Research/FluxGradient"

# CSVs shipped in this repo (relative to project root)
REPO_ROOT <- here::here()
CSV_SUMMARY    <- file.path(REPO_ROOT,
  "workflows/check/OUTPUT/NEON_total_flux_sign_summary_by_site.csv")
CSV_DIAGNOSTIC <- file.path(REPO_ROOT,
  "workflows/check/OUTPUT/NEON_total_flux_gapfill_sign_bias_diagnostic.csv")

# ERA5 gapfill outputs (written by NEON.ERA5.HalfHourlyGapfill.R)
ERA5_COMPARISON <- file.path(LOCALDIR_CH4, "OUTPUT",
  "NEON_ERA5_vs_monthly_bin_budget_comparison.csv")
ERA5_ANNUAL     <- file.path(LOCALDIR_CH4, "OUTPUT",
  "NEON_ERA5_gapfilled_annual_budget_by_year.csv")

OUTDIR <- file.path(REPO_ROOT, "FIGURES_IMPROVED")
dir.create(OUTDIR, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 1. CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════
BEHAVIOR_COLORS <- c(
  "Weak-sink"   = "#2166AC",
  "Fluctuating"       = "#4D4D4D",
  "Weak-source" = "#B2182B"
)
SEASON_COLORS <- c(
  Winter = "#5E81AC", Spring = "#4C9F70",
  Summer = "#D99000", Autumn = "#A65628"
)
DRIVER_GROUP_COLORS <- c(
  "Moisture"              = "#1B9E77",
  "Soil texture/physical" = "#D95F02",
  "Soil chemistry/depth"  = "#7570B3",
  "Climate"               = "#E7298A",
  "Vegetation/canopy"     = "#66A61E"
)
QUALITY_COLORS <- c(
  "moderate observed coverage"     = "#1a9641",
  "low observed coverage"          = "#fdae61",
  "very low observed coverage"     = "#d7191c"
)
SOURCE_THRESHOLD <- 1

classify_behavior <- function(p) {
  case_when(
    p >= SOURCE_THRESHOLD ~ "Weak-source",
    p <= (1 - SOURCE_THRESHOLD) ~ "Weak-sink",
    !is.na(p) ~ "Fluctuating",
    TRUE ~ NA_character_
  )
}

# Helper: rename first matching candidate column to new_name (no-op if none found)
safe_rename <- function(df, new_name, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0 || hit[1] == new_name) return(df)
  rename(df, !!sym(new_name) := !!sym(hit[1]))
}

theme_neon <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey92", colour = NA),
      strip.text        = element_text(face = "bold"),
      legend.background = element_rect(fill = alpha("white", 0.92),
                                       colour = "grey75", linewidth = 0.3),
      legend.key        = element_blank(),
      plot.title        = element_text(face = "bold", size = 11),
      plot.subtitle     = element_text(size = 8.5, colour = "grey40"),
      axis.title        = element_text(size = 10)
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. LOAD SITE SUMMARY (CSVs — always available)
# ══════════════════════════════════════════════════════════════════════════════
site_summary <- read_csv(CSV_SUMMARY, na = "NA", show_col_types = FALSE) %>%
  rename(
    prop_pos     = prop_positive_total_flux,
    mean_total   = mean_total_flux,
    median_total = median_total_flux,
    median_grad  = median_gradient_flux,
    median_stor  = median_storage_flux
  ) %>%
  mutate(
    behavior      = classify_behavior(prop_pos),
    changed_class = !is.na(mean_total) & !is.na(median_grad) &
                    (sign(mean_total) != sign(median_grad))
  )

annual_diag <- read_csv(CSV_DIAGNOSTIC, na = "NA", show_col_types = FALSE) %>%
  left_join(select(site_summary, SITE_ID, prop_pos, behavior), by = "SITE_ID") %>%
  mutate(
    quality_flag = factor(quality_flag,
      levels = c("moderate observed coverage",
                 "low observed coverage",
                 "very low observed coverage"))
  )

cat(sprintf("Site summary: %d sites\n", nrow(site_summary)))
cat(sprintf("Annual diagnostic: %d site-years across %d sites\n",
            nrow(annual_diag), n_distinct(annual_diag$SITE_ID)))
cat(sprintf("Behaviors present: %s\n",
            paste(na.omit(unique(site_summary$behavior)), collapse = ", ")))

# ERA5 outputs (optional — used in Figs 5 & 6)
has_era5 <- file.exists(ERA5_COMPARISON) && file.exists(ERA5_ANNUAL)
if (has_era5) {
  era5_comparison <- read_csv(ERA5_COMPARISON, show_col_types = FALSE) %>%
    left_join(select(site_summary, SITE_ID, behavior), by = "SITE_ID") %>%
    rename(
      era5_budget   = mean_era5_gapfilled_annual_budget_gC_m2_yr,
      era5_sd       = sd_era5_gapfilled_annual_budget_gC_m2_yr,
      obs_budget    = monthly_bin_annual_budget_gC_m2_yr,
      obs_lwr       = monthly_bin_lwr_gC_m2_yr,
      obs_upr       = monthly_bin_upr_gC_m2_yr,
      prop_source   = monthly_bin_prop_source_months,
      abs_diff      = absolute_difference_gC_m2_yr,
      budget_diff   = budget_difference_era5_minus_monthly_bin_gC_m2_yr
    ) %>%
    mutate(sign_agree = sign_agreement == "TRUE" | sign_agreement == TRUE)

  era5_annual <- read_csv(ERA5_ANNUAL, show_col_types = FALSE) %>%
    left_join(select(site_summary, SITE_ID, behavior), by = "SITE_ID") %>%
    rename(era5_budget_yr = annual_budget_gC_m2_yr) %>%
    mutate(Year = as.character(Year))

  cat(sprintf("ERA5 comparison: %d sites, annual by year: %d rows\n",
              nrow(era5_comparison), nrow(era5_annual)))
} else {
  message("ERA5 output not found — Figs 5 & 6 will use model-gapfilled only.")
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. LOAD 30-MIN DATA (MaloneLab drive — needed for Figs 1, 2, 7)
# ══════════════════════════════════════════════════════════════════════════════
RDATA_PATH <- file.path(LOCALDIR_CH4,
  "SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata")
SOILDATA_PATH <- file.path(LOCALDIR_CH4, "Soildata_YearMon.Rdata")
METADATA_PATH <- file.path(LOCALDIR, "Ameriflux_NEON field-sites.csv")

has_30min <- file.exists(RDATA_PATH)
has_soil  <- file.exists(SOILDATA_PATH)
has_meta  <- file.exists(METADATA_PATH)

if (has_30min) {
  .loaded <- load(RDATA_PATH)               # load() returns names of loaded objects
  # Score each object: data.frame → nrow; list of data.frames → total rows; else 0
  .scores <- sapply(.loaded, function(x) {
    obj <- get(x)
    if (is.data.frame(obj))                   nrow(obj)
    else if (is.list(obj) && !is.data.frame(obj))
      sum(sapply(obj, function(o) if (is.data.frame(o)) nrow(o) else 0L))
    else 0L
  })
  .best <- get(.loaded[which.max(.scores)])
  # If it's a named list of data frames, bind into one — use .id to recover SITE_ID
  if (is.list(.best) && !is.data.frame(.best)) {
    id_col <- if (!is.null(names(.best))) "SITE_ID" else NULL
    dat30  <- bind_rows(.best, .id = id_col)
  } else {
    dat30 <- .best
  }
  cat(sprintf("  Columns in loaded data: %s\n",
              paste(names(dat30)[seq_len(min(20, ncol(dat30)))], collapse = ", ")))
  # Standardise key column names
  dat30 <- dat30 %>%
    safe_rename("flux_total",    c("totalFlux",           "flux_total",    "Total_Flux")) %>%
    safe_rename("flux_gradient", c("FG_ENSEMBLE_RSHP",    "flux_gradient", "Gradient_Flux")) %>%
    safe_rename("storage_flux",  c("storage_flux_filled",  "storage_flux",  "Storage_Flux")) %>%
    safe_rename("datetime",      c("startDateTime",        "time.rounded",  "datetime",
                                   "time",                 "Time")) %>%
    safe_rename("Tair_f",        c("Tair_C",              "Tair_f",        "airTemp"))
  # Derive month/season/hour from datetime only if not already in the data
  if (!"month" %in% names(dat30) && "datetime" %in% names(dat30)) {
    dat30 <- dat30 %>%
      mutate(
        datetime = ymd_hms(datetime),
        hour     = hour(datetime) + minute(datetime)/60,
        month    = month(datetime),
        season   = case_when(
          month %in% c(12, 1, 2) ~ "Winter",
          month %in% c( 3, 4, 5) ~ "Spring",
          month %in% c( 6, 7, 8) ~ "Summer",
          month %in% c( 9,10,11) ~ "Autumn"
        )
      )
  }
  dat30 <- dat30 %>%
    mutate(
      hour   = as.numeric(hour),
      month  = as.integer(month),
      season = factor(as.character(season),
                      levels = c("Winter","Spring","Summer","Autumn"))
    ) %>%
    left_join(select(site_summary, SITE_ID, behavior), by = "SITE_ID") %>%
    filter(!is.na(behavior), !is.na(flux_total))
  cat(sprintf("30-min data: %d rows, %d sites\n",
              nrow(dat30), n_distinct(dat30$SITE_ID)))
} else {
  warning(sprintf("30-min Rdata not found at:\n  %s\n  Figs 1, 2, 7 will be skipped.",
                  RDATA_PATH))
}

if (has_soil) {
  .loaded2 <- load(SOILDATA_PATH)
  .scores2 <- sapply(.loaded2, function(x) {
    obj <- get(x)
    if (is.data.frame(obj))                   nrow(obj)
    else if (is.list(obj) && !is.data.frame(obj))
      sum(sapply(obj, function(o) if (is.data.frame(o)) nrow(o) else 0L))
    else 0L
  })
  .best2 <- get(.loaded2[which.max(.scores2)])
  if (is.list(.best2) && !is.data.frame(.best2)) {
    id_col2  <- if (!is.null(names(.best2))) "SITE_ID" else NULL
    soil_dat <- bind_rows(.best2, .id = id_col2)
  } else {
    soil_dat <- .best2
  }
  # Standardise SITE_ID column name in case it differs
  soil_dat <- soil_dat %>%
    safe_rename("SITE_ID", c("siteID", "site_id", "Site", "site", "SITE_ID"))
  cat(sprintf("Soil data: %d rows, columns: %s\n", nrow(soil_dat),
              paste(names(soil_dat)[seq_len(min(15, ncol(soil_dat)))], collapse = ", ")))
}

# Site coordinates from metadata or fallback lookup
if (has_meta) {
  .meta_raw <- read_csv(METADATA_PATH, show_col_types = FALSE) %>%
    rename_with(tolower) %>%
    rename_with(~ str_replace_all(., "[^a-z0-9]", "_"))

  cat(sprintf("Metadata columns: %s\n",
              paste(names(.meta_raw)[seq_len(min(20, ncol(.meta_raw)))], collapse = ", ")))

  # Find the right columns by name patterns rather than position
  .id_col  <- intersect(c("field_site_id","site_id","neon_site_id","siteID","site"),
                        names(.meta_raw))[1]
  .lat_col <- names(.meta_raw)[str_detect(names(.meta_raw), "^lat")][1]
  .lon_col <- names(.meta_raw)[str_detect(names(.meta_raw), "^lon|^long")][1]

  if (!is.na(.id_col) && !is.na(.lat_col) && !is.na(.lon_col)) {
    site_coords <- .meta_raw %>%
      select(SITE_ID = !!.id_col, lat = !!.lat_col, lon = !!.lon_col) %>%
      filter(!is.na(lat), !is.na(lon))
    cat(sprintf("Coords loaded from metadata: %d sites  (id=%s, lat=%s, lon=%s)\n",
                nrow(site_coords), .id_col, .lat_col, .lon_col))
  } else {
    message(sprintf(
      "Could not identify id/lat/lon in metadata (found: id=%s, lat=%s, lon=%s) — using hardcoded coords.",
      .id_col, .lat_col, .lon_col))
    has_meta <- FALSE   # fall through to hardcoded block below
  }
}

if (!has_meta) {
  # Hardcoded NEON site coordinates as fallback
  site_coords <- tribble(
    ~SITE_ID,  ~lat,     ~lon,
    "NOGP", 46.77,-100.91, "MOAB", 38.25,-109.39, "LENO", 31.85, -88.16,
    "OAES", 35.41, -99.06, "RMNP", 40.28,-105.55, "JORN", 32.59,-106.84,
    "STER", 40.46,-103.03, "KONZ", 39.10, -96.56, "SCBI", 38.89, -78.14,
    "GRSM", 35.69, -83.50, "TALL", 32.95, -87.39, "BART", 44.06, -71.29,
    "SRER", 31.91,-110.84, "ORNL", 35.96, -84.29, "MLBS", 37.38, -80.52,
    "YELL", 44.95,-110.54, "HEAL", 63.88,-149.21, "KONA", 39.11, -96.61,
    "SJER", 37.11,-119.73, "JERC", 31.19, -84.47, "WOOD", 47.13, -99.24,
    "UNDE", 46.23, -89.54, "DSNY", 28.13, -81.44, "BARR", 71.28,-156.62,
    "LAJA", 18.02, -67.08, "SERC", 38.89, -76.56, "CLBJ", 33.40, -97.57,
    "SOAP", 37.03,-119.26, "BLAN", 39.04, -78.07, "TREE", 45.49, -89.59,
    "UKFS", 39.04, -95.19, "CPER", 40.82,-104.75, "STEI", 45.51, -89.59,
    "WREF", 45.82,-121.95, "DCFS", 47.16, -99.11, "NIWO", 40.05,-105.58,
    "DEJU", 63.88,-145.75, "ONAQ", 40.18,-112.45, "BONA", 65.15,-147.50,
    "ABBY", 45.76,-122.33, "DELA", 32.54, -87.80, "GUAN", 18.11, -65.98,
    "HARV", 42.54, -72.17, "OSBS", 29.69, -81.99, "PUUM", 19.55,-155.32,
    "TEAK", 37.01,-119.01, "TOOL", 68.66,-149.37
  )
}

site_summary <- site_summary %>% left_join(site_coords, by = "SITE_ID")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Diel CH4 flux by season and behavior class
# Requires: 30-min Rdata (MaloneLab drive)
# ══════════════════════════════════════════════════════════════════════════════
if (!has_30min) {
  message("Skipping Fig 1 — 30-min data not found.")
} else {
  # Bin to integer hour then compute mean ± SE per behavior × season × hour
  diel_df <- dat30 %>%
    mutate(hour_bin = floor(hour)) %>%
    group_by(behavior, season, hour_bin) %>%
    summarise(
      n        = n(),
      flux_mean = mean(flux_total, na.rm = TRUE),
      flux_se   = sd(flux_total,   na.rm = TRUE) / sqrt(n()),
      .groups  = "drop"
    ) %>%
    filter(n >= 5) %>%
    mutate(behavior = factor(behavior, levels = names(BEHAVIOR_COLORS)))

  fig1 <- ggplot(diel_df, aes(hour_bin, flux_mean,
                               colour = behavior, fill = behavior)) +
    # night shading
    annotate("rect", xmin = 0,  xmax = 6,  ymin = -Inf, ymax = Inf,
             fill = "navy", alpha = 0.05) +
    annotate("rect", xmin = 20, xmax = 24, ymin = -Inf, ymax = Inf,
             fill = "navy", alpha = 0.05) +
    geom_hline(yintercept = 0, colour = "grey50",
               linewidth = 0.5, linetype = "dashed") +
    geom_ribbon(aes(ymin = flux_mean - flux_se,
                    ymax = flux_mean + flux_se),
                alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1.6) +
    scale_colour_manual(values = BEHAVIOR_COLORS, name = "Behavior class") +
    scale_fill_manual(  values = BEHAVIOR_COLORS, name = "Behavior class") +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24),
                       limits = c(0, 24)) +
    facet_wrap(~season, ncol = 2) +
    labs(
      title    = "Figure 1 — Diel CH₄ Flux Structure by Season and Behavior Class",
      subtitle = "Mean ± SE across all 30-min observations (shaded) by behavior class",
      x        = "Hour of day",
      y        = "Total CH₄ flux (mg C m⁻² 30 min⁻¹)"
    ) +
    theme_neon() +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal"
    ) +
    guides(colour = guide_legend(override.aes = list(linewidth = 2)),
           fill   = "none")

  ggsave(file.path(OUTDIR, "fig1_diel_combined.png"), fig1,
         width = 11, height = 8, dpi = 180)
  cat("  fig1 done\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — GAM driver effects
# Requires: 30-min Rdata (MaloneLab drive)
# ══════════════════════════════════════════════════════════════════════════════
if (!has_30min) {
  message("Skipping Fig 2 — 30-min data not found.")
} else {
  # Candidate names for each driver (first match wins)
  driver_candidates <- list(
    temp    = c("Tair_f","Tair_C","airTemp","TA_F","ta"),
    moisture = c("VSWC","vswc","SWC","soilMoisture","VSWCMean","swc_mean"),
    par     = c("PAR","par","SW_IN","PPFD_IN","ppfd")
  )
  driver_labels <- c(
    temp     = "Air temperature (°C)",
    moisture = "Volumetric soil water content",
    par      = "log(PAR + 1)  (μmol m⁻² s⁻¹)"
  )

  # Resolve actual column names present in dat30
  driver_cols <- sapply(driver_candidates, function(cands) {
    hit <- intersect(cands, names(dat30))
    if (length(hit) == 0) NA_character_ else hit[1]
  })
  driver_cols <- driver_cols[!is.na(driver_cols)]

  if (length(driver_cols) < 1) {
    message("Skipping Fig 2 — no driver columns found in 30-min data.")
  } else {
    cat(sprintf("  GAM drivers resolved: %s\n",
                paste(names(driver_cols), "=", driver_cols, collapse = ", ")))

    gam_dat <- dat30 %>%
      select(flux_total, SITE_ID, all_of(unname(driver_cols))) %>%
      filter(complete.cases(.)) %>%
      mutate(
        log_par = if ("par" %in% names(driver_cols))
                    log1p(pmax(.data[[driver_cols["par"]]], 0))
                  else NA_real_,
        SITE_ID = factor(SITE_ID)
      ) %>%
      slice_sample(n = min(200000, nrow(.)))

    # Build GAM formula dynamically from resolved columns
    smooth_terms <- c(
      if ("temp"     %in% names(driver_cols)) sprintf("s(%s, k=8)", driver_cols["temp"]),
      if ("moisture" %in% names(driver_cols)) sprintf("s(%s, k=8)", driver_cols["moisture"]),
      if ("par"      %in% names(driver_cols)) "s(log_par, k=8)"
    )
    gam_formula <- as.formula(
      paste("flux_total ~", paste(c(smooth_terms, "s(SITE_ID, bs='re')"), collapse = " + "))
    )

    cat("  Fitting GAM for Fig 2...\n")
    m_flux <- bam(gam_formula, data = gam_dat, method = "fREML", discrete = TRUE)

    # Prediction: vary one driver at a time, hold others at median
    pred_vars <- c(
      if ("temp"     %in% names(driver_cols)) setNames(driver_cols["temp"],     "temp"),
      if ("moisture" %in% names(driver_cols)) setNames(driver_cols["moisture"], "moisture"),
      if ("par"      %in% names(driver_cols)) c(par = "log_par")
    )
    all_pred_cols <- c(unname(driver_cols[names(driver_cols) != "par"]), "log_par")
    all_pred_cols <- intersect(all_pred_cols, names(gam_dat))

    pred_base <- gam_dat %>%
      summarise(across(all_of(all_pred_cols), ~ median(.x, na.rm = TRUE))) %>%
      mutate(SITE_ID = levels(gam_dat$SITE_ID)[1])

    make_pred <- function(col, vals, label) {
      df        <- pred_base[rep(1, length(vals)), ]
      df[[col]] <- vals
      pred      <- predict(m_flux, newdata = df, se.fit = TRUE,
                           exclude = "s(SITE_ID)")
      tibble(driver = label, x = vals,
             fit = pred$fit,
             lo  = pred$fit - 1.64 * pred$se.fit,
             hi  = pred$fit + 1.64 * pred$se.fit)
    }

    pred_df <- bind_rows(lapply(names(pred_vars), function(key) {
      col <- pred_vars[[key]]
      src <- if (key == "par") "log_par" else driver_cols[key]
      make_pred(col, seq(quantile(gam_dat[[col]], 0.02, na.rm = TRUE),
                         quantile(gam_dat[[col]], 0.98, na.rm = TRUE),
                         length.out = 80),
                label = driver_labels[key])
    }))

    fig2 <- ggplot(pred_df, aes(x, fit)) +
      geom_hline(yintercept = 0, colour = "grey55",
                 linewidth = 0.5, linetype = "dashed") +
      geom_ribbon(aes(ymin = lo, ymax = hi),
                  fill = "grey50", alpha = 0.20) +
      geom_line(linewidth = 1.8, colour = "grey15") +
      facet_wrap(~driver, scales = "free_x", ncol = length(pred_vars)) +
      labs(
        title    = "Figure 2 — Population-Level GAM: Driver Effects on CH₄ Flux",
        subtitle = "Site random effect excluded · 90% CI shaded · fitted from 30-min observations",
        x        = NULL,
        y        = "CH₄ flux (mg C m⁻² 30 min⁻¹)"
      ) +
      theme_neon()

    ggsave(file.path(OUTDIR, "fig2_driver_effects.png"), fig2,
           width = 13, height = 5, dpi = 180)
    cat("  fig2 done\n")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Phase space: prop_positive vs mean total flux (REAL)
# ══════════════════════════════════════════════════════════════════════════════
sm3 <- site_summary %>% filter(!is.na(prop_pos), !is.na(mean_total))

fig3 <- ggplot(sm3, aes(prop_pos, mean_total, colour = behavior, fill = behavior)) +
  annotate("rect", xmin = 0,    xmax = 0, ymin = -Inf, ymax = Inf,
           fill = BEHAVIOR_COLORS["Weak-sink"],   alpha = 0.07) +
  annotate("rect", xmin = 0,    xmax = 1, ymin = -Inf, ymax = Inf,
           fill = BEHAVIOR_COLORS["Fluctuating"],       alpha = 0.05) +
  annotate("rect", xmin = 1,    xmax = 1, ymin = -Inf, ymax = Inf,
           fill = BEHAVIOR_COLORS["Weak-source"], alpha = 0.07) +
  geom_vline(xintercept = c(0, 1),
             colour = "grey65", linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0,
             colour = "grey55", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = 0.125, y = Inf,
           label = "Consistent\nsink", vjust = 1.4, size = 3.0,
           fontface = "bold", colour = BEHAVIOR_COLORS["Weak-sink"]) +
  annotate("text", x = 0.500, y = Inf,
           label = "Fluctuating", vjust = 1.4, size = 3.0,
           fontface = "bold", colour = BEHAVIOR_COLORS["Fluctuating"]) +
  annotate("text", x = 0.875, y = Inf,
           label = "Consistent\nsource", vjust = 1.4, size = 3.0,
           fontface = "bold", colour = BEHAVIOR_COLORS["Weak-source"]) +
  geom_point(shape = 21, size = 3.0, stroke = 0.4,
             colour = "white", alpha = 0.92) +
  geom_text_repel(
    aes(label = SITE_ID), size = 2.6, family = "mono",
    max.overlaps = 50,
    segment.colour = "grey55", segment.size = 0.3,
    box.padding = 0.35, point.padding = 0.15,
    min.segment.length = 0.2, show.legend = FALSE
  ) +
  scale_colour_manual(values = BEHAVIOR_COLORS, guide = "none") +
  scale_fill_manual(  values = BEHAVIOR_COLORS, name = "Behavior class") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1), expand = expansion(add = 0.01)) +
  labs(
    title    = "Figure 3 — Site Behavior Phase Space",
    subtitle = "Fraction of months with positive mean CH₄ flux vs mean half-hourly flux",
    x        = "Fraction of months with positive mean CH₄ flux",
    y        = "Mean total CH₄ flux (nmol m⁻² s⁻¹)"
  ) +
  theme_neon() +
  # legend anchored top-left — sink zone is empty, giving clear whitespace there
  theme(
    legend.position     = c(0.03, 0.97),
    legend.justification = c(0, 1)
  ) +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 4, colour = "white"), nrow = 3
  ))

ggsave(file.path(OUTDIR, "fig3_phase_space.png"), fig3,
       width = 9.5, height = 7, dpi = 180)
cat("  fig3 done\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Total vs gradient flux, split by whether storage changed class
# ══════════════════════════════════════════════════════════════════════════════
sm4  <- site_summary %>% filter(!is.na(mean_total), !is.na(median_grad))
lim4 <- max(abs(sm4$mean_total), abs(sm4$median_grad), na.rm = TRUE) * 1.20

panel4 <- function(data, changed, title_txt) {
  sub <- filter(data, changed_class == changed)
  ggplot(sub, aes(median_grad, mean_total,
                  colour = behavior, fill = behavior)) +
    geom_hline(yintercept = 0, colour = "grey60",
               linewidth = 0.4, linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "grey60",
               linewidth = 0.4, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0,
                colour = "grey80", linewidth = 0.8, linetype = "dotted") +
    geom_segment(aes(xend = median_grad, yend = median_grad),
                 alpha = 0.35, linewidth = 0.7) +
    geom_point(shape = 21, size = 2.8, stroke = 0.4,
               colour = "white", alpha = 0.92) +
    {if (changed)
       geom_text_repel(aes(label = SITE_ID), size = 2.8, family = "mono",
                       max.overlaps = 20, show.legend = FALSE,
                       segment.colour = "grey60", segment.size = 0.3,
                       box.padding = 0.45)} +
    scale_colour_manual(values = BEHAVIOR_COLORS, guide = "none") +
    scale_fill_manual(  values = BEHAVIOR_COLORS, guide = "none") +
    coord_fixed(xlim = c(-lim4, lim4), ylim = c(-lim4, lim4)) +
    labs(title = title_txt,
         x = "Median gradient CH₄ flux (nmol m⁻² s⁻¹)",
         y = "Mean total CH₄ flux (nmol m⁻² s⁻¹)") +
    theme_neon()
}

p4a <- panel4(sm4, FALSE, "Class unchanged by storage")
p4b <- panel4(sm4, TRUE,  "Class changed by storage flux")

# Build shared legend from a dummy plot
legend4_df <- distinct(sm4, behavior) %>% filter(!is.na(behavior)) %>%
  mutate(behavior = factor(behavior, levels = names(BEHAVIOR_COLORS)))
legend4 <- get_legend(
  ggplot(legend4_df, aes(1, 1, fill = behavior)) +
    geom_point(shape = 21, size = 4, colour = "white") +
    scale_fill_manual(values = BEHAVIOR_COLORS, name = "Behavior class") +
    theme_void() +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_text(face = "bold", size = 9),
      legend.text      = element_text(size = 9)
    ) +
    guides(fill = guide_legend(
      override.aes = list(size = 4, colour = "white"), nrow = 1
    ))
)

fig4 <- (p4a | p4b) / wrap_elements(full = legend4) +
  plot_layout(heights = c(10, 1)) +
  plot_annotation(
    title    = "Figure 4 — Storage Flux Contribution: Total vs Gradient CH₄",
    subtitle = paste0("Vertical segments show how storage shifts each site relative to 1:1 line",
                      " · labeled sites changed class"),
    theme = theme(plot.title    = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 8.5, colour = "grey40"))
  )

ggsave(file.path(OUTDIR, "fig4_total_vs_gradient.png"), fig4,
       width = 13, height = 7, dpi = 180)
cat("  fig4 done\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Annual budgets by site (all years shown per site, grouped)
# Sites sorted by their mean annual budget; years dodged within each site
# ══════════════════════════════════════════════════════════════════════════════
df5 <- annual_diag %>%
  filter(!is.na(annual_budget_gC_m2_yr), !is.na(behavior)) %>%
  mutate(
    budget = if_else(interpretable_budget %in% c(TRUE, "TRUE") &
                       !is.na(annual_budget_interpretable_gC_m2_yr),
                     annual_budget_interpretable_gC_m2_yr,
                     annual_budget_gC_m2_yr),
    lwr    = if_else(interpretable_budget %in% c(TRUE, "TRUE") &
                       !is.na(annual_budget_interpretable_lower95_gC_m2_yr),
                     annual_budget_interpretable_lower95_gC_m2_yr,
                     annual_budget_lower95_gC_m2_yr),
    upr    = if_else(interpretable_budget %in% c(TRUE, "TRUE") &
                       !is.na(annual_budget_interpretable_upper95_gC_m2_yr),
                     annual_budget_interpretable_upper95_gC_m2_yr,
                     annual_budget_upper95_gC_m2_yr),
    Year = as.character(Year)
  ) %>%
  group_by(SITE_ID) %>%
  mutate(site_mean = mean(budget, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(SITE_ID = fct_reorder(SITE_ID, site_mean))

# Dodge offset per year (up to 4 years: ±0.3 spacing)
year_levels <- sort(unique(df5$Year))
n_yr        <- length(year_levels)
dodge_offsets <- setNames(
  seq(-0.3, 0.3, length.out = n_yr),
  year_levels
)
df5 <- df5 %>%
  mutate(y_pos = as.numeric(SITE_ID) + dodge_offsets[Year])

year_shapes <- c("2021" = 21, "2022" = 22, "2023" = 23, "2024" = 24)

# Join ERA5 year-by-year budgets if available
if (has_era5) {
  df5_era5 <- era5_annual %>%
    select(SITE_ID, Year, era5_budget_yr) %>%
    right_join(df5 %>% select(SITE_ID, Year, y_pos, behavior), by = c("SITE_ID","Year"))
}

fig5 <- ggplot(df5) +
  # thin range line per site (min to max model budget across years)
  geom_linerange(
    data = df5 %>%
      group_by(SITE_ID, behavior) %>%
      summarise(xmin = min(budget, na.rm = TRUE),
                xmax = max(budget, na.rm = TRUE),
                y    = mean(as.numeric(SITE_ID)),
                .groups = "drop"),
    aes(y = y, xmin = xmin, xmax = xmax, colour = behavior),
    linewidth = 0.35, alpha = 0.40
  ) +
  geom_vline(xintercept = 0, colour = "grey35",
             linewidth = 0.8, linetype = "dashed") +
  # 95% CI bars (model-gapfilled)
  geom_errorbarh(aes(y = y_pos, xmin = lwr, xmax = upr, colour = behavior),
                 height = 0, linewidth = 0.7, alpha = 0.50) +
  # ERA5 budget as a small diamond (if available)
  {if (has_era5)
    geom_point(data = df5_era5,
               aes(x = era5_budget_yr, y = y_pos + 0.12),
               shape = 23, size = 1.8, fill = "#6A0DAD", colour = "white",
               stroke = 0.4, alpha = 0.85, na.rm = TRUE)
  } +
  # Model-gapfilled points, shape = year
  geom_point(aes(x = budget, y = y_pos, fill = behavior, shape = Year),
             colour = "white", stroke = 0.5, size = 2.6, alpha = 0.95) +
  scale_colour_manual(values = BEHAVIOR_COLORS, name = "Behavior class") +
  scale_fill_manual(  values = BEHAVIOR_COLORS, name = "Behavior class") +
  scale_shape_manual( values = year_shapes,     name = "Year") +
  scale_y_continuous(
    breaks = seq_along(levels(df5$SITE_ID)),
    labels = levels(df5$SITE_ID),
    expand = expansion(add = 0.6)
  ) +
  scale_x_continuous(labels = label_number(accuracy = 0.01)) +
  labs(
    title    = "Figure 5 — Annual CH₄ Budgets by Site",
    subtitle = if (has_era5)
      "Circles = model-gapfilled (95% CI) · purple diamonds = ERA5-gapfilled · shape = year · sorted by mean budget"
    else
      "Each symbol = one site-year (95% CI shown) · sites sorted by mean annual budget",
    x        = "Annual CH₄ budget (g C m⁻² yr⁻¹)",
    y        = NULL
  ) +
  theme_neon() +
  theme(
    axis.text.y      = element_text(family = "mono", size = 8),
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.title     = element_text(face = "bold", size = 9),
    legend.text      = element_text(size = 8.5)
  ) +
  guides(
    colour = guide_legend(override.aes = list(shape = 21, size = 4),
                          order = 1, nrow = 1),
    fill   = "none",
    shape  = guide_legend(
               override.aes = list(fill = "grey50", colour = "white"),
               order = 2, nrow = 1)
  )

ggsave(file.path(OUTDIR, "fig5_annual_budgets.png"), fig5,
       width = 9, height = 10, dpi = 180)
cat("  fig5 done\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — Two-panel gapfill validation
#   Left:  Model-gapfilled vs observed-only (site-year level, 120 points)
#   Right: ERA5-gapfilled vs observed monthly-bin (site mean level, 39 points)
# ══════════════════════════════════════════════════════════════════════════════

# ── Left panel data ──────────────────────────────────────────────────────────
df6 <- annual_diag %>%
  filter(!is.na(observed_scaled_annual_gC_m2_yr),
         !is.na(annual_budget_gC_m2_yr), !is.na(behavior)) %>%
  mutate(
    gapfill_budget = annual_budget_gC_m2_yr,
    obs_budget     = observed_scaled_annual_gC_m2_yr,
    sign_flip      = sign_flip %in% c(TRUE, "TRUE"),
    yr_label       = paste0("'", substr(as.character(Year), 3, 4))
  )

quant99 <- quantile(c(abs(df6$obs_budget), abs(df6$gapfill_budget)),
                    0.99, na.rm = TRUE)
lim6L   <- max(quant99, max(abs(df6$gapfill_budget), na.rm = TRUE)) * 1.15

df6_clip <- df6 %>% filter(abs(obs_budget) <= lim6L)
rho6L  <- cor(df6_clip$obs_budget, df6_clip$gapfill_budget,
              method = "spearman", use = "complete.obs")
rmse6L <- sqrt(mean((df6_clip$gapfill_budget - df6_clip$obs_budget)^2, na.rm = TRUE))
nflip  <- sum(df6$sign_flip, na.rm = TRUE)

multi_site <- df6 %>% group_by(SITE_ID) %>%
  filter(n() > 1) %>% arrange(SITE_ID, Year) %>% ungroup()

# Helper that builds a scatter panel (used for both left and right)
scatter_panel <- function(df, x_col, y_col, lim,
                          x_lab, y_lab, panel_title,
                          rho, rmse, n_pts,
                          label_col = NULL, flip_col = NULL,
                          path_df = NULL, sign_agree_col = NULL) {
  p <- ggplot(df, aes(.data[[x_col]], .data[[y_col]])) +
    annotate("rect", xmin=-lim, xmax=0,  ymin=0,   ymax= lim,
             fill="#ff7043", alpha=0.07) +
    annotate("rect", xmin=0,   xmax=lim, ymin=-lim, ymax=0,
             fill="#ff7043", alpha=0.07) +
    geom_hline(yintercept=0, colour="grey55", linewidth=0.5, linetype="dashed") +
    geom_vline(xintercept=0, colour="grey55", linewidth=0.5, linetype="dashed") +
    geom_abline(slope=1, intercept=0, colour="grey25", linewidth=1.0)
  if (!is.null(path_df))
    p <- p + geom_path(data=path_df,
                       aes(group=SITE_ID),
                       colour="grey55", linewidth=0.4, alpha=0.55,
                       arrow=arrow(length=unit(0.10,"cm"), type="open", ends="last"))
  p <- p +
    geom_point(aes(fill=behavior), shape=21,
               colour="white", stroke=0.5, size=2.5, alpha=0.90)
  if (!is.null(label_col))
    p <- p + geom_text_repel(aes(label=.data[[label_col]]),
                              size=2.0, colour="grey35", family="mono",
                              segment.colour="grey70", segment.size=0.2,
                              box.padding=0.22, max.overlaps=40, show.legend=FALSE)
  if (!is.null(flip_col))
    p <- p + geom_text_repel(
      data = df %>% filter(.data[[flip_col]]),
      aes(label=SITE_ID),
      size=2.3, colour="#c62828", family="mono", fontface="bold",
      nudge_y=0.06, segment.size=0.3, segment.colour="#c62828",
      box.padding=0.4, max.overlaps=20, show.legend=FALSE)
  p +
    annotate("label", x=-lim*0.97, y=lim*0.97,
             label=sprintf("Spearman ρ = %.2f\nRMSE = %.3f g C m⁻² yr⁻¹\nn = %d",
                           rho, rmse, n_pts),
             hjust=0, vjust=1, size=2.9,
             fill=alpha("white",0.92), label.size=0.3, colour="grey30") +
    scale_fill_manual(values=BEHAVIOR_COLORS, name="Behavior class") +
    coord_fixed(xlim=c(-lim,lim), ylim=c(-lim,lim)) +
    labs(title=panel_title, x=x_lab, y=y_lab) +
    theme_neon() +
    theme(legend.position="none",
          plot.title=element_text(size=10, face="bold"))
}

p6L <- scatter_panel(
  df      = df6,
  x_col   = "obs_budget", y_col = "gapfill_budget", lim = lim6L,
  x_lab   = "Observed-only scaled (g C m⁻² yr⁻¹)",
  y_lab   = "Model-gapfilled (g C m⁻² yr⁻¹)",
  panel_title = "(a) Model-gapfilled vs observed  [site-year]",
  rho = rho6L, rmse = rmse6L, n_pts = nrow(df6_clip),
  label_col = "yr_label", flip_col = "sign_flip",
  path_df = multi_site
)

# ── Right panel: ERA5 (if available) ─────────────────────────────────────────
if (has_era5) {
  df6R  <- era5_comparison %>% filter(!is.na(behavior))
  lim6R <- max(abs(df6R$obs_budget), abs(df6R$era5_budget), na.rm=TRUE) * 1.20
  rho6R  <- cor(df6R$obs_budget, df6R$era5_budget,
                method="spearman", use="complete.obs")
  rmse6R <- sqrt(mean((df6R$era5_budget - df6R$obs_budget)^2, na.rm=TRUE))

  p6R <- scatter_panel(
    df      = df6R,
    x_col   = "obs_budget", y_col = "era5_budget", lim = lim6R,
    x_lab   = "Monthly-bin observed (g C m⁻² yr⁻¹)",
    y_lab   = "ERA5-gapfilled (g C m⁻² yr⁻¹)",
    panel_title = "(b) ERA5-gapfilled vs observed  [site mean]",
    rho = rho6R, rmse = rmse6R, n_pts = nrow(df6R),
    label_col = "SITE_ID", flip_col = "sign_agree"  # sign_agree FALSE = disagreement
  )
  # flip_col labels disagreements; invert logic
  p6R <- scatter_panel(
    df      = df6R %>% mutate(.flip = !sign_agree),
    x_col   = "obs_budget", y_col = "era5_budget", lim = lim6R,
    x_lab   = "Monthly-bin observed (g C m⁻² yr⁻¹)",
    y_lab   = "ERA5-gapfilled (g C m⁻² yr⁻¹)",
    panel_title = "(b) ERA5-gapfilled vs observed  [site mean]",
    rho = rho6R, rmse = rmse6R, n_pts = nrow(df6R),
    label_col = "SITE_ID", flip_col = ".flip"
  )
} else {
  p6R <- ggplot() +
    annotate("text", x=0.5, y=0.5, label="ERA5 output not found\nRun NEON.ERA5.HalfHourlyGapfill.R",
             hjust=0.5, vjust=0.5, size=4, colour="grey50") +
    theme_void()
}

# ── Shared legend ─────────────────────────────────────────────────────────────
leg6_df <- distinct(df6, behavior) %>% filter(!is.na(behavior)) %>%
  mutate(behavior = factor(behavior, levels = names(BEHAVIOR_COLORS)))
leg6 <- get_legend(
  ggplot(leg6_df, aes(1, 1, fill = behavior)) +
    geom_point(shape = 21, size = 4, colour = "white") +
    scale_fill_manual(values = BEHAVIOR_COLORS, name = "Behavior class") +
    theme_void() +
    theme(legend.position  = "bottom",
          legend.direction = "horizontal",
          legend.title     = element_text(face = "bold", size = 9),
          legend.text      = element_text(size = 9)) +
    guides(fill = guide_legend(override.aes = list(size=4, colour="white"), nrow=1))
)

fig6 <- (p6L | p6R) / wrap_elements(full = leg6) +
  plot_layout(heights = c(10, 0.6)) +
  plot_annotation(
    title    = "Figure 6 — Gapfill Validation: Model vs ERA5",
    subtitle = paste0("(a) Model-gapfilled by site-year · arrows = same site across years · ",
                      "red = sign flip\n",
                      "(b) ERA5-gapfilled site means · red labels = ERA5 sign disagrees with observed"),
    theme = theme(plot.title    = element_text(face="bold", size=11),
                  plot.subtitle = element_text(size=8.5, colour="grey40"))
  )

ggsave(file.path(OUTDIR, "fig6_observed_vs_gapfilled.png"), fig6,
       width = 14, height = 8, dpi = 180)
cat("  fig6 done\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 7 — Cliff's delta: environmental driver differences between classes
# Requires: Soildata_YearMon.Rdata (MaloneLab drive)
# ══════════════════════════════════════════════════════════════════════════════
if (!has_soil) {
  message("Skipping Fig 7 — Soildata_YearMon.Rdata not found.")
} else {
  # Join site behaviors to soil/climate driver data
  cat(sprintf("  Fig 7: soil_dat has %d rows; SITE_ID present: %s\n",
              nrow(soil_dat), "SITE_ID" %in% names(soil_dat)))
  drivers_wide <- soil_dat %>%
    left_join(select(site_summary, SITE_ID, behavior), by = "SITE_ID") %>%
    filter(!is.na(behavior))

  # Print ALL column names so we can update the candidate lists if needed
  cat(sprintf("  Fig 7: drivers_wide columns: %s\n",
              paste(names(drivers_wide), collapse = ", ")))

  # Driver columns to test — candidates cover common naming conventions
  driver_cols <- list(
    "Mean VSWC"       = c("VSWCMean_site","VSWC_mean","vswc","VSWCMean","meanVSWC",
                           "soilMoisture","SWC","swc"),
    "VSWC variance"   = c("VSWCVar_site","VSWC_var","vswc_var","VSWCVar","varVSWC"),
    "Total C"         = c("carbonTot","carbon_tot","totalC","total_carbon","C_tot",
                           "organicC","org_carbon"),
    "C:N ratio"       = c("ctonRatio","cton","cn_ratio","CN","cnRatio","c_n_ratio"),
    "Bulk density"    = c("bulkDens","bulk_density","bulkdensity","bulkDensity",
                           "BD","bd"),
    "Clay"            = c("clayTotal","clay_total","clay","Clay","claypct",
                           "clay_pct"),
    "MAP"             = c("MAP","map","mean_annual_precip","meanAnnualPrecip",
                           "precipitation","annualPrecip"),
    "MAT"             = c("MAT","mat","mean_annual_temp","meanAnnualTemp",
                           "temperature","annualTemp"),
    "LAI"             = c("LAI","lai","leafAreaIndex","leaf_area_index"),
    "Canopy height"   = c("canopyHeight","canopy_height","height","Height",
                           "canopy_ht","treeHeight")
  )

  cliff_delta <- function(x, y) {
    # Non-parametric Cliff's delta (x vs y)
    pairs <- outer(x, y, `-`)
    (sum(pairs > 0) - sum(pairs < 0)) / length(pairs)
  }

  # Compare "Weak-source" vs "Fluctuating" for each driver
  src <- drivers_wide %>% filter(behavior == "Weak-source")
  flu <- drivers_wide %>% filter(behavior == "Fluctuating")

  cliff_results <- imap_dfr(driver_cols, function(candidates, label) {
    col <- intersect(candidates, names(drivers_wide))[1]
    if (is.na(col)) {
      cat(sprintf("    [skip] '%s' — no matching column found\n", label))
      return(NULL)
    }
    x <- na.omit(src[[col]])
    y <- na.omit(flu[[col]])
    if (length(x) < 3 || length(y) < 3) {
      cat(sprintf("    [skip] '%s' (col=%s) — too few observations (src=%d, flu=%d)\n",
                  label, col, length(x), length(y)))
      return(NULL)
    }
    cat(sprintf("    [ok]   '%s' using column '%s'\n", label, col))
    d  <- cliff_delta(x, y)
    wt <- wilcox.test(x, y, exact = FALSE)
    tibble(
      label = label,
      group = case_when(
        str_detect(label, "VSWC|MAP|MAT") ~ "Moisture/Climate",
        str_detect(label, "C:|clay|Bulk")  ~ "Soil",
        TRUE ~ "Vegetation/canopy"
      ),
      delta = d,
      pval  = wt$p.value,
      n_src = length(x),
      n_flu = length(y)
    )
  })

  if (nrow(cliff_results) == 0) {
    message(sprintf(
      paste0("Skipping Fig 7 — no driver columns matched in soil data.\n",
             "  Actual columns available: %s\n",
             "  Update the driver_cols candidate lists above to match."),
      paste(names(drivers_wide), collapse = ", ")
    ))
  } else {

  cliff_results <- cliff_results %>%
    arrange(delta) %>%
    mutate(
      label    = factor(label, levels = label),
      y_num    = as.integer(label),          # numeric position for continuous y scale
      evidence = cut(pval, breaks = c(-Inf, 0.05, 0.10, Inf),
                     labels = c("p < 0.05", "p < 0.10", "Weak"))
    )

  n7 <- nrow(cliff_results)

  # alpha must live in the data so each rect gets its own value via aes()
  bands7 <- tribble(
    ~xmin, ~xmax, ~fill,      ~label_x, ~band_label,         ~alpha,
    -1.0,  -0.5,  "#c62828",  -0.75,   "large neg",          0.18,
    -0.5,  -0.3,  "#ef9a9a",  -0.40,   "medium",             0.22,
    -0.3,  -0.1,  "#ffccbc",  -0.20,   "small",              0.25,
    -0.1,   0.1,  "#eeeeee",   0.00,   "negligible",         0.35,
     0.1,   0.3,  "#c8e6c9",   0.20,   "small",              0.25,
     0.3,   0.5,  "#a5d6a7",   0.40,   "medium",             0.22,
     0.5,   1.0,  "#2e7d32",   0.75,   "large pos",          0.18
  )

  # Use numeric y throughout so geom_segment works on a continuous scale;
  # scale_y_continuous maps the integers back to driver labels.
  fig7 <- ggplot(cliff_results, aes(y = y_num)) +
    geom_rect(data = bands7,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = fill, alpha = alpha),
              inherit.aes = FALSE) +
    scale_fill_identity() +
    scale_alpha_identity() +
    annotate("text",
             x     = bands7$label_x,
             y     = rep(n7 + 0.55, nrow(bands7)),
             label = bands7$band_label,
             size = 2.4, colour = "grey55", fontface = "italic", hjust = 0.5) +
    geom_vline(xintercept = 0, colour = "grey35",
               linewidth = 0.8, linetype = "dashed") +
    geom_segment(aes(x = 0, xend = delta, y = y_num, yend = y_num, colour = group),
                 linewidth = 1.6, alpha = 0.70) +
    geom_point(aes(x = delta, y = y_num, colour = group, shape = evidence),
               size = 3.8, fill = "white", stroke = 1.2) +
    scale_colour_manual(values = DRIVER_GROUP_COLORS, name = "Driver group") +
    scale_shape_manual(
      values = c("p < 0.05" = 23, "p < 0.10" = 24, "Weak" = 21),
      name   = "Evidence"
    ) +
    scale_x_continuous(limits = c(-0.9, 0.9)) +
    scale_y_continuous(
      breaks = cliff_results$y_num,
      labels = as.character(cliff_results$label),
      expand = expansion(add = 0.7)
    ) +
    labs(
      title    = "Figure 7 — Environmental Driver Differences: Consistent Source vs Fluctuating",
      subtitle = "Cliff's δ (positive = higher in Weak-source) · bands show effect-size thresholds",
      x        = "Cliff's δ effect size",
      y        = NULL
    ) +
    theme_neon() +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position    = "bottom",
      legend.box         = "horizontal",
      legend.title       = element_text(face = "bold", size = 9),
      legend.text        = element_text(size = 8.5)
    ) +
    guides(
      colour = guide_legend(order = 1, nrow = 2,
                            override.aes = list(shape = NA, linewidth = 2)),
      shape  = guide_legend(order = 2, nrow = 1,
                            override.aes = list(colour = "grey40", fill = "white"))
    )

  ggsave(file.path(OUTDIR, "fig7_cliffs_delta.png"), fig7,
         width = 9, height = 7, dpi = 180)
  cat("  fig7 done\n")
  } # end else (cliff_results non-empty)
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 8 — NEON site map: behavior class + flux magnitude
# Uses ggplot2::borders() — no extra packages needed beyond what's loaded
# ══════════════════════════════════════════════════════════════════════════════
{
  # Embedded NEON coordinates — no file dependency
  .coords8 <- tribble(
    ~SITE_ID,  ~lat,     ~lon,
    "ABBY",  45.762, -122.330, "BARR",  71.283, -156.620,
    "BART",  44.064,  -71.287, "BLAN",  39.036,  -78.072,
    "BONA",  65.154, -147.503, "CLBJ",  33.401,  -97.570,
    "CPER",  40.816, -104.746, "DCFS",  47.162,  -99.106,
    "DEJU",  63.882, -145.751, "DELA",  32.542,  -87.804,
    "DSNY",  28.125,  -81.436, "GRSM",  35.689,  -83.502,
    "GUAN",  18.113,  -65.979, "HARV",  42.537,  -72.172,
    "HEAL",  63.876, -149.213, "JERC",  31.195,  -84.469,
    "JORN",  32.590, -106.843, "KONZ",  39.101,  -96.563,
    "LAJA",  18.023,  -67.076, "LENO",  31.853,  -88.161,
    "MLBS",  37.379,  -80.524, "MOAB",  38.248, -109.388,
    "NIWO",  40.054, -105.582, "NOGP",  46.770, -100.913,
    "OAES",  35.411,  -99.058, "ONAQ",  40.178, -112.452,
    "ORNL",  35.964,  -84.282, "OSBS",  29.689,  -81.993,
    "PUUM",  19.553, -155.318, "RMNP",  40.276, -105.546,
    "SCBI",  38.893,  -78.140, "SERC",  38.890,  -76.560,
    "SJER",  37.109, -119.731, "SOAP",  37.033, -119.264,
    "SRER",  31.911, -110.836, "STEI",  45.509,  -89.586,
    "STER",  40.462, -103.029, "TALL",  32.950,  -87.393,
    "TEAK",  37.006, -119.006, "TOOL",  68.661, -149.373,
    "TREE",  45.494,  -89.586, "UKFS",  39.040,  -95.192,
    "UNDE",  46.234,  -89.538, "WREF",  45.821, -121.952,
    "WOOD",  47.128,  -99.241, "YELL",  44.954, -110.539
  )

  .f8 <- site_summary %>%
    select(SITE_ID, behavior, median_total) %>%
    mutate(behavior = replace_na(behavior, "Fluctuating")) %>%
    inner_join(.coords8, by = "SITE_ID") %>%
    mutate(
      flux_abs = abs(coalesce(median_total, 0)),
      pt_size  = 2 + (flux_abs / max(flux_abs, na.rm = TRUE)) * 6
    )

  cat(sprintf("  Fig 8: %d sites matched\n", nrow(.f8)))

  # Shared theme — no aspect-ratio constraint so plot_grid can resize freely
  .mt8 <- theme_bw(base_size = 10) +
    theme(
      panel.background  = element_rect(fill = "#cfe2f3"),
      panel.grid        = element_blank(),
      axis.title        = element_blank(),
      axis.text         = element_text(size = 6),
      axis.ticks        = element_line(linewidth = 0.25),
      panel.border      = element_rect(colour = "grey40", linewidth = 0.4),
      legend.background = element_rect(fill = alpha("white", 0.9),
                                       colour = "grey60", linewidth = 0.3),
      legend.key        = element_blank()
    )

  # Single panel builder using borders() — clips via xlim/ylim on scales
  .panel8 <- function(dat, xlo, xhi, ylo, yhi, fsize = 2.1,
                      title = NULL, show_leg = FALSE) {

    p <- ggplot(dat, aes(x = lon, y = lat)) +
      borders("world", colour = "grey55", fill = "#ede8e0", linewidth = 0.22) +
      borders("state", colour = "grey70", fill = NA,        linewidth = 0.12) +
      geom_point(aes(fill = behavior, size = pt_size),
                 shape = 21, colour = "white", stroke = 0.55, alpha = 0.93) +
      geom_text_repel(aes(label = SITE_ID),
                      size = fsize, family = "mono", max.overlaps = 80,
                      segment.colour = "grey35", segment.size = 0.20,
                      box.padding = 0.25, point.padding = 0.10,
                      show.legend = FALSE) +
      scale_fill_manual(values = BEHAVIOR_COLORS, name = "Behavior") +
      scale_size_identity() +
      scale_x_continuous(limits = c(xlo, xhi), expand = c(0, 0)) +
      scale_y_continuous(limits = c(ylo, yhi), expand = c(0, 0)) +
      .mt8

    if (!is.null(title))
      p <- p + ggtitle(title) +
        theme(plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5))

    if (show_leg)
      p <- p + theme(
        legend.position      = c(0.01, 0.02),
        legend.justification = c(0, 0),
        legend.title         = element_text(size = 7.5, face = "bold"),
        legend.text          = element_text(size = 7),
        legend.key.size      = unit(0.38, "cm")
      ) + guides(fill = guide_legend(
        override.aes = list(size = 3.5, colour = "white"), nrow = 3
      ))
    else
      p <- p + theme(legend.position = "none")

    p
  }

  p8_conus <- .panel8(.f8 %>% filter(lat < 55, lon > -128),
                      xlo=-126, xhi=-63, ylo=22, yhi=50, show_leg=TRUE)
  p8_ak    <- .panel8(.f8 %>% filter(lat >= 58),
                      xlo=-170, xhi=-140, ylo=58, yhi=72, title="Alaska")
  p8_hi    <- .panel8(.f8 %>% filter(lat>=18, lat<=23, lon< -154),
                      xlo=-161, xhi=-154, ylo=18.8, yhi=22.3, title="Hawaii")
  p8_pr    <- .panel8(.f8 %>% filter(lat>=17, lat<=19, lon> -68),
                      xlo=-67.5, xhi=-65.3, ylo=17.8, yhi=18.65, title="Puerto Rico")

  # Stack insets, then place beside CONUS — no ggdraw, no aspect-ratio conflict
  p8_insets <- plot_grid(p8_ak, p8_hi, p8_pr,
                         ncol = 1, rel_heights = c(2.5, 1.8, 1.0))

  fig8 <- plot_grid(p8_conus, p8_insets,
                    ncol = 2, rel_widths = c(2, 1)) +
    theme(plot.background = element_rect(fill = "white", colour = NA))

  fig8 <- ggdraw(fig8) +
    draw_label("Figure 8 — NEON CH₄ Behavior by Site",
               x = 0.02, y = 0.985, hjust = 0, vjust = 1,
               size = 11, fontface = "bold", colour = "grey10") +
    draw_label("Color = behavior class  ·  Point size ∝ |median flux|",
               x = 0.02, y = 0.958, hjust = 0, vjust = 1,
               size = 8.5, colour = "grey35")

  ggsave(file.path(OUTDIR, "fig8_site_map.png"), fig8,
         width = 15, height = 8, dpi = 180)
  cat("  fig8 done\n")
}

cat(sprintf("\nDone. Figures saved to: %s\n", OUTDIR))

