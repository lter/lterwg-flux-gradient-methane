# NEON CH4 Flux — Plotting Script 
# Author (SLM)
#
# Run AFTER the three analysis scripts in this order:
#   1. flow.30min.analysis.R
#   2. NEON.30min.Gapfill.r
#   3. NEON.ERA5.HalfHourlyGapfill.R
#
# Reads CSVs from OUTPUT/ and writes all figures to FIGURES/.

library(tidyverse)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(ggrepel)
library(ggpubr)

localdir.ch4 <- Sys.getenv(
  "LOCALDIR_CH4",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/Methane")
setwd(localdir.ch4)
dir.create("FIGURES", showWarnings = FALSE, recursive = TRUE)

# ── Figure size constants (manuscript layout, 7" text block width) ─────────────
# FIG_BASE      — simple / single-panel figures (scatter, map, bar)
# FIG_BASE_CPLEX — complex faceted / many-site figures; smaller so labels fit
FIG_BASE        <- 9         # base font size for simple figures
FIG_BASE_CPLEX  <- 7         # base font size for complex multi-panel / many-site figures
FIG_W_FULL      <- 7.0       # full text-block width (in)
FIG_W_HALF      <- 3.4       # single-column width (in)
FIG_H_TALL      <- 10.0      # tall panels (many sites)
FIG_H_MED       <- 8.0       # medium multi-panel
FIG_H_MAP       <- 5.5       # map figures

# ── Constants ──────────────────────────────────────────────────────────────────

behavior_levels <- c("Weak-sink", "Fluctuating", "Weak-source")
behavior_colors <- c(
  "Weak-sink"   = "#2166AC",
  "Fluctuating" = "#4D4D4D",
  "Weak-source" = "#B2182B"
)
ecotype_colors <- c(
  "Cropland"  = "#E69F00",
  "Forest"    = "#009E73",
  "Grassland" = "#F0E442",
  "Shrubland" = "#CC79A7",
  "Wetland"   = "#56B4E9"
)
scale_levels <- c("30 min", "Daily", "Annual")

fill_source_levels <- c(
  "observed",
  "site_month_hour_lookup",
  "site_season_hour_lookup",
  "site_biseason_hour_lookup",
  "site_annual_hour_lookup",
  "site_annual_mean_lookup",
  "not_filled"
)
fill_source_labels <- c(
  observed                  = "Observed",
  site_month_hour_lookup    = "Site-month-hour",
  site_season_hour_lookup   = "Site-season-hour",
  site_biseason_hour_lookup = "Site-biseason-hour",
  site_annual_hour_lookup   = "Site-annual-hour",
  site_annual_mean_lookup   = "Site annual mean",
  not_filled                = "Not filled"
)
fill_source_colors <- c(
  observed                  = "#222222",
  site_month_hour_lookup    = "#1B9E77",
  site_season_hour_lookup   = "#66A61E",
  site_biseason_hour_lookup = "#A6D854",
  site_annual_hour_lookup   = "#B2DF8A",
  site_annual_mean_lookup   = "#E6AB02",
  not_filled                = "#CCCCCC"
)

# ── Base theme — applied to every figure; ensures labels are never clipped ─────
# Wraps theme_bw with consistent margins and subtitle sizing.
theme_neon <- function(base_size = FIG_BASE) {
  theme_bw(base_size = base_size) +
    theme(
      plot.margin   = margin(t = 5, r = 14, b = 5, l = 5, unit = "pt"),
      plot.subtitle = element_text(size = pmax(6.5, base_size - 2),
                                   lineheight = 1.15,
                                   margin = margin(b = 3, unit = "pt")),
      strip.text    = element_text(size = pmax(6.5, base_size - 1.5),
                                   lineheight = 0.95)
    )
}

# ── Helper: duplicate the top column-strip labels at the bottom of a facet_grid ──
# facet_grid() only draws column strips once, at the top. This clones the already-
# rendered top-strip grobs into a new row inserted just above the legend (guide-box)
# — i.e. still within the plot panel area, not stuck below the legend at the very
# bottom of the image. Falls back to appending at the very bottom if no legend is
# found. Aligned to the same panel columns as the top strip.
duplicate_top_strip_to_bottom <- function(gg) {
  g <- ggplotGrob(gg)
  strip_idx <- grep("^strip-t", g$layout$name)
  if (length(strip_idx) == 0) {
    warning("duplicate_top_strip_to_bottom(): no top strips found; returning plot unchanged.")
    return(g)
  }
  strip_layout <- g$layout[strip_idx, ]
  strip_height <- g$heights[unique(strip_layout$t)]

  # Insert directly above the legend (guide-box) if present; otherwise at the bottom.
  # ggplotGrob() always reserves a "guide-box-<side>" slot for all four legend
  # positions, filling the unused ones with an empty zeroGrob placeholder near the
  # top of the gtable -- so we must exclude those and keep only the real legend.
  guide_candidates <- grep("^guide-box", g$layout$name)
  guide_idx <- guide_candidates[
    vapply(guide_candidates, function(i) !inherits(g$grobs[[i]], "zeroGrob"), logical(1))
  ]
  insert_pos <- if (length(guide_idx) > 0) min(g$layout$t[guide_idx]) - 1 else nrow(g)

  g <- gtable::gtable_add_rows(g, strip_height, pos = insert_pos)
  new_row <- insert_pos + 1

  for (i in seq_len(nrow(strip_layout))) {
    g <- gtable::gtable_add_grob(
      g,
      g$grobs[[strip_idx[i]]],
      t = new_row, b = new_row,
      l = strip_layout$l[i], r = strip_layout$r[i],
      name = paste0("strip-b-dup-", i)
    )
  }
  g
}

# ── Helper: render "-2"/"-1" unit exponents as proper superscripts ─────────────
# e.g. "nmol C m-2 s-1" -> "nmol C m⁻² s⁻¹" (m⁻² s⁻¹)
superscript_units <- function(x) {
  x <- gsub("m-2", "m⁻²", x, fixed = TRUE)
  x <- gsub("s-1", "s⁻¹", x, fixed = TRUE)
  x <- gsub("d-1", "d⁻¹", x, fixed = TRUE)
  x <- gsub("yr-1", "yr⁻¹", x, fixed = TRUE)
  x
}

# ── Load CSVs ─────────────────────────────────────────────────────────────────

# ── ERA5 behavior lookup (loaded first; used to reclassify all figures) ────────
# Authoritative ERA5 annual budget classification (sink/source/fluctuating).
# Every figure reclassifies sites through this lookup so all labels are ERA5-based.
# Falls back to lookup-fill behavior if the ERA5 budget script has not been run.
era5_behavior_file <- "OUTPUT/NEON_scale_ERA5_annual_budget_summary.csv"
if (file.exists(era5_behavior_file)) {
  era5_behavior_lookup <- read.csv(era5_behavior_file) %>%
    transmute(
      SITE_ID,
      annual_behavior = factor(era5_annual_behavior, levels = behavior_levels)
    )
  message("Using ERA5 annual behavior classification for all figures (",
          nrow(era5_behavior_lookup), " sites).")
} else {
  era5_behavior_lookup <- NULL
  message("ERA5 behavior file not found — using lookup-fill behavior.")
}

# Replace behavior column(s) in a site-level data frame with ERA5 classification.
# drop_cols: character vector of existing behavior column names to remove first.
# new_col:   name to give the joined ERA5 annual_behavior column.
apply_era5_behavior <- function(df, drop_cols, new_col = "annual_behavior") {
  if (is.null(era5_behavior_lookup)) {
    # Fallback: keep the first existing behavior column, factor it, rename to new_col
    existing <- intersect(drop_cols, names(df))[1]
    return(df %>%
      mutate(!!new_col := factor(.data[[existing]], levels = behavior_levels)) %>%
      select(-any_of(setdiff(drop_cols, new_col))))
  }
  df %>%
    select(-any_of(drop_cols)) %>%
    left_join(era5_behavior_lookup, by = "SITE_ID") %>%
    rename(!!new_col := annual_behavior)
}

# ── Core datasets ──────────────────────────────────────────────────────────────

annual_budget_summary <- read.csv("OUTPUT/NEON_scale_annual_budget_summary.csv") %>%
  mutate(EcoType = as.character(EcoType)) %>%
  apply_era5_behavior(drop_cols = "behavior_annual_scaled", new_col = "behavior_annual_scaled")

# Site order: sorted by ERA5 behavior then annual flux magnitude
annual_site_order <- annual_budget_summary %>%
  arrange(behavior_annual_scaled, flux_scaled_daily_annual_gC_m2_yr, SITE_ID) %>%
  pull(SITE_ID)

scale_long_summary <- read.csv("OUTPUT/NEON_scale_long_flux_budget_summary.csv") %>%
  mutate(
    scale   = recode(scale, "Annual scaled" = "Annual"),
    scale   = factor(scale, levels = scale_levels),
    EcoType = as.character(EcoType)
  ) %>%
  apply_era5_behavior(
    drop_cols = c("behavior", "behavior_30min", "behavior_daily", "behavior_annual_scaled"),
    new_col   = "behavior"
  ) %>%
  mutate(
    behavior_30min         = behavior,
    behavior_daily         = behavior,
    behavior_annual_scaled = behavior
  )

all_site_flux_magnitude_summary <- read.csv("OUTPUT/NEON_all_site_flux_magnitude_summary.csv") %>%
  mutate(
    EcoType = as.character(EcoType),
    scale   = recode(scale, "Annual scaled" = "Annual"),
    # 30-min flux: convert umol -> nmol C m-2 s-1 (x1000)
    across(c(flux_native, flux_lower_native, flux_upper_native),
           ~ if_else(scale == "30 min", .x * 1000, .x)),
    flux_unit   = if_else(scale == "30 min", "nmol C m-2 s-1", flux_unit),
    flux_unit   = superscript_units(flux_unit),
    scale       = factor(scale, levels = scale_levels),
    scale_label = factor(
      paste0(scale, "\n", flux_unit),
      levels = paste0(scale_levels, "\n", superscript_units(c("nmol C m-2 s-1", "mg C m-2 d-1", "g C m-2 yr-1")))
    )
  ) %>%
  apply_era5_behavior(drop_cols = "annual_behavior") %>%
  mutate(SITE_ID_plot = factor(SITE_ID, levels = rev(annual_site_order)))

# Diel summary: re-aggregate from raw site-hour data using ERA5 behavior
site_diel_30min_raw <- read.csv("OUTPUT/NEON_site_diel_30min_by_behavior.csv") %>%
  select(-any_of("behavior_annual_scaled")) %>%
  { if (!is.null(era5_behavior_lookup))
      left_join(., era5_behavior_lookup %>% rename(behavior_annual_scaled = annual_behavior),
                by = "SITE_ID")
    else mutate(., behavior_annual_scaled = NA_character_) } %>%
  filter(!is.na(behavior_annual_scaled), is.finite(hour_num))

diel_behavior_summary <- site_diel_30min_raw %>%
  reframe(
    .by = c(behavior_annual_scaled, hour_num),
    n_sites                = n_distinct(SITE_ID),
    mean_flux_umolC_m2_s   = mean(flux_umolC_m2_s,   na.rm = TRUE),
    sd_flux_umolC_m2_s     = sd(flux_umolC_m2_s,     na.rm = TRUE),
    se_flux_umolC_m2_s     = sd_flux_umolC_m2_s / sqrt(n_sites),
    mean_source_probability = mean(source_probability, na.rm = TRUE),
    sd_source_probability   = sd(source_probability,  na.rm = TRUE),
    se_source_probability   = sd_source_probability / sqrt(n_sites)
  ) %>%
  mutate(
    se_flux_umolC_m2_s      = replace_na(se_flux_umolC_m2_s, 0),
    se_source_probability   = replace_na(se_source_probability, 0),
    behavior_annual_scaled  = factor(behavior_annual_scaled, levels = behavior_levels)
  )

# Behavior counts: recompute from ERA5 lookup
if (!is.null(era5_behavior_lookup)) {
  annual_behavior_counts <- era5_behavior_lookup %>%
    count(annual_behavior, name = "n_sites") %>%
    complete(annual_behavior = factor(behavior_levels, levels = behavior_levels),
             fill = list(n_sites = 0L)) %>%
    rename(behavior = annual_behavior)

  annual_behavior_ecotype_counts <- era5_behavior_lookup %>%
    left_join(annual_budget_summary %>% transmute(SITE_ID, EcoType), by = "SITE_ID") %>%
    filter(!is.na(EcoType)) %>%
    count(annual_behavior, EcoType, name = "n_sites") %>%
    complete(annual_behavior = factor(behavior_levels, levels = behavior_levels),
             EcoType, fill = list(n_sites = 0L)) %>%
    rename(behavior = annual_behavior)
} else {
  annual_behavior_counts <- read.csv("OUTPUT/NEON_annual_behavior_site_counts.csv") %>%
    mutate(behavior = factor(behavior, levels = behavior_levels))
  annual_behavior_ecotype_counts <- read.csv("OUTPUT/NEON_annual_behavior_ecotype_counts.csv") %>%
    mutate(behavior = factor(behavior, levels = behavior_levels))
}

annual_site_map_data <- read.csv("OUTPUT/NEON_annual_site_map_data.csv") %>%
  apply_era5_behavior(drop_cols = "annual_behavior")

daily_fill_source_summary <- read.csv("OUTPUT/NEON_daily_fill_source_summary.csv") %>%
  mutate(fill_source = factor(fill_source, levels = fill_source_levels))

daily_fill_source_by_site <- read.csv("OUTPUT/NEON_daily_fill_source_by_site.csv") %>%
  mutate(
    fill_source       = factor(fill_source, levels = fill_source_levels),
    fill_source_label = factor(fill_source_label, levels = unname(fill_source_labels))
  )

daily_flux <- read.csv("OUTPUT/NEON_scale_daily_flux_all_sites.csv") %>%
  mutate(Date = as.Date(Date))

# ── ERA5 full outputs (model diagnostics + comparison figures) ─────────────────
era5_files_present <- all(file.exists(c(
  "OUTPUT/NEON_ERA5_halfhour_gapfill_fit_plot_data.csv",
  "OUTPUT/NEON_ERA5_halfhour_gapfill_model_effects.csv",
  "OUTPUT/NEON_ERA5_halfhour_gapfill_fit_metrics.csv",
  era5_behavior_file,
  "OUTPUT/NEON_ERA5_all_site_flux_magnitude_summary.csv",
  "OUTPUT/NEON_ERA5_annual_site_map_data.csv",
  "OUTPUT/NEON_ERA5_diel_behavior_summary.csv",
  "OUTPUT/NEON_ERA5_annual_behavior_site_counts.csv",
  "OUTPUT/NEON_ERA5_annual_behavior_ecotype_counts.csv",
  "OUTPUT/NEON_ERA5_reference_annual_class_changes.csv",
  "OUTPUT/NEON_ERA5_scaled_vs_era5_annual_flux_by_site.csv"
)))

if (era5_files_present) {
  fit_plot_data <- read.csv("OUTPUT/NEON_ERA5_halfhour_gapfill_fit_plot_data.csv")
  fit_metrics   <- read.csv("OUTPUT/NEON_ERA5_halfhour_gapfill_fit_metrics.csv")
  effect_grid   <- read.csv("OUTPUT/NEON_ERA5_halfhour_gapfill_model_effects.csv")

  budget_comparison <- read.csv("OUTPUT/NEON_ERA5_vs_model_standardized_budget_comparison.csv") %>%
    mutate(
      reference_annual_behavior = factor(reference_annual_behavior, levels = behavior_levels),
      sign_agree = sign(mean_era5_gapfilled_annual_budget_gC_m2_yr) ==
                   sign(model_standardized_annual_budget_gC_m2_yr)
    )

  era5_annual_class_change_summary <- read.csv("OUTPUT/NEON_ERA5_reference_annual_class_changes.csv") %>%
    mutate(
      reference_annual_behavior = factor(reference_annual_behavior, levels = behavior_levels),
      era5_annual_behavior      = factor(era5_annual_behavior, levels = behavior_levels)
    )

  era5_annual_class_change_matrix <- expand_grid(
    reference_annual_behavior = factor(behavior_levels, levels = behavior_levels),
    era5_annual_behavior      = factor(behavior_levels, levels = behavior_levels)
  ) %>%
    left_join(era5_annual_class_change_summary,
              by = c("reference_annual_behavior", "era5_annual_behavior")) %>%
    mutate(
      n_sites = replace_na(n_sites, 0L),
      changed = reference_annual_behavior != era5_annual_behavior,
      label   = if_else(n_sites > 0, as.character(n_sites), "")
    )

  era5_mag_raw <- read.csv("OUTPUT/NEON_ERA5_all_site_flux_magnitude_summary.csv") %>%
    mutate(
      scale = recode(as.character(scale), "Annual scaled" = "Annual"),
      annual_behavior = factor(annual_behavior, levels = behavior_levels),
      scale = factor(scale, levels = c("30 min", "Daily", "Annual", "Annual ERA5")),
      scale_label = factor(
        paste0(scale, "\n", flux_unit),
        levels = paste0(
          c("30 min", "Daily", "Annual", "Annual ERA5"),
          "\n",
          c("umol C m-2 s-1", "mg C m-2 d-1", "g C m-2 yr-1", "g C m-2 yr-1")
        )
      )
    )

  era5_site_order <- era5_mag_raw %>%
    filter(scale == "Annual ERA5") %>%
    arrange(annual_behavior, flux_native, SITE_ID) %>%
    pull(SITE_ID)

  era5_all_site_flux_magnitude_summary <- era5_mag_raw %>%
    mutate(SITE_ID_plot = factor(SITE_ID, levels = rev(era5_site_order)))

  era5_annual_site_map_data <- read.csv("OUTPUT/NEON_ERA5_annual_site_map_data.csv") %>%
    mutate(annual_behavior = factor(annual_behavior, levels = behavior_levels))

  era5_annual_method_flux_summary <- read.csv("OUTPUT/NEON_ERA5_scaled_vs_era5_annual_flux_by_site.csv") %>%
    mutate(
      annual_behavior = factor(annual_behavior, levels = behavior_levels),
      annual_method   = factor(annual_method, levels = c("Scaled annual", "ERA5 annual"))
    )

  seasonal_file <- "OUTPUT/NEON_ERA5_seasonal_behavior_summary.csv"
  if (file.exists(seasonal_file)) {
    era5_seasonal_behavior_summary <- read.csv(seasonal_file) %>%
      mutate(annual_behavior = factor(annual_behavior, levels = behavior_levels))
  } else {
    era5_seasonal_behavior_summary <- NULL
    message("Seasonal behavior summary not found — re-run NEON.ERA5.HalfHourlyGapfill.R to generate Panel E.")
  }

  era5_diel_behavior_summary <- read.csv("OUTPUT/NEON_ERA5_diel_behavior_summary.csv") %>%
    mutate(annual_behavior = factor(annual_behavior, levels = behavior_levels))

  era5_annual_behavior_counts <- read.csv("OUTPUT/NEON_ERA5_annual_behavior_site_counts.csv") %>%
    mutate(annual_behavior = factor(annual_behavior, levels = behavior_levels))

  era5_annual_behavior_ecotype_counts <- read.csv("OUTPUT/NEON_ERA5_annual_behavior_ecotype_counts.csv") %>%
    mutate(annual_behavior = factor(annual_behavior, levels = behavior_levels))
} else {
  message("ERA5 output files not found — skipping ERA5 plots. Run NEON.ERA5.HalfHourlyGapfill.R first.")
}

# ── flow.30min.analysis.R Figures ─────────────────────────────────────────────

plot_scale_flux <- scale_long_summary %>%
  filter(is.finite(flux_native)) %>%
  mutate(SITE_ID_plot = fct_reorder(SITE_ID, flux_native, .na_rm = TRUE)) %>%
  ggplot(aes(x = flux_native, y = SITE_ID_plot, color = behavior)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_errorbar(
    aes(xmin = flux_lower_native, xmax = flux_upper_native),
    orientation = "y", width = 0.25, linewidth = 0.45, alpha = 0.6
  ) +
  geom_point(size = 2.5, alpha = 0.8) +
  facet_wrap(~ paste0(scale, "\n", flux_unit), scales = "free_x", nrow = 1) +
  scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  theme_neon(FIG_BASE_CPLEX) +
  labs(
    x = "Flux in native units",
    y = NULL,
    color = "State class",
    title = "NEON CH4 Native-Scale Flux Products",
    subtitle = "30-min values are standardized by balanced site-month-hour bins; daily values use site-only lookup fills. Bars show +/- 1 SD."
  ) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text    = element_text(face = "bold", size = 9)
  )

ggsave("FIGURES/NEON_scale_native_flux_by_site.png", plot_scale_flux, width = FIG_W_FULL, height = FIG_H_TALL, units = "in", dpi = 300)

plot_flux_by_behavior <- scale_long_summary %>%
  filter(is.finite(flux_native), !is.na(behavior)) %>%
  ggplot(aes(x = behavior, y = flux_native, color = behavior)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_boxplot(outlier.shape = NA, width = 0.55, color = "grey35", fill = "grey92") +
  geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.75) +
  facet_wrap(~ scale + flux_unit, scales = "free_y", nrow = 1) +
  scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  coord_cartesian(clip = "off") +
  theme_neon(FIG_BASE) +
  labs(x = NULL, y = "Flux in native units", title = "A. Native-scale flux distributions") +
  guides(color = "none") +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", size = 10),
    axis.text.x     = element_text(angle = 25, hjust = 1),
    strip.text      = element_text(face = "bold", size = 9),
    plot.margin     = margin(t = 5, r = 30, b = 5, l = 5, unit = "pt")
  )

plot_diel_flux <- diel_behavior_summary %>%
  ggplot(aes(x = hour_num, y = mean_flux_umolC_m2_s,
             color = behavior_annual_scaled, fill = behavior_annual_scaled)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_ribbon(
    aes(ymin = mean_flux_umolC_m2_s - se_flux_umolC_m2_s,
        ymax = mean_flux_umolC_m2_s + se_flux_umolC_m2_s),
    color = NA, alpha = 0.16
  ) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(0, 23.5)) +
  theme_neon(FIG_BASE) +
  labs(
    x = "Hour of day",
    y = expression(paste("30-min CH"[4], " flux (", mu, "mol C ", m^-2, " ", s^-1, ")")),
    color = "State class",
    title = "B. Mean diel flux pattern"
  ) +
  guides(fill = "none") +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10))

plot_diel_source_probability <- diel_behavior_summary %>%
  ggplot(aes(x = hour_num, y = mean_source_probability,
             color = behavior_annual_scaled, fill = behavior_annual_scaled)) +
  geom_ribbon(
    aes(ymin = pmax(0, mean_source_probability - se_source_probability),
        ymax = pmin(1, mean_source_probability + se_source_probability)),
    color = NA, alpha = 0.16
  ) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(0, 23.5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_neon(FIG_BASE) +
  labs(
    x = "Hour of day", y = "Probability of positive flux",
    color = "State class", title = "C. Diel source probability"
  ) +
  guides(fill = "none") +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10))

plot_annual_behavior_counts <- annual_behavior_counts %>%
  ggplot(aes(x = n_sites, y = fct_rev(behavior), fill = behavior)) +
  geom_col(width = 0.68, color = "grey30", linewidth = 0.25) +
  geom_text(aes(label = n_sites), hjust = -0.25, size = 3.1) +
  scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16)), breaks = scales::breaks_width(5)) +
  theme_neon(FIG_BASE) +
  labs(x = "Number of sites", y = NULL, title = "D. Sites per state class") +
  guides(fill = "none") +
  theme(plot.title = element_text(face = "bold", size = 10), axis.text.y = element_text(size = 9))

plot_annual_behavior_ecotypes <- annual_behavior_ecotype_counts %>%
  ggplot(aes(x = n_sites, y = fct_rev(behavior), fill = EcoType)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = ecotype_colors, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04)), breaks = scales::breaks_width(5)) +
  theme_neon(FIG_BASE) +
  labs(x = "Number of sites", y = NULL, fill = "Ecosystem type", title = "By ecosystem type") +
  theme(
    plot.title    = element_text(face = "bold", size = 10),
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    legend.position = "bottom",
    legend.title  = element_text(size = 9),
    legend.text   = element_text(size = 9)
  )

plot_flux_pattern_diel_panel <- (
  plot_flux_by_behavior /
    ((plot_diel_flux | plot_diel_source_probability) /
       (plot_annual_behavior_counts | plot_annual_behavior_ecotypes))
) +
  plot_layout(heights = c(0.95, 1.7), guides = "collect") +
  plot_annotation(
    title    = "NEON CH4 Flux Patterns By State Class",
    subtitle = paste(strwrap(
      "Categories are defined from annual fluxes; finer-scale panels show 30-min and daily behavior within those state classes. Diel ribbons show +/- 1 SE among sites.",
      width = 90), collapse = "\n")
  ) &
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave("FIGURES/NEON_flux_pattern_diel_behavior_panel.png", plot_flux_pattern_diel_panel, width = FIG_W_FULL, height = FIG_H_TALL, units = "in", dpi = 300)

# plot_all_site_flux_magnitude is built below:
# - with ERA5 as a 4th column  (inside if (era5_files_present))
# - with 3 columns only        (fallback, outside the block)

make_flux_magnitude_plot <- function(dat, subtitle_txt) {
  dat <- dat %>% filter(is.finite(flux_native), !is.na(annual_behavior))

  # Abbreviate "Fluctuating" row label to "F" when strip space is tight
  # (≤ 4 sites → strip height may be too small for the full word)
  n_fluct <- dat %>% distinct(SITE_ID, annual_behavior) %>%
    filter(annual_behavior == "Fluctuating") %>% nrow()
  row_labeller <- labeller(
    annual_behavior = function(x) ifelse(x == "Fluctuating" & n_fluct <= 4, "F", x)
  )
  # Wrap long subtitle so it doesn't overflow the figure width
  subtitle_wrapped <- paste(strwrap(subtitle_txt, width = 80), collapse = "\n")

  dat %>%
    ggplot(aes(x = flux_native, y = SITE_ID_plot, color = annual_behavior, shape = EcoType)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_errorbar(
      aes(xmin = flux_lower_native, xmax = flux_upper_native),
      orientation = "y", width = 0.18, linewidth = 0.45,
      color = "grey35", alpha = 0.65, na.rm = TRUE
    ) +
    geom_point(size = 3.35, alpha = 0.55, stroke = 0.55) +
    facet_grid(rows = vars(annual_behavior), cols = vars(scale_label),
               scales = "free", space = "free_y",
               labeller = row_labeller) +
    scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
    theme_neon(FIG_BASE_CPLEX) +
    labs(
      x = "Flux magnitude in native units", y = NULL, shape = "Ecosystem type",
      #title    = expression(paste("NEON CH"[4], " Site Categories And Flux Magnitudes")),
      subtitle = subtitle_wrapped
    ) +
    guides(color = "none", shape = guide_legend(override.aes = list(size = 3.7, alpha = 1))) +
    theme(
      plot.title      = element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.y     = element_text(size = 7),
      axis.text.x     = element_text(size = 7),
      strip.background = element_rect(fill = "black", color = NA),
      strip.text.x    = element_text(face = "bold", size = 7, lineheight = 0.95, color = "white"),
      strip.text.y    = element_text(face = "bold", size = 7, color = "white"),
      panel.spacing.x = unit(0.45, "lines"),
      panel.spacing.y = unit(0.30, "lines")
    )
}

north_america_map <- rnaturalearth::ne_countries(continent = "North America", returnclass = "sf")

plot_annual_site_map <- ggplot() +
  geom_sf(data = north_america_map, fill = "grey94", color = "white", linewidth = 0.25) +
  geom_point(
    data = annual_site_map_data,
    aes(x = longitude, y = latitude, color = annual_behavior,
        shape = EcoType, size = annual_flux_magnitude_gC_m2_yr),
    alpha = 0.55, stroke = 0.65
  ) +
  scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
  scale_size_continuous(
    range = c(2.4, 8), breaks = scales::breaks_pretty(n = 4),
    name = expression(paste("|Annual flux| (g C ", m^-2, " ", yr^-1, ")"))
  ) +
  coord_sf(xlim = c(-170, -60), ylim = c(15, 72), expand = FALSE) +
  theme_neon(FIG_BASE) +
  labs(
    x = NULL, y = NULL,
    color    = "State class",
    shape    = "Ecosystem type",
    title    = "NEON CH4 Annual Flux Categories Across Sites",
    subtitle = "Color shows state class; symbol size shows absolute annual flux magnitude."
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4.2, alpha = 1)),
    shape = guide_legend(override.aes = list(size = 4.2, alpha = 1))
  ) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box    = "vertical",
    panel.grid.major = element_line(color = "grey88", linewidth = 0.2)
  )

ggsave("FIGURES/NEON_annual_site_category_map.png", plot_annual_site_map, width = FIG_W_FULL, height = FIG_H_MAP, units = "in", dpi = 300)

plot_daily_fill_source_frequency <- daily_fill_source_summary %>%
  mutate(fill_source_label = fct_reorder(fill_source_label, pct_halfhour_slots)) %>%
  ggplot(aes(x = pct_halfhour_slots, y = fill_source_label, fill = fill_source)) +
  geom_col(width = 0.72, color = "grey30", linewidth = 0.2) +
  geom_text(
    aes(label = if_else(pct_halfhour_slots >= 0.1, sprintf("%.1f%%", pct_halfhour_slots), "<0.1%")),
    hjust = -0.12, size = 3
  ) +
  scale_fill_manual(values = fill_source_colors, drop = FALSE, na.translate = FALSE) +
  scale_x_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.16))) +
  theme_neon(FIG_BASE) +
  labs(
    x = "Half-hour slots used in daily flux calculation", y = NULL,
    #title    = "Frequency Of Data Sources Used For Daily CH4 Gap Filling",
    #subtitle = "Observed half-hour fluxes are used first; missing slots are filled by increasingly broad site-specific lookup means."
  ) +
  guides(fill = "none") +
  theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 9))

plot_daily_fill_source_by_site <- daily_fill_source_by_site %>%
  left_join(
    annual_budget_summary %>%
      transmute(SITE_ID, annual_behavior = behavior_annual_scaled),
    by = "SITE_ID"
  ) %>%
  mutate(
    SITE_ID_plot      = fct_reorder(SITE_ID, prop_halfhour_slots, .fun = function(x) sum(x, na.rm = TRUE)),
    fill_source_label = factor(fill_source_label, levels = unname(fill_source_labels))
  ) %>%
  ggplot(aes(x = SITE_ID_plot, y = pct_halfhour_slots, fill = fill_source)) +
  geom_col(width = 0.78, color = NA) +
  facet_grid(
    ~ annual_behavior, scales = "free_x", space = "free_x",
    labeller = labeller(annual_behavior = c(
      "Weak-sink" = "Weak-sink", "Fluctuating" = "F", "Weak-source" = "Weak-source"
    ))
  ) +
  scale_fill_manual(
    values = fill_source_colors, breaks = fill_source_levels,
    labels = fill_source_labels, drop = FALSE, na.translate = FALSE
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.02))) +
  theme_neon(FIG_BASE_CPLEX) +
  labs(
    x = NULL, y = "Half-hour slots", fill = "Data source",
    #title    = "Daily Gap-Fill Data Sources By Site",
    #subtitle = "Bars show the fraction of site-date half-hour slots that were observed or filled from each lookup level."
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title    = element_text(size = 9),
    legend.text     = element_text(size = 9),
    axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
    strip.text      = element_text(face = "bold", size = 9)
  )

plot_daily_fill_source_panel <- plot_daily_fill_source_frequency / plot_daily_fill_source_by_site +
  plot_layout(heights = c(0.85, 1.15))

ggsave("FIGURES/NEON_daily_gapfill_source_frequency.png", plot_daily_fill_source_panel, width = 9.0, height = 11.0, units = "in", dpi = 300)

plot_daily_distribution <- daily_flux %>%
  ggplot(aes(x = daily_mgC_m2_day, y = fct_reorder(SITE_ID, daily_mgC_m2_day, median, .na_rm = TRUE))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_boxplot(outlier.alpha = 0.12, width = 0.55, fill = "grey90", color = "grey30") +
  theme_neon(FIG_BASE_CPLEX) +
  labs(
    x = "Gap-filled daily CH4 flux (mg C m-2 d-1)", y = NULL,
    title    = "Daily Total CH4 Flux By Site",
    subtitle = "Missing half-hour slots are filled from site-only month, season, biseason, and annual-hour means."
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave("FIGURES/NEON_daily_flux_by_site.png", plot_daily_distribution, width = FIG_W_HALF, height = FIG_H_TALL, units = "in", dpi = 300)

# ── NEON.ERA5.HalfHourlyGapfill.R Figures ─────────────────────────────────────

if (era5_files_present) {

  # ── Augment all_site_flux_magnitude_summary with ERA5 annual column ──────────
  # era5_behavior_lookup already loaded at top of script; annual_behavior in
  # all_site_flux_magnitude_summary is already ERA5-based.
  scale_levels_4 <- c("30 min", "Daily", "Annual", "Annual ERA5")
  flux_units_4   <- superscript_units(c("nmol C m-2 s-1", "mg C m-2 d-1", "g C m-2 yr-1", "g C m-2 yr-1"))

  # Site order: ERA5 behavior → ERA5 annual flux magnitude
  era5_site_order_4 <- budget_comparison %>%
    mutate(annual_behavior = factor(era5_annual_behavior, levels = behavior_levels)) %>%
    arrange(annual_behavior, mean_era5_gapfilled_annual_budget_gC_m2_yr, SITE_ID) %>%
    pull(SITE_ID)

  era5_flux_rows <- budget_comparison %>%
    left_join(annual_budget_summary %>% transmute(SITE_ID, EcoType), by = "SITE_ID") %>%
    transmute(
      SITE_ID,
      EcoType,
      annual_behavior   = factor(era5_annual_behavior, levels = behavior_levels),
      scale             = "Annual ERA5",
      flux_unit         = superscript_units("g C m-2 yr-1"),
      flux_native       = mean_era5_gapfilled_annual_budget_gC_m2_yr,
      flux_sd_native    = sd_era5_gapfilled_annual_budget_gC_m2_yr,
      flux_lower_native = flux_native - flux_sd_native,
      flux_upper_native = flux_native + flux_sd_native
    )

  all_site_flux_magnitude_4 <- all_site_flux_magnitude_summary %>%
    select(-scale_label) %>%
    mutate(scale = as.character(scale)) %>%
    bind_rows(era5_flux_rows) %>%
    mutate(
      scale       = factor(scale, levels = scale_levels_4),
      scale_label = factor(
        paste0(scale, "\n", flux_unit),
        levels = paste0(scale_levels_4, "\n", flux_units_4)
      ),
      SITE_ID_plot = factor(SITE_ID, levels = rev(era5_site_order_4))
    )

  plot_all_site_flux_magnitude <- make_flux_magnitude_plot(
    all_site_flux_magnitude_4,
    ""
  )

  ggsave("FIGURES/NEON_all_site_category_flux_magnitudes.png",
         duplicate_top_strip_to_bottom(plot_all_site_flux_magnitude),
         width = FIG_W_FULL, height = FIG_H_TALL, units = "in", dpi = 300)
  fit_axis_lim <- quantile(
    abs(c(fit_plot_data$CH4_mgC_30min, fit_plot_data$fitted_CH4_mgC_30min)),
    0.995, na.rm = TRUE
  )

  plot_fit_observed <- fit_plot_data %>%
    ggplot(aes(x = fitted_CH4_mgC_30min, y = CH4_mgC_30min)) +
    geom_bin2d(bins = 70) +
    geom_abline(slope = 1, intercept = 0, color = "white", linewidth = 0.9) +
    geom_hline(yintercept = 0, color = "grey75", linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "grey75", linetype = "dashed", linewidth = 0.5) +
    coord_cartesian(xlim = c(-fit_axis_lim, fit_axis_lim), ylim = c(-fit_axis_lim, fit_axis_lim)) +
    scale_fill_viridis_c(option = "magma", trans = "sqrt", name = "Count") +
    labs(
      title    = "A. Observed vs fitted",
      subtitle = paste0("RMSE = ", signif(fit_metrics$rmse_mgC_m2_30min, 3),
                        "; r = ", signif(fit_metrics$correlation_observed_fitted, 3)),
      x = "Fitted CH4 flux (mg C m-2 30 min-1)",
      y = "Observed CH4 flux (mg C m-2 30 min-1)"
    ) +
    theme_neon(FIG_BASE) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

  plot_fit_residuals <- fit_plot_data %>%
    ggplot(aes(x = fitted_CH4_mgC_30min, y = residual_CH4_mgC_30min)) +
    geom_bin2d(bins = 70) +
    geom_hline(yintercept = 0, color = "white", linewidth = 0.9) +
    coord_cartesian(
      xlim = c(-fit_axis_lim, fit_axis_lim),
      ylim = quantile(fit_plot_data$residual_CH4_mgC_30min, c(0.005, 0.995), na.rm = TRUE)
    ) +
    scale_fill_viridis_c(option = "magma", trans = "sqrt", name = "Count") +
    labs(
      title    = "B. Residuals vs fitted",
      subtitle = paste0("Bias = ", signif(fit_metrics$bias_mgC_m2_30min, 3),
                        "; MAE = ", signif(fit_metrics$mae_mgC_m2_30min, 3)),
      x = "Fitted CH4 flux (mg C m-2 30 min-1)",
      y = "Residual (observed - fitted)"
    ) +
    theme_neon(FIG_BASE) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

  plot_model_effects <- effect_grid %>%
    ggplot(aes(x = driver_value, y = pred_CH4_mgC_30min)) +
    geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.5) +
    geom_ribbon(aes(ymin = lower_CH4_mgC_30min, ymax = upper_CH4_mgC_30min),
                fill = "grey35", alpha = 0.18) +
    geom_line(color = "grey10", linewidth = 1.0) +
    facet_wrap(~ driver_label, scales = "free_x", ncol = 2) +
    labs(
      title    = "C. Population-level RF effects",
      subtitle = "Partial dependence: other covariates held at reference values; population RF (no SITE_ID).",
      x = NULL, y = "Predicted CH4 flux (mg C m-2 30 min-1)"
    ) +
    theme_neon(FIG_BASE) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey35"),
      strip.background = element_rect(fill = "grey94", color = "grey40"),
      strip.text    = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    )

  plot_residual_distribution <- fit_plot_data %>%
    ggplot(aes(x = residual_CH4_mgC_30min)) +
    geom_vline(xintercept = 0, color = "grey35", linetype = "dashed") +
    geom_histogram(bins = 80, fill = "#5E81AC", color = "white", linewidth = 0.15) +
    coord_cartesian(xlim = quantile(fit_plot_data$residual_CH4_mgC_30min, c(0.005, 0.995), na.rm = TRUE)) +
    labs(
      title = "D. Residual distribution",
      x = "Residual (mg C m-2 30 min-1)", y = "Training observations"
    ) +
    theme_neon(FIG_BASE) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

  plot_era5_model_fit_summary <- (plot_fit_observed | plot_fit_residuals) /
    (plot_model_effects | plot_residual_distribution) +
    plot_layout(heights = c(1, 1.25), guides = "collect") +
    plot_annotation(
      title    = "ERA5 Half-Hourly Gapfill RF Fit Summary (OOB)",
      subtitle = paste0("Training rows = ", fit_metrics$n_training,
                        "; Site RF OOB R² = ", round(fit_metrics$oob_r2_site, 3),
                        "; Pop RF OOB R² = ", round(fit_metrics$oob_r2_pop, 3)),
      caption  = "Observed/fitted diagnostics use OOB predictions from the site RF (includes SITE_ID). Model-effect panels use population RF (no SITE_ID)."
    ) &
    theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, color = "grey35"),
      plot.caption  = element_text(size = 9, color = "grey35"),
      legend.position = "bottom"
    )

  ggsave("FIGURES/NEON_ERA5_halfhour_gapfill_model_fit_summary.png",
         plot_era5_model_fit_summary, width = FIG_W_FULL, height = FIG_H_TALL, units = "in", dpi = 300)

  axis_lim <- max(
    abs(c(budget_comparison$model_standardized_annual_budget_gC_m2_yr,
          budget_comparison$mean_era5_gapfilled_annual_budget_gC_m2_yr,
          budget_comparison$model_standardized_lwr_gC_m2_yr,
          budget_comparison$model_standardized_upr_gC_m2_yr)),
    na.rm = TRUE
  ) * 1.12

  comparison_cor  <- suppressWarnings(cor(
    budget_comparison$mean_era5_gapfilled_annual_budget_gC_m2_yr,
    budget_comparison$model_standardized_annual_budget_gC_m2_yr,
    method = "spearman", use = "complete.obs"
  ))
  comparison_rmse <- sqrt(mean(budget_comparison$budget_difference_era5_minus_model_standardized_gC_m2_yr^2, na.rm = TRUE))

  plot_budget_comparison <- budget_comparison %>%
    ggplot(aes(x = model_standardized_annual_budget_gC_m2_yr,
               y = mean_era5_gapfilled_annual_budget_gC_m2_yr,
               color = era5_annual_behavior)) +
    annotate("rect", xmin = -axis_lim, xmax = 0, ymin = 0, ymax =  axis_lim, fill = "#ff7043", alpha = 0.08) +
    annotate("rect", xmin =  0, xmax =  axis_lim, ymin = -axis_lim, ymax = 0, fill = "#ff7043", alpha = 0.08) +
    geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, color = "grey20", linewidth = 1.2) +
    geom_errorbar(aes(xmin = model_standardized_lwr_gC_m2_yr, xmax = model_standardized_upr_gC_m2_yr),
                  orientation = "y", alpha = 0.30, width = 0, linewidth = 0.8) +
    geom_point(aes(shape = sign_agree, size = sign_agree), alpha = 0.87) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                       labels = c(`TRUE` = "Sign agrees", `FALSE` = "Sign disagrees"), name = NULL) +
    scale_size_manual(values = c(`TRUE` = 2.2, `FALSE` = 3.5), guide = "none") +
    ggrepel::geom_text_repel(
      data = budget_comparison %>% filter(!sign_agree),
      aes(label = SITE_ID), size = 2.7, max.overlaps = 30, show.legend = FALSE
    ) +
    scale_color_manual(values = behavior_colors, na.translate = FALSE) +
    coord_cartesian(xlim = c(-axis_lim, axis_lim), ylim = c(-axis_lim, axis_lim)) +
    labs(
    #  title    = paste(strwrap(
    #    "",
    #    width = 45), collapse = "\n"),
     # subtitle = paste0("Spearman rho = ", signif(comparison_cor, 3),
     #                   "; RMSE = ", signif(comparison_rmse, 3), " g C m⁻² yr⁻¹",
     #                   "\nSign agrees: ",
     #                   sum(budget_comparison$sign_agree, na.rm = TRUE), " of ",
    #                    nrow(budget_comparison), " sites"),
      x       = expression(paste("Balanced (g C ", m^-2, " yr"^-1, ")")),
      y       = expression(paste("ERA5 gapfilled (g C ", m^-2, " yr"^-1, ")")),
      color   = "State class",
      #caption = paste(strwrap(
      #  "Balanced = site-month-hour lookup from NEON.30min.Gapfill.R. Orange quadrants: methods disagree on sign. Bars: 95% simulation CI.",
     #   width = 65), collapse = "\n")
    ) +
    theme_neon(FIG_BASE) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 7.5, color = "grey35"),
      plot.caption  = element_text(size = 7, color = "grey40", hjust = 0),
      legend.position  = "bottom",
      legend.box       = "vertical",
      panel.grid.minor = element_blank()
    )

  ggsave("FIGURES/NEON_ERA5_vs_lut_budget_scatter.png",
         plot_budget_comparison, width = FIG_W_FULL, height = FIG_H_MAP + 1, units = "in", dpi = 300)

  plot_budget_difference <- budget_comparison %>%
    mutate(SITE_ID = fct_reorder(SITE_ID, budget_difference_era5_minus_model_standardized_gC_m2_yr)) %>%
    ggplot(aes(x = budget_difference_era5_minus_model_standardized_gC_m2_yr,
               y = SITE_ID, fill = era5_annual_behavior)) +
    geom_vline(xintercept = 0, color = "grey45", linetype = "dashed") +
    geom_col(width = 0.7, alpha = 0.9) +
    scale_fill_manual(values = behavior_colors, na.translate = FALSE) +
    labs(
      title = paste(strwrap("ERA5 vs. Balanced Annual Budget Difference", width = 35),
                    collapse = "\n"),
      subtitle = "ERA5 minus Balanced budget; positive = ERA5 higher",
      x     = expression(paste(Delta, " Annual budget (g C ", m^-2, " yr"^-1, ")")),
      y     = NULL, fill = "State class"
    ) +
    theme_neon(FIG_BASE_CPLEX) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
          axis.text.y = element_text(size = 7))

  ggsave("FIGURES/NEON_ERA5_vs_lut_budget_differences.png",
         plot_budget_difference, width = FIG_W_HALF, height = FIG_H_TALL, units = "in", dpi = 300)

  plot_era5_annual_class_changes <- era5_annual_class_change_matrix %>%
    mutate(
      tile_fill = case_when(
        !changed & n_sites > 0 ~ "agreement",
        changed  & n_sites > 0 ~ "changed",
        TRUE                    ~ "empty"
      )
    ) %>%
    ggplot(aes(x = era5_annual_behavior, y = reference_annual_behavior)) +
    geom_tile(aes(fill = tile_fill), color = "white", linewidth = 1.4) +
    geom_text(aes(label = label), fontface = "bold", size = 8, color = "grey15") +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_discrete(limits = rev(behavior_levels), drop = FALSE) +
    scale_fill_manual(
      values = c("agreement" = "#c8e6c9", "changed" = "#ffe0b2", "empty" = "grey97"),
      labels = c("agreement" = "Unchanged class", "changed" = "Class changed", "empty" = "No sites"),
      name = NULL, na.translate = FALSE
    ) +
    labs(
      title    = "Reference vs ERA5 Annual-Budget Class",
      subtitle = paste(strwrap(
        "Rows: lookup-filled daily state classes. Columns: ERA5 state classes. Weak-source: 100% of years positive; Weak-sink: 0%; Fluctuating: in between.",
        width = 48), collapse = "\n"),
      x = "ERA5 state class", y = "Reference state class"
    ) +
    theme_neon(FIG_BASE) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 7, color = "grey35"),
      legend.position  = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x  = element_text(face = "bold", angle = 25, hjust = 0),
      axis.text.y  = element_text(face = "bold")
    )

  ggsave("FIGURES/NEON_ERA5_reference_annual_class_changes.png",
         plot_era5_annual_class_changes, width = FIG_W_HALF + 0.5, height = FIG_H_MAP - 1, units = "in", dpi = 300)

  era5_n_fluct <- era5_all_site_flux_magnitude_summary %>%
    filter(!is.na(annual_behavior)) %>%
    distinct(SITE_ID, annual_behavior) %>%
    filter(annual_behavior == "Fluctuating") %>% nrow()
  era5_row_labeller <- labeller(
    annual_behavior = function(x) ifelse(x == "Fluctuating" & era5_n_fluct <= 4, "F", x)
  )

  plot_era5_all_site_flux_magnitude <- era5_all_site_flux_magnitude_summary %>%
    filter(is.finite(flux_native), !is.na(annual_behavior)) %>%
    ggplot(aes(x = flux_native, y = SITE_ID_plot, color = annual_behavior, shape = EcoType)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_errorbar(
      aes(xmin = flux_lower_native, xmax = flux_upper_native),
      orientation = "y", width = 0.18, linewidth = 0.45,
      color = "grey35", alpha = 0.65, na.rm = TRUE
    ) +
    geom_point(size = 3.35, alpha = 0.55, stroke = 0.55) +
    facet_grid(rows = vars(annual_behavior), cols = vars(scale_label),
               scales = "free", space = "free_y",
               labeller = era5_row_labeller) +
    scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
    theme_neon(FIG_BASE_CPLEX) +
    labs(
      x = "Flux magnitude in native units", y = NULL, shape = "Ecosystem type",
      #title = expression(paste("NEON CH"[4], " ERA5-Gapfilled Flux Magnitudes"))
    ) +
    guides(color = "none", shape = guide_legend(override.aes = list(size = 3.7, alpha = 1))) +
    theme(
      plot.title    = element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.y   = element_text(size = 7),
      axis.text.x   = element_text(size = 7),
      strip.text.x  = element_text(face = "bold", size = 7, lineheight = 0.95),
      strip.text.y  = element_text(face = "bold", size = 7),
      panel.spacing.x = unit(0.45, "lines"),
      panel.spacing.y = unit(0.30, "lines")
    )

  ggsave("FIGURES/NEON_ERA5_all_site_category_flux_magnitudes.png",
         plot_era5_all_site_flux_magnitude, width = FIG_W_FULL, height = FIG_H_TALL, units = "in", dpi = 300)

  plot_era5_annual_site_map <- ggplot() +
    geom_sf(data = north_america_map, fill = "grey94", color = "white", linewidth = 0.25) +
    geom_point(
      data = era5_annual_site_map_data,
      aes(x = longitude, y = latitude, color = annual_behavior,
          shape = EcoType, size = annual_flux_magnitude_gC_m2_yr),
      alpha = 0.55, stroke = 0.65
    ) +
    scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
    scale_size_continuous(
      range = c(2.4, 8), breaks = scales::breaks_pretty(n = 4),
      name = expression(paste("ERA5 Annual Flux (g C ", m^-2, " ", yr^-1, ")"))
    ) +
    coord_sf(xlim = c(-170, -60), ylim = c(15, 72), expand = FALSE) +
    theme_neon(FIG_BASE) +
    labs(
      x = NULL, y = NULL,
      color    = "State class",
      shape    = "Ecosystem type",
      title    = "ERA5 Annual Flux Behavior Across Sites",
      subtitle = "Color shows state class; symbol size shows absolute ERA5 annual flux magnitude."
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 4.2, alpha = 1)),
      shape = guide_legend(override.aes = list(size = 4.2, alpha = 1))
    ) +
    theme(
      plot.title    = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box    = "vertical",
      panel.grid.major = element_line(color = "grey88", linewidth = 0.2)
    )

  ggsave("FIGURES/NEON_ERA5_annual_site_category_map.png",
         plot_era5_annual_site_map, width = FIG_W_FULL, height = FIG_H_MAP, units = "in", dpi = 300)

  plot_era5_annual_method_boxplots <- era5_annual_method_flux_summary %>%
    filter(annual_method == "ERA5 annual") %>%
    ggplot(aes(x = annual_behavior, y = annual_flux_gC_m2_yr,
               fill = annual_behavior, col = annual_behavior)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.5,
                 position = position_dodge(width = 0.72)) +
    geom_point(aes(group = annual_method, col = annual_behavior),
               position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72, seed = 20260522),
               size = 2.1, alpha = 0.55, stroke = 0.35) +
    scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
    scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
    theme_neon(FIG_BASE) +
    labs(
      x = "State class",
      y = expression(paste("Annual CH"[4], " flux (g C ", m^-2, " ", yr^-1, ")")),
      title = "A. "
    ) +
    guides(fill = guide_legend(nrow = 1, order = 1, override.aes = list(alpha = 1))) +
    theme(
      plot.title   = element_text(face = "bold", size = 12),
      axis.text.x  = element_text(angle = 35, hjust = 1),
      legend.position = "none"
    )

  # Smoothness check for diel panels ─────────────────────────────────────────
  # Fluctuating sites fluctuate between source and sink by definition, so their
  # mean diel cycle is often erratic (low signal, high inter-site variance).
  # Compute roughness (sum of squared consecutive differences in source
  # probability) for each behavior class.  If Fluctuating is > 2× the median
  # roughness of the other two classes, drop it from both diel panels.
  diel_roughness <- era5_diel_behavior_summary %>%
    arrange(annual_behavior, hour_num) %>%
    group_by(annual_behavior) %>%
    summarise(
      roughness = sum(diff(mean_source_probability)^2, na.rm = TRUE),
      .groups   = "drop"
    )

  fluct_rough  <- diel_roughness %>%
    filter(annual_behavior == "Fluctuating") %>% pull(roughness)
  other_median <- diel_roughness %>%
    filter(annual_behavior != "Fluctuating") %>%
    summarise(m = median(roughness, na.rm = TRUE)) %>% pull(m)

  smooth_threshold <- 2   # Fluctuating must be ≤ this × other median to be shown
  show_fluctuating_diel <- !(
    length(fluct_rough)  > 0 && !is.na(fluct_rough)  &&
    length(other_median) > 0 && !is.na(other_median) &&
    fluct_rough > smooth_threshold * other_median
  )

  if (!show_fluctuating_diel) {
    message(sprintf(
      "Diel panels: removing Fluctuating (roughness %.4f > %.0f× other median %.4f).",
      fluct_rough, smooth_threshold, other_median))
  }

  era5_diel_smooth <- if (show_fluctuating_diel) {
    era5_diel_behavior_summary
  } else {
    era5_diel_behavior_summary %>% filter(annual_behavior != "Fluctuating")
  }

  plot_era5_diel_flux <- era5_diel_smooth %>%
    mutate(
      mean_flux_nmolC_m2_s = mean_flux_umolC_m2_s * 1000,
      se_flux_nmolC_m2_s    = se_flux_umolC_m2_s * 1000
    ) %>%
    ggplot(aes(x = hour_num, y = mean_flux_nmolC_m2_s,
               color = annual_behavior, fill = annual_behavior)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_ribbon(
      aes(ymin = mean_flux_nmolC_m2_s - se_flux_nmolC_m2_s,
          ymax = mean_flux_nmolC_m2_s + se_flux_nmolC_m2_s),
      color = NA, alpha = 0.16
    ) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
    scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
    scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(0, 23.5)) +
    theme_neon(FIG_BASE) +
    labs(
      x = "Hour of day",
      y = expression(paste("30-min CH"[4], " flux (nmol C ", m^-2, " ", s^-1, ")")),
      color = "State Class", title = "C. "
    ) +
    guides(fill = "none", color = guide_legend(nrow = 1, order = 2, override.aes = list(linewidth = 1.2))) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 12))

  plot_era5_diel_source_probability <- era5_diel_smooth %>%
    ggplot(aes(x = hour_num, y = mean_source_probability,
               color = annual_behavior, fill = annual_behavior)) +
    geom_ribbon(
      aes(ymin = pmax(0, mean_source_probability - se_source_probability),
          ymax = pmin(1, mean_source_probability + se_source_probability)),
      color = NA, alpha = 0.16
    ) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
    scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
    scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(0, 23.5)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    theme_neon(FIG_BASE) +
    labs(
      x = "Hour of day", y = "Probability of positive flux",
      color = "State Class", title = "D."
    ) +
    guides(fill = "none", color = guide_legend(nrow = 1, order = 2, override.aes = list(linewidth = 1.2))) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 12))

  plot_era5_annual_behavior_counts <- era5_annual_behavior_counts %>%
    ggplot(aes(x = n_sites, y = fct_rev(annual_behavior), fill = annual_behavior)) +
    geom_col(width = 0.68, color = "grey30", linewidth = 0.25) +
    geom_text(aes(label = n_sites), hjust = -0.25, size = 3.1) +
    scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16)), breaks = scales::breaks_width(5)) +
    theme_neon(FIG_BASE) +
    labs(x = "Number of sites", y = NULL, title = "B. ") +
    guides(fill = "none") +
    theme(plot.title = element_text(face = "bold", size = 12), axis.text.y = element_text(size = 9))

  # Panel E — seasonal (monthly) source probability by behavior class (if data available)
  if (!is.null(era5_seasonal_behavior_summary)) {
    plot_era5_seasonal_source_probability <- era5_seasonal_behavior_summary %>%
      ggplot(aes(x = month, y = mean_source_probability,
                 color = annual_behavior, fill = annual_behavior)) +
      geom_ribbon(
        aes(ymin = pmax(0, mean_source_probability - se_source_probability),
            ymax = pmin(1, mean_source_probability + se_source_probability)),
        color = NA, alpha = 0.16
      ) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.8) +
      scale_color_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE,
                         name = "State Class") +
      scale_fill_manual(values = behavior_colors, drop = FALSE, na.translate = FALSE) +
      scale_x_continuous(breaks = 1:12, labels = month.abb, limits = c(1, 12),
                         expand = expansion(add = 0.3)) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      theme_neon(FIG_BASE) +
      labs(
        x = "Month", y = "Probability of positive flux",
        title = "E."
      ) +
      guides(
        fill  = "none",
        color = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2, size = 2.5))
      ) +
      theme(
        plot.title      = element_text(face = "bold", size = 12),
        legend.position = "bottom",
        axis.text.x     = element_text(angle = 0, hjust = 0.5)
      )

    # 5-panel layout: 2×2 top + E spanning full width on bottom
    top_row    <- ggarrange(plot_era5_annual_method_boxplots, plot_era5_annual_behavior_counts,
                            ncol = 2, nrow = 1)
    middle_row <- ggarrange(plot_era5_diel_flux, plot_era5_diel_source_probability,
                            ncol = 2, nrow = 1)
    plot_era5_flux_pattern_panel <- ggarrange(
      top_row, middle_row, plot_era5_seasonal_source_probability,
      ncol = 1, nrow = 3, heights = c(1, 1, 1.1)
    ) %>%
      annotate_figure(top = text_grob("Flux Patterns by State Class",
                                      color = "black", face = "bold", size = 13,
                                      hjust = 0.5, x = 0.5))
    panel_height <- FIG_H_TALL
  } else {
    # Fallback: 4-panel layout without seasonal panel
    plot_era5_flux_pattern_panel <- ggarrange(
      plot_era5_annual_method_boxplots, plot_era5_annual_behavior_counts,
      plot_era5_diel_flux, plot_era5_diel_source_probability,
      ncol = 2, nrow = 2
    ) %>%
      annotate_figure(top = text_grob("Flux Patterns by State Class",
                                      color = "black", face = "bold", size = 13,
                                      hjust = 0.5, x = 0.5))
    panel_height <- FIG_H_MED
  }

  ggsave("FIGURES/Figure2_NEON_ERA5_flux_pattern_diel_behavior_panel.png",
         plot_era5_flux_pattern_panel,
         width = FIG_W_FULL, height = panel_height, units = "in", dpi = 300)

} else {
  # ERA5 not available — save 3-scale version of flux magnitude figure
  plot_all_site_flux_magnitude <- make_flux_magnitude_plot(
    all_site_flux_magnitude_summary,
    "")
  ggsave("FIGURES/Figure1_NEON_all_site_category_flux_magnitudes.png",
         duplicate_top_strip_to_bottom(plot_all_site_flux_magnitude),
         width = FIG_W_FULL, height = FIG_H_TALL, units = "in", dpi = 300)
  message("ERA5 output files not found — skipping ERA5 plots. Run NEON.ERA5.HalfHourlyGapfill.R first.")
}

message("Wrote all NEON CH4 figures to FIGURES/.")

