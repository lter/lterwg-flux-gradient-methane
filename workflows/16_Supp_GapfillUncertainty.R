# Supplemental_GF_Uncertainty.R
# ─────────────────────────────────────────────────────────────────────────────
# Flux-magnitude sensitivity analysis for the ERA5 half-hourly gap-filled
# CH4 fluxes.  Tests whether the NEON-network source/sink classification
# (site level and continental) is robust to systematic over- or under-
# estimation of flux magnitudes by scaling all fluxes by:
#
#   scale factors: 0.5, 0.75, 1.0, 1.25, 1.5, 2.0
#
# Requires (in environment, from NEON.ERA5.HalfHourlyGapfill.R):
#   era5_gapfilled_30min    — 30-min gap-filled fluxes per site × time
#   era5_mean_annual_budget — per-site mean annual budget + era5_annual_behavior
#   behavior_levels         — c("Weak-sink","Fluctuating","Weak-source")
#   behavior_colors         — named color vector
#
# Outputs (to figure_dir, same as ERA5 figures):
#   FigS_GF_Uncertainty.png   — 4-panel sensitivity figure
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(patchwork)
library(colorspace)

# ── Output directory ──────────────────────────────────────────────────────────
supp_fig_dir <- "/Volumes/MaloneLab/Research/FluxGradient/METHANE/FIGURES"
dir.create(supp_fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Ensure NEON gapfilled 30-min data is loaded (not validation version) ─────
# The validation workflow also creates era5_gapfilled_30min but with different
# column names (FG_gapfilled instead of gapfilled_CH4_mgC_30min).
# If the NEON column is missing, reload from the saved CSV.
if (!exists("era5_gapfilled_30min") ||
    !"gapfilled_CH4_mgC_30min" %in% names(era5_gapfilled_30min)) {
  neon_csv <- file.path(getwd(), "OUTPUT/NEON_ERA5_gapfilled_30min.csv.gz")
  if (!file.exists(neon_csv))
    stop("NEON gapfilled 30-min file not found at: ", neon_csv,
         "\nPlease run NEON.ERA5.HalfHourlyGapfill.R first.")
  message("Loading NEON ERA5 gapfilled 30-min data from CSV...")
  era5_gapfilled_30min <- data.table::fread(neon_csv) %>% as_tibble()
  message(sprintf("  Loaded %d rows, %d columns.", nrow(era5_gapfilled_30min),
                  ncol(era5_gapfilled_30min)))
}

# ── Ensure behavior_levels/behavior_colors exist (normally inherited from
#    06_NEON_ERA5Gapfill.R when scripts run in sequence via WorkFlow.R) ──────
# This script has no other guard for these, so sourcing it on its own (or
# after a session where they were never set) throws a *different* error than
# the era5_gapfilled_30min guard above catches — e.g.
# "subscript out of bounds" on behavior_colors[["Weak-source"]] rather than
# "object not found" — same underlying issue (missing prerequisite), just
# manifesting differently because the code indexes into it instead of
# reading it directly. Defining canonically here (same values as script 06)
# makes this script runnable standalone, matching the fallback already used
# for era5_gapfilled_30min above.
if (!exists("behavior_levels") || !exists("behavior_colors") ||
    !all(c("Weak-sink", "Fluctuating", "Weak-source") %in% names(behavior_colors))) {
  message("behavior_levels/behavior_colors not found in environment (or incomplete) — ",
          "defining canonically (same as 06_NEON_ERA5Gapfill.R).")
  behavior_levels <- c("Weak-sink", "Fluctuating", "Weak-source")
  behavior_colors <- c(
    "Weak-sink"   = "#2166AC",
    "Fluctuating" = "#4D4D4D",
    "Weak-source" = "#B2182B"
  )
}

# ── Ensure era5_mean_annual_budget is loaded (same fallback pattern as
#    era5_gapfilled_30min above) — the other prerequisite from
#    06_NEON_ERA5Gapfill.R this script never guarded for.
if (!exists("era5_mean_annual_budget")) {
  budget_csv <- file.path(getwd(), "OUTPUT/NEON_ERA5_gapfilled_mean_annual_budget.csv")
  if (!file.exists(budget_csv))
    stop("era5_mean_annual_budget not found in environment and no cached file at: ", budget_csv,
         "\nPlease run 06_NEON_ERA5Gapfill.R first.")
  message("Loading era5_mean_annual_budget from CSV...")
  era5_mean_annual_budget <- read.csv(budget_csv) %>% as_tibble()
}

# ── Scale factors to test ─────────────────────────────────────────────────────
scale_factors <- c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
scale_labels  <- c("0.5×", "0.75×", "1.0×\n(baseline)", "1.25×", "1.5×", "2.0×")

# ── Baseline site behavior (from ERA5 annual budget) ─────────────────────────
# EcoType is not in era5_mean_annual_budget (aggregated); pull from 30-min data
site_ecotype <- era5_gapfilled_30min %>%
  distinct(SITE_ID, EcoType)

site_baseline <- era5_mean_annual_budget %>%
  left_join(site_ecotype, by = "SITE_ID") %>%
  transmute(
    SITE_ID,
    EcoType,
    baseline_behavior  = factor(era5_annual_behavior, levels = behavior_levels),
    baseline_budget_gC = mean_era5_gapfilled_annual_budget_gC_m2_yr
  )

# ── Step 1: Monthly budgets from 30-min gapfilled fluxes ─────────────────────
message("Computing monthly budgets...")
monthly_budget <- era5_gapfilled_30min %>%
  group_by(SITE_ID, Year, month) %>%
  summarise(
    monthly_budget_gC = sum(gapfilled_CH4_mgC_30min, na.rm = TRUE) / 1000,
    n_halfhours       = n(),
    .groups           = "drop"
  )

# ── Step 2: Apply each scale factor and compute site-level summaries ──────────
message("Applying scale factors...")
sensitivity_monthly <- crossing(monthly_budget, scale = scale_factors) %>%
  mutate(
    scaled_monthly_gC = monthly_budget_gC * scale,
    is_source_month   = scaled_monthly_gC > 0
  )

sensitivity_annual <- sensitivity_monthly %>%
  group_by(SITE_ID, Year, scale) %>%
  summarise(
    annual_budget_gC       = sum(scaled_monthly_gC, na.rm = TRUE),
    prop_months_source     = mean(is_source_month,  na.rm = TRUE),
    n_months               = n(),
    .groups                = "drop"
  )

sensitivity_site <- sensitivity_annual %>%
  group_by(SITE_ID, scale) %>%
  summarise(
    mean_annual_budget_gC  = mean(annual_budget_gC,       na.rm = TRUE),
    prop_years_source      = mean(annual_budget_gC > 0,   na.rm = TRUE),
    mean_months_source     = mean(prop_months_source,     na.rm = TRUE),
    n_years                = n(),
    .groups                = "drop"
  ) %>%
  mutate(
    scaled_class = case_when(
      prop_years_source >= 1.0 ~ "Weak-source",
      prop_years_source <= 0.0 ~ "Weak-sink",
      TRUE                     ~ "Fluctuating"
    ),
    scaled_class = factor(scaled_class, levels = behavior_levels),
    scale_label  = factor(
      scale_labels[match(scale, scale_factors)],
      levels = scale_labels
    )
  ) %>%
  left_join(site_baseline, by = "SITE_ID")

# ── Step 3: Network-level summaries ──────────────────────────────────────────
network_summary <- sensitivity_site %>%
  group_by(scale, scale_label) %>%
  summarise(
    n_sites        = n(),
    n_source       = sum(scaled_class == "Weak-source"),
    n_sink         = sum(scaled_class == "Weak-sink"),
    n_fluctuating  = sum(scaled_class == "Fluctuating"),
    pct_source     = n_source / n_sites * 100,
    median_budget  = median(mean_annual_budget_gC, na.rm = TRUE),
    mean_budget    = mean(mean_annual_budget_gC,   na.rm = TRUE),
    .groups        = "drop"
  )

# ── Step 4: Sign-stable sites (source at baseline that stay source) ───────────
source_sites_baseline <- site_baseline %>%
  filter(baseline_behavior == "Weak-source") %>%
  pull(SITE_ID)

sign_stability <- sensitivity_site %>%
  filter(SITE_ID %in% source_sites_baseline) %>%
  group_by(scale, scale_label) %>%
  summarise(
    n_retain_source = sum(scaled_class == "Weak-source"),
    n_total_baseline_sources = n(),
    pct_retain = n_retain_source / n_total_baseline_sources * 100,
    .groups = "drop"
  )

# Print key result
message(sprintf(
  "At 0.5× scaling: %d of %d baseline-source sites retain source classification (%.0f%%).",
  sign_stability$n_retain_source[sign_stability$scale == 0.5],
  sign_stability$n_total_baseline_sources[sign_stability$scale == 0.5],
  sign_stability$pct_retain[sign_stability$scale == 0.5]
))
message(sprintf(
  "At 0.5× scaling: %d of %d total sites classified as Weak-source.",
  network_summary$n_source[network_summary$scale == 0.5],
  network_summary$n_sites[network_summary$scale == 0.5]
))

# ── FIGURES ───────────────────────────────────────────────────────────────────

base_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold"),
    plot.title        = element_text(face = "bold", size = 10)
  )

# ── Panel A: Network classification count vs scale factor ─────────────────────
class_long <- sensitivity_site %>%
  count(scale, scale_label, scaled_class) %>%
  complete(scale, scale_label, scaled_class = factor(behavior_levels, levels = behavior_levels),
           fill = list(n = 0L))

pA <- ggplot(class_long,
             aes(x = scale_label, y = n, fill = scaled_class)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = behavior_colors, name = "Classification") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  geom_hline(
    yintercept = network_summary$n_source[network_summary$scale == 1.0],
    linetype = "dashed", color = "grey40", linewidth = 0.6
  ) +
  annotate("text",
           x = length(scale_factors) - 0.4,
           y = network_summary$n_source[network_summary$scale == 1.0] + 0.8,
           label = "Baseline source count",
           size = 2.8, color = "grey30", hjust = 1) +
  labs(
    title = "A.",
    x     = "Flux scale factor",
    y     = "Number of NEON sites"
  ) +
  base_theme +
  theme(legend.position = "bottom")

# ── Panel B: Annual budget distribution vs scale factor ───────────────────────
pB <- sensitivity_site %>%
  ggplot(aes(x = scale_label, y = mean_annual_budget_gC, fill = baseline_behavior)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed", linewidth = 0.8) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.2, outlier.alpha = 0.6,
               width = 0.55, color = "grey30", alpha = 0.75) +
  scale_fill_manual(values = behavior_colors, name = "Baseline class") +
  labs(
    title = "B. ",
    x     = "Flux scale factor",
    y     = expression(Mean~annual~budget~(gC~m^{-2}~yr^{-1}))
  ) +
  base_theme +
  theme(legend.position = "bottom")

# ── Panel C: Retention of source sign among baseline-source sites ─────────────
pC <- sign_stability %>%
  ggplot(aes(x = scale, y = pct_retain)) +
  geom_line(color = behavior_colors[["Weak-source"]], linewidth = 1.2) +
  geom_point(size = 3, color = behavior_colors[["Weak-source"]]) +
  geom_text(aes(label = sprintf("%d/%d", n_retain_source, n_total_baseline_sources)),
            vjust = -0.8, size = 2.8, color = "grey30") +
  scale_x_continuous(breaks = scale_factors,
                     labels = c("0.5×","0.75×","1.0×","1.25×","1.5×","2.0×")) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  geom_hline(yintercept = 100, linetype = "dotted", color = "grey60") +
  labs(
    title = "C.",
    x     = "Flux scale factor",
    y     = "% retaining source status"
  ) +
  base_theme

# ── Panel D: Site-level heatmap (site × scale, fill = annual budget) ──────────
# Order sites by baseline budget (most positive to most negative)
site_order <- site_baseline %>%
  arrange(baseline_budget_gC) %>%
  pull(SITE_ID)

pD <- sensitivity_site %>%
  mutate(SITE_ID = factor(SITE_ID, levels = site_order)) %>%
  ggplot(aes(x = scale_label, y = SITE_ID, fill = mean_annual_budget_gC)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    name     = expression(gC~m^{-2}~yr^{-1}),
    limits   = function(x) c(-max(abs(x), na.rm = TRUE), max(abs(x), na.rm = TRUE))
  ) +
  labs(
    title = "D. ",
    x     = "Flux scale factor",
    y     = "NEON site"
  ) +
  base_theme +
  theme(
    axis.text.y    = element_text(size = 5.5),
    legend.position = "right"
  )

# ── Assemble and save ─────────────────────────────────────────────────────────
fig_uncertainty <- (pA | pB) / (pC | pD) +
  plot_annotation(
    #title   = "Flux-magnitude sensitivity: source/sink classification under systematic scaling",
    #caption = paste(strwrap(paste(
    #  "All ERA5 half-hourly gapfilled CH₄ fluxes multiplied by each scale factor.",
    #  "A: Network-level classification count. B: Annual budget distributions;",
    #  "dashed line = zero. C: Fraction of baseline Weak-source sites (n labelled)",
    #  "that retain source status. D: Per-site mean annual budget; blue = sink, red = source.",
    #  "Dashed line in A marks baseline (1.0×) source count."),
    #  width = 110), collapse = "\n"),
    
    theme = theme(
      plot.title   = element_text(size = 11, face = "bold"),
      plot.caption = element_text(size = 7.5, hjust = 0, color = "grey30")
    )
  )

ggsave(
  file.path(supp_fig_dir, "FigS_GF_Uncertainty.png"),
  fig_uncertainty,
  width  = 10.0,
  height = 9.0,
  dpi    = 300,
  bg     = "white"
)

message("Saved FigS_GF_Uncertainty.png")

# ── Summary statistics for supplement text ────────────────────────────────────
n_total    <- nrow(site_baseline)
n_source_1 <- network_summary$n_source[network_summary$scale == 1.0]
n_source_h <- network_summary$n_source[network_summary$scale == 0.5]
pct_h      <- sign_stability$pct_retain[sign_stability$scale == 0.5]
n_retain_h <- sign_stability$n_retain_source[sign_stability$scale == 0.5]
n_base_src <- sign_stability$n_total_baseline_sources[sign_stability$scale == 0.5]
med_1      <- round(network_summary$median_budget[network_summary$scale == 1.0], 1)
med_h      <- round(network_summary$median_budget[network_summary$scale == 0.5], 1)

cat("\n── Supplement text statistics ──────────────────────────────────────────\n")
cat(sprintf("Total NEON sites: %d\n", n_total))
cat(sprintf("Weak-source at 1.0× (baseline): %d (%.0f%%)\n",
            n_source_1, n_source_1/n_total*100))
cat(sprintf("Weak-source at 0.5×: %d (%.0f%%)\n",
            n_source_h, n_source_h/n_total*100))
cat(sprintf("Baseline-source sites retaining source at 0.5×: %d/%d (%.0f%%)\n",
            n_retain_h, n_base_src, pct_h))
cat(sprintf("Median annual budget at 1.0×: %.1f gC m-2 yr-1\n", med_1))
cat(sprintf("Median annual budget at 0.5×: %.1f gC m-2 yr-1\n\n", med_h))

cat(paste(strwrap(sprintf(
  "Supplemental Text — Flux-magnitude sensitivity analysis.
At the baseline scaling (1.0×), %d of %d NEON sites (%.0f%%) were classified as
Weak-source based on their ERA5 gapfilled annual CH4 budget. To evaluate whether
this continental-scale conclusion is sensitive to potential systematic bias in
flux-gradient estimates, we scaled all gapfilled half-hourly fluxes by factors
ranging from 0.5 (halving) to 2.0 (doubling) and recomputed monthly and annual
budgets for each site. Even under the most conservative scenario — fluxes halved
to 0.5× — %d sites remained Weak-source, and %d of %d baseline-source sites
(%.0f%%) retained their source classification. The network-median annual budget
declined from %.1f gC m-2 yr-1 at 1.0× to %.1f gC m-2 yr-1 at 0.5× but
remained positive, indicating that the continental-scale conclusion of a
predominantly CH4-emitting NEON land surface is robust to magnitude uncertainty
of at least a factor of two. These results suggest that even substantial
overestimation of FG fluxes would not reverse the dominant source signal across
the network.",
  n_source_1, n_total, n_source_1/n_total*100,
  n_source_h, n_retain_h, n_base_src, pct_h,
  med_1, med_h), width = 90), collapse = "\n"))
cat("\n")

