# VAL_Justification_Supplement.R
# Generates publication figures for the FG validation supplement.
# FG figures use FG_total (gradient + storage flux, from 02b_VAL_TotalFlux.R)
# rather than FG_mean alone, since EC_mean inherently includes storage.
#
# Requires (in localdir.ch4):
#   SITEval_DATA_FILTERED_RSHP_EnSEMBLE_TotalFlux.Rdata
#   CCC_CH4.Rdata
#   VAL_RSHP_pairs.Rdata
#
# Outputs (saved to localdir.ch4):
#   ValSupp_Fig1_TowerProfiles.png   — measurement-height diagrams
#   ValSupp_Fig2_CCC.png             — RSHP concordance by season/approach
#   ValSupp_Fig3_FGvsEC.png          — 30-min FG vs EC scatter
#   ValSupp_Fig4_Annual.png          — annual/seasonal performance summary

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(patchwork)

# ── Output directory ───────────────────────────────────────────────────────────
figure_dir <- "/Volumes/MaloneLab/Research/FluxGradient/Validation_Sites/FIGURES"
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# ── Unit conversion ────────────────────────────────────────────────────────────
unit <- 2 * 0.0000288872 * 1000   # µmol CH4 m-2 s-1 → mgC m-2 per 30 min

# ── Site metadata (from attr files and Val_canopy.csv) ─────────────────────────
site_meta <- list(
  "SE-Sto" = list(
    label     = "SE-Sto\n(Subarctic mire)",
    heights   = c(0.43, 1.00, 1.63, 2.00, 3.00),
    canopy_hc = 0.25,
    ec_ht     = 2.00,
    IGBP      = "WET",
    col       = "#1B7837"
  ),
  "SE-Svb" = list(
    label     = "SE-Svb\n(Boreal tall tower)",
    heights   = c(1, 4, 10, 15, 20, 25, 30, 35, 42, 50, 60, 70, 85, 100, 125, 150),
    canopy_hc = 23.5,
    ec_ht     = 85,
    IGBP      = "ENF",
    col       = "#762A83"
  ),
  "US-Uaf" = list(
    label     = "US-Uaf\n(Boreal black spruce)",
    heights   = c(1, 2, 4, 8),
    canopy_hc = 5.0,
    ec_ht     = 8,
    IGBP      = "ENF",
    col       = "#E08214"
  )
)

# ── RSHP pairs with their selected dLevelsAminusB (canopy-structure filtered) ──
rshp_selected <- list(
  "SE-Sto" = c("3_1","3_2","4_1","4_2","4_3","5_1","5_2","5_3","5_4"),
  "SE-Svb" = c("13_8"),
  "US-Uaf" = c("4_1","4_2","4_3")
)

# ── Load data ──────────────────────────────────────────────────────────────────
load(fs::path(localdir.ch4, "SITEval_DATA_FILTERED_RSHP_EnSEMBLE_TotalFlux.Rdata"))
load(fs::path(localdir.ch4, "CCC_CH4.Rdata"))
load(fs::path(localdir.ch4, "VAL_RSHP_pairs.Rdata"))

season_cols <- c(Winter = "#7FBFEF", Spring = "#A8D5A2",
                 Summer = "#F4A261", Autumn = "#C77DFF")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE S1 — Tower measurement-height profiles
# ══════════════════════════════════════════════════════════════════════════════
# Each panel: site schematic with canopy zone, sensor levels, EC height, and
# the selected RSHP candidate pairs highlighted.

make_tower_panel <- function(site, meta, rshp_pairs_site) {

  hts   <- meta$heights
  hc    <- meta$canopy_hc
  ec    <- meta$ec_ht
  ymax  <- max(hts) * 1.05

  # Assign canopy zone label to each level
  lev_df <- tibble(
    height = hts,
    level  = seq_along(hts),
    zone   = ifelse(hts > hc, "Above canopy", "Within canopy")
  )

  # Pairs in the RSHP candidate set — decode A and B level numbers
  pair_segs <- tibble(dLevelsAminusB = rshp_pairs_site) %>%
    separate(dLevelsAminusB, into = c("A","B"), sep = "_", convert = TRUE) %>%
    mutate(zA = hts[A], zB = hts[B],
           x  = 0.6, xend = 0.6)

  ggplot() +
    # Canopy zone shading
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = hc,
             fill = "#D5E8D4", alpha = 0.6) +
    # Canopy-top line
    geom_hline(yintercept = hc, linetype = "dashed", color = "#4A7C59", linewidth = 0.6) +
    # RSHP candidate pairs
    geom_segment(data = pair_segs,
                 aes(x = x, xend = xend, y = zB, yend = zA),
                 color = "grey60", linewidth = 0.4, alpha = 0.7) +
    # Sensor levels
    geom_point(data = lev_df,
               aes(x = 0.5, y = height, fill = zone),
               shape = 21, size = 3, color = "white", stroke = 0.4) +
    # EC height
    geom_point(data = tibble(h = ec), aes(x = 0.5, y = h),
               shape = 18, size = 4, color = meta$col) +
    # Canopy-height label — offset 3% of ymax above the dashed canopy line
    annotate("text", x = 0.85, y = hc + ymax * 0.03,
             label = paste0("h[c]==", hc, "~m"),
             parse = TRUE, size = 3, hjust = 0, vjust = 0, color = "#4A7C59") +
    # EC-height label — offset 3% above EC marker, capped at 96% of ymax
    annotate("text", x = 0.85, y = pmin(ec + ymax * 0.03, ymax * 0.96),
             label = paste0("EC: ", ec, " m"),
             size = 3, hjust = 0, vjust = 0, color = meta$col) +
    scale_fill_manual(values = c("Within canopy" = "#A8D5A2",
                                  "Above canopy"  = "#AEC6E8"),
                      name = NULL) +
    scale_x_continuous(limits = c(0.3, 1.2), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-ymax * 0.10, ymax),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title = meta$label, y = "Height (m)", x = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(size = 9, face = "bold"))
}

tower_plots <- imap(site_meta, ~ make_tower_panel(.y, .x, rshp_selected[[.y]]))

fig1 <- wrap_plots(tower_plots, nrow = 1) +
  plot_annotation(
    title   = "Measurement height configurations at validation towers",
    caption = paste(strwrap(paste(
      "Green shading = canopy zone (0 to h_c). Open circles = concentration sensor levels",
      "(green = within canopy, blue = above canopy). Diamond = EC flux measurement height.",
      "Grey segments = candidate RSHP height pairs after canopy-structure filtering."),
      width = 100), collapse = "\n"),
    theme   = theme(plot.title   = element_text(size = 11, face = "bold"),
                    plot.caption = element_text(size = 8, hjust = 0))
  )

ggsave(fs::path(figure_dir, "ValSupp_Fig1_TowerProfiles.png"),
       fig1, width = 7.5, height = 5.5, dpi = 300)

message("Saved ValSupp_Fig1_TowerProfiles.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE S2 — CCC concordance by site × season × approach
# ══════════════════════════════════════════════════════════════════════════════

# Keep only CH4 and only canopy-structure-allowed pairs
ccc_allowed <- CCC_VAL %>%
  filter(gas == "CH4") %>%
  rowwise() %>%
  filter(dLevelsAminusB %in% rshp_selected[[Site]]) %>%
  ungroup()

# Join tier labels from rshp_pairs (NA = not selected in that season)
ccc_plot <- ccc_allowed %>%
  left_join(rshp_pairs %>% select(Site, gas, Season, Approach, dLevelsAminusB, RSHP_tier),
            by = c("Site","gas","Season","Approach","dLevelsAminusB")) %>%
  mutate(
    RSHP_tier = replace_na(RSHP_tier, "not selected"),
    RSHP_tier = factor(RSHP_tier,
                       levels = c("standard","best_avail","unreliable","not selected")),
    Season    = factor(Season, levels = c("Winter","Spring","Summer","Autumn")),
    Site      = factor(Site, levels = c("SE-Sto","SE-Svb","US-Uaf"),
                       labels = c("SE-Sto\n(mire)","SE-Svb\n(tall forest)",
                                  "US-Uaf\n(black spruce)"))
  )

tier_cols <- c(standard       = "#1B7837",
               best_avail     = "#E08214",
               unreliable     = "#C0392B",
               `not selected` = "grey75")

fig2 <- ggplot(ccc_plot,
               aes(x = CCC, y = dLevelsAminusB,
                   color = RSHP_tier, shape = Approach)) +
  geom_vline(xintercept = 0,   linetype = "solid",  color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_point(size = 2.5, alpha = 0.85) +
  facet_grid(Site ~ Season, scales = "free_y", space = "free_y") +
  scale_color_manual(values = tier_cols, name = "RSHP tier") +
  scale_shape_manual(values = c(AE = 16, MBR = 17, WP = 15), name = "Approach") +
  scale_x_continuous(limits = c(-0.2, 0.45), breaks = c(-0.1, 0, 0.1, 0.3)) +
  labs(
    title   = "Concordance correlation coefficients for candidate height pairs",
    x       = "Concordance Correlation Coefficient (CCC)",
    y       = "Height pair (upper_lower level)",
    caption = paste(strwrap(paste(
      "Dashed line: CCC = 0.1 (Tier-1 'standard' threshold).",
      "Solid line: CCC = 0. Pairs not reaching Tier-1 in a given season",
      "fall to Tier-2 (best_avail) or Tier-3 (unreliable)."),
      width = 100), collapse = "\n")
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "right",
        strip.text.y     = element_text(angle = 0, size = 8),
        strip.text.x     = element_text(size = 9),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 11, face = "bold"),
        plot.caption     = element_text(size = 8, hjust = 0))

ggsave(fs::path(figure_dir, "ValSupp_Fig2_CCC.png"),
       fig2, width = 7.5, height = 6.5, dpi = 300)

message("Saved ValSupp_Fig2_CCC.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE S3 — 30-min FG vs EC scatter (one panel per site)
# ══════════════════════════════════════════════════════════════════════════════

scatter_list <- imap(site_meta, function(meta, site) {

  df <- SITEval_DATA_FILTERED_RSHPc_H_total[[site]] %>%
    filter(gas == "CH4") %>%
    mutate(FG = FG_total * unit,
           EC = EC_mean * unit) %>%
    drop_na(FG, EC) %>%
    filter(is.finite(FG), is.finite(EC)) %>%
    mutate(Season = factor(Season, levels = c("Winter","Spring","Summer","Autumn")))

  # Stats
  fit  <- lm(FG ~ EC, data = df)
  r2   <- summary(fit)$r.squared
  slp  <- coef(fit)[2]
  n    <- nrow(df)

  # Axis limits — symmetric, clipped at 99th percentile
  lim  <- quantile(c(df$FG, df$EC), 0.995, na.rm = TRUE)
  lim  <- c(-lim * 0.05, lim)

  lbl  <- sprintf("italic(r)^2==%.2f~';'~slope==%.2f~';'~italic(n)==%d",
                  r2, slp, n)

  ggplot(df, aes(x = EC, y = FG, color = Season)) +
    geom_point(size = 0.5, alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                color = "black", linewidth = 0.8) +
    annotate("text", x = lim[1] + diff(lim) * 0.05,
             y = lim[2] * 0.92, label = lbl, parse = TRUE,
             size = 2.8, hjust = 0) +
    scale_color_manual(values = season_cols, name = NULL) +
    coord_cartesian(xlim = lim, ylim = lim) +
    labs(title = meta$label,
         x = expression(EC~flux~(mgC~m^{-2}~per~30~min)),
         y = expression(FG~flux~(mgC~m^{-2}~per~30~min))) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          plot.title      = element_text(size = 9, face = "bold"))
})

fig3 <- wrap_plots(scatter_list, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")

fig3 <- fig3 +
  plot_annotation(
    title   = "Half-hourly flux-gradient vs eddy covariance CH₄ fluxes",
    caption = paste(strwrap(paste(
      "Each point = one 30-min median FG (across RSHP-selected pairs) vs concurrent EC flux.",
      "Dashed line = 1:1. Solid line = OLS regression. Units: mgC m⁻² per 30 min."),
      width = 100), collapse = "\n"),
    theme   = theme(plot.title   = element_text(size = 11, face = "bold"),
                    plot.caption = element_text(size = 8, hjust = 0))
  )

ggsave(fs::path(figure_dir, "ValSupp_Fig3_FGvsEC.png"),
       fig3, width = 7.5, height = 4.5, dpi = 300)

message("Saved ValSupp_Fig3_FGvsEC.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE S4 — Seasonal FG/EC ratio and pseudo-annual budget comparison
# ══════════════════════════════════════════════════════════════════════════════

# Seasonal mean fluxes
seasonal_summary <- imap_dfr(site_meta, function(meta, site) {
  SITEval_DATA_FILTERED_RSHPc_H_total[[site]] %>%
    filter(gas == "CH4") %>%
    mutate(FG = FG_total * unit,
           EC = EC_mean * unit) %>%
    drop_na(FG, EC) %>%
    filter(is.finite(FG), is.finite(EC)) %>%
    group_by(Season) %>%
    summarise(FG_mean_s = mean(FG, na.rm = TRUE),
              EC_mean_s = mean(EC, na.rm = TRUE),
              ratio     = FG_mean_s / EC_mean_s,
              n         = n(),
              .groups   = "drop") %>%
    mutate(Site = site)
}) %>%
  filter(!is.nan(ratio), is.finite(ratio)) %>%
  mutate(Season = factor(Season, levels = c("Winter","Spring","Summer","Autumn")),
         Site   = factor(Site, levels = c("SE-Sto","SE-Svb","US-Uaf"),
                         labels = c("SE-Sto\n(mire)",
                                    "SE-Svb\n(tall forest)",
                                    "US-Uaf\n(black spruce)")))

# Ratio plot
p_ratio <- ggplot(seasonal_summary,
                  aes(x = Season, y = ratio, fill = Season)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.1f×", ratio),
                y = pmin(ratio + 0.3, max(seasonal_summary$ratio) * 0.98)),
            size = 2.5, vjust = 0) +
  facet_wrap(~ Site, nrow = 1) +
  scale_fill_manual(values = season_cols, guide = "none") +
  scale_y_continuous(trans  = "log10",
                     breaks = c(0.5, 1, 2, 5, 10, 20),
                     labels = c("0.5","1","2","5","10","20")) +
  labs(y = "FG / EC ratio (log scale)", x = NULL,
       title = "FG over-estimation relative to EC by season") +
  theme_bw(base_size = 10) +
  theme(strip.text   = element_text(size = 9, face = "bold"),
        plot.title   = element_text(size = 10, face = "bold"),
        axis.text.x  = element_text(angle = 30, hjust = 1))

# Annual pseudo-budget bar chart (FG and EC side by side)
annual_summary <- imap_dfr(site_meta, function(meta, site) {
  df <- SITEval_DATA_FILTERED_RSHPc_H_total[[site]] %>%
    filter(gas == "CH4") %>%
    mutate(FG = FG_total * unit,
           EC = EC_mean * unit) %>%
    drop_na(FG, EC) %>%
    filter(is.finite(FG), is.finite(EC))
  tibble(Site  = site,
         FG_ann = mean(df$FG, na.rm = TRUE) * 48 * 365,
         EC_ann = mean(df$EC, na.rm = TRUE) * 48 * 365)
}) %>%
  pivot_longer(c(FG_ann, EC_ann), names_to = "Method",
               values_to = "Annual_mgC") %>%
  mutate(Method = recode(Method, FG_ann = "Flux-gradient", EC_ann = "Eddy covariance"),
         Site   = factor(Site, levels = c("SE-Sto","SE-Svb","US-Uaf"),
                         labels = c("SE-Sto\n(mire)",
                                    "SE-Svb\n(tall forest)",
                                    "US-Uaf\n(black spruce)")))

p_annual <- ggplot(annual_summary,
                   aes(x = Method, y = Annual_mgC, fill = Method)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = round(Annual_mgC, 0)),
            vjust = -0.4, size = 2.8) +
  facet_wrap(~ Site, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("Flux-gradient"    = "#4292C6",
                                "Eddy covariance" = "#41AB5D"),
                    guide = "none") +
  labs(y = expression(Annual~CH[4]~budget~(mgC~m^{-2}~yr^{-1})),
       x = NULL,
       title = "Scaled-annual CH₄ budget: flux-gradient vs eddy covariance") +
  theme_bw(base_size = 10) +
  theme(strip.text  = element_text(size = 9, face = "bold"),
        plot.title  = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))

fig4 <- p_ratio / p_annual +
  plot_annotation(
    title   = "Seasonal FG/EC ratio and scaled-annual CH₄ budget comparison",
    caption = paste(strwrap(paste(
      "Top: ratio of mean FG to mean EC flux by season (log scale); dashed line = 1:1.",
      "Bottom: scaled-annual budget derived from mean 30-min flux × 48 × 365.",
      "Budgets are representative of sampled periods only (not gap-filled)."),
      width = 100), collapse = "\n"),
    theme   = theme(plot.title   = element_text(size = 11, face = "bold"),
                    plot.caption = element_text(size = 8, hjust = 0))
  )

ggsave(fs::path(figure_dir, "ValSupp_Fig4_Annual.png"),
       fig4, width = 7.5, height = 7.0, dpi = 300)

message("Saved ValSupp_Fig4_Annual.png")
message("All validation supplement figures complete.")
