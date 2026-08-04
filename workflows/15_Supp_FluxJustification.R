# Flux_Justification_Supplement.R
#
# Addresses three reviewer concerns about the NEON flux-gradient CH4 dataset:
#
#   (a) Is the flux-gradient method positively biased in forest canopies due to
#       within-canopy CH4 storage and lateral mixing inflating apparent surface fluxes?
#
#   (b) Do NEON sites integrate CH4 from nearby wet microhabitats (riparian margins,
#       wet depressions, standing water) that upland chambers would miss?
#
#   (c) Does the high ERA5-RF gap-fill rate dampen extreme sink values while
#       preserving modest source values?
#
# Outputs:
#   FIGURES/Flux_Justification_Supplement.png
#   OUTPUT/Flux_Justification_Supplement_stats.csv
#
# Requires outputs from:
#   NEON.ERA5.HalfHourlyGapfill.R       -> NEON_ERA5_halfhour_gapfill_fit_metrics.csv
#                                           NEON_ERA5_gapfilled_30min.csv.gz
#                                           NEON_ERA5_rf_oob_sign_accuracy_by_site.csv
#                                           NEON_ERA5_rf_sign_sensitivity_by_site.csv
#   flow.30min.analysis.R / flow.TotalFlux.R
#                                        -> 30min_site_behavior.csv
#                                           30min_total_vs_gradient_behavior_comparison.csv
#   NEON.ERA5.AnnualBudget.R             -> NEON_30min_gapfill_annual_budgets.csv

library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
library(ggrepel)

# -- Paths --------------------------------------------------------------------

localdir.ch4 <- "/Volumes/MaloneLab/Research/FluxGradient/METHANE"
setwd(localdir.ch4)

dir.create("FIGURES", showWarnings = FALSE)
dir.create("OUTPUT",  showWarnings = FALSE)

# -- Colors / theme -----------------------------------------------------------

behavior_colors <- c(
  "Consistent sink"   = "#2166AC",
  "Fluctuating"       = "#4D4D4D",
  "Consistent source" = "#B2182B"
)

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey92"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    plot.title        = element_text(face = "bold", size = 10),
    plot.subtitle     = element_text(size = 8, color = "grey40")
  )

# -- Load data ----------------------------------------------------------------

# Site-level 30-min behavior (storage fractions, fluxes, canopy height)
site_beh <- read.csv("OUTPUT/30min_site_behavior.csv", stringsAsFactors = FALSE)

# Annual ERA5 gap-filled behavior classification
annual_beh <- read.csv("OUTPUT/NEON_30min_gapfill_annual_budgets.csv",
                       stringsAsFactors = FALSE) %>%
  select(SITE_ID, annual_behavior, standardized_behavior, prob_annual_source,
         annual_budget_mean_gC_m2_yr)

# Merge so each site has its annual behavior label
site_meta <- site_beh %>%
  left_join(annual_beh, by = "SITE_ID") %>%
  mutate(
    behavior = factor(annual_behavior,
                      levels = c("Consistent sink", "Fluctuating", "Consistent source"))
  )

# Total-flux vs gradient-only comparison (storage term analysis)
tvg <- read.csv("OUTPUT/30min_total_vs_gradient_behavior_comparison.csv",
                stringsAsFactors = FALSE)

# ERA5 gap-fill overall OOB metrics
fit_metrics <- read.csv("OUTPUT/NEON_ERA5_halfhour_gapfill_fit_metrics.csv",
                        stringsAsFactors = FALSE)

# Per-site OOB sign accuracy
sign_acc <- read.csv("OUTPUT/NEON_ERA5_rf_oob_sign_accuracy_by_site.csv",
                     stringsAsFactors = FALSE) %>%
  left_join(annual_beh %>% select(SITE_ID, annual_behavior), by = "SITE_ID") %>%
  mutate(
    behavior = factor(annual_behavior,
                      levels = c("Consistent sink", "Fluctuating", "Consistent source")),
    sign_pct = 100 * sign_accuracy_30min
  )

# Annual budget sign stability (bootstrap sensitivity analysis)
sign_stab <- read.csv("OUTPUT/NEON_ERA5_rf_sign_sensitivity_by_site.csv",
                      stringsAsFactors = FALSE) %>%
  left_join(annual_beh %>% select(SITE_ID, annual_behavior), by = "SITE_ID") %>%
  mutate(
    behavior = factor(annual_behavior,
                      levels = c("Consistent sink", "Fluctuating", "Consistent source"))
  )

# -- Pre-computed OOB QQ quantiles --------------------------------------------
# Extracted from NEON_ERA5_gapfilled_30min.csv.gz (n = 136,497 observed pairs).
# Key finding: variance is compressed symmetrically -- both extreme sinks AND
# extreme sources are pulled toward zero; mean bias is +0.000225 mg C m-2 30min-1.

oob_qq <- data.frame(
  pct    = seq(0, 100, 5),
  obs_q  = c(-10.363733, -0.201337, -0.070289, -0.021098, -0.004160,
             -0.000725, -0.000064,  0.000362,  0.001281,  0.004285,
              0.008401,  0.013360,  0.019863,  0.028382,  0.039695,
              0.054806,  0.074845,  0.104401,  0.150950,  0.246175,
             19.956567),
  pred_q = c(-4.420469, -0.142850, -0.049911, -0.017009, -0.004692,
             -0.000269,  0.002260,  0.004966,  0.007872,  0.011110,
              0.014855,  0.019303,  0.024677,  0.031356,  0.039696,
              0.050350,  0.064350,  0.083961,  0.114603,  0.174959,
              7.066143)
)

oob_n    <- fit_metrics$n_training[1]
oob_bias <- fit_metrics$bias_site_mgC_m2_30min[1]
oob_sign <- fit_metrics$sign_accuracy_site_30min[1]
oob_r2   <- fit_metrics$oob_r2_site[1]


# =============================================================================
# PANEL A -- Concern (a): Within-canopy storage vs. canopy height
#
# If storage inflates source fluxes in tall forests, storage fraction should
# scale with canopy height and be higher for source sites.
#
# Finding: no significant height effect (R2 near zero); network median < 3%;
# max < 16% even at 53-m canopy (WREF). Sink/source classification unchanged.
# =============================================================================

pA_data <- site_meta %>%
  filter(!is.na(canopyHeight_m), !is.na(median_storage_abs_fraction)) %>%
  mutate(
    storage_pct = 100 * median_storage_abs_fraction,
    canopy_m    = as.numeric(canopyHeight_m)
  )

lm_a  <- lm(storage_pct ~ canopy_m, data = pA_data)
r2_a  <- summary(lm_a)$r.squared
p_a   <- coef(summary(lm_a))["canopy_m", "Pr(>|t|)"]
net_median_storage_pct <- median(pA_data$storage_pct, na.rm = TRUE)
max_storage_pct        <- max(pA_data$storage_pct, na.rm = TRUE)

label_sites_A <- pA_data %>%
  filter(storage_pct > quantile(storage_pct, 0.80, na.rm = TRUE) | canopy_m > 25)

pA <- ggplot(pA_data, aes(x = canopy_m, y = storage_pct)) +
  geom_hline(yintercept = net_median_storage_pct,
             linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_smooth(method = "lm", se = TRUE,
              color = "grey55", fill = "grey82", linewidth = 0.7, alpha = 0.25) +
  geom_point(aes(color = behavior, shape = EcoType), size = 2.5, alpha = 0.85) +
  geom_text_repel(data = label_sites_A,
                  aes(label = SITE_ID, color = behavior),
                  size = 2.2, max.overlaps = 12, segment.size = 0.3,
                  show.legend = FALSE) +
  scale_color_manual(values = behavior_colors, name = "Annual exchange",
                     na.value = "grey70") +
  scale_shape_manual(values = c(Forest = 16, Grassland = 17, Shrubland = 15,
                                Cropland = 18, Wetland = 8),
                     name = "Ecosystem") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  annotate("text",
           x     = 1,
           y     = max_storage_pct * 0.94,
           label = sprintf("R² = %.2f, p = %.2f\nNetwork median = %.1f%%",
                           r2_a, p_a, net_median_storage_pct),
           hjust = 0, size = 2.6, color = "grey35") +
  labs(
    title    = "A.",#  Within-canopy storage vs. canopy height",
    #subtitle = "Storage fraction = |storage term| / |total flux| (30-min median per site)",
    x        = "Canopy height (m)",
    y        = "Median |storage| / |total flux| (%)"
  ) +
  base_theme +
  theme(legend.position  = c(0.70, 0.68),
        legend.key.size  = unit(0.35, "cm"),
        legend.text      = element_text(size = 7),
        legend.title     = element_text(size = 7.5))


# =============================================================================
# PANEL B -- Concern (b): Total flux vs. gradient-only flux
#
# If lateral CH4 from wet microhabitats contaminates the FG signal,
# total flux should exceed gradient flux systematically.
#
# Finding: all 46 sites on 1:1 line (r > 0.99); zero sites change behavior
# class when storage is added/removed.
# =============================================================================

pB_data <- tvg %>%
  left_join(annual_beh %>% select(SITE_ID, annual_behavior), by = "SITE_ID") %>%
  left_join(site_beh  %>% select(SITE_ID, EcoType),         by = "SITE_ID") %>%
  mutate(
    behavior         = factor(annual_behavior,
                              levels = c("Consistent sink","Fluctuating","Consistent source")),
    behavior_changed = factor(as.character(behavior_changed_from_gradient),
                              levels = c("FALSE","TRUE"),
                              labels = c("Unchanged","Changed"))
  )

r_b       <- cor(pB_data$mean_total_mgC_30min, pB_data$mean_gradient_mgC_30min,
                 use = "complete.obs")
n_changed <- sum(as.character(pB_data$behavior_changed_from_gradient) == "TRUE",
                 na.rm = TRUE)
n_B  <- nrow(pB_data)
lim_b <- max(abs(c(pB_data$mean_total_mgC_30min,
                   pB_data$mean_gradient_mgC_30min)), na.rm = TRUE) * 1.1

pB <- ggplot(pB_data,
             aes(x = mean_gradient_mgC_30min, y = mean_total_mgC_30min)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey55", linewidth = 0.7) +
  geom_point(aes(color = behavior, shape = behavior_changed),
             size = 2.5, alpha = 0.85) +
  geom_text_repel(
    data = pB_data %>% filter(behavior_changed == "Changed"),
    aes(label = SITE_ID, color = behavior),
    size = 2.2, max.overlaps = 10, segment.size = 0.3, show.legend = FALSE
  ) +
  scale_color_manual(values = behavior_colors, name = "Annual exchange",
                     na.value = "grey70") +
  scale_shape_manual(values = c(Unchanged = 16, Changed = 4),
                     name   = "Behavior changed\nby storage?") +
  coord_equal(xlim = c(-lim_b, lim_b), ylim = c(-lim_b, lim_b)) +
  annotate("text",
           x     = -lim_b * 0.97,
           y     = lim_b  * 0.90,
           label = sprintf("r = %.3f\n%d / %d sites unchanged\nby storage correction",
                           r_b, n_B - n_changed, n_B),
           hjust = 0, size = 2.6, color = "grey35") +
  labs(
    title    = "B.", #  Total vs. gradient-only site-mean flux",
    #subtitle = "All 46 sites on 1:1 line; storage term does not change exchange class",
    x        = expression("Mean gradient flux (mg C m"^{-2}~"30min"^{-1}~")"),
    y        = expression("Mean total flux (mg C m"^{-2}~"30min"^{-1}~")")
  ) +
  base_theme +
  theme(legend.position  = c(0.76, 0.20),
        legend.key.size  = unit(0.35, "cm"),
        legend.text      = element_text(size = 7),
        legend.title     = element_text(size = 7.5))


# =============================================================================
# PANEL C -- Concern (c): RF variance compression -- QQ plot
#
# OOB predicted vs. observed quantiles. The S-shaped deviation from 1:1
# shows variance compression at both tails: extreme sinks AND extreme sources
# are pulled toward zero equally. Near-zero mean bias is preserved.
# =============================================================================

qq_inner <- oob_qq %>% filter(pct >= 5, pct <= 95)
lim_c    <- max(abs(c(qq_inner$obs_q, qq_inner$pred_q))) * 1.12

pC <- ggplot(qq_inner, aes(x = obs_q, y = pred_q)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey55", linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_line(color = "#1F78B4", linewidth = 0.8) +
  geom_point(size = 2.2, color = "#1F78B4") +
  # Mark p5, p50, p95 explicitly
  geom_point(data = qq_inner %>% filter(pct %in% c(5, 50, 95)),
             size = 3.5, shape = 21, fill = "white", color = "#1F78B4", stroke = 1.2) +
  geom_text_repel(
    data = qq_inner %>% filter(pct %in% c(5, 50, 95)),
    aes(label = paste0("p", pct)),
    size = 2.5, nudge_x = 0.01, segment.size = 0.3
  ) +
  coord_equal(xlim = c(-lim_c, lim_c), ylim = c(-lim_c, lim_c)) +
  annotate("text",
           x     = -lim_c * 0.97,
           y     = lim_c  * 0.82,
           label = sprintf(
             "n = %s half-hours\nMean bias = %+.5f mg C m⁻² 30min⁻¹\nSign acc = %.1f%%  |  OOB R² = %.3f\n(p5–p95 shown; outliers excluded)",
             formatC(oob_n, big.mark = ",", format = "d"),
             oob_bias, 100 * oob_sign, oob_r2),
           hjust = 0, size = 2.4, color = "grey30") +
  labs(
    title    = "C.", #  OOB predicted vs. observed flux quantiles (p5 - p95)",
   # subtitle = "S-shape = symmetric variance compression; mean bias ~ 0",
    x        = expression("Observed quantile (mg C m"^{-2}~"30min"^{-1}~")"),
    y        = expression("OOB predicted quantile (mg C m"^{-2}~"30min"^{-1}~")")
  ) +
  base_theme


# =============================================================================
# PANEL D -- Concern (c) continued: Annual budget sign stability
#
# Bootstrapped annual budgets confirm that variance compression does NOT flip
# annual source/sink classification for the vast majority of sites.
# =============================================================================

stab_order <- c("Stable sink",
                "Probable sink (sign-uncertain)",
                "Uncertain",
                "Probable source (sign-uncertain)",
                "Stable source")

stab_colors <- c(
  "Stable sink"                      = "#2166AC",
  "Probable sink (sign-uncertain)"   = "#74ADD1",
  "Uncertain"                        = "#D9D9D9",
  "Probable source (sign-uncertain)" = "#F4A582",
  "Stable source"                    = "#B2182B"
)

stab_summary <- sign_stab %>%
  count(sign_sensitivity_class) %>%
  mutate(sign_sensitivity_class = factor(sign_sensitivity_class, levels = stab_order))

pD <- ggplot(stab_summary,
             aes(x = "", y = n, fill = sign_sensitivity_class)) +
  geom_bar(stat = "identity", width = 0.55, alpha = 0.85) +
  geom_text(aes(label = paste0(n, " sites")),
            position = position_stack(vjust = 0.5),
            size = 2.9, color = "white", fontface = "bold") +
  scale_fill_manual(values = stab_colors,
                    name   = "Annual sign stability",
                    drop   = FALSE,
                    limits = stab_order) +
  scale_y_continuous(breaks = seq(0, 50, 10)) +
  labs(
    title    = "D.",#  Annual budget sign stability (bootstrap)",
  #  subtitle = "Gap-fill variance rarely flips annual source/sink classification",
    x        = NULL,
    y        = "Number of NEON sites"
  ) +
  base_theme +
  theme(legend.key.size  = unit(0.38, "cm"),
        legend.text      = element_text(size = 7),
        legend.title     = element_text(size = 7.5),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank())


# =============================================================================
# Assemble and save figure
# =============================================================================

top_row    <- plot_grid(pA, pB, ncol = 2, align = "hv")
bottom_row <- plot_grid(pC, pD, ncol = 2, align = "hv")

fig_supp <- plot_grid(
  top_row, bottom_row,
  ncol        = 1,
  rel_heights = c(1, 1)
)

ggsave(
  "FIGURES/Flux_Justification_Supplement.png",
  fig_supp,
  width = 11, height = 9.5, dpi = 300, bg = "white"
)
message("Saved: FIGURES/Flux_Justification_Supplement.png")


# =============================================================================
# Summary statistics CSV
# =============================================================================

stats_out <- data.frame(
  concern = c(
    rep("(a) Storage bias", 5),
    rep("(b) Lateral mixing", 2),
    rep("(c) Gap-fill damping", 7)
  ),
  metric = c(
    "Network median |storage|/|total| (%)",
    "Max |storage|/|total| across all sites (%)",
    "R2: storage_pct ~ canopy height",
    "p-value: storage ~ canopy height",
    "Sites with behavior class unchanged by storage",
    "Pearson r: total vs. gradient site-mean flux",
    "Sites behavior unchanged by storage (%)",
    "n OOB observed-predicted pairs",
    "Mean bias, site RF (mg C m-2 30min-1)",
    "Sign accuracy, site RF (%)",
    "OOB R2, site RF",
    "p5: obs / pred (mg C m-2 30min-1)",
    "p95: obs / pred (mg C m-2 30min-1)",
    "Sites with stable annual sign (n / 46)"
  ),
  value = as.character(c(
    round(net_median_storage_pct, 2),
    round(max_storage_pct, 2),
    round(r2_a, 3),
    round(p_a, 3),
    n_B - n_changed,
    round(r_b, 4),
    round(100 * (n_B - n_changed) / n_B, 1),
    oob_n,
    round(oob_bias, 6),
    round(100 * oob_sign, 1),
    round(oob_r2, 3),
    paste0(round(oob_qq$obs_q[oob_qq$pct == 5],  3), " / ",
           round(oob_qq$pred_q[oob_qq$pct == 5], 3)),
    paste0(round(oob_qq$obs_q[oob_qq$pct == 95],  3), " / ",
           round(oob_qq$pred_q[oob_qq$pct == 95], 3)),
    sum(sign_stab$sign_sensitivity_class %in% c("Stable sink", "Stable source"))
  )),
  stringsAsFactors = FALSE
)

write.csv(stats_out, "OUTPUT/Flux_Justification_Supplement_stats.csv", row.names = FALSE)
message("Saved: OUTPUT/Flux_Justification_Supplement_stats.csv")

# Print key stats to console
cat("\n")
cat(strrep("=", 64), "\n")
cat("KEY STATS FOR REVIEWER RESPONSE\n")
cat(strrep("=", 64), "\n")
cat(sprintf(
"(a) Storage bias
    Median storage/total flux:  %.1f%%
    Max storage/total flux:     %.1f%%
    R2 (storage ~ height):      %.3f  (p = %.2f)
    Behavior unchanged (n/46):  %d

(b) Lateral mixing
    r(total, gradient):         %.4f
    Behavior unchanged (n/46):  %d

(c) Gap-fill damping
    n pairs:                    %s
    Mean bias:                  %+.6f mg C m-2 30min-1
    Sign accuracy:              %.1f%%
    OOB R2:                     %.3f
    p5  obs / pred:             %.3f / %.3f
    p95 obs / pred:             %.3f / %.3f
    Stable annual sign (n/46):  %d\n",
  net_median_storage_pct, max_storage_pct,
  r2_a, p_a,
  n_B - n_changed,
  r_b,
  n_B - n_changed,
  formatC(oob_n, big.mark = ",", format = "d"),
  oob_bias,
  100 * oob_sign,
  oob_r2,
  oob_qq$obs_q[oob_qq$pct == 5],  oob_qq$pred_q[oob_qq$pct == 5],
  oob_qq$obs_q[oob_qq$pct == 95], oob_qq$pred_q[oob_qq$pct == 95],
  sum(sign_stab$sign_sensitivity_class %in% c("Stable sink","Stable source"))
))
cat(strrep("=", 64), "\n")

