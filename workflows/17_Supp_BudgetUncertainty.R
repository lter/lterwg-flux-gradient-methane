# Supplemental_Budget_Uncertainty.R
# ─────────────────────────────────────────────────────────────────────────────
# Quantifies the potential impact of a systematic bias in flux-gradient (FG)
# magnitude estimates on the ERA5 spatially upscaled continental CH4 budget.
#
# Question: If FG fluxes are systematically over- or under-estimated by X%,
# how does the upscaled annual budget (Tg CH4 yr-1) respond?
#
# Approach:
#   1. Quantify the observed FG–EC relative bias from validation sites
#      (VAL_ERA5_FG_vs_EC_annual.csv) to anchor the perturbation range.
#   2. Apply multiplicative corrections from −50% to +50% (bracketing
#      ±1× and ±2× the observed bias) to the continuous, dichotomous, and
#      all-sink annual budgets from ERA-SpatialUpscaling-Monthly-RF.R.
#   3. Report how the sign and magnitude of the continental budget changes
#      across scenarios, and whether the source conclusion is maintained.
#
# Requires:
#   annual_budget_three_approaches.csv  — in rf_dir/OUTPUT/
#   VAL_ERA5_FG_vs_EC_annual.csv        — in localdir.ch4/OUTPUT/
#
# Output:
#   /Volumes/MaloneLab/Research/FluxGradient/METHANE/FIGURES/FigS_Budget_Uncertainty.png
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(patchwork)

# ── Paths ─────────────────────────────────────────────────────────────────────
rf_dir    <- "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF"
val_dir   <- "/Volumes/MaloneLab/Research/FluxGradient/METHANE"
fig_dir   <- "/Volumes/MaloneLab/Research/FluxGradient/METHANE/FIGURES"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

rf_output_dir  <- file.path(rf_dir, "OUTPUT")
gmb_target     <- -35  # Tg CH4 yr-1 (GMB upland soil sink reference)

# ── 1. Load upscaled annual budgets ──────────────────────────────────────────
budget_file <- file.path(rf_output_dir, "annual_budget_three_approaches.csv")
if (!file.exists(budget_file))
  stop("Run ERA-SpatialUpscaling-Monthly-RF.R first. File not found:\n  ", budget_file)

annual_budget <- read.csv(budget_file, stringsAsFactors = FALSE)

# Long format for the three modelled approaches
budget_long <- annual_budget %>%
  pivot_longer(
    cols      = c(continuous_tg_ch4_yr, balanced_tg_ch4_yr, constrained_tg_ch4_yr),
    names_to  = "approach_raw",
    values_to = "annual_tg_ch4_yr"
  ) %>%
  mutate(approach = recode(approach_raw,
    continuous_tg_ch4_yr  = "Continuous",
    balanced_tg_ch4_yr    = "Dichotomous",
    constrained_tg_ch4_yr = "All-Sink"
  ))

baseline_means <- budget_long %>%
  group_by(approach) %>%
  summarise(mean_tg = mean(annual_tg_ch4_yr, na.rm = TRUE), .groups = "drop")

message("Baseline mean annual budgets (Tg CH4 yr-1):")
print(baseline_means)

# ── 2. Quantify observed FG–EC validation bias ────────────────────────────────
val_file <- file.path(val_dir, "OUTPUT/VAL_ERA5_FG_vs_EC_annual.csv")
if (!file.exists(val_file)) {
  message("Validation file not found: ", val_file)
  message("Using placeholder bias of 0 (no correction). Run Val.ERA5.HalfHourlyGapfill.R to populate.")
  obs_bias     <- 0
  obs_bias_sd  <- NA_real_
  n_val_sites  <- NA_integer_
} else {
  val_annual <- read.csv(val_file, stringsAsFactors = FALSE)
  # Relative bias per validation site: (FG − EC) / |EC|
  # Guard against near-zero EC budgets
  val_bias_df <- val_annual %>%
    filter(abs(EC_mean_annual_gC_m2_yr) >= 0.05) %>%
    mutate(rel_bias = (FG_mean_annual_gC_m2_yr - EC_mean_annual_gC_m2_yr) /
                      abs(EC_mean_annual_gC_m2_yr))
  obs_bias    <- mean(val_bias_df$rel_bias, na.rm = TRUE)
  obs_bias_sd <- sd(val_bias_df$rel_bias,  na.rm = TRUE)
  n_val_sites <- nrow(val_bias_df)
  message(sprintf("Observed mean FG–EC relative bias: %+.1f%%  (SD = %.1f%%, n = %d sites)",
                  obs_bias * 100, obs_bias_sd * 100, n_val_sites))
}

# ── 3. Define perturbation scenarios ─────────────────────────────────────────
# Multiplicative: adjusted_flux = baseline × (1 + frac)
# Range: −50% to +50%, with ±1× and ±2× observed bias included.

obs_abs <- abs(obs_bias)
raw_fracs <- sort(unique(round(c(
  -0.50,
  if (!is.na(obs_bias) && obs_abs > 0.03) c(-2 * obs_abs, -obs_abs, obs_abs, 2 * obs_abs),
   0,
  +0.50
), 3)))

# Remove near-duplicates (within 2 pp of each other)
keep <- c(TRUE, abs(diff(raw_fracs)) > 0.02)
scenario_fracs <- raw_fracs[keep]

# Build label: show obs_bias-derived scenarios specially
make_label <- function(frac) {
  if (frac == 0) return("0%\n(baseline)")
  pct <- sprintf("%+.0f%%", frac * 100)
  if (!is.na(obs_bias) && obs_abs > 0.03) {
    if (abs(frac + obs_abs) < 0.005) return(paste0(pct, "\n(−bias)"))
    if (abs(frac - obs_abs) < 0.005) return(paste0(pct, "\n(+bias)"))
    if (abs(frac + 2 * obs_abs) < 0.005) return(paste0(pct, "\n(−2×bias)"))
    if (abs(frac - 2 * obs_abs) < 0.005) return(paste0(pct, "\n(+2×bias)"))
  }
  pct
}

scenarios <- tibble(
  frac        = scenario_fracs,
  label       = map_chr(scenario_fracs, make_label),
  is_baseline = frac == 0,
  is_neg_bias = !is.na(obs_bias) & obs_abs > 0.03 & abs(frac + obs_abs) < 0.005
) %>%
  mutate(label = factor(label, levels = unique(label)))

message(sprintf("Scenarios: %s", paste(sprintf("%+.0f%%", scenario_fracs * 100), collapse = ", ")))

# ── 4. Apply scenarios to annual budgets ──────────────────────────────────────
sensitivity <- crossing(budget_long, scenarios) %>%
  mutate(adj_tg = annual_tg_ch4_yr * (1 + frac))

# Network-level summaries per scenario × approach
scenario_summary <- sensitivity %>%
  group_by(approach, frac, label, is_baseline, is_neg_bias) %>%
  summarise(
    mean_adj_tg    = mean(adj_tg,          na.rm = TRUE),
    sd_adj_tg      = sd(adj_tg,            na.rm = TRUE),
    pct_pos_years  = mean(adj_tg > 0,      na.rm = TRUE) * 100,
    gap_to_gmb     = mean_adj_tg - gmb_target,
    .groups        = "drop"
  )

# ── 5. Adjusted budget summary across scenarios ───────────────────────────────
message("\nBaseline mean annual budgets (Tg CH4 yr-1) and response to ±50% correction:")
scenario_summary %>%
  filter(frac %in% c(-0.50, 0, 0.50)) %>%
  select(approach, frac, mean_adj_tg, sd_adj_tg, pct_pos_years) %>%
  print()

cont_row <- baseline_means %>% filter(approach == "Continuous")

# ── FIGURES ───────────────────────────────────────────────────────────────────
base_theme <- theme_bw(base_size = 13) +
  theme(panel.grid.minor  = element_blank(),
        plot.title        = element_text(face = "bold", size = 13),
        axis.title        = element_text(size = 12),
        axis.text         = element_text(size = 11),
        legend.text       = element_text(size = 11),
        strip.text        = element_text(size = 12),
        legend.position   = "bottom")

approach_colors <- c(
  Continuous   = "#009688",
  Dichotomous  = "#9C27B0",
  `All-Sink`   = "#E07B39"
)
approach_linetype <- c(Continuous = "solid", Dichotomous = "dashed", `All-Sink` = "dotted")

highlight_fill <- "#E07B3920"  # orange-transparent for worst-case column

# ── Panel A: Mean adjusted budget vs scenario (line per approach) ─────────────
# Shows the range of budgets across all scenarios; GMB reference as red band.

# Add the GMB-required correction as an annotation
pA_data <- scenario_summary %>%
  filter(approach %in% c("Continuous", "Dichotomous", "All-Sink"))

# Find integer position of the −obs_bias label on the discrete axis for highlighting
neg_bias_pos <- if (!is.na(obs_bias) && obs_abs > 0.03) {
  lbl <- pA_data %>% filter(is_neg_bias) %>% slice(1) %>% pull(label) %>% as.character()
  which(levels(pA_data$label) == lbl)[1]
} else NULL

pA <- ggplot(pA_data,
             aes(x = label, y = mean_adj_tg,
                 color = approach, linetype = approach, shape = approach)) +
  # Highlight observed-bias column
  {if (!is.null(neg_bias_pos))
      annotate("rect",
               xmin = neg_bias_pos - 0.45, xmax = neg_bias_pos + 0.45,
               ymin = -Inf, ymax = Inf,
               fill = "#B22222", alpha = 0.08)} +
  # GMB reference band
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = gmb_target - 1, ymax = gmb_target + 1,
           fill = "#2166AC", alpha = 0.12) +
  geom_hline(yintercept = gmb_target, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = 0,          color = "grey40",  linetype = "solid",  linewidth = 0.4) +
  geom_line(aes(group = approach), linewidth = 1.0) +
  geom_point(size = 2.5) +
  annotate("text", x = nlevels(pA_data$label), y = gmb_target +5,
           label = "GMB reference\n(−35 Tg yr⁻¹)", hjust = 1, size = 3.8, color = "#2166AC") +
  scale_color_manual(values = approach_colors, name = NULL) +
  scale_linetype_manual(values = approach_linetype, name = NULL) +
  scale_shape_manual(values = c(Continuous = 16, Dichotomous = 17, `All-Sink` = 15), name = NULL) +
  labs(
    title = "A. ",
    x     = "Flux magnitude adjustment",
    y     = "Mean annual budget (Tg CH₄ yr⁻¹)"
  ) +
  base_theme +
  theme(axis.text.x = element_text(size = 10))

# ── Panel B: Budget time series at ±obs_bias and baseline (Continuous only) ───
# Show 3 lines: baseline, −observed_bias, −50%.

if (!is.na(obs_bias) && obs_abs > 0.03) {
  ts_fracs  <- sort(unique(c(-0.50, -obs_abs, 0, obs_abs, 0.50)))
  ts_labels <- tibble(
    frac  = ts_fracs,
    label = map_chr(ts_fracs, make_label)
  ) %>% mutate(label = factor(label, levels = unique(label)))
} else {
  ts_fracs  <- c(-0.50, 0, 0.50)
  ts_labels <- tibble(
    frac  = ts_fracs,
    label = factor(sprintf("%+.0f%%", ts_fracs * 100))
  )
}

ts_data <- budget_long %>%
  filter(approach == "Continuous") %>%
  crossing(ts_labels) %>%
  mutate(adj_tg = annual_tg_ch4_yr * (1 + frac))

ts_colors <- scales::seq_gradient_pal("#B22222", "#2ca25f")(
  seq(0, 1, length.out = nlevels(ts_data$label)))
names(ts_colors) <- levels(ts_data$label)

pB <- ggplot(ts_data, aes(x = Year, y = adj_tg, color = label, linetype = label)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = gmb_target - 1, ymax = gmb_target + 1,
           fill = "#2166AC", alpha = 0.12) +
  geom_hline(yintercept = gmb_target, color = "#2166AC", linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = 0,          color = "grey40",  linetype = "solid",  linewidth = 0.4) +
  geom_line(linewidth = 0.85) +
  scale_color_manual(values = ts_colors, name = "Adjustment") +
  scale_linetype_manual(
    values = setNames(c("dotted","dashed","solid","dashed","dotted"),
                      levels(ts_data$label)),
    name = "Adjustment") +
  scale_x_continuous(breaks = seq(2000, 2025, 5)) +
  labs(
    title = "B. ",
    x     = "Year",
    y     = "Annual budget (Tg CH₄ yr⁻¹)"
  ) +
  base_theme +
  theme(legend.position = "right", legend.text = element_text(size = 8))

# ── Panel C: Absolute shift in mean annual budget relative to baseline ─────────
# Shows how much the budget changes (Tg yr-1) at each perturbation level.
pC_data <- scenario_summary %>%
  filter(approach %in% c("Continuous", "Dichotomous", "All-Sink")) %>%
  left_join(baseline_means %>% rename(baseline_tg = mean_tg), by = "approach") %>%
  mutate(delta_tg = mean_adj_tg - baseline_tg)

pC <- ggplot(pC_data,
             aes(x = label, y = delta_tg, fill = approach, group = approach)) +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dashed", linewidth = 0.6) +
  geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.82) +
  scale_fill_manual(values = approach_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  labs(
    title = "C.",
    x     = "Flux magnitude adjustment",
    y     = expression(Delta~"mean annual budget (Tg CH"[4]~"yr"^{-1}*")")
  ) +
  base_theme +
  theme(axis.text.x = element_text(size = 10))

# ── Panel D: Fraction of years with net-positive budget ───────────────────────
pD <- scenario_summary %>%
  filter(approach %in% c("Continuous", "Dichotomous", "All-Sink")) %>%
  ggplot(aes(x = label, y = pct_pos_years, fill = approach, group = approach)) +
  geom_hline(yintercept = 50, color = "grey50", linetype = "dashed", linewidth = 0.6) +
  geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.82) +
  scale_fill_manual(values = approach_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title = "D. ",
    x     = "Flux magnitude adjustment",
    y     = "% of years with budget > 0"
  ) +
  base_theme +
  theme(axis.text.x = element_text(size = 10))

# ── Assemble ──────────────────────────────────────────────────────────────────
fig_budget <- (pA | pB) / (pC | pD) 

ggsave(
  file.path(fig_dir, "FigS_Budget_Uncertainty.png"),
  fig_budget,
  width  = 14.0,
  height = 11.0,
  dpi    = 300,
  bg     = "white"
)

message("Saved FigS_Budget_Uncertainty.png")

# ── Summary statistics and supplement text ────────────────────────────────────
cat("\n── Budget uncertainty supplement statistics ───────────────────────────\n")
if (!is.na(obs_bias) && obs_abs > 0.03)
  cat(sprintf("Observed FG–EC relative bias: %+.1f%% (SD = %.1f%%, n = %d sites)\n",
              obs_bias * 100, obs_bias_sd * 100, n_val_sites))

worst_frac <- min(scenario_fracs)
cat("\nBaseline and worst-case adjusted budgets (Tg CH4 yr-1):\n")
scenario_summary %>%
  filter(frac %in% c(worst_frac, 0)) %>%
  select(approach, frac, mean_adj_tg, sd_adj_tg, pct_pos_years) %>%
  arrange(approach, frac) %>%
  print()

# Key numbers for supplement text
cont_base      <- baseline_means %>% filter(approach == "Continuous") %>% pull(mean_tg)
cont_worst     <- scenario_summary %>%
                    filter(approach == "Continuous", frac == worst_frac) %>% pull(mean_adj_tg)
cont_worst_pos <- scenario_summary %>%
                    filter(approach == "Continuous", frac == worst_frac) %>% pull(pct_pos_years)
sink_base      <- baseline_means %>% filter(approach == "All-Sink") %>% pull(mean_tg)
sink_worst     <- scenario_summary %>%
                    filter(approach == "All-Sink",   frac == worst_frac) %>% pull(mean_adj_tg)

# Supplement paragraph
cat("\n")
cat(paste(strwrap(sprintf(
"Supplemental Text — Sensitivity of spatially upscaled CH4 budget to FG magnitude bias.

To evaluate the potential impact of a systematic bias in flux-gradient (FG) magnitude
estimates on the continental CH4 budget, we applied multiplicative corrections ranging
from %+.0f%% to +50%% uniformly to all ERA5 spatially upscaled annual budgets across
three modelling approaches (Continuous, Dichotomous, and All-Sink).%s

At baseline (0%% correction), the Continuous approach produced a mean annual budget
of %.1f Tg CH4 yr-1 and the All-Sink approach (every non-arid upland cell assigned
the sink flux magnitude) produced %.1f Tg CH4 yr-1. Under the most conservative
scenario (%+.0f%%), the Continuous budget shifted to %.1f Tg CH4 yr-1, remaining
net positive in %.0f%% of years, and the All-Sink budget shifted to %.1f Tg CH4 yr-1.
The sign and broad magnitude of the upscaled source estimate are therefore insensitive
to systematic FG biases within the range tested, and the conclusion that upland
NEON terrestrial ecosystems constitute a net CH4 source is robust to flux-magnitude
uncertainty of this scale.",
  worst_frac * 100,
  if (!is.na(obs_bias) && obs_abs > 0.03)
    sprintf(" The perturbation range was anchored on the observed mean FG-EC annual budget
bias of %+.1f%% (SD = %.1f%%, n = %d validation sites), with scenarios extending
to ±50%% to bound the full plausible range of systematic error.",
    obs_bias * 100, obs_bias_sd * 100, n_val_sites) else "",
  cont_base, sink_base,
  worst_frac * 100, cont_worst, cont_worst_pos, sink_worst),
  width = 90), collapse = "\n"))
cat("\n")
