# Descriptive threshold sensitivity for the existing source-probability model.
# This does not recalibrate global prevalence from NEON. It compares
# area-weighted global annual classifications with leave-one-site-out NEON
# site-years, which is a transfer diagnostic rather than a global truth set.

library(data.table)
library(ggplot2)

rfdir <- Sys.getenv("MONTHLY_RF_DIR", unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF")
out <- file.path(rfdir, "OUTPUT")
figdir <- file.path(rfdir, "FIGURES")
thresholds <- seq(0.1, 0.9, by = 0.05)

loso <- fread(file.path(out, "source_probability_loso_predictions.csv"))
siteyear <- loso[, .(
  n_months = .N,
  observed_month_source_fraction = mean(observed_source),
  mean_probability = mean(p_balanced)
), by = .(SITE_ID, EcoType, Year)][n_months >= 10]
siteyear[, observed_source_year := observed_month_source_fraction >= 0.5]

global_file <- file.path(out, "monthly_cell_predictions_2000_2025.rds")
message("Reading global monthly predictions for annual cell summaries...")
x <- readRDS(global_file)
setDT(x)
annual_cell <- x[, .(
  area_mha = first(area_mha),
  mean_probability = mean(source_prob_A, na.rm = TRUE)
), by = .(Year, cell, EcoType, is_arid)]
annual_cell <- annual_cell[is_arid == 0]

global_long <- annual_cell[, {
  area <- sum(area_mha)
  .(threshold = thresholds,
    source_area_share = vapply(thresholds, function(t) sum(area_mha[mean_probability >= t]) / area, numeric(1)))
}, by = .(Year, EcoType)]
global_summary <- global_long[, .(
  estimate = mean(source_area_share),
  sd_annual = sd(source_area_share)
), by = .(EcoType, threshold)]
global_area <- annual_cell[, .(area_mha = sum(area_mha)), by = .(Year, EcoType)][,
  .(area_mha = mean(area_mha)), by = EcoType]

loso_summary <- siteyear[, {
  .(threshold = thresholds,
    predicted_source_year_fraction = vapply(thresholds, function(t) mean(mean_probability >= t), numeric(1)),
    observed_source_year_fraction = mean(observed_source_year),
    n_siteyears = .N,
    n_sites = uniqueN(SITE_ID))
}, by = EcoType]

result <- merge(global_summary, loso_summary, by = c("EcoType", "threshold"), all = TRUE)
result <- merge(result, global_area, by = "EcoType", all.x = TRUE)
setorder(result, EcoType, threshold)
fwrite(result, file.path(out, "source_probability_threshold_sensitivity.csv"))

plotdat <- melt(result,
  id.vars = c("EcoType", "threshold", "observed_source_year_fraction", "n_siteyears", "n_sites"),
  measure.vars = c("estimate", "predicted_source_year_fraction", "observed_source_year_fraction"),
  variable.name = "series", value.name = "source_fraction")
plotdat[, series := factor(series,
  levels = c("estimate", "predicted_source_year_fraction", "observed_source_year_fraction"),
  labels = c("Global area classified source", "LOSO predicted source site-years", "LOSO observed source site-years"))]
g <- ggplot(plotdat, aes(threshold, source_fraction, color = series, linetype = series)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(data = plotdat[series != "LOSO observed source site-years"], size = 1.5, na.rm = TRUE) +
  geom_hline(data = unique(plotdat[, .(EcoType, observed_source_year_fraction)]),
             aes(yintercept = observed_source_year_fraction), color = "grey35", linetype = "dotted", inherit.aes = FALSE) +
  facet_wrap(~EcoType, nrow = 1) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1), limits = c(0.1, 0.9)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom") +
  labs(title = "Source-area sensitivity to probability threshold",
       subtitle = "Global area-weighted annual cell probabilities vs leave-one-site-out NEON site-years",
       x = "Threshold for classifying a site-year / grid cell as source",
       y = "Fraction classified as source", color = NULL, linetype = NULL,
       caption = "Dotted line: observed NEON site-year source fraction. NEON is a transfer test, not a global prevalence sample.")
ggsave(file.path(figdir, "FIGURE_source_probability_threshold_sensitivity.png"), g,
       width = 11, height = 4.8, dpi = 250, bg = "white")

message("Saved threshold sensitivity CSV and figure; no production predictions changed.")
