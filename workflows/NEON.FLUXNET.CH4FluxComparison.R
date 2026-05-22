# Compare NEON daily CH4 fluxes with published FLUXNET-CH4 ecosystem-class rates.
#
# Published ecosystem-class values are annual CH4 fluxes from Delwiche et al. (2021),
# Table 1, converted from g C m-2 yr-1 to mg C m-2 d-1.
# NEON values are site-level medians of lookup-filled daily total fluxes from
# OUTPUT/NEON_scale_daily_flux_all_sites.csv. Daily rates are built on a
# complete 48 half-hour grid and reported as mg C m-2 d-1.

library(tidyverse)
library(ggplot2)
library(patchwork)
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

model_data_file <- "OUTPUT/30min_ch4_model_data.csv"
site_behavior_file <- "OUTPUT/30min_site_behavior.csv"
adjusted_daily_flux_file <- "OUTPUT/NEON_scale_daily_flux_all_sites.csv"

if (!file.exists(site_behavior_file)) {
  stop("Missing OUTPUT/30min_site_behavior.csv. Run flow.30min.analysis.R first.")
}

if (!file.exists(adjusted_daily_flux_file)) {
  stop(
    "Missing ", adjusted_daily_flux_file, ". Run flow.30min.analysis.R first ",
    "to create lookup-filled daily fluxes."
  )
}

behavior_levels <- c("Consistent sink", "Fluctuating", "Consistent source")
behavior_colors <- c(
  "Consistent sink" = "red3",
  "Fluctuating" = "grey35",
  "Consistent source" = "blue4"
)

base_plot_size <- 13
panel_title_size <- 15
axis_title_size <- 13
axis_text_size <- 11.5
legend_text_size <- 12
legend_title_size <- 13

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
    daily_mgC_m2_day = -uptake_mgCH4_m2_yr * 12 / 16 / 365,
    daily_sd_mgC_m2_day = uptake_sd_mgCH4_m2_yr * 12 / 16 / 365,
    daily_low_mgC_m2_day = daily_mgC_m2_day - daily_sd_mgC_m2_day,
    daily_high_mgC_m2_day = daily_mgC_m2_day + daily_sd_mgC_m2_day
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

site_behavior <- read.csv(site_behavior_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    CH4_behavior = factor(CH4_behavior, levels = behavior_levels)
  )

if (file.exists(adjusted_daily_flux_file)) {
  daily_neon <- read.csv(adjusted_daily_flux_file) %>%
    mutate(
      SITE_ID = as.character(SITE_ID),
      Date = as.Date(Date)
    )
}

daily_neon <- daily_neon %>%
  left_join(
    site_behavior %>% dplyr::select(SITE_ID, CH4_behavior),
    by = "SITE_ID"
  ) %>%
  filter(
    is.finite(daily_mgC_m2_day),
    !is.na(Date),
    !is.na(CH4_behavior)
  )

write.csv(daily_neon, "OUTPUT/NEON_daily_CH4_flux_rates_for_FLUXNET_comparison.csv", row.names = FALSE)

neon_site_summary <- daily_neon %>%
  group_by(SITE_ID, CH4_behavior) %>%
  summarise(
    n_days = n(),
    median_n_30min_per_day = median(n_30min, na.rm = TRUE),
    min_n_30min_per_day = min(n_30min, na.rm = TRUE),
    daily_mgC_m2_day = median(daily_mgC_m2_day, na.rm = TRUE),
    daily_q25_mgC_m2_day = quantile(daily_mgC_m2_day, 0.25, na.rm = TRUE),
    daily_q75_mgC_m2_day = quantile(daily_mgC_m2_day, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    neon_class = recode(
      as.character(CH4_behavior),
      "Consistent sink" = "NEON strong sink",
      "Fluctuating" = "NEON fluctuating",
      "Consistent source" = "NEON source"
    )
  )

neon_class_summary <- neon_site_summary %>%
  group_by(neon_class, CH4_behavior) %>%
  summarise(
    comparison_group = "NEON half-hour total flux",
    ecosystem_class = first(neon_class),
    n_sites = n_distinct(SITE_ID),
    n_site_years = NA_integer_,
    daily_median_mgC_m2_day = median(daily_mgC_m2_day, na.rm = TRUE),
    daily_low_mgC_m2_day = quantile(daily_mgC_m2_day, 0.25, na.rm = TRUE),
    daily_high_mgC_m2_day = quantile(daily_mgC_m2_day, 0.75, na.rm = TRUE),
    daily_sd_mgC_m2_day = sd(daily_mgC_m2_day, na.rm = TRUE),
    reference_id = "This analysis: median lookup-filled daily total fluxes from OUTPUT/NEON_scale_daily_flux_all_sites.csv",
    source_type = "NEON behavior class",
    .groups = "drop"
  ) %>%
  mutate(
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
  neon_class_summary %>%
    dplyr::select(
      comparison_group, ecosystem_class, CH4_behavior, annual_gC_m2_yr,
      annual_sd_gC_m2_yr, daily_mgC_m2_day, daily_low_mgC_m2_day,
      daily_high_mgC_m2_day, n_sites, n_site_years, reference_id, source_type
    )
) %>%
  arrange(desc(daily_mgC_m2_day)) %>%
  mutate(
    plot_rank = rev(seq_len(dplyr::n())),
    source_type = factor(source_type, levels = c("Published ecosystem class", "Soil chamber literature", "NEON behavior class"))
  )

write.csv(comparison_table, "OUTPUT/CH4_flux_FLUXNET_NEON_comparison_values.csv", row.names = FALSE)
write.csv(neon_site_summary, "OUTPUT/NEON_site_median_daily_CH4_flux_classes.csv", row.names = FALSE)
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
  "## NEON values",
  "NEON values are from this repository's lookup-filled daily total flux output: OUTPUT/NEON_scale_daily_flux_all_sites.csv and OUTPUT/30min_site_behavior.csv.",
  "Daily NEON fluxes were computed on a complete 48 half-hour grid. Observed half-hours were retained; missing half-hours were filled from site month-hour, season-hour, biseason-hour, then annual-hour mean rates before summing to mg C m-2 d-1.",
  "NEON points in the figure are site-level medians of those lookup-filled daily total fluxes. NEON behavior classes are from OUTPUT/30min_site_behavior.csv: Consistent sink, Fluctuating, and Consistent source.",
  "",
  "## Output tables",
  "- OUTPUT/CH4_flux_FLUXNET_NEON_comparison_values.csv",
  "- OUTPUT/NEON_daily_CH4_flux_rates_for_FLUXNET_comparison.csv",
  "- OUTPUT/NEON_site_median_daily_CH4_flux_classes.csv",
  "- OUTPUT/process_model_upland_CH4_uptake_values.csv",
  "- OUTPUT/soil_chamber_CH4_flux_reference_values.csv"
)
writeLines(reference_lines, "OUTPUT/CH4_flux_FLUXNET_NEON_comparison_references.md")

plot_full <- comparison_table %>%
  ggplot(aes(x = daily_mgC_m2_day, y = plot_rank)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  geom_errorbar(
    data = comparison_table %>% filter(source_type == "Published ecosystem class"),
    aes(xmin = daily_low_mgC_m2_day, xmax = daily_high_mgC_m2_day),
    orientation = "y",
    width = 0.18,
    color = "grey35",
    linewidth = 0.65
  ) +
  geom_point(
    data = comparison_table %>% filter(source_type == "Published ecosystem class"),
    shape = 21,
    fill = "white",
    color = "black",
    size = 3.2,
    stroke = 0.8
  ) +
  geom_errorbar(
    data = comparison_table %>% filter(source_type == "Soil chamber literature"),
    aes(xmin = daily_low_mgC_m2_day, xmax = daily_high_mgC_m2_day),
    orientation = "y",
    width = 0.18,
    color = "#238B45",
    linewidth = 0.75
  ) +
  geom_point(
    data = comparison_table %>% filter(source_type == "Soil chamber literature"),
    shape = 24,
    fill = "#A1D99B",
    color = "#238B45",
    size = 3.4,
    stroke = 0.9
  ) +
  geom_errorbar(
    data = comparison_table %>% filter(source_type == "NEON behavior class"),
    aes(xmin = daily_low_mgC_m2_day, xmax = daily_high_mgC_m2_day, color = CH4_behavior),
    orientation = "y",
    width = 0.22,
    linewidth = 0.75
  ) +
  geom_point(
    data = neon_site_summary %>%
      left_join(
        comparison_table %>% dplyr::select(ecosystem_class, plot_rank),
        by = c("neon_class" = "ecosystem_class")
      ),
    aes(x = daily_mgC_m2_day, y = plot_rank, color = CH4_behavior),
    alpha = 0.78,
    size = 2.2,
    position = position_jitter(height = 0.12, width = 0)
  ) +
  geom_point(
    data = comparison_table %>% filter(source_type == "NEON behavior class"),
    aes(color = CH4_behavior),
    size = 3.8
  ) +
  scale_color_manual(values = behavior_colors, name = "NEON class", na.translate = FALSE) +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 0.01),
    breaks = c(-20, -5, -1, -0.1, 0, 0.1, 1, 10, 75, 100),
    labels = c("-20", "-5", "-1", "-0.1", "0", "0.1", "1", "10", "75", "100")
  ) +
  scale_y_continuous(
    breaks = comparison_table$plot_rank,
    labels = comparison_table$ecosystem_class,
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  labs(
    title = "A. Published ecosystem and soil-chamber CH4 fluxes ranked high to low, with NEON classes overlaid",
    x = expression("Daily CH"[4] * " flux (mg C m"^-2 * " d"^-1 * "; pseudo-log scale)"),
    y = NULL
  ) +
  theme_bw(base_size = base_plot_size) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = panel_title_size),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    legend.title = element_text(size = legend_title_size),
    legend.text = element_text(size = legend_text_size),
    panel.grid.minor = element_blank()
  )

plot_process_model <- process_model_uptake %>%
  ggplot(aes(x = daily_mgC_m2_day, y = plot_rank)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  geom_errorbar(
    aes(xmin = daily_low_mgC_m2_day, xmax = daily_high_mgC_m2_day),
    orientation = "y",
    width = 0.18,
    linewidth = 0.65,
    color = "grey35"
  ) +
  geom_point(shape = 21, fill = "white", color = "black", size = 3.1, stroke = 0.8) +
  geom_point(
    data = neon_site_summary,
    aes(x = daily_mgC_m2_day, y = 0, color = CH4_behavior),
    inherit.aes = FALSE,
    show.legend = FALSE,
    alpha = 0.75,
    size = 2.1,
    position = position_jitter(height = 0.13, width = 0)
  ) +
  scale_y_continuous(
    breaks = c(process_model_uptake$plot_rank, 0),
    labels = c(process_model_uptake$ecosystem_class, "NEON site medians"),
    expand = expansion(mult = c(0.05, 0.04))
  ) +
  scale_color_manual(values = behavior_colors, guide = "none", na.translate = FALSE) +
  labs(
    title = "B. NEON compared with process-based upland soil CH4 uptake estimates",
    x = expression("Daily CH"[4] * " flux (mg C m"^-2 * " d"^-1 * "; negative = uptake)"),
    y = NULL
  ) +
  theme_bw(base_size = base_plot_size) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = panel_title_size),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    panel.grid.minor = element_blank()
  )

comparison_figure <- plot_full / plot_process_model +
  plot_layout(heights = c(1.45, 1.2), guides = "collect") +
  plot_annotation(
    title = "NEON CH4 fluxes are weak relative to published FLUXNET-CH4 ecosystem classes",
    subtitle = paste(
      "Published FLUXNET-CH4 values are ecosystem-class annual means +/- SD converted to daily units;",
      "soil-chamber values are literature benchmarks; process-model values are MeMo v1.0 uptake means +/- SD by ecosystem type."
    ),
    caption = paste(
      "Negative values indicate CH4 uptake. NEON absolute magnitudes depend on the repository conversion to CH4_mgC_30min" ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 12.5),
      plot.caption = element_text(size = 10, color = "grey25")
    )
  )

comparison_figure <- comparison_figure &
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.title = element_text(size = legend_title_size),
    legend.text = element_text(size = legend_text_size)
  )

ggsave(
  "FIGURES/NEON_FLUXNET_CH4_flux_comparison.png",
  plot = comparison_figure,
  width = 16,
  height = 10.5,
  units = "in",
  dpi = 300
)

ggsave(
  "FIGURES/NEON_FLUXNET_CH4_flux_comparison.pdf",
  plot = comparison_figure,
  width = 16,
  height = 10.5,
  units = "in"
)

message("Wrote FIGURES/NEON_FLUXNET_CH4_flux_comparison.png and .pdf")
