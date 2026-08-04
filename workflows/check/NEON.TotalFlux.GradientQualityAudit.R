# Audit whether retained CH4 gradients are physically trustworthy.

library(tidyverse)
library(ggplot2)

script_arg <- commandArgs(FALSE) %>%
  str_subset("^--file=") %>%
  str_remove("^--file=") %>%
  first()

script_dir <- if (is.na(script_arg)) {
  getwd()
} else {
  dirname(normalizePath(script_arg, mustWork = FALSE))
}

output_dir <- file.path(script_dir, "OUTPUT")
figure_dir <- file.path(script_dir, "FIGURES")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

localdir_ch4_candidates <- c(
  "/Volumes/MaloneLab/Research/FluxGradient/Methane",
  "/Volumes/MaloneLab/Research/FluxGradient/METHANE"
)

filtered_file <- file.path(
  localdir_ch4_candidates,
  "NEON_GradientFlux_Data_Filter/SITE_DATA_FILTERED_CH4.Rdata"
) %>%
  .[file.exists(.)] %>%
  first()

total_flux_file <- file.path(
  localdir_ch4_candidates,
  "SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata"
) %>%
  .[file.exists(.)] %>%
  first()

if (is.na(filtered_file)) {
  stop(
    "Could not find SITE_DATA_FILTERED_CH4.Rdata in: ",
    paste(localdir_ch4_candidates, collapse = ", ")
  )
}

if (is.na(total_flux_file)) {
  stop(
    "Could not find SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata in: ",
    paste(localdir_ch4_candidates, collapse = ", ")
  )
}

load(filtered_file)
load(total_flux_file)

if (!exists("SITE_DATA_FILTERED")) {
  stop("Expected object SITE_DATA_FILTERED was not found in filtered CH4 RData.")
}

if (!exists("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage")) {
  stop("Expected object SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage was not found in total-flux RData.")
}

as_flagged <- function(x) {
  ifelse(is.na(x), NA, as.character(x) != "0" & tolower(as.character(x)) != "false")
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

make_gradient_key <- function(df, site, approach) {
  time_col <- intersect(c("timeEnd_A", "timeEndA", "timeEndA.local", "timeEnd.local"), names(df)) %>%
    first()

  if (is.na(time_col)) {
    time_value <- seq_len(nrow(df))
  } else {
    time_value <- as.character(df[[time_col]])
  }

  paste(
    site,
    approach,
    time_value,
    as.character(df$TowerPosition_A),
    as.character(df$TowerPosition_B),
    as.character(df$dLevelsAminusB),
    as.character(df$gas),
    sep = "|"
  )
}

safe_load_rdata <- function(file, env) {
  tryCatch(
    {
      load(file, envir = env)
      NULL
    },
    error = function(e) {
      conditionMessage(e)
    }
  )
}

prep_gradient_rows <- function(df, site, approach, stage) {
  needed <- c("gas", "dConc", "FG_mean", "TowerPosition_A", "TowerPosition_B", "dLevelsAminusB")
  missing <- setdiff(needed, names(df))

  if (length(missing) > 0) {
    return(tibble(
      SITE_ID = site,
      Approach = approach,
      stage = stage,
      missing_cols = paste(missing, collapse = ", ")
    ))
  }

  optional_numeric <- c(
    "dConcSNR", "dConc_sd", "dConc_pvalue", "mean_A", "mean_B",
    "TowerHeight_A", "TowerHeight_B", "cross_grad_flag"
  )

  for (col in optional_numeric) {
    if (!col %in% names(df)) {
      df[[col]] <- NA_real_
    }
  }

  df %>%
    filter(gas == "CH4", is.finite(dConc), dConc != 0) %>%
    mutate(
      SITE_ID = site,
      Approach = approach,
      stage = stage,
      gradient_key = make_gradient_key(pick(everything()), site, approach),
      TowerPosition_A_num = safe_numeric(TowerPosition_A),
      TowerPosition_B_num = safe_numeric(TowerPosition_B),
      tower_order = case_when(
        TowerPosition_A_num > TowerPosition_B_num ~ "A above B",
        TowerPosition_A_num < TowerPosition_B_num ~ "A below B",
        TRUE ~ "unknown"
      ),
      pair_span = abs(TowerPosition_A_num - TowerPosition_B_num),
      dConc_direction = case_when(
        dConc < 0 ~ "negative dConc/source-like",
        dConc > 0 ~ "positive dConc/sink-like",
        TRUE ~ "zero"
      ),
      FG_direction = case_when(
        is.finite(FG_mean) & FG_mean > 0 ~ "positive FG/source",
        is.finite(FG_mean) & FG_mean < 0 ~ "negative FG/sink",
        TRUE ~ "missing or zero FG"
      ),
      cross_grad_flagged = if ("cross_grad_flag" %in% names(.)) as_flagged(cross_grad_flag) else NA,
      dConcSNR_calc = case_when(
        is.finite(dConcSNR) ~ dConcSNR,
        is.finite(dConc_sd) & dConc_sd > 0 ~ abs(dConc) / dConc_sd,
        TRUE ~ NA_real_
      ),
      missing_cols = NA_character_
    ) %>%
    select(
      SITE_ID, Approach, stage, gradient_key, dLevelsAminusB, tower_order, pair_span,
      TowerPosition_A, TowerPosition_B, TowerHeight_A, TowerHeight_B,
      dConc, dConc_sd, dConcSNR_calc, dConc_pvalue, mean_A, mean_B, FG_mean,
      dConc_direction, FG_direction, cross_grad_flagged, missing_cols
    )
}

retained_compact <- purrr::imap_dfr(SITE_DATA_FILTERED, function(site_df, site) {
  needed <- c("FG_mean", "gas", "TowerPosition_A", "TowerPosition_B", "dLevelsAminusB", "Approach")
  missing <- setdiff(needed, names(site_df))

  if (length(missing) > 0) {
    return(tibble(SITE_ID = site, missing_cols = paste(missing, collapse = ", ")))
  }

  site_df %>%
    filter(gas == "CH4", is.finite(FG_mean), FG_mean != 0) %>%
    mutate(
      SITE_ID = site,
      TowerPosition_A_num = safe_numeric(TowerPosition_A),
      TowerPosition_B_num = safe_numeric(TowerPosition_B),
      tower_order = case_when(
        TowerPosition_A_num > TowerPosition_B_num ~ "A above B",
        TowerPosition_A_num < TowerPosition_B_num ~ "A below B",
        TRUE ~ "unknown"
      ),
      pair_span = abs(TowerPosition_A_num - TowerPosition_B_num),
      FG_direction = case_when(
        FG_mean > 0 ~ "positive FG/source",
        FG_mean < 0 ~ "negative FG/sink",
        TRUE ~ "missing or zero FG"
      ),
      cross_grad_flagged = if ("cross_grad_flag" %in% names(.)) as_flagged(cross_grad_flag) else NA,
      cross_grad_class = case_when(
        is.na(cross_grad_flagged) ~ "missing",
        cross_grad_flagged ~ "flagged",
        TRUE ~ "not flagged"
      ),
      missing_cols = NA_character_
    ) %>%
    select(
      SITE_ID, Approach, Canopy_L1, dLevelsAminusB, tower_order, pair_span,
      TowerPosition_A, TowerPosition_B, FG_mean, FG_direction,
      cross_grad_flagged, cross_grad_class, missing_cols
    )
})

retained_quality_by_site <- retained_compact %>%
  filter(is.na(missing_cols)) %>%
  reframe(
    .by = SITE_ID,
    n_retained = n(),
    n_pairs = n_distinct(paste(Approach, dLevelsAminusB)),
    prop_positive_FG = mean(FG_mean > 0, na.rm = TRUE),
    prop_cross_grad_flagged = mean(cross_grad_flagged %in% TRUE, na.rm = TRUE),
    prop_cross_grad_missing = mean(is.na(cross_grad_flagged)),
    median_pair_span = median(pair_span, na.rm = TRUE),
    dominant_pair_fraction = {
      pair_counts <- table(paste(Approach, dLevelsAminusB))
      max(pair_counts) / sum(pair_counts)
    }
  ) %>%
  arrange(desc(prop_positive_FG))

cross_gradient_summary <- retained_compact %>%
  filter(is.na(missing_cols)) %>%
  reframe(
    .by = c(Approach, cross_grad_class),
    n = n(),
    prop_positive_FG = mean(FG_mean > 0, na.rm = TRUE),
    median_FG = median(FG_mean, na.rm = TRUE),
    median_abs_FG = median(abs(FG_mean), na.rm = TRUE),
    n_sites = n_distinct(SITE_ID)
  ) %>%
  arrange(Approach, cross_grad_class)

sensor_pair_summary <- retained_compact %>%
  filter(is.na(missing_cols)) %>%
  reframe(
    .by = c(SITE_ID, Approach, Canopy_L1, dLevelsAminusB, tower_order, pair_span),
    n = n(),
    prop_positive_FG = mean(FG_mean > 0, na.rm = TRUE),
    prop_cross_grad_flagged = mean(cross_grad_flagged %in% TRUE, na.rm = TRUE),
    median_FG = median(FG_mean, na.rm = TRUE),
    median_abs_FG = median(abs(FG_mean), na.rm = TRUE)
  ) %>%
  arrange(desc(prop_positive_FG), desc(n))

total_flux_summary <- purrr::imap_dfr(
  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage,
  function(site_df, site) {
    site_df %>%
      filter(is.finite(flux_total)) %>%
      reframe(
        SITE_ID = site,
        n_total_flux = n(),
        prop_positive_total_flux = mean(flux_total > 0, na.rm = TRUE),
        median_total_flux = median(flux_total, na.rm = TRUE)
      )
  }
) %>%
  arrange(desc(prop_positive_total_flux))

# Raw evaluation files can be hundreds of MB on a network volume. Keep the
# default small; set NEON_GRADIENT_AUDIT_MAX_RAW_SITES for a broader audit.
max_raw_sites <- Sys.getenv("NEON_GRADIENT_AUDIT_MAX_RAW_SITES", "1") %>%
  as.integer()

sample_sites <- total_flux_summary %>%
  filter(n_total_flux > 0) %>%
  summarise(
    sites = list(unique(c(
      "KONZ",
      head(SITE_ID, 12),
      tail(SITE_ID, 8),
      retained_quality_by_site %>% slice_max(prop_cross_grad_flagged, n = 6, with_ties = FALSE) %>% pull(SITE_ID),
      "BARR"
    )))
  ) %>%
  pull(sites) %>%
  .[[1]] %>%
  na.omit() %>%
  unique() %>%
  head(max_raw_sites)

load_site_filter_rows <- function(site) {
  message("Loading retained filter file for ", site)

  filter_file_candidates <- file.path(
    c(
      "/Volumes/MaloneLab/Research/FluxGradient/NEON_GradientFlux_Data_Filter",
      "/Volumes/MaloneLab/Research/FluxGradient/METHANE/NEON_GradientFlux_Data_Filter",
      "/Volumes/MaloneLab/Research/FluxGradient/Methane/NEON_GradientFlux_Data_Filter"
    ),
    site,
    paste0(site, c("_FILTER_AA_AW.Rdata", "_FILTER.Rdata"))
  )

  filter_file <- filter_file_candidates[file.exists(filter_file_candidates)] %>% first()

  if (is.na(filter_file)) {
    return(tibble(SITE_ID = site, stage = "retained", missing_cols = "site filter file not found"))
  }

  env <- new.env(parent = emptyenv())
  load_error <- safe_load_rdata(filter_file, env)

  if (!is.null(load_error)) {
    return(tibble(
      SITE_ID = site,
      stage = "retained",
      missing_cols = paste("could not load site filter file:", load_error)
    ))
  }

  purrr::map_dfr(
    c("MBR_9min_FILTER", "AE_9min_FILTER", "WP_9min_FILTER"),
    function(object_name) {
      if (!exists(object_name, envir = env)) {
        return(tibble())
      }

      prep_gradient_rows(
        get(object_name, envir = env),
        site,
        str_remove(object_name, "_9min_FILTER"),
        "retained"
      )
    }
  )
}

load_site_evaluation_rows <- function(site) {
  message("Loading raw evaluation file for ", site)

  eval_file_candidates <- file.path(
    c(
      "/Volumes/MaloneLab/Research/FluxGradient/NEON_GradientFlux_Data",
      "/Volumes/MaloneLab/Research/FluxGradient/METHANE/NEON_GradientFlux_Data",
      "/Volumes/MaloneLab/Research/FluxGradient/Methane/NEON_GradientFlux_Data"
    ),
    site,
    paste0(site, "_Evaluation.Rdata")
  )

  eval_file <- eval_file_candidates[file.exists(eval_file_candidates)] %>% first()

  if (is.na(eval_file)) {
    return(tibble(SITE_ID = site, stage = "all_evaluable", missing_cols = "site evaluation file not found"))
  }

  env <- new.env(parent = emptyenv())
  load_error <- safe_load_rdata(eval_file, env)

  if (!is.null(load_error)) {
    return(tibble(
      SITE_ID = site,
      stage = "all_evaluable",
      missing_cols = paste("could not load site evaluation file:", load_error)
    ))
  }

  purrr::map_dfr(
    c("MBR_9min.df.final", "AE_9min.df.final", "WP_9min.df.final"),
    function(object_name) {
      if (!exists(object_name, envir = env)) {
        return(tibble())
      }

      prep_gradient_rows(
        get(object_name, envir = env),
        site,
        str_remove(object_name, "_9min.df.final"),
        "all_evaluable"
      )
    }
  )
}

message("Loading sampled raw evaluation and retained filter files for ", length(sample_sites), " sites.")
all_evaluable_raw <- purrr::map_dfr(sample_sites, load_site_evaluation_rows)
retained_raw <- purrr::map_dfr(sample_sites, load_site_filter_rows)

raw_load_issues <- bind_rows(
  all_evaluable_raw %>%
    filter(!is.na(missing_cols)) %>%
    distinct(SITE_ID, stage, missing_cols),
  retained_raw %>%
    filter(!is.na(missing_cols)) %>%
    distinct(SITE_ID, stage, missing_cols)
)

retained_keys <- retained_raw %>%
  filter(is.na(missing_cols)) %>%
  distinct(gradient_key) %>%
  pull(gradient_key)

filter_bias_rows <- all_evaluable_raw %>%
  filter(is.na(missing_cols)) %>%
  mutate(retained_by_filter = gradient_key %in% retained_keys)

dconc_filter_bias_summary <- filter_bias_rows %>%
  reframe(
    .by = c(SITE_ID, Approach, dLevelsAminusB),
    n_evaluable = n(),
    n_retained = sum(retained_by_filter, na.rm = TRUE),
    prop_negative_dConc_evaluable = mean(dConc < 0, na.rm = TRUE),
    prop_negative_dConc_retained = ifelse(
      sum(retained_by_filter, na.rm = TRUE) > 0,
      mean(dConc[retained_by_filter] < 0, na.rm = TRUE),
      NA_real_
    ),
    retention_rate_negative_dConc = mean(retained_by_filter[dConc < 0], na.rm = TRUE),
    retention_rate_positive_dConc = mean(retained_by_filter[dConc > 0], na.rm = TRUE),
    retention_rate_ratio_negative_to_positive = retention_rate_negative_dConc / retention_rate_positive_dConc,
    median_abs_dConc_evaluable = median(abs(dConc), na.rm = TRUE),
    median_abs_dConc_retained = ifelse(
      sum(retained_by_filter, na.rm = TRUE) > 0,
      median(abs(dConc[retained_by_filter]), na.rm = TRUE),
      NA_real_
    ),
    median_dConcSNR_evaluable = median(dConcSNR_calc, na.rm = TRUE),
    median_dConcSNR_retained = ifelse(
      sum(retained_by_filter, na.rm = TRUE) > 0,
      median(dConcSNR_calc[retained_by_filter], na.rm = TRUE),
      NA_real_
    )
  ) %>%
  filter(n_evaluable > 0) %>%
  arrange(desc(retention_rate_ratio_negative_to_positive), desc(n_retained))

dconc_filter_bias_overall <- filter_bias_rows %>%
  reframe(
    .by = Approach,
    n_evaluable = n(),
    n_retained = sum(retained_by_filter, na.rm = TRUE),
    prop_negative_dConc_evaluable = mean(dConc < 0, na.rm = TRUE),
    prop_negative_dConc_retained = mean(dConc[retained_by_filter] < 0, na.rm = TRUE),
    retention_rate_negative_dConc = mean(retained_by_filter[dConc < 0], na.rm = TRUE),
    retention_rate_positive_dConc = mean(retained_by_filter[dConc > 0], na.rm = TRUE),
    retention_rate_ratio_negative_to_positive = retention_rate_negative_dConc / retention_rate_positive_dConc,
    median_abs_dConc_evaluable = median(abs(dConc), na.rm = TRUE),
    median_abs_dConc_retained = median(abs(dConc[retained_by_filter]), na.rm = TRUE),
    median_dConcSNR_evaluable = median(dConcSNR_calc, na.rm = TRUE),
    median_dConcSNR_retained = median(dConcSNR_calc[retained_by_filter], na.rm = TRUE)
  ) %>%
  arrange(Approach)

concentration_offset_summary <- retained_raw %>%
  filter(is.na(missing_cols)) %>%
  reframe(
    .by = c(SITE_ID, Approach, dLevelsAminusB, tower_order, pair_span),
    n_retained = n(),
    prop_negative_dConc = mean(dConc < 0, na.rm = TRUE),
    prop_positive_FG = mean(FG_mean > 0, na.rm = TRUE),
    median_mean_A = median(mean_A, na.rm = TRUE),
    median_mean_B = median(mean_B, na.rm = TRUE),
    median_dConc = median(dConc, na.rm = TRUE),
    median_abs_dConc = median(abs(dConc), na.rm = TRUE),
    median_dConc_sd = median(dConc_sd, na.rm = TRUE),
    median_dConcSNR = median(dConcSNR_calc, na.rm = TRUE),
    prop_p_lt_0_05 = mean(dConc_pvalue < 0.05, na.rm = TRUE),
    prop_cross_grad_flagged = mean(cross_grad_flagged %in% TRUE, na.rm = TRUE)
  ) %>%
  arrange(desc(prop_negative_dConc), desc(n_retained))

readr::write_csv(retained_quality_by_site, file.path(output_dir, "NEON_CH4_retained_gradient_quality_by_site.csv"))
readr::write_csv(cross_gradient_summary, file.path(output_dir, "NEON_CH4_cross_gradient_flag_summary.csv"))
readr::write_csv(sensor_pair_summary, file.path(output_dir, "NEON_CH4_sensor_pair_behavior_summary.csv"))
readr::write_csv(dconc_filter_bias_summary, file.path(output_dir, "NEON_CH4_dConc_filter_bias_by_pair.csv"))
readr::write_csv(dconc_filter_bias_overall, file.path(output_dir, "NEON_CH4_dConc_filter_bias_overall.csv"))
readr::write_csv(concentration_offset_summary, file.path(output_dir, "NEON_CH4_concentration_offset_summary.csv"))
readr::write_csv(raw_load_issues, file.path(output_dir, "NEON_CH4_gradient_quality_raw_load_issues.csv"))

cross_flag_plot <- cross_gradient_summary %>%
  mutate(cross_grad_class = factor(cross_grad_class, levels = c("not flagged", "flagged", "missing"))) %>%
  ggplot(aes(x = Approach, y = prop_positive_FG, fill = cross_grad_class)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set2", na.translate = FALSE) +
  theme_bw(base_size = 11) +
  labs(
    x = NULL,
    y = "Retained rows with positive GF",
    fill = "Cross-gradient flag",
    title = "Retained CH4 source fraction by cross-gradient flag"
  ) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(
  file.path(figure_dir, "NEON_CH4_cross_gradient_flag_audit.png"),
  cross_flag_plot,
  width = 7,
  height = 4.5,
  units = "in",
  dpi = 300
)

filter_bias_plot_data <- dconc_filter_bias_overall %>%
  select(Approach, prop_negative_dConc_evaluable, prop_negative_dConc_retained) %>%
  pivot_longer(
    cols = starts_with("prop_negative"),
    names_to = "stage",
    values_to = "prop_negative_dConc"
  ) %>%
  mutate(
    stage = recode(
      stage,
      prop_negative_dConc_evaluable = "All evaluable",
      prop_negative_dConc_retained = "Retained"
    )
  )

filter_bias_plot <- filter_bias_plot_data %>%
  ggplot(aes(x = Approach, y = prop_negative_dConc, fill = stage)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw(base_size = 11) +
  labs(
    x = NULL,
    y = "Rows with negative dConc",
    fill = NULL,
    title = "Does dConc filtering preferentially retain source-like gradients?",
    subtitle = "Negative dConc corresponds to upward-positive/source-like GF for A-above-B pairs."
  ) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(
  file.path(figure_dir, "NEON_CH4_dConc_filter_bias_audit.png"),
  filter_bias_plot,
  width = 7,
  height = 4.5,
  units = "in",
  dpi = 300
)

top_pair_plot <- sensor_pair_summary %>%
  filter(n >= 50) %>%
  slice_max(order_by = n, n = 40, with_ties = FALSE) %>%
  mutate(pair_label = paste(SITE_ID, Approach, dLevelsAminusB, sep = " / ")) %>%
  ggplot(aes(x = reorder(pair_label, prop_positive_FG), y = prop_positive_FG, fill = prop_cross_grad_flagged)) +
  geom_col(color = "black", linewidth = 0.15) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_viridis_c(labels = scales::percent_format(accuracy = 1), option = "C", end = 0.9) +
  theme_bw(base_size = 9) +
  labs(
    x = NULL,
    y = "Retained rows with positive GF",
    fill = "Flagged",
    title = "Dominant retained sensor-pair behavior"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(
  file.path(figure_dir, "NEON_CH4_sensor_pair_behavior_audit.png"),
  top_pair_plot,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

overall_lines <- dconc_filter_bias_overall %>%
  mutate(
    line = paste0(
      "- ", Approach, ": evaluable n = ", n_evaluable,
      ", retained n = ", n_retained,
      ", negative dConc evaluable = ", round(100 * prop_negative_dConc_evaluable, 1), "%",
      ", retained = ", round(100 * prop_negative_dConc_retained, 1), "%",
      ", retention negative/positive ratio = ", signif(retention_rate_ratio_negative_to_positive, 3),
      ", retained median |dConc| = ", signif(median_abs_dConc_retained, 3),
      ", retained median dConcSNR = ", signif(median_dConcSNR_retained, 3)
    )
  ) %>%
  pull(line)

cross_lines <- cross_gradient_summary %>%
  mutate(
    line = paste0(
      "- ", Approach, " / ", cross_grad_class, ": n = ", n,
      ", positive GF = ", round(100 * prop_positive_FG, 1), "%",
      ", median GF = ", signif(median_FG, 3),
      ", sites = ", n_sites
    )
  ) %>%
  pull(line)

top_offset_lines <- concentration_offset_summary %>%
  filter(n_retained >= 50) %>%
  slice_max(prop_negative_dConc, n = 10, with_ties = FALSE) %>%
  mutate(
    line = paste0(
      "- ", SITE_ID, " / ", Approach, " / ", dLevelsAminusB,
      ": n = ", n_retained,
      ", negative dConc = ", round(100 * prop_negative_dConc, 1), "%",
      ", median dConc = ", signif(median_dConc, 3),
      ", median |dConc| = ", signif(median_abs_dConc, 3),
      ", median SNR = ", signif(median_dConcSNR, 3)
    )
  ) %>%
  pull(line)

writeLines(
  c(
    "# NEON CH4 Gradient Quality Audit",
    "",
    "## What This Checks",
    "This audit separates broad retained-gradient behavior from sampled raw filter behavior. The all-site retained checks use `SITE_DATA_FILTERED_CH4.Rdata`; the raw-versus-retained checks load sampled site `*_Evaluation.Rdata` and `*_FILTER_AA_AW.Rdata` files so `dConc`, `mean_A`, `mean_B`, and dConc SNR are available.",
    "",
    paste0("Raw dConc filtering diagnostics were run for sampled sites: ", paste(sample_sites, collapse = ", "), ". Set `NEON_GRADIENT_AUDIT_MAX_RAW_SITES` to broaden this sample, but note that the raw evaluation files are large."),
    "",
    "## Cross-Gradient Flags",
    cross_lines,
    "",
    "## dConc Filtering Bias",
    overall_lines,
    "",
    "A retention negative/positive ratio above 1 means the filter retained negative/source-like CH4 gradients more readily than positive/sink-like gradients in the sampled raw files. A value below 1 means the opposite.",
    "",
    "## Largest Retained Negative Concentration Offsets",
    top_offset_lines,
    "",
    "## Outputs",
    "- `OUTPUT/NEON_CH4_retained_gradient_quality_by_site.csv`",
    "- `OUTPUT/NEON_CH4_cross_gradient_flag_summary.csv`",
    "- `OUTPUT/NEON_CH4_sensor_pair_behavior_summary.csv`",
    "- `OUTPUT/NEON_CH4_dConc_filter_bias_by_pair.csv`",
    "- `OUTPUT/NEON_CH4_dConc_filter_bias_overall.csv`",
    "- `OUTPUT/NEON_CH4_concentration_offset_summary.csv`",
    "- `OUTPUT/NEON_CH4_gradient_quality_raw_load_issues.csv`",
    "- `FIGURES/NEON_CH4_cross_gradient_flag_audit.png`",
    "- `FIGURES/NEON_CH4_dConc_filter_bias_audit.png`",
    "- `FIGURES/NEON_CH4_sensor_pair_behavior_audit.png`"
  ),
  file.path(output_dir, "NEON_CH4_gradient_quality_audit.md")
)

message("Wrote NEON CH4 gradient quality audit outputs.")
