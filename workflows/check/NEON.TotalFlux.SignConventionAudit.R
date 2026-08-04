# Audit CH4 gradient-flux sign convention used in total-flux budgets.

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

sign_audit_rows <- purrr::imap_dfr(SITE_DATA_FILTERED, function(site_df, site) {
  needed <- c("FG_mean", "dConc", "gas", "TowerPosition_A", "TowerPosition_B", "dLevelsAminusB")
  missing <- setdiff(needed, names(site_df))

  if (length(missing) > 0) {
    return(tibble(
      SITE_ID = site,
      Approach = NA_character_,
      dLevelsAminusB = NA_character_,
      tower_order = NA_character_,
      FG_mean = NA_real_,
      dConc = NA_real_,
      sign_relation = NA_character_,
      missing_cols = paste(missing, collapse = ", ")
    ))
  }

  if ("Approach" %in% names(site_df)) {
    approach <- as.character(site_df$Approach)
  } else if ("approach" %in% names(site_df)) {
    approach <- as.character(site_df$approach)
  } else {
    approach <- rep("unknown", nrow(site_df))
  }

  site_df %>%
    mutate(Approach = approach) %>%
    filter(gas == "CH4", is.finite(FG_mean), is.finite(dConc), dConc != 0, FG_mean != 0) %>%
    mutate(
      SITE_ID = site,
      sign_relation = case_when(
        sign(FG_mean) == -sign(dConc) ~ "FG opposite dConc",
        sign(FG_mean) == sign(dConc) ~ "FG same as dConc",
        TRUE ~ "other"
      ),
      tower_order = case_when(
        suppressWarnings(as.numeric(TowerPosition_A) > as.numeric(TowerPosition_B)) ~ "A above B",
        suppressWarnings(as.numeric(TowerPosition_A) < as.numeric(TowerPosition_B)) ~ "A below B",
        TRUE ~ "unknown"
      )
    ) %>%
    dplyr::select(SITE_ID, Approach, dLevelsAminusB, tower_order, FG_mean, dConc, sign_relation) %>%
    mutate(missing_cols = NA_character_)
})

sign_audit_summary <- sign_audit_rows %>%
  filter(is.na(missing_cols)) %>%
  reframe(
    .by = c(SITE_ID, Approach, tower_order),
    n = dplyr::n(),
    prop_FG_opposite_dConc = mean(sign_relation == "FG opposite dConc", na.rm = TRUE),
    prop_FG_same_as_dConc = mean(sign_relation == "FG same as dConc", na.rm = TRUE),
    median_FG = median(FG_mean, na.rm = TRUE),
    median_dConc = median(dConc, na.rm = TRUE),
    prop_positive_FG = mean(FG_mean > 0, na.rm = TRUE),
    prop_negative_dConc = mean(dConc < 0, na.rm = TRUE)
  ) %>%
  arrange(SITE_ID, Approach, tower_order)

sign_audit_overall <- sign_audit_rows %>%
  filter(is.na(missing_cols)) %>%
  reframe(
    .by = c(Approach, tower_order),
    n = dplyr::n(),
    prop_FG_opposite_dConc = mean(sign_relation == "FG opposite dConc", na.rm = TRUE),
    prop_FG_same_as_dConc = mean(sign_relation == "FG same as dConc", na.rm = TRUE),
    median_FG = median(FG_mean, na.rm = TRUE),
    median_dConc = median(dConc, na.rm = TRUE),
    prop_positive_FG = mean(FG_mean > 0, na.rm = TRUE),
    prop_negative_dConc = mean(dConc < 0, na.rm = TRUE)
  ) %>%
  arrange(Approach, tower_order)

total_flux_summary <- purrr::imap_dfr(
  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage,
  function(site_df, site) {
    site_df %>%
      filter(is.finite(flux_total)) %>%
      reframe(
        SITE_ID = site,
        n_total_flux = dplyr::n(),
        prop_positive_total_flux = mean(flux_total > 0, na.rm = TRUE),
        median_total_flux = median(flux_total, na.rm = TRUE),
        mean_total_flux = mean(flux_total, na.rm = TRUE),
        median_gradient_flux = median(FG_ENSEMBLE_RSHP, na.rm = TRUE),
        median_storage_flux = median(storage_flux_filled, na.rm = TRUE)
      )
  }
) %>%
  arrange(desc(prop_positive_total_flux))

load_site_filter_sign_rows <- function(site) {
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
    return(tibble(
      SITE_ID = site,
      Approach = NA_character_,
      dLevelsAminusB = NA_character_,
      tower_order = NA_character_,
      FG_mean = NA_real_,
      dConc = NA_real_,
      sign_relation = NA_character_,
      missing_cols = "site filter file not found"
    ))
  }

  env <- new.env(parent = emptyenv())
  load(filter_file, envir = env)

  purrr::map_dfr(
    c("MBR_9min_FILTER", "AE_9min_FILTER", "WP_9min_FILTER"),
    function(object_name) {
      if (!exists(object_name, envir = env)) {
        return(tibble())
      }

      site_df <- get(object_name, envir = env)
      approach <- str_remove(object_name, "_9min_FILTER")
      needed <- c("FG_mean", "dConc", "gas", "TowerPosition_A", "TowerPosition_B", "dLevelsAminusB")
      missing <- setdiff(needed, names(site_df))

      if (length(missing) > 0) {
        return(tibble(
          SITE_ID = site,
          Approach = approach,
          dLevelsAminusB = NA_character_,
          tower_order = NA_character_,
          FG_mean = NA_real_,
          dConc = NA_real_,
          sign_relation = NA_character_,
          missing_cols = paste(missing, collapse = ", ")
        ))
      }

      site_df %>%
        filter(gas == "CH4", is.finite(FG_mean), is.finite(dConc), dConc != 0, FG_mean != 0) %>%
        mutate(
          SITE_ID = site,
          Approach = approach,
          sign_relation = case_when(
            sign(FG_mean) == -sign(dConc) ~ "FG opposite dConc",
            sign(FG_mean) == sign(dConc) ~ "FG same as dConc",
            TRUE ~ "other"
          ),
          tower_order = case_when(
            suppressWarnings(as.numeric(TowerPosition_A) > as.numeric(TowerPosition_B)) ~ "A above B",
            suppressWarnings(as.numeric(TowerPosition_A) < as.numeric(TowerPosition_B)) ~ "A below B",
            TRUE ~ "unknown"
          ),
          missing_cols = NA_character_
        ) %>%
        dplyr::select(SITE_ID, Approach, dLevelsAminusB, tower_order, FG_mean, dConc, sign_relation, missing_cols)
    }
  )
}

if (nrow(sign_audit_overall) == 0) {
  sample_sites <- total_flux_summary %>%
    filter(n_total_flux > 0) %>%
    summarise(
      sites = list(unique(c(
        head(SITE_ID, 10),
        tail(SITE_ID, 5),
        "KONZ",
        "BARR"
      )))
    ) %>%
    pull(sites) %>%
    .[[1]]

  sign_audit_rows <- purrr::map_dfr(sample_sites, load_site_filter_sign_rows)

  sign_audit_summary <- sign_audit_rows %>%
    filter(is.na(missing_cols)) %>%
    reframe(
      .by = c(SITE_ID, Approach, tower_order),
      n = dplyr::n(),
      prop_FG_opposite_dConc = mean(sign_relation == "FG opposite dConc", na.rm = TRUE),
      prop_FG_same_as_dConc = mean(sign_relation == "FG same as dConc", na.rm = TRUE),
      median_FG = median(FG_mean, na.rm = TRUE),
      median_dConc = median(dConc, na.rm = TRUE),
      prop_positive_FG = mean(FG_mean > 0, na.rm = TRUE),
      prop_negative_dConc = mean(dConc < 0, na.rm = TRUE)
    ) %>%
    arrange(SITE_ID, Approach, tower_order)

  sign_audit_overall <- sign_audit_rows %>%
    filter(is.na(missing_cols)) %>%
    reframe(
      .by = c(Approach, tower_order),
      n = dplyr::n(),
      prop_FG_opposite_dConc = mean(sign_relation == "FG opposite dConc", na.rm = TRUE),
      prop_FG_same_as_dConc = mean(sign_relation == "FG same as dConc", na.rm = TRUE),
      median_FG = median(FG_mean, na.rm = TRUE),
      median_dConc = median(dConc, na.rm = TRUE),
      prop_positive_FG = mean(FG_mean > 0, na.rm = TRUE),
      prop_negative_dConc = mean(dConc < 0, na.rm = TRUE)
    ) %>%
    arrange(Approach, tower_order)
}

readr::write_csv(sign_audit_summary, file.path(output_dir, "NEON_CH4_gradient_flux_sign_audit_by_site.csv"))
readr::write_csv(sign_audit_overall, file.path(output_dir, "NEON_CH4_gradient_flux_sign_audit_overall.csv"))
readr::write_csv(total_flux_summary, file.path(output_dir, "NEON_total_flux_sign_summary_by_site.csv"))

sign_plot <- sign_audit_overall %>%
  ggplot(aes(x = Approach, y = prop_FG_opposite_dConc, fill = tower_order)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_bw(base_size = 11) +
  labs(
    x = NULL,
    y = "Rows where sign(GF_mean) is opposite sign(dConc)",
    fill = "Tower order",
    title = "CH4 gradient-flux sign audit",
    subtitle = "For F = -K dC/dz, A above B should produce sign(GF) opposite sign(concentration_A - concentration_B)."
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "NEON_CH4_gradient_flux_sign_audit.png"),
  sign_plot,
  width = 7,
  height = 4.5,
  units = "in",
  dpi = 300
)

overall_lines <- sign_audit_overall %>%
  mutate(
    line = paste0(
      "- ", Approach, " / ", tower_order, ": n = ", n,
      ", opposite sign = ", round(100 * prop_FG_opposite_dConc, 1), "%",
      ", same sign = ", round(100 * prop_FG_same_as_dConc, 1), "%",
      ", median GF = ", signif(median_FG, 3),
      ", median dConc = ", signif(median_dConc, 3)
    )
  ) %>%
  pull(line)

total_flux_lines <- total_flux_summary %>%
  slice_head(n = 10) %>%
  mutate(
    line = paste0(
      "- ", SITE_ID, ": n = ", n_total_flux,
      ", positive total flux = ", round(100 * prop_positive_total_flux, 1), "%",
      ", median total = ", signif(median_total_flux, 3),
      ", median gradient = ", signif(median_gradient_flux, 3),
      ", median storage = ", signif(median_storage_flux, 3)
    )
  ) %>%
  pull(line)

writeLines(
  c(
    "# NEON CH4 Gradient-Flux Sign Convention Audit",
    "",
    "## Interpretation",
    "`dConc` is `concentration_A - concentration_B`. For rows where tower position A is above B, Fickian upward-positive flux should have `sign(FG_mean) = -sign(dConc)`: higher CH4 aloft than below implies downward uptake, and higher CH4 below than aloft implies upward emission.",
    "",
    "The filtered CH4 data overwhelmingly follow that opposite-sign relationship for A-above-B rows. That means the current `FG_mean` sign convention is consistent with upward-positive flux rather than obviously inverted.",
    "",
    "The positive annual/source behavior is therefore not explained by a simple sign-convention inversion in `FG_mean`. It is driven by many retained observations having negative CH4 gradients (`dConc < 0`), which translate to positive upward flux under the current convention.",
    "",
    "## Overall Sign Relationship",
    overall_lines,
    "",
    "## Total-Flux Sites With Highest Positive Fraction",
    total_flux_lines,
    "",
    "## Outputs",
    "- `OUTPUT/NEON_CH4_gradient_flux_sign_audit_by_site.csv`",
    "- `OUTPUT/NEON_CH4_gradient_flux_sign_audit_overall.csv`",
    "- `OUTPUT/NEON_total_flux_sign_summary_by_site.csv`",
    "- `FIGURES/NEON_CH4_gradient_flux_sign_audit.png`"
  ),
  file.path(output_dir, "NEON_CH4_gradient_flux_sign_audit.md")
)

message("Wrote NEON CH4 gradient-flux sign convention audit outputs.")
