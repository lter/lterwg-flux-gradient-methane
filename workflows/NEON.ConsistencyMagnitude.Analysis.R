# When/where/magnitude interpretation for NEON 30-minute CH4 behavior.
#
# This script turns the model-standardized gapfill outputs into interpretable
# site and seasonal summaries:
# - when sites are expected to be sources/sinks through the year,
# - where consistent versus fluctuating behavior appears in site-attribute space,
# - and the expected annual CH4 magnitude with uncertainty.

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(rpart)
library(rpart.plot)

localdir.ch4 <- "/Volumes/MaloneLab/Research/FluxGradient/Methane"
setwd(localdir.ch4)

dir.create("OUTPUT", showWarnings = FALSE)
dir.create("FIGURES", showWarnings = FALSE)

monthly_file <- "OUTPUT/NON_30min_gapfill_monthly_budgets.csv"
annual_file <- "OUTPUT/NON_30min_gapfill_annual_budgets.csv"
site_driver_file <- "OUTPUT/NEON_site_driver_values_for_sink_comparison.csv"
multivariate_matrix_file <- "OUTPUT/NEON_site_attribute_multivariate_matrix.csv"
variable_quality_file <- "OUTPUT/NEON_site_attribute_variables_used.csv"

required_files <- c(monthly_file, annual_file, site_driver_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required inputs: ", paste(missing_files, collapse = ", "),
       ". Run NON.30min.Gapfill.r and NEON.StrongSink.DriverComparison.R first.")
}

behavior_levels <- c("Consistent sink", "Fluctuating", "Consistent source")
behavior_colors <- c(
  "Consistent sink" = "#2166AC",
  "Fluctuating" = "#4D4D4D",
  "Consistent source" = "#B2182B"
)
consistency_colors <- c("Consistent" = "#2166AC", "Fluctuating" = "#4D4D4D")

monthly_budgets <- read.csv(monthly_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    month = as.integer(month),
    month_name = factor(month_name, levels = month.abb),
    monthly_budget_gC_m2 = budget_mean_mgC_m2_month / 1000
  )

annual_budgets <- read.csv(annual_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    standardized_behavior = factor(standardized_behavior, levels = behavior_levels),
    observed_behavior = factor(observed_behavior, levels = behavior_levels),
    consistency_group = factor(
      if_else(standardized_behavior == "Fluctuating", "Fluctuating", "Consistent"),
      levels = c("Fluctuating", "Consistent")
    ),
    annual_uncertainty_width_gC_m2_yr = annual_budget_upr_gC_m2_yr - annual_budget_lwr_gC_m2_yr
  )

site_drivers <- read.csv(site_driver_file) %>%
  mutate(SITE_ID = as.character(SITE_ID))

if (file.exists(multivariate_matrix_file)) {
  site_attributes <- read.csv(multivariate_matrix_file) %>%
    mutate(SITE_ID = as.character(SITE_ID)) %>%
    dplyr::select(-any_of(c("CH4_behavior", "standardized_behavior")))
} else {
  flux_summary_cols <- c(
    "n_months", "prop_source_months", "prop_sink_months", "prop_source_gradient_months",
    "mean_CH4_mgC_30min", "mean_gradient_mgC_30min", "mean_storage_mgC_30min",
    "median_CH4_mgC_30min", "median_gradient_mgC_30min", "median_storage_mgC_30min",
    "median_storage_abs_fraction", "sd_monthly_CH4", "sign_changes",
    "CH4_behavior", "CH4_gradient_behavior", "behavior_changed_from_gradient",
    "behavior_comparison", "n_30min", "n_months_with_vswc"
  )
  site_attributes <- site_drivers %>%
    dplyr::select(-any_of(flux_summary_cols))
}

analysis_data_raw <- annual_budgets %>%
  left_join(site_attributes, by = "SITE_ID")

numeric_predictors <- analysis_data_raw %>%
  dplyr::select(where(is.numeric)) %>%
  names() %>%
  setdiff(c(
    "annual_budget_mean_mgC_m2_yr", "annual_budget_median_mgC_m2_yr",
    "annual_budget_lwr_mgC_m2_yr", "annual_budget_upr_mgC_m2_yr",
    "annual_budget_sd_mgC_m2_yr", "annual_budget_mean_gC_m2_yr",
    "annual_budget_median_gC_m2_yr", "annual_budget_lwr_gC_m2_yr",
    "annual_budget_upr_gC_m2_yr", "annual_budget_sd_gC_m2_yr",
    "prob_annual_source", "standardized_source_months",
    "standardized_prop_source_months", "annual_uncertainty_width_gC_m2_yr"
  ))

predictor_quality <- analysis_data_raw %>%
  dplyr::select(all_of(numeric_predictors)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  reframe(
    .by = variable,
    n_sites = sum(is.finite(value)),
    prop_missing = mean(!is.finite(value)),
    sd_value = sd(value, na.rm = TRUE),
    median_value = median(value, na.rm = TRUE),
    used = n_sites >= 20 & prop_missing <= 0.4 & is.finite(sd_value) & sd_value > 0
  )

usable_predictors <- predictor_quality %>%
  filter(used) %>%
  pull(variable)

analysis_data <- analysis_data_raw
for (variable in usable_predictors) {
  fill_value <- median(analysis_data[[variable]], na.rm = TRUE)
  analysis_data[[variable]][!is.finite(analysis_data[[variable]])] <- fill_value
}

predictor_scores <- purrr::map_dfr(usable_predictors, function(variable) {
  x <- analysis_data[[variable]]
  budget_rho <- suppressWarnings(cor(x, analysis_data$annual_budget_mean_gC_m2_yr, method = "spearman"))
  source_rho <- suppressWarnings(cor(x, analysis_data$standardized_prop_source_months, method = "spearman"))

  class_test <- tryCatch(
    kruskal.test(x ~ analysis_data$consistency_group),
    error = function(e) NULL
  )

  tibble(
    variable = variable,
    budget_spearman = budget_rho,
    source_month_spearman = source_rho,
    consistency_p = if (is.null(class_test)) NA_real_ else class_test$p.value,
    score = max(abs(budget_rho), abs(source_rho), na.rm = TRUE)
  )
}) %>%
  left_join(predictor_quality, by = "variable") %>%
  arrange(desc(score), consistency_p, variable)

tree_predictors <- predictor_scores %>%
  slice_head(n = min(12, nrow(predictor_scores))) %>%
  pull(variable)

if (length(tree_predictors) < 3) {
  stop("Need at least three usable predictors for consistency/magnitude interpretation.")
}

tree_formula_class <- as.formula(paste("consistency_group ~", paste(tree_predictors, collapse = " + ")))
tree_formula_budget <- as.formula(paste("annual_budget_mean_gC_m2_yr ~", paste(tree_predictors, collapse = " + ")))

classification_tree <- rpart(
  tree_formula_class,
  data = analysis_data,
  method = "class",
  control = rpart.control(cp = 0.01, minsplit = 6, minbucket = 3, xval = 10)
)

magnitude_tree <- rpart(
  tree_formula_budget,
  data = analysis_data,
  method = "anova",
  control = rpart.control(cp = 0.01, minsplit = 6, minbucket = 3, xval = 10)
)

loocv_results <- purrr::map_dfr(seq_len(nrow(analysis_data)), function(i) {
  train <- analysis_data[-i, , drop = FALSE]
  test <- analysis_data[i, , drop = FALSE]

  class_fit <- rpart(
    tree_formula_class,
    data = train,
    method = "class",
    control = rpart.control(cp = 0.01, minsplit = 6, minbucket = 3, xval = 0)
  )

  budget_fit <- rpart(
    tree_formula_budget,
    data = train,
    method = "anova",
    control = rpart.control(cp = 0.01, minsplit = 6, minbucket = 3, xval = 0)
  )

  tibble(
    SITE_ID = test$SITE_ID,
    observed_consistency = test$consistency_group,
    predicted_consistency = factor(
      predict(class_fit, newdata = test, type = "class"),
      levels = levels(analysis_data$consistency_group)
    ),
    observed_budget_gC_m2_yr = test$annual_budget_mean_gC_m2_yr,
    predicted_budget_gC_m2_yr = as.numeric(predict(budget_fit, newdata = test))
  )
})

classification_accuracy <- mean(loocv_results$observed_consistency == loocv_results$predicted_consistency)
budget_rmse <- sqrt(mean((loocv_results$observed_budget_gC_m2_yr - loocv_results$predicted_budget_gC_m2_yr)^2))
budget_mae <- mean(abs(loocv_results$observed_budget_gC_m2_yr - loocv_results$predicted_budget_gC_m2_yr))

variable_importance <- bind_rows(
  enframe(classification_tree$variable.importance, name = "variable", value = "importance") %>%
    mutate(model = "Consistency tree"),
  enframe(magnitude_tree$variable.importance, name = "variable", value = "importance") %>%
    mutate(model = "Magnitude tree")
) %>%
  group_by(model) %>%
  mutate(relative_importance = importance / max(importance, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(model, desc(relative_importance))

site_interpretation <- analysis_data %>%
  dplyr::select(
    SITE_ID, observed_behavior, standardized_behavior, consistency_group,
    standardized_prop_source_months, annual_budget_mean_gC_m2_yr,
    annual_budget_lwr_gC_m2_yr, annual_budget_upr_gC_m2_yr,
    annual_uncertainty_width_gC_m2_yr, prob_annual_source,
    all_of(tree_predictors)
  ) %>%
  left_join(loocv_results %>% dplyr::select(SITE_ID, predicted_consistency, predicted_budget_gC_m2_yr), by = "SITE_ID") %>%
  arrange(standardized_behavior, desc(abs(annual_budget_mean_gC_m2_yr)))

write.csv(site_interpretation, "OUTPUT/NEON_consistency_magnitude_site_summary.csv", row.names = FALSE)
write.csv(predictor_scores, "OUTPUT/NEON_consistency_magnitude_predictor_scores.csv", row.names = FALSE)
write.csv(variable_importance, "OUTPUT/NEON_consistency_magnitude_tree_importance.csv", row.names = FALSE)
write.csv(loocv_results, "OUTPUT/NEON_consistency_magnitude_loocv.csv", row.names = FALSE)

save(
  classification_tree, magnitude_tree, tree_predictors,
  file = "OUTPUT/NEON_consistency_magnitude_trees.Rdata"
)

monthly_plot_data <- monthly_budgets %>%
  left_join(annual_budgets %>% dplyr::select(SITE_ID, standardized_behavior, annual_budget_mean_gC_m2_yr), by = "SITE_ID") %>%
  mutate(
    SITE_ID = fct_reorder(SITE_ID, annual_budget_mean_gC_m2_yr),
    standardized_behavior = factor(standardized_behavior, levels = behavior_levels)
  )

plot_monthly_budget_heatmap <- monthly_plot_data %>%
  ggplot(aes(x = month_name, y = SITE_ID, fill = monthly_budget_gC_m2)) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(standardized_behavior ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
  scale_fill_gradient2(
    low = behavior_colors[["Consistent sink"]],
    mid = "white",
    high = behavior_colors[["Consistent source"]],
    midpoint = 0
  ) +
  labs(
    title = "A. When: Model-Standardized Monthly CH4 Budgets",
    x = NULL,
    y = NULL,
    fill = expression(paste("g C ", m^-2, " month"^-1))
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(face = "bold", angle = 0),
    panel.grid = element_blank()
  )

plot_source_probability_heatmap <- monthly_plot_data %>%
  ggplot(aes(x = month_name, y = SITE_ID, fill = prob_source_month)) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(standardized_behavior ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
  scale_fill_gradient(
    low = behavior_colors[["Consistent sink"]],
    high = behavior_colors[["Consistent source"]],
    limits = c(0, 1)
  ) +
  labs(
    title = "B. When: Monthly Source Probability",
    x = NULL,
    y = NULL,
    fill = "P(source)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(face = "bold", angle = 0),
    panel.grid = element_blank()
  )

monthly_profile <- monthly_plot_data %>%
  reframe(
    .by = c(standardized_behavior, month, month_name),
    mean_budget = mean(monthly_budget_gC_m2, na.rm = TRUE),
    se_budget = sd(monthly_budget_gC_m2, na.rm = TRUE) / sqrt(dplyr::n()),
    mean_prob_source = mean(prob_source_month, na.rm = TRUE),
    se_prob_source = sd(prob_source_month, na.rm = TRUE) / sqrt(dplyr::n())
  ) %>%
  mutate(
    se_budget = replace_na(se_budget, 0),
    se_prob_source = replace_na(se_prob_source, 0)
  )

plot_monthly_profiles <- monthly_profile %>%
  ggplot(aes(x = month, y = mean_budget, color = standardized_behavior, fill = standardized_behavior)) +
  geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
  geom_ribbon(aes(ymin = mean_budget - se_budget, ymax = mean_budget + se_budget), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_color_manual(values = behavior_colors, na.translate = FALSE) +
  scale_fill_manual(values = behavior_colors, na.translate = FALSE) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title = "C. Expected Seasonal Magnitude by Standardized Behavior",
    x = NULL,
    y = expression(paste("Mean monthly CH"[4], " budget (g C ", m^-2, ")")),
    color = "Standardized behavior",
    fill = "Standardized behavior"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

when_figure <- (plot_monthly_budget_heatmap / plot_source_probability_heatmap / plot_monthly_profiles) +
  plot_layout(heights = c(1.2, 1.2, 0.8)) +
  plot_annotation(
    title = "When Sites Are Expected to Be Sources or Sinks",
    subtitle = "Model-standardized monthly budgets from a balanced month x half-hour grid",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

ggsave(
  "FIGURES/NEON_when_monthly_consistency_magnitude.png",
  plot = when_figure,
  width = 12,
  height = 14,
  units = "in",
  dpi = 300
)

ggsave(
  "FIGURES/NEON_when_monthly_consistency_magnitude.pdf",
  plot = when_figure,
  width = 12,
  height = 14,
  units = "in"
)

plot_annual_magnitude <- annual_budgets %>%
  mutate(SITE_ID = fct_reorder(SITE_ID, annual_budget_mean_gC_m2_yr)) %>%
  ggplot(aes(x = annual_budget_mean_gC_m2_yr, y = SITE_ID, color = standardized_behavior)) +
  geom_vline(xintercept = 0, color = "grey45", linetype = "dashed") +
  geom_errorbar(
    aes(xmin = annual_budget_lwr_gC_m2_yr, xmax = annual_budget_upr_gC_m2_yr),
    orientation = "y",
    width = 0,
    alpha = 0.7
  ) +
  geom_point(size = 2.4) +
  facet_grid(standardized_behavior ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
  scale_color_manual(values = behavior_colors, na.translate = FALSE) +
  labs(
    title = "A. Expected Annual Magnitude",
    x = expression(paste("Annual CH"[4], " budget (g C ", m^-2, " yr"^-1, ")")),
    y = NULL,
    color = "Standardized behavior"
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.y = element_text(size = 7),
    strip.text.y = element_text(face = "bold", angle = 0)
  )

plot_prop_vs_magnitude <- annual_budgets %>%
  ggplot(aes(x = standardized_prop_source_months, y = annual_budget_mean_gC_m2_yr, color = standardized_behavior)) +
  geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
  geom_vline(xintercept = c(0.25, 0.75), color = "grey70", linetype = "dotted") +
  geom_errorbar(aes(ymin = annual_budget_lwr_gC_m2_yr, ymax = annual_budget_upr_gC_m2_yr), width = 0.015, alpha = 0.45) +
  geom_point(size = 2.6) +
  ggrepel::geom_text_repel(aes(label = SITE_ID), size = 3, max.overlaps = 50, show.legend = FALSE) +
  scale_color_manual(values = behavior_colors, na.translate = FALSE) +
  labs(
    title = "B. Consistency Versus Magnitude",
    x = "Fraction of standardized source months",
    y = expression(paste("Annual CH"[4], " budget (g C ", m^-2, " yr"^-1, ")")),
    color = "Standardized behavior"
  ) +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

magnitude_figure <- plot_annual_magnitude / plot_prop_vs_magnitude +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(
    title = "Expected CH4 Magnitude and Consistency by Site",
    subtitle = "Annual budgets are from the model-standardized balanced grid",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

ggsave(
  "FIGURES/NEON_annual_magnitude_consistency.png",
  plot = magnitude_figure,
  width = 11,
  height = 13,
  units = "in",
  dpi = 300
)

ggsave(
  "FIGURES/NEON_annual_magnitude_consistency.pdf",
  plot = magnitude_figure,
  width = 11,
  height = 13,
  units = "in"
)

tree_png <- "FIGURES/NEON_where_consistency_tree.png"
png(tree_png, width = 1800, height = 1100, res = 180)
rpart.plot::rpart.plot(
  classification_tree,
  main = "Where: Site Attributes Associated with Consistent vs Fluctuating Behavior",
  type = 2,
  extra = 104,
  under = TRUE,
  fallen.leaves = TRUE,
  box.palette = "GnBu"
)
dev.off()

tree_budget_png <- "FIGURES/NEON_where_magnitude_tree.png"
png(tree_budget_png, width = 1800, height = 1100, res = 180)
rpart.plot::rpart.plot(
  magnitude_tree,
  main = "Where: Site Attributes Associated with Annual CH4 Magnitude",
  type = 2,
  extra = 101,
  under = TRUE,
  fallen.leaves = TRUE,
  box.palette = "Blues"
)
dev.off()

pdf("FIGURES/NEON_where_consistency_tree.pdf", width = 10, height = 7)
rpart.plot::rpart.plot(
  classification_tree,
  main = "Where: Site Attributes Associated with Consistent vs Fluctuating Behavior",
  type = 2,
  extra = 104,
  under = TRUE,
  fallen.leaves = TRUE,
  box.palette = "GnBu"
)
dev.off()

pdf("FIGURES/NEON_where_magnitude_tree.pdf", width = 10, height = 7)
rpart.plot::rpart.plot(
  magnitude_tree,
  main = "Where: Site Attributes Associated with Annual CH4 Magnitude",
  type = 2,
  extra = 101,
  under = TRUE,
  fallen.leaves = TRUE,
  box.palette = "Blues"
)
dev.off()

top_site_lines <- site_interpretation %>%
  arrange(desc(abs(annual_budget_mean_gC_m2_yr))) %>%
  slice_head(n = 10) %>%
  mutate(
    line = paste0(
      "- ", SITE_ID, " (", standardized_behavior, "): ",
      signif(annual_budget_mean_gC_m2_yr, 3), " g C m-2 yr-1 [",
      signif(annual_budget_lwr_gC_m2_yr, 3), ", ",
      signif(annual_budget_upr_gC_m2_yr, 3), "], source-month fraction = ",
      signif(standardized_prop_source_months, 3)
    )
  ) %>%
  pull(line)

top_predictor_lines <- predictor_scores %>%
  slice_head(n = 10) %>%
  mutate(
    line = paste0(
      "- ", variable,
      ": |score| = ", signif(score, 3),
      ", budget rho = ", signif(budget_spearman, 3),
      ", source-month rho = ", signif(source_month_spearman, 3),
      ", consistency p = ", signif(consistency_p, 3)
    )
  ) %>%
  pull(line)

importance_lines <- variable_importance %>%
  group_by(model) %>%
  slice_head(n = 6) %>%
  ungroup() %>%
  mutate(line = paste0("- ", model, " / ", variable, ": ", signif(relative_importance, 3))) %>%
  pull(line)

writeLines(
  c(
    "# NEON Consistency and Magnitude Analysis",
    "",
    "## Purpose",
    "This analysis is organized around when sites are likely to switch or remain stable, where consistent/fluctuating behavior occurs in site-attribute space, and what annual CH4 magnitude is expected.",
    "",
    "## Inputs",
    "- Model-standardized monthly and annual budgets from `NON.30min.Gapfill.r`.",
    "- Site attributes from the NEON driver comparison and, when available, the expanded multivariate site-attribute matrix.",
    "",
    "## Predictive Summaries",
    paste0("- Tree predictors screened from ", length(usable_predictors), " usable site attributes; trees use the top ", length(tree_predictors), " by univariate association with consistency or magnitude."),
    paste0("- Leave-one-site-out consistency tree accuracy: ", scales::percent(classification_accuracy, accuracy = 1), "."),
    paste0("- Leave-one-site-out annual-budget RMSE: ", signif(budget_rmse, 3), " g C m-2 yr-1; MAE: ", signif(budget_mae, 3), " g C m-2 yr-1."),
    "",
    "## Largest Absolute Annual Budgets",
    top_site_lines,
    "",
    "## Top Predictor Screens",
    top_predictor_lines,
    "",
    "## Tree Variable Importance",
    importance_lines,
    "",
    "## Interpretation Notes",
    "- These trees are small-sample exploratory summaries, not final predictive models.",
    "- Use the monthly heatmaps for timing, the annual budget figure for expected magnitude, and the tree summaries for candidate site attributes associated with consistency/magnitude.",
    "",
    "## Outputs",
    "- `OUTPUT/NEON_consistency_magnitude_site_summary.csv`",
    "- `OUTPUT/NEON_consistency_magnitude_predictor_scores.csv`",
    "- `OUTPUT/NEON_consistency_magnitude_tree_importance.csv`",
    "- `OUTPUT/NEON_consistency_magnitude_loocv.csv`",
    "- `OUTPUT/NEON_consistency_magnitude_trees.Rdata`",
    "- `FIGURES/NEON_when_monthly_consistency_magnitude.png`",
    "- `FIGURES/NEON_when_monthly_consistency_magnitude.pdf`",
    "- `FIGURES/NEON_annual_magnitude_consistency.png`",
    "- `FIGURES/NEON_annual_magnitude_consistency.pdf`",
    "- `FIGURES/NEON_where_consistency_tree.png`",
    "- `FIGURES/NEON_where_consistency_tree.pdf`",
    "- `FIGURES/NEON_where_magnitude_tree.png`",
    "- `FIGURES/NEON_where_magnitude_tree.pdf`"
  ),
  "OUTPUT/NEON_consistency_magnitude_results.md"
)

message("Wrote NEON consistency, timing, and magnitude analysis outputs.")
