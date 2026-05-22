# Multivariate site-attribute structure for NEON CH4 behavior categories.
#
# This is an exploratory companion to NEON.StrongSink.DriverComparison.R.
# It asks whether consistent/fluctuating CH4 behavior classes occupy different
# regions of multivariate site-attribute space. It should not be interpreted as
# causal driver attribution.

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cluster)
library(vegan)

localdir.ch4 <- "/Volumes/MaloneLab/Research/FluxGradient/Methane"
localdir <- "/Volumes/MaloneLab/Research/FluxGradient"
setwd(localdir.ch4)

dir.create("OUTPUT", showWarnings = FALSE)
dir.create("FIGURES", showWarnings = FALSE)

site_driver_file <- "OUTPUT/NEON_site_driver_values_for_sink_comparison.csv"
standardized_budget_file <- "OUTPUT/NON_30min_gapfill_annual_budgets.csv"

if (!file.exists(site_driver_file)) {
  stop("Missing OUTPUT/NEON_site_driver_values_for_sink_comparison.csv. Run NEON.StrongSink.DriverComparison.R first.")
}

site_drivers <- read.csv(site_driver_file) %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    CH4_behavior = factor(CH4_behavior, levels = c("Consistent sink", "Fluctuating", "Consistent source"))
  )

canopy_file <- file.path(localdir, "canopy_commbined.csv")

if (file.exists(canopy_file)) {
  canopy_exclude <- c(
    "X", "Year", "dist.m", "wedge",
    "MeasurementHeight_m_A", "MeasurementHeight_m_B",
    "TowerPosition_A", "TowerPosition_B",
    "Max1_Tower_Position", "Max2_Tower_Position",
    "EVI.years", "NDVI.years", "LAI.years", "PRI.years", "CHM.years", "SAVI.years"
  )

  canopy_site <- read.csv(canopy_file) %>%
    rename(SITE_ID = Site) %>%
    mutate(SITE_ID = as.character(SITE_ID)) %>%
    dplyr::select(
      SITE_ID,
      where(is.numeric),
      -any_of(canopy_exclude)
    ) %>%
    reframe(
      .by = SITE_ID,
      across(where(is.numeric), ~ mean(.x, na.rm = TRUE))
    ) %>%
    mutate(across(where(is.numeric), ~ if_else(is.nan(.x), NA_real_, .x))) %>%
    rename_with(~ paste0("canopy_", .x), -SITE_ID)

  site_drivers <- site_drivers %>%
    left_join(canopy_site, by = "SITE_ID")
} else {
  canopy_site <- tibble(SITE_ID = character())
}

if (file.exists(standardized_budget_file)) {
  standardized_behavior <- read.csv(standardized_budget_file) %>%
    transmute(
      SITE_ID = as.character(SITE_ID),
      standardized_behavior = factor(
        standardized_behavior,
        levels = c("Consistent sink", "Fluctuating", "Consistent source")
      )
    )
} else {
  standardized_behavior <- tibble(
    SITE_ID = character(),
    standardized_behavior = factor(levels = c("Consistent sink", "Fluctuating", "Consistent source"))
  )
}

candidate_variables <- c(
  "VSWCMean_site", "VSWCVar_site", "VSWCMean_obs_sd",
  "VSWC_monthly_mean", "VSWC_monthly_sd", "VSWC_monthly_range",
  "sulfurTot", "biogeoTopDepth", "biogeoBottomDepth", "carbonTot",
  "nitrogenTot", "ctonRatio", "estimatedOC", "acidity",
  "bulkDensOvenDry", "sandTotal", "siltTotal", "clayTotal", "dryMass",
  "MAP", "MAT", "canopyHeight_m", "LAI.mean", "CHM.mean",
  names(site_drivers)[startsWith(names(site_drivers), "canopy_")]
)

variable_labels <- c(
  VSWCMean_site = "Mean VSWC",
  VSWCVar_site = "Mean VSWC variance",
  VSWCMean_obs_sd = "VSWC observation SD",
  VSWC_monthly_mean = "Monthly VSWC mean",
  VSWC_monthly_sd = "Monthly VSWC SD",
  VSWC_monthly_range = "Monthly VSWC range",
  sulfurTot = "Sulfur total",
  biogeoTopDepth = "Biogeochem top depth",
  biogeoBottomDepth = "Biogeochem bottom depth",
  carbonTot = "Total C",
  nitrogenTot = "Total N",
  ctonRatio = "C:N",
  estimatedOC = "Estimated OC",
  acidity = "pH/acidity",
  bulkDensOvenDry = "Bulk density",
  sandTotal = "Sand",
  siltTotal = "Silt",
  clayTotal = "Clay",
  dryMass = "Dry mass",
  MAP = "MAP",
  MAT = "MAT",
  canopyHeight_m = "Canopy height",
  LAI.mean = "LAI",
  CHM.mean = "CHM"
)

available_variables <- unique(candidate_variables[candidate_variables %in% names(site_drivers)])

variable_quality <- site_drivers %>%
  dplyr::select(all_of(available_variables)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  reframe(
    .by = variable,
    n_sites = sum(is.finite(value)),
    prop_missing = mean(!is.finite(value)),
    sd_value = sd(value, na.rm = TRUE),
    median_value = median(value, na.rm = TRUE)
  ) %>%
  mutate(
    label = recode(variable, !!!variable_labels, .default = str_remove(variable, "^canopy_")),
    used = n_sites >= 20 & prop_missing <= 0.4 & is.finite(sd_value) & sd_value > 0
  ) %>%
  arrange(desc(used), prop_missing, variable)

driver_variables <- variable_quality %>%
  filter(used) %>%
  pull(variable)

if (length(driver_variables) < 3) {
  stop("Need at least three usable site-attribute variables for MDS/clustering.")
}

site_matrix_raw <- site_drivers %>%
  dplyr::select(SITE_ID, CH4_behavior, all_of(driver_variables)) %>%
  left_join(standardized_behavior, by = "SITE_ID")

site_matrix_imputed <- site_matrix_raw
for (variable in driver_variables) {
  fill_value <- median(site_matrix_imputed[[variable]], na.rm = TRUE)
  site_matrix_imputed[[variable]][!is.finite(site_matrix_imputed[[variable]])] <- fill_value
}

scaled_matrix <- site_matrix_imputed %>%
  dplyr::select(all_of(driver_variables)) %>%
  scale() %>%
  as.matrix()

rownames(scaled_matrix) <- site_matrix_imputed$SITE_ID

distance_matrix <- dist(scaled_matrix, method = "euclidean")
mds_fit <- cmdscale(distance_matrix, k = 2, eig = TRUE)
variance_explained <- mds_fit$eig[seq_len(2)] / sum(abs(mds_fit$eig))

cluster_tree <- hclust(distance_matrix, method = "ward.D2")
cluster_range <- 2:min(6, nrow(site_matrix_imputed) - 1)
silhouette_summary <- purrr::map_dfr(cluster_range, function(k) {
  cluster_id <- cutree(cluster_tree, k = k)
  sil <- silhouette(cluster_id, distance_matrix)
  tibble(k = k, mean_silhouette = mean(sil[, "sil_width"]))
})

best_k <- silhouette_summary %>%
  arrange(desc(mean_silhouette), k) %>%
  slice(1) %>%
  pull(k)

cluster_id <- cutree(cluster_tree, k = best_k)

site_scores <- tibble(
  SITE_ID = rownames(mds_fit$points),
  MDS1 = mds_fit$points[, 1],
  MDS2 = mds_fit$points[, 2],
  cluster = factor(cluster_id[rownames(mds_fit$points)])
) %>%
  left_join(site_matrix_imputed %>% dplyr::select(SITE_ID, CH4_behavior, standardized_behavior), by = "SITE_ID") %>%
  mutate(
    standardized_behavior = fct_na_value_to_level(standardized_behavior, level = "Not available")
  )

pca_fit <- prcomp(scaled_matrix, center = FALSE, scale. = FALSE)
pca_loadings <- as_tibble(pca_fit$rotation[, 1:2], rownames = "variable") %>%
  mutate(
    label = recode(variable, !!!variable_labels, .default = variable),
    loading_strength = sqrt(PC1^2 + PC2^2)
  ) %>%
  arrange(desc(loading_strength))

permanova_observed <- vegan::adonis2(
  distance_matrix ~ CH4_behavior,
  data = site_scores,
  permutations = 999
) %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  as_tibble() %>%
  mutate(classification = "Observed balanced half-hour behavior")

if (n_distinct(na.omit(site_scores$standardized_behavior)) > 1) {
  permanova_standardized <- vegan::adonis2(
    distance_matrix ~ standardized_behavior,
    data = site_scores,
    permutations = 999
  ) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    mutate(classification = "Model-standardized behavior")
} else {
  permanova_standardized <- tibble()
}

permanova_results <- bind_rows(permanova_observed, permanova_standardized)

cluster_behavior_table <- site_scores %>%
  count(cluster, CH4_behavior, standardized_behavior, name = "n_sites")

write.csv(site_matrix_imputed, "OUTPUT/NEON_site_attribute_multivariate_matrix.csv", row.names = FALSE)
write.csv(variable_quality, "OUTPUT/NEON_site_attribute_variables_used.csv", row.names = FALSE)
write.csv(site_scores, "OUTPUT/NEON_site_attribute_mds_scores.csv", row.names = FALSE)
write.csv(pca_loadings, "OUTPUT/NEON_site_attribute_pca_loadings.csv", row.names = FALSE)
write.csv(silhouette_summary, "OUTPUT/NEON_site_attribute_cluster_silhouette.csv", row.names = FALSE)
write.csv(cluster_behavior_table, "OUTPUT/NEON_site_attribute_cluster_behavior_table.csv", row.names = FALSE)
write.csv(permanova_results, "OUTPUT/NEON_site_attribute_permanova.csv", row.names = FALSE)

behavior_colors <- c("Consistent sink" = "red3", "Fluctuating" = "grey35", "Consistent source" = "blue4")
standardized_colors <- c(behavior_colors, "Not available" = "grey80")

top_loading_labels <- pca_loadings %>%
  slice_head(n = 6) %>%
  mutate(
    PC1 = PC1 * max(abs(site_scores$MDS1), na.rm = TRUE) * 1.8,
    PC2 = PC2 * max(abs(site_scores$MDS2), na.rm = TRUE) * 1.8
  )

plot_mds_observed <- site_scores %>%
  ggplot(aes(x = MDS1, y = MDS2, color = CH4_behavior, shape = cluster)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_point(size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = SITE_ID), size = 3, max.overlaps = 50, show.legend = FALSE) +
  geom_segment(
    data = top_loading_labels,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.12, "inches")),
    color = "grey35",
    alpha = 0.75
  ) +
  ggrepel::geom_text_repel(
    data = top_loading_labels,
    aes(x = PC1, y = PC2, label = label),
    inherit.aes = FALSE,
    size = 3,
    color = "grey20"
  ) +
  scale_color_manual(values = behavior_colors, na.translate = FALSE) +
  labs(
    title = "A. Site Attribute MDS: Observed Behavior Classes",
    subtitle = paste0("Euclidean distance on scaled site attributes; Ward clusters, k = ", best_k),
    x = paste0("MDS1 (", scales::percent(variance_explained[1], accuracy = 1), " eig. share)"),
    y = paste0("MDS2 (", scales::percent(variance_explained[2], accuracy = 1), " eig. share)"),
    color = "Observed behavior",
    shape = "Cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

plot_mds_standardized <- site_scores %>%
  ggplot(aes(x = MDS1, y = MDS2, color = standardized_behavior, shape = cluster)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_point(size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = SITE_ID), size = 3, max.overlaps = 50, show.legend = FALSE) +
  scale_color_manual(values = standardized_colors, na.translate = FALSE) +
  labs(
    title = "B. Same Attribute Space: Model-Standardized Classes",
    x = paste0("MDS1 (", scales::percent(variance_explained[1], accuracy = 1), " eig. share)"),
    y = paste0("MDS2 (", scales::percent(variance_explained[2], accuracy = 1), " eig. share)"),
    color = "Standardized behavior",
    shape = "Cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

plot_cluster_heatmap <- cluster_behavior_table %>%
  reframe(.by = c(cluster, CH4_behavior), n_sites = sum(n_sites)) %>%
  complete(cluster, CH4_behavior, fill = list(n_sites = 0)) %>%
  ggplot(aes(x = cluster, y = CH4_behavior, fill = n_sites)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = if_else(n_sites > 0, as.character(n_sites), "")), fontface = "bold", size = 5) +
  scale_fill_gradient(low = "grey96", high = "blue4") +
  labs(
    title = "C. Cluster Membership by Observed Behavior",
    x = "Ward cluster",
    y = "Observed behavior",
    fill = "Sites"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid = element_blank())

site_attribute_multivariate_figure <- (plot_mds_observed / plot_mds_standardized / plot_cluster_heatmap) +
  plot_layout(heights = c(1, 1, 0.8)) +
  plot_annotation(
    title = "NEON Site Attribute Multivariate Structure",
    subtitle = "Exploratory MDS/clustering of site attributes used to describe CH4 behavior categories",
    caption = "Variables are site-level attributes only; flux-derived behavior summaries are excluded from the distance matrix.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      plot.caption = element_text(size = 9, color = "grey30")
    )
  )

ggsave(
  "FIGURES/NEON_site_attribute_mds_clusters.png",
  plot = site_attribute_multivariate_figure,
  width = 11,
  height = 15,
  units = "in",
  dpi = 300
)

ggsave(
  "FIGURES/NEON_site_attribute_mds_clusters.pdf",
  plot = site_attribute_multivariate_figure,
  width = 11,
  height = 15,
  units = "in"
)

png("FIGURES/NEON_site_attribute_cluster_dendrogram.png", width = 1800, height = 1100, res = 180)
plot(
  cluster_tree,
  labels = site_scores$SITE_ID[match(cluster_tree$labels, site_scores$SITE_ID)],
  main = paste0("Hierarchical Clustering of NEON Site Attributes (k = ", best_k, ")"),
  xlab = "",
  sub = "Ward.D2 clustering on Euclidean distances among scaled site attributes"
)
rect.hclust(cluster_tree, k = best_k, border = "blue4")
dev.off()

pdf("FIGURES/NEON_site_attribute_cluster_dendrogram.pdf", width = 10, height = 7)
plot(
  cluster_tree,
  labels = site_scores$SITE_ID[match(cluster_tree$labels, site_scores$SITE_ID)],
  main = paste0("Hierarchical Clustering of NEON Site Attributes (k = ", best_k, ")"),
  xlab = "",
  sub = "Ward.D2 clustering on Euclidean distances among scaled site attributes"
)
rect.hclust(cluster_tree, k = best_k, border = "blue4")
dev.off()

permanova_lines <- permanova_results %>%
  filter(term != "Total") %>%
  mutate(
    line = paste0(
      "- ", classification, " / ", term,
      ": R2 = ", signif(R2, 3),
      ", F = ", signif(F, 3),
      ", p = ", signif(`Pr(>F)`, 3)
    )
  ) %>%
  pull(line)

silhouette_lines <- silhouette_summary %>%
  mutate(line = paste0("- k = ", k, ": mean silhouette = ", signif(mean_silhouette, 3))) %>%
  pull(line)

loading_lines <- pca_loadings %>%
  slice_head(n = 8) %>%
  mutate(line = paste0("- ", label, ": loading strength = ", signif(loading_strength, 3))) %>%
  pull(line)

writeLines(
  c(
    "# NEON Site Attribute MDS and Cluster Analysis",
    "",
    "## Purpose",
    "This analysis describes whether CH4 behavior categories occupy different regions of multivariate site-attribute space. It is exploratory and does not establish causal drivers.",
    "",
    "## Design",
    paste0("- Sites: ", nrow(site_scores)),
    paste0("- Variables used: ", length(driver_variables)),
    "- Distance: Euclidean distance on centered/scaled site attributes.",
    paste0("- Clustering: Ward.D2 hierarchical clustering; selected k = ", best_k, " by maximum mean silhouette across k = ", min(cluster_range), "-", max(cluster_range), "."),
    "- Ordination: classical multidimensional scaling on the same distance matrix.",
    "- Separation test: PERMANOVA using the same distance matrix.",
    "",
    "## PERMANOVA",
    permanova_lines,
    "",
    "## Cluster Silhouette",
    silhouette_lines,
    "",
    "## Strongest PCA Loadings",
    loading_lines,
    "",
    "## Outputs",
    "- `OUTPUT/NEON_site_attribute_multivariate_matrix.csv`",
    "- `OUTPUT/NEON_site_attribute_variables_used.csv`",
    "- `OUTPUT/NEON_site_attribute_mds_scores.csv`",
    "- `OUTPUT/NEON_site_attribute_pca_loadings.csv`",
    "- `OUTPUT/NEON_site_attribute_cluster_silhouette.csv`",
    "- `OUTPUT/NEON_site_attribute_cluster_behavior_table.csv`",
    "- `OUTPUT/NEON_site_attribute_permanova.csv`",
    "- `FIGURES/NEON_site_attribute_mds_clusters.png`",
    "- `FIGURES/NEON_site_attribute_mds_clusters.pdf`",
    "- `FIGURES/NEON_site_attribute_cluster_dendrogram.png`",
    "- `FIGURES/NEON_site_attribute_cluster_dendrogram.pdf`"
  ),
  "OUTPUT/NEON_site_attribute_mds_cluster_results.md"
)

message("Wrote NEON site-attribute MDS and cluster analysis outputs.")
