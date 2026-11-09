# Empirical uptake analogs for the deliberately forced-arid sink scenario.
# Inputs are monthly gapfilled site budgets, not independent chamber observations.
# The dry donor cutoff uses the pipeline's MAP/(MAT+10) index. It is wider
# than the forced-Arid cutoff (15) because those sites have no sink months.
fit_dry_sink_analogs <- function(training, dry_index_max = 20, k = 12L,
                                min_rows = 6L, min_sites = 2L,
                                max_climate_rms = 2) {
  features <- c("mean_ERA5_Tair_C", "mean_ERA5_VSWC", "MAT", "MAP")
  ok <- training$weak_source == 0 & training$monthly_flux_gC_m2_month < 0 &
    is.finite(training$aridity_index) & training$aridity_index <= dry_index_max
  ok <- ok & Reduce(`&`, lapply(training[features], is.finite))
  if ("n_finite_halfhours" %in% names(training)) {
    days <- as.integer(format(as.Date(paste0(training$Year + (training$month == 12), "-",
      sprintf("%02d", training$month %% 12 + 1), "-01")) - 1, "%d"))
    ok <- ok & training$n_finite_halfhours == days * 48L
  }
  donors <- as.data.frame(training[which(ok), c("SITE_ID", "Year", "month",
    "monthly_flux_gC_m2_month", "aridity_index", features,
    intersect("n_finite_halfhours", names(training)))])
  donors$SITE_ID <- as.character(donors$SITE_ID)
  if (nrow(donors) < min_rows || length(unique(donors$SITE_ID)) < min_sites)
    stop("Insufficient dry weak-sink donors; do not substitute wetter-system uptake.")
  days <- as.integer(format(as.Date(paste0(donors$Year + (donors$month == 12), "-",
    sprintf("%02d", donors$month %% 12 + 1), "-01")) - 1, "%d"))
  donors$daily_flux_gC_m2 <- donors$monthly_flux_gC_m2_month / days
  center <- vapply(donors[features], mean, numeric(1))
  scale <- pmax(vapply(donors[features], sd, numeric(1)), c(5, 0.05, 5, 300))
  names(scale) <- features
  list(donors = donors, features = features, center = center, scale = scale,
       dry_index_max = dry_index_max, k = k, min_rows = min_rows,
       min_sites = min_sites, max_climate_rms = max_climate_rms,
       method = "dry_weak_sink_analogs_v1")
}

predict_dry_sink_analogs <- function(fit, newdata) {
  if (!requireNamespace("FNN", quietly = TRUE)) stop("Dry analog matching requires FNN.")
  n <- nrow(newdata)
  out <- data.frame(flux = rep(NA_real_, n), climate_supported = rep(FALSE, n),
    season_matched = rep(FALSE, n), nearest_climate_rms = rep(NA_real_, n),
    donor_sites = integer(n), method = rep("", n))
  if (!n) return(out)
  features <- fit$features
  query <- as.matrix(newdata[features])
  if (any(!is.finite(query))) stop("Nonfinite climate supplied to dry analog matcher.")
  target_month <- as.integer(newdata$month)
  if (any(!target_month %in% 1:12)) stop("Invalid target month.")
  # Seasonal phase follows hemisphere; donor sites are Northern Hemisphere.
  phase <- target_month
  if ("y" %in% names(newdata)) {
    south <- is.finite(newdata$y) & newdata$y < 0
    phase[south] <- (phase[south] + 5L) %% 12L + 1L
  }
  target_days <- as.integer(format(as.Date(paste0(newdata$Year + (target_month == 12),
    "-", sprintf("%02d", target_month %% 12 + 1), "-01")) - 1, "%d"))
  standardize <- function(x) sweep(sweep(as.matrix(x), 2, fit$center, "-"), 2, fit$scale, "/")
  query <- standardize(query)
  donors <- fit$donors
  donors$SITE_ID <- as.character(donors$SITE_ID)
  for (m in unique(phase)) {
    rows <- which(phase == m)
    delta <- abs(donors$month - m)
    pool <- donors[pmin(delta, 12 - delta) <= 1, , drop = FALSE]
    seasonal <- nrow(pool) >= fit$min_rows && length(unique(pool$SITE_ID)) >= fit$min_sites
    if (!seasonal) pool <- donors
    knn <- FNN::get.knnx(standardize(pool[features]), query[rows, , drop = FALSE],
                        k = min(fit$k, nrow(pool)))
    ix <- knn$nn.index
    distance <- knn$nn.dist / sqrt(length(features))
    # Equal total influence per donor site before distance weighting.
    rates <- matrix(pool$daily_flux_gC_m2[ix], nrow = nrow(ix))
    sites <- matrix(as.character(pool$SITE_ID[ix]), nrow = nrow(ix))
    weights <- 1 / (1 + distance^2)
    nsites <- integer(nrow(ix))
    for (site in unique(pool$SITE_ID)) {
      hit <- sites == site
      count <- rowSums(hit)
      nsites <- nsites + (count > 0)
      weights[hit] <- weights[hit] / rep(pmax(count, 1), ncol(ix))[hit]
    }
    supported <- distance[, 1] <= fit$max_climate_rms & nsites >= fit$min_sites
    rate <- rowSums(rates * weights) / rowSums(weights)
    fallback <- mean(vapply(split(pool$daily_flux_gC_m2, pool$SITE_ID, drop = TRUE), mean, numeric(1)))
    rate[!supported] <- fallback
    out$flux[rows] <- rate * target_days[rows]
    out$climate_supported[rows] <- supported
    out$season_matched[rows] <- seasonal
    out$nearest_climate_rms[rows] <- distance[, 1]
    out$donor_sites[rows] <- ifelse(supported, nsites, length(unique(pool$SITE_ID)))
    out$method[rows] <- ifelse(supported,
      if (seasonal) "climate_and_season_analogs" else "climate_analogs_all_seasons",
      if (seasonal) "dry_seasonal_mean_fallback" else "dry_mean_fallback")
  }
  stopifnot(all(is.finite(out$flux)), all(out$flux < 0))
  out
}
