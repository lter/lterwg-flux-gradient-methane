# Retrospective missing-block validation. All flux-derived predictors must be
# rebuilt after removing the complete outer test block. Weather is unchanged.
audit_fit_climatology <- function(d) {
  keys <- list(site_month_hour = c("SITE_ID", "month", "hour_num"),
    site_season_hour = c("SITE_ID", "season", "hour_num"),
    site_hour = c("SITE_ID", "hour_num"), global_month_hour = c("month", "hour_num"))
  tables <- lapply(keys, function(k) d[, .(clim_sum = sum(CH4_mgC_30min),
                                           clim_n = .N), by = k])
  list(keys = keys, tables = tables, n_training = nrow(d))
}

audit_apply_climatology <- function(d, fit, leave_one_out = FALSE) {
  value <- rep(NA_real_, nrow(d)); level <- rep(NA_character_, nrow(d))
  donors <- integer(nrow(d))
  for (name in names(fit$keys)) {
    z <- fit$tables[[name]][d, on = fit$keys[[name]], .(clim_sum, clim_n)]
    if (leave_one_out) {
      # Only valid when every row in d contributed exactly once to fit.
      contributed <- rep(TRUE, nrow(d))
      if (name == "site_season_hour" && "clim_donor_season" %in% names(d))
        contributed <- d$season == d$clim_donor_season
      z[, `:=`(clim_sum = clim_sum - ifelse(contributed, d$CH4_mgC_30min, 0),
                clim_n = clim_n - as.integer(contributed))]
    }
    next_value <- z$clim_sum / z$clim_n
    use <- !is.finite(value) & z$clim_n > 0 & is.finite(next_value)
    use[is.na(use)] <- FALSE
    value[use] <- next_value[use]; level[use] <- name; donors[use] <- z$clim_n[use]
  }
  data.table::data.table(climatological_CH4_mgC_30min = value,
                         climatology_level = level, climatology_donors = donors)
}

audit_make_design <- function(d, scenarios, repeats = 3L, seed = 20261021L,
                              min_test_obs = 8L) {
  records <- list(); eligibility <- list(); nr <- 0L; ne <- 0L
  origin <- as.POSIXct("2021-01-01", tz = "UTC")
  first <- min(d$time.rounded); last <- max(d$time.rounded)
  for (j in seq_along(scenarios)) {
    scenario <- scenarios[j]; set.seed(seed + j)
    for (site in sort(unique(d$SITE_ID))) {
      s <- d[SITE_ID == site]
      if (scenario == "random_rows") {
        n_test <- floor(nrow(s) * 0.2)
        ne <- ne + 1L
        eligibility[[ne]] <- data.table::data.table(scenario, SITE_ID = site,
          n_observed = nrow(s), eligible_blocks = if (n_test >= min_test_obs) repeats else 0L)
        if (n_test < min_test_obs) next
        for (rep in seq_len(repeats)) {
          nr <- nr + 1L
          records[[nr]] <- data.table::data.table(scenario, repetition = rep, SITE_ID = site,
            block_id = paste0("random_", rep), start = as.POSIXct(NA, tz = "UTC"),
            end = as.POSIXct(NA, tz = "UTC"), duration_days = NA_real_,
            row_ids = list(sample(s$row_id, n_test)))
        }
      } else {
        if (scenario == "site_year") {
          s[, block_id := format(time.rounded, "%Y")]
          s[, start := as.POSIXct(paste0(block_id, "-01-01"), tz = "UTC")]
          s[, end := as.POSIXct(paste0(as.integer(block_id) + 1L, "-01-01"), tz = "UTC")]
        } else {
          days <- as.numeric(sub("d$", "", scenario))
          stopifnot(is.finite(days), days > 0)
          s[, block_id := as.character(floor(as.numeric(difftime(time.rounded, origin, units = "days")) / days))]
          s[, start := origin + as.numeric(block_id) * days * 86400]
          s[, end := start + days * 86400]
        }
        candidates <- s[, .(start = start[1], end = end[1], n_test = .N,
                             row_ids = list(row_id)), by = block_id]
        candidates <- candidates[n_test >= min_test_obs & start >= first & end <= last]
        ne <- ne + 1L
        eligibility[[ne]] <- data.table::data.table(scenario, SITE_ID = site,
          n_observed = nrow(s), eligible_blocks = nrow(candidates))
        if (!nrow(candidates)) next
        chosen <- sample.int(nrow(candidates), min(repeats, nrow(candidates)))
        for (rep in seq_along(chosen)) {
          c <- candidates[chosen[rep]]; nr <- nr + 1L
          records[[nr]] <- data.table::data.table(scenario, repetition = rep, SITE_ID = site,
            block_id = c$block_id, start = c$start, end = c$end,
            duration_days = as.numeric(difftime(c$end, c$start, units = "days")), row_ids = c$row_ids)
        }
      }
    }
  }
  list(blocks = data.table::rbindlist(records), eligibility = data.table::rbindlist(eligibility))
}

audit_score <- function(y, p) {
  ok <- is.finite(y) & is.finite(p); total <- length(y); y <- y[ok]; p <- p[ok]
  n <- length(y)
  if (!n) return(list(n_scored = 0L, n_test = total, rmse = NA_real_, mae = NA_real_,
    bias = NA_real_, r = NA_real_, r2 = NA_real_, sign_accuracy = NA_real_,
    observed_sum_gC_m2 = NA_real_, predicted_sum_gC_m2 = NA_real_))
  list(n_scored = n, n_test = total, rmse = sqrt(mean((p-y)^2)), mae = mean(abs(p-y)),
    bias = mean(p-y), r = if (n > 2 && sd(y)>0 && sd(p)>0) cor(y,p) else NA_real_,
    r2 = if (n > 1 && sum((y-mean(y))^2)>0) 1-sum((p-y)^2)/sum((y-mean(y))^2) else NA_real_,
    sign_accuracy = mean(sign(y)==sign(p)), observed_sum_gC_m2 = sum(y)/1000,
    predicted_sum_gC_m2 = sum(p)/1000)
}
