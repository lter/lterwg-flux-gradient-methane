# 03_NEON_EnsembleNoEC.R ------------------------------------------------------------------------
# Nov 2026: reads the v3b fluxes (NEON_GradientFlux_EDI_2021_2026_v3b). With the CH4 pair rule of 28 Oct 2026, the
# ensemble written here is replaced by steps 02-03; this script supplies the filtered + despiked CH4 data
# (SITE_DATA_FILTERED_CH4.Rdata) and the no-EC diagnostics (NEON_CH4_noEC_*.csv).
# Rebuilds SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE for CH4 from the 2021-2026 gradient fluxes using the
# EC-free pair selection of Supplement S2 (lterwg-flux-gradient-eval, flow.evaluation.noEC_pair_selection.R).
# Replaces steps 01-03 for the CH4 ensemble. (The previous 02/03 selected pairs on CO2/H2O CCC against EC
# and averaged FG_mean over ALL gases in SITE_DATA_FILTERED - CO2, H2O and CH4 together.)
#
# Per site, CH4 only:
#  1. Filter (same rule as filter_fluxes(): dConc SNR >= 3, tracer SNR >= 3 for MBR, u* >= threshold),
#     drop WW pairs (canopy_commbined.csv), then despike (despike_fluxes(): median +/- 6 MAD per
#     series x day/night (07-18 h local) x season; MBR series are pair x tracer).
#  2. 30-min values per series (series = approach x pair, x tracer for MBR); time.rounded as in the old
#     step 03 (round_date to 30 min, local time).
#  3. Reference (REF) = equal-weight 30-min mean of above-canopy (AA) MBR series. Which MBR series:
#       - per tracer, check the AA MBR pairs for a sign split (pair medians of opposite sign);
#       - if neither tracer is split, use both tracers (tracer_CCC between the two is reported only);
#       - if a tracer is split, candidates are each tracer, each sign group of a split tracer, and both; the candidate
#         whose REF agrees best with the independent AE/WP pairs (mean of the 3 highest correlations,
#         >= min.windows shared windows each) is chosen.
#     No AA MBR at a site -> AW MBR is used as the reference and the site is flagged.
#  4. AE/WP series (AA or AW) with >= min.windows shared windows and CCC vs REF >= ccc.cut are kept and
#     rescaled to REF by reduced-major-axis:  b = sign(r) sd(F)/sd(REF), a = mean(F) - b mean(REF),
#     F' = (F - a)/b.  ccc.cut = 0.45 from the CO2/H2O threshold scan (Spike_Eval/noEC_threshold_scan_*).
#  5. FG_ENSEMBLE_RSHP = equal-weight 30-min mean of the REF series and the scaled AE/WP series.
#
# Outputs (localdir.ch4 = METHANE/):
#   SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata   same list/columns as before + FG_REF_MBR, n_series, site_flag
#   NEON_GradientFlux_Data_Filter/SITES_<MBR|AE|WP>_9min_FILTER_CH4.Rdata, SITE_DATA_FILTERED_CH4.Rdata
#       (CH4 only, filtered + despiked, 2021-2026; same columns as before + tracer)
#   NEON_CH4_noEC_reference_by_site.csv, NEON_CH4_noEC_pairs.csv, NEON_CH4_noEC_reference_candidates.csv
# Files that are overwritten are first copied to METHANE/archive_<archive.tag>/ (once).
# -------------------------------------------------------------------------------------------------
library(data.table)
library(lutz)
library(sf)

if (!exists("localdir"))     localdir     <- "/Volumes/MaloneLab/Research/FluxGradient"
if (!exists("localdir.ch4")) localdir.ch4 <- file.path(localdir, "METHANE")
if (!exists("gf.dir"))       gf.dir       <- file.path(localdir, "NEON_GradientFlux_EDI_2021_2026_v3b")  # per-site <SITE>_<APP>_9min.df.final.csv (v3b, Nov 2026)
if (!exists("ccc.cut"))      ccc.cut      <- 0.45
if (!exists("min.windows"))  min.windows  <- 30
if (!exists("snr.min"))      snr.min      <- 3
if (!exists("despike.k"))    despike.k    <- 6
if (!exists("archive.tag"))  archive.tag  <- "2026-11_pre_v3b"   # Oct 2026 (pre-v3b) CH4 files are copied here before overwrite
if (!exists("sites.run"))    sites.run    <- NULL   # e.g. c("KONZ","TALL") to test; NULL = all sites
if (!exists("checkpoint.dir")) checkpoint.dir <- NULL  # optional local per-site restart cache
if (!exists("noec.checkpoint.only")) noec.checkpoint.only <- FALSE  # TRUE: only write per-site checkpoints (cluster array jobs); a later
                                                                  # run over all sites restores them and writes the combined files

# ---- helpers -------------------------------------------------------------------------------------
lccc <- function(x, y) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(NA_real_); mx <- mean(x); my <- mean(y)
  2 * mean((x - mx) * (y - my)) / (mean((x - mx)^2) + mean((y - my)^2) + (mx - my)^2) }
rcor <- function(x, y) { ok <- is.finite(x) & is.finite(y); if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok])) }
season_of <- function(m) fcase(m %in% c(12, 1, 2), "Winter", m %in% 3:5, "Spring", m %in% 6:8, "Summer", default = "Autumn")

gf.cols <- c("timeEnd_A", "dLevelsAminusB", "FG_mean", "gas", "tracer", "dConc", "dConc_sd", "dConc.tracer",
             "dConc.tracer_sd", "ustar_interp", "ustar_threshold", "TowerPosition_A", "TowerPosition_B",
             "PAR", "Tair_K", "RH", "cross_grad_flag")

# 1. read + filter one approach for one site (CH4 rows only)
read_ch4 <- function(site, app, tz) {
  f <- file.path(gf.dir, site, sprintf("%s_%s_9min.df.final.csv", site, app))
  if (!file.exists(f)) return(NULL)
  message("  reading ", app, " (", round(file.info(f)$size / 1024^3, 2), " GiB)")
  hdr <- names(fread(f, nrows = 0))
  d <- fread(f, select = intersect(gf.cols, hdr), showProgress = FALSE)[gas == "CH4"]
  for (cc in setdiff(gf.cols, names(d))) set(d, j = cc, value = if (cc == "tracer") NA_character_ else NA_real_)
  if (app != "MBR") { set(d, j = c("tracer", "dConc.tracer", "dConc.tracer_sd"), value = NULL)
                      set(d, j = "tracer", value = ""); set(d, j = "dConc.tracer", value = as.numeric(snr.min))
                      set(d, j = "dConc.tracer_sd", value = 1) }
  d <- d[TowerPosition_A != TowerPosition_B &
         abs(dConc) / dConc_sd >= snr.min &
         abs(as.numeric(dConc.tracer)) / as.numeric(dConc.tracer_sd) >= snr.min &
         ustar_interp >= ustar_threshold & is.finite(FG_mean)]
  t.utc <- if (inherits(d$timeEnd_A, "POSIXt")) as.POSIXct(d$timeEnd_A, tz = "UTC") else
    as.POSIXct(d$timeEnd_A, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
  # Preserve the instant while changing its display zone; character round trips
  # can discard subdaily times when the first observation is at midnight.
  d[, timeEndA.local := as.POSIXct(as.numeric(t.utc), origin = "1970-01-01", tz = tz)]
  d[, `:=`(Approach = app, tracer = as.character(tracer), Tair_C = Tair_K - 273.15)]
  message("  ", app, ": ", nrow(d), " CH4 rows after quality filters")
  d[]
}

# despike (same rule as despike_fluxes() in lterwg-flux-gradient-eval/functions/calc.filter_FG.R)
despike <- function(d) {
  d[, `:=`(.dn = as.integer(format(timeEndA.local, "%H")) %in% 7:18,
           .ss = season_of(as.integer(format(timeEndA.local, "%m"))))]
  d[, `:=`(.med = median(FG_mean), .mad = mad(FG_mean), .n = .N), by = .(Approach, dLevelsAminusB, tracer, .dn, .ss)]
  n0 <- nrow(d)
  d <- d[.n < 20 | abs(FG_mean - .med) <= despike.k * pmax(.mad, 1e-6)]
  d[, c(".dn", ".ss", ".med", ".mad", ".n") := NULL]
  list(d = d, spike_pct = 100 * (1 - nrow(d) / n0))
}

# 3. choose the MBR reference
choose_ref <- function(P) {
  M <- P[Approach == "MBR" & Canopy_L1 == "AA"]; ref.level <- "AA"
  if (!nrow(M)) { M <- P[Approach == "MBR" & Canopy_L1 == "AW"]; ref.level <- "AW" }
  if (!nrow(M)) return(list(pick = NA_character_, series = character(0), ref.level = NA_character_,
                            split = NA, tracer_ccc = NA_real_, cands = data.table()))
  O  <- P[Approach != "MBR"]
  pm <- M[, .(med = median(F), n = .N), by = .(tracer, dLevelsAminusB, series)][, sgn := sign(med)]
  split.by.tracer <- pm[, .(split = uniqueN(sgn[sgn != 0]) > 1), by = tracer]
  cands <- list()
  for (tr in sort(unique(pm$tracer))) {
    cands[[tr]] <- pm[tracer == tr, series]
    if (split.by.tracer[tracer == tr, split])
      for (g in c(1, -1)) if (any(pm[tracer == tr, sgn] == g))
        cands[[paste0(tr, ifelse(g > 0, "+", "-"))]] <- pm[tracer == tr & sgn == g, series]
  }
  if (uniqueN(pm$tracer) > 1) cands[["both"]] <- pm$series
  RT <- dcast(M[, .(R = mean(F)), by = .(tracer, w)], w ~ tracer, value.var = "R")
  tracer_ccc <- if (all(c("CO2", "H2O") %in% names(RT))) lccc(RT$CO2, RT$H2O) else NA_real_
  sc <- rbindlist(lapply(names(cands), function(k) {
    R <- M[series %in% cands[[k]], .(REF = mean(F)), by = w]
    a <- O[, { m <- merge(.SD[, .(w, F)], R, by = "w")
               if (nrow(m) >= min.windows) .(r = rcor(m$F, m$REF), ccc = lccc(m$F, m$REF)) else .(r = NA_real_, ccc = NA_real_) },
           by = series][is.finite(r)]
    data.table(candidate = k, n_series = length(cands[[k]]), n_windows = nrow(R), median_REF = median(R$REF),
               n_AEWP_compared = nrow(a),
               top3_r = if (nrow(a)) mean(head(sort(a$r, decreasing = TRUE), 3)) else NA_real_,
               max_ccc = if (nrow(a)) max(a$ccc) else NA_real_)
  }))
  any.split <- any(split.by.tracer$split)
  pick <- if (!any.split) {
    if (uniqueN(pm$tracer) > 1) "both" else unique(pm$tracer)
  } else if (any(is.finite(sc$top3_r))) sc[which.max(top3_r), candidate] else if ("both" %in% names(cands)) "both" else names(cands)[1]
  list(pick = pick, series = cands[[pick]], ref.level = ref.level, split = any.split,
       split_tracers = paste(split.by.tracer[split == TRUE, tracer], collapse = ","),
       tracer_ccc = tracer_ccc, cands = sc, pair_medians = pm)
}

# 4-5. S2 selection, scaling, ensemble
s2_ensemble <- function(P, ref.series) {
  REF <- P[series %in% ref.series, .(REF = mean(F)), by = w]
  cand <- P[Approach != "MBR"]
  sel <- cand[, { m <- merge(.SD[, .(w, F)], REF, by = "w"); n <- nrow(m)
      if (n >= 3) { r <- rcor(m$F, m$REF); b <- sign(r) * sd(m$F) / sd(m$REF)
        .(n_shared = n, CCC_ref = lccc(m$F, m$REF), r_ref = r, scale_b = b, scale_a = mean(m$F) - b * mean(m$REF)) }
      else .(n_shared = n, CCC_ref = NA_real_, r_ref = NA_real_, scale_b = NA_real_, scale_a = NA_real_) },
    by = .(series, Approach, dLevelsAminusB, Canopy_L1)]
  sel[, selected := n_shared >= min.windows & is.finite(CCC_ref) & CCC_ref >= ccc.cut & is.finite(scale_b) & scale_b != 0]
  K <- merge(cand, sel[selected == TRUE, .(series, scale_a, scale_b)], by = "series")[, F := (F - scale_a) / scale_b]
  E <- rbind(P[series %in% ref.series, .(w, F)], K[, .(w, F)])[, .(FG_ENSEMBLE_RSHP = mean(F), n_series = .N), by = w]
  E <- merge(E, REF, by = "w", all = TRUE)
  list(E = E, sel = sel)
}

# ---- site loop -----------------------------------------------------------------------------------
metadata  <- read.csv(file.path(localdir, "Ameriflux_NEON field-sites.csv"))
site.list <- if (is.null(sites.run)) unique(metadata$Site_Id.NEON) else sites.run
sites.location <- st_as_sf(metadata, coords = c("Longitude..degrees.", "Latitude..degrees."), crs = "+proj=longlat +datum=WGS84")
sites.location$TZ <- tz_lookup(sites.location, method = "accurate")
canopy <- unique(fread(file.path(localdir, "canopy_commbined.csv"))[, .(Site, dLevelsAminusB, Canopy_L1)])

SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE <- list()
SITES_MBR_9min_FILTER <- list(); SITES_AE_9min_FILTER <- list(); SITES_WP_9min_FILTER <- list()
SITE_DATA_FILTERED <- list()
ref.tab <- list(); pair.tab <- list(); cand.tab <- list()
if (!is.null(checkpoint.dir)) dir.create(checkpoint.dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_helpers <- lapply(list(read_ch4, despike, choose_ref, s2_ensemble, lccc, rcor, season_of),
  function(f) deparse(body(f)))

for (site in site.list) {
  message("CH4 no-EC ensemble: ", site)
  tz <- sites.location$TZ[sites.location$Site_Id.NEON == site][1]
  checkpoint.file <- if (!is.null(checkpoint.dir)) file.path(checkpoint.dir, paste0(site, ".rds")) else NULL
  signature <- NULL
  if (!is.null(checkpoint.file)) {
    inputs <- file.path(gf.dir, site, sprintf("%s_%s_9min.df.final.csv", site, c("MBR", "AE", "WP")))
    signature <- list(inputs = file.info(inputs)[, c("size", "mtime")],
      timezone = tz, canopy = as.data.frame(canopy[Site == site]),
      config = list(ccc.cut, min.windows, snr.min, despike.k, gf.cols), helpers = checkpoint_helpers)
    cached <- if (file.exists(checkpoint.file))
      tryCatch(readRDS(checkpoint.file), error = function(e) NULL) else NULL
    if (!is.null(cached) && identical(cached$signature, signature)) {
      SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]] <- cached$ensemble
      X9 <- cached$filtered
      SITE_DATA_FILTERED[[site]] <- X9
      SITES_MBR_9min_FILTER[[site]] <- X9[X9$Approach == "MBR", ]
      SITES_AE_9min_FILTER[[site]] <- X9[X9$Approach == "AE", ]
      SITES_WP_9min_FILTER[[site]] <- X9[X9$Approach == "WP", ]
      ref.tab[[site]] <- cached$reference
      pair.tab[[site]] <- cached$pairs
      cand.tab[[site]] <- cached$candidates
      message("  restored site checkpoint")
      if (exists("cleanup_noec_staged_inputs", mode = "function")) cleanup_noec_staged_inputs(site)
      rm(cached, X9); invisible(gc())
      next
    }
    rm(cached)
  }
  D <- rbindlist(lapply(c("MBR", "AE", "WP"), function(a) read_ch4(site, a, tz)), fill = TRUE)
  if (!nrow(D)) {
    message("  no CH4 data")
    if (exists("cleanup_noec_staged_inputs", mode = "function")) cleanup_noec_staged_inputs(site)
    next
  }
  D <- merge(D, canopy[Site == site, .(dLevelsAminusB, Canopy_L1)], by = "dLevelsAminusB")[Canopy_L1 != "WW"]
  ds <- despike(D); D <- ds$d; spike.pct <- ds$spike_pct; rm(ds)
  D[, series := trimws(paste(Approach, dLevelsAminusB, tracer))]
  D[, time.rounded := lubridate::round_date(timeEndA.local, unit = "30 minutes")]
  D[, w := as.numeric(time.rounded)]

  # filtered 9-min CH4 objects (same columns as the old step 01 output, + tracer)
  keep9 <- c("timeEndA.local", "FG_mean", "Approach", "Canopy_L1", "gas", "TowerPosition_A", "TowerPosition_B",
             "cross_grad_flag", "dLevelsAminusB", "RH", "PAR", "Tair_K", "Tair_C", "tracer")
  X9 <- as.data.frame(D[, ..keep9][, EC_mean := NA_real_])
  SITES_MBR_9min_FILTER[[site]] <- X9[X9$Approach == "MBR", ]
  SITES_AE_9min_FILTER[[site]]  <- X9[X9$Approach == "AE", ]
  SITES_WP_9min_FILTER[[site]]  <- X9[X9$Approach == "WP", ]
  SITE_DATA_FILTERED[[site]]    <- X9

  P <- D[, .(F = mean(FG_mean)), by = .(series, Approach, dLevelsAminusB, tracer, Canopy_L1, w)]
  ch <- choose_ref(P)
  met <- D[, .(PAR = mean(PAR, na.rm = TRUE), Tair_C = mean(Tair_C, na.rm = TRUE), n_9min = .N), by = .(w, time.rounded)]

  if (length(ch$series)) {
    res <- s2_ensemble(P, ch$series)
    E <- merge(met, res$E, by = "w", all.x = TRUE)
    pair.tab[[site]] <- cbind(Site = site, res$sel)
  } else {
    E <- met[, `:=`(FG_ENSEMBLE_RSHP = NA_real_, n_series = 0L, REF = NA_real_)]
    res <- list(sel = data.table(selected = logical(0), series = character(0)))
  }
  if (nrow(ch$cands)) cand.tab[[site]] <- cbind(Site = site, ch$cands, chosen = ch$cands$candidate == ch$pick)

  flag <- fcase(is.na(ch$pick), "no_MBR_reference",
                ch$ref.level == "AW", "reference_AW_MBR",
                isTRUE(ch$split), "reference_chosen_by_AEWP_agreement",
                default = "ok")
  if (flag == "reference_chosen_by_AEWP_agreement" && sum(res$sel$selected) == 0) flag <- "reference_unconfirmed_no_AEWP_agreement"

  ref.tab[[site]] <- data.table(Site = site, tz = tz, reference = ch$pick, reference_level = ch$ref.level,
    reference_series = paste(ch$series, collapse = "; "), tracer_CCC = ch$tracer_ccc,
    sign_split = ch$split, split_tracers = if (is.null(ch$split_tracers)) NA_character_ else ch$split_tracers,
    n_AEWP_selected = sum(res$sel$selected), AEWP_selected = paste(res$sel[selected == TRUE, series], collapse = "; "),
    windows_REF = sum(is.finite(E$REF)), windows_ensemble = sum(is.finite(E$FG_ENSEMBLE_RSHP)),
    median_REF = median(E$REF, na.rm = TRUE), median_ensemble = median(E$FG_ENSEMBLE_RSHP, na.rm = TRUE),
    spike_pct = spike.pct, site_flag = flag)

  setorder(E, time.rounded)
  E[, `:=`(month = as.numeric(format(time.rounded, "%m")), hour = format(time.rounded, "%H"))]
  E[, season := season_of(month)]
  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]] <- as.data.frame(E[, .(time.rounded, FG_ENSEMBLE_RSHP, count = n_9min,
    PAR, Tair_C, month, season, hour, FG_REF_MBR = REF, n_series = fifelse(is.na(n_series), 0L, n_series), site_flag = flag)])
  if (!is.null(checkpoint.file)) {
    cached <- list(signature = signature, ensemble = SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]],
      filtered = X9, reference = ref.tab[[site]], pairs = pair.tab[[site]], candidates = cand.tab[[site]])
    tmp <- tempfile(tmpdir = checkpoint.dir, fileext = ".rds")
    saveRDS(cached, tmp)
    if (!file.rename(tmp, checkpoint.file)) stop("Cannot finalize site checkpoint: ", site)
    message("  saved site checkpoint")
    if (exists("cleanup_noec_staged_inputs", mode = "function")) cleanup_noec_staged_inputs(site)
    rm(cached)
  }
  rm(D, P, X9, E); invisible(gc())
}

if (isTRUE(noec.checkpoint.only)) {
  if (is.null(checkpoint.dir)) stop("noec.checkpoint.only needs checkpoint.dir")
  message("Checkpoints written for: ", paste(site.list, collapse = ", "))
} else {

ref.tab  <- rbindlist(ref.tab, fill = TRUE)
pair.tab <- rbindlist(pair.tab, fill = TRUE)
cand.tab <- rbindlist(cand.tab, fill = TRUE)
print(ref.tab[, .(Site, reference, tracer_CCC = round(tracer_CCC, 2), sign_split, n_AEWP_selected, windows_ensemble, site_flag)])
print(ref.tab[, .N, by = site_flag])

# ---- save (archive the files being replaced, once) -----------------------------------------------
if (is.null(sites.run)) {
  arch <- file.path(localdir.ch4, paste0("archive_", archive.tag)); dir.create(arch, showWarnings = FALSE)
  old.files <- c(file.path(localdir.ch4, "SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata"),
                 file.path(localdir.ch4, "SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata"),
                 file.path(localdir.ch4, "NEON_RSHP_CCC.Rdata"),
                 file.path(localdir.ch4, "NEON_GradientFlux_Data_Filter",
                           c("SITE_DATA_FILTERED_CH4.Rdata", paste0("SITES_", c("MBR", "AE", "WP"), "_9min_FILTER_CH4.Rdata"))))
  for (f in old.files) if (file.exists(f) && !file.exists(file.path(arch, basename(f)))) file.copy(f, file.path(arch, basename(f)))
  out.ens <- file.path(localdir.ch4, "SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata")
  fdir <- file.path(localdir.ch4, "NEON_GradientFlux_Data_Filter")
} else {
  out.ens <- file.path(localdir.ch4, "SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TEST.Rdata")
  fdir <- file.path(localdir.ch4, "NEON_GradientFlux_Data_Filter_TEST"); dir.create(fdir, showWarnings = FALSE)
}
save(SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE, file = out.ens)
save(SITE_DATA_FILTERED,    file = file.path(fdir, "SITE_DATA_FILTERED_CH4.Rdata"))
save(SITES_MBR_9min_FILTER, file = file.path(fdir, "SITES_MBR_9min_FILTER_CH4.Rdata"))
save(SITES_AE_9min_FILTER,  file = file.path(fdir, "SITES_AE_9min_FILTER_CH4.Rdata"))
save(SITES_WP_9min_FILTER,  file = file.path(fdir, "SITES_WP_9min_FILTER_CH4.Rdata"))
sfx <- if (is.null(sites.run)) "" else "_TEST"
fwrite(ref.tab,  file.path(localdir.ch4, paste0("NEON_CH4_noEC_reference_by_site", sfx, ".csv")))
fwrite(pair.tab, file.path(localdir.ch4, paste0("NEON_CH4_noEC_pairs", sfx, ".csv")))
fwrite(cand.tab, file.path(localdir.ch4, paste0("NEON_CH4_noEC_reference_candidates", sfx, ".csv")))
message("Saved ", out.ens)
}
