# Compare Stage-1 classifier predictor support between NEON and the global
# projection. NEON prevalence is not used as a global prior or calibration.

library(data.table)
library(terra)
library(FNN)

root <- Sys.getenv("LOCALDIR_CH4", unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE")
spatial_dir <- Sys.getenv("MONTHLY_UPSCALING_DIR", unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly")
rfdir <- Sys.getenv("MONTHLY_RF_DIR", unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF")
out <- file.path(rfdir, "OUTPUT")
pred_file <- file.path(out, "monthly_cell_predictions_2000_2025.rds")
cache_file <- file.path(out, "mat_map_climatology.tif")
era_dir <- file.path(spatial_dir, "DATA/era5_land_monthly")
set.seed(20261021)
numeric_predictors <- c("mean_ERA5_Tair_C", "mean_ERA5_VSWC", "MAP", "MAT", "aridity_index")
wquantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

# NEON site-month predictors and outcomes, using the same aggregation and
# aridity rule as Stage 1. The outcomes are retained only for later validation.
attrs <- fread(file.path(root, "OUTPUT/NEON_scale_annual_budget_summary.csv"))
attrs[, SITE_ID := as.character(SITE_ID)]
attrs <- attrs[!is.na(EcoType) &
  !grepl("wetland|inundat|flood|marsh|swamp|bog|fen|lake|rice|crop|agri", EcoType, ignore.case=TRUE),
  .(SITE_ID, EcoType, MAP, MAT)]
attrs <- unique(attrs)
neon_raw <- fread(file.path(root, "OUTPUT/NEON_ERA5_gapfilled_30min.csv.gz"),
  select=c("SITE_ID","EcoType","Year","month","ERA5_Tair_C","ERA5_VSWC","gapfilled_CH4_mgC_30min"))
neon_raw[, SITE_ID := as.character(SITE_ID)]
neon_raw <- merge(neon_raw, attrs, by=c("SITE_ID","EcoType"))
neon <- neon_raw[is.finite(ERA5_Tair_C) & is.finite(ERA5_VSWC), .(
  budget=sum(gapfilled_CH4_mgC_30min, na.rm=TRUE),
  mean_ERA5_Tair_C=mean(ERA5_Tair_C, na.rm=TRUE),
  mean_ERA5_VSWC=mean(ERA5_VSWC, na.rm=TRUE),
  MAP=first(MAP), MAT=first(MAT)
), by=.(SITE_ID,EcoType,Year,month)]
neon[, aridity_index := fifelse(MAT > -9, MAP/(MAT+10), NA_real_)]
neon[, is_arid := as.integer(is.finite(aridity_index) & aridity_index < 15)]
neon <- neon[is_arid == 0 & is.finite(aridity_index)]
neon[, observed_source := as.integer(budget > 0)]

# Draw an area-proportional sample within each non-Arid ecosystem from all
# global monthly predictions, then recover the exact predictors used to make
# those predictions from the underlying ERA5-Land and climatology rasters.
message("Reading spatial prediction index...")
pred <- readRDS(pred_file)
setDT(pred)
support_years <- as.integer(strsplit(Sys.getenv("SUPPORT_YEARS", unset="2000,2010,2020,2025"), ",", fixed=TRUE)[[1]])
support_months <- as.integer(strsplit(Sys.getenv("SUPPORT_MONTHS", unset=paste(1:12,collapse=",")), ",", fixed=TRUE)[[1]])
support_tag <- paste0("years_",paste(support_years,collapse="-"),"_months_",paste(support_months,collapse="-"))
pred <- pred[Year %in% support_years & month %in% support_months & is_arid == 0 & area_mha > 0 & EcoType %chin% c("Forest","Grassland","Shrubland")]
global_sample_index <- pred[, {
  n <- min(.N, 10000L)
  .SD[sample.int(.N, n, replace=FALSE)]
}, by=EcoType]
rm(pred); invisible(gc())

matmap <- rast(cache_file)
names(matmap) <- c("MAT","MAP")
years <- sort(unique(global_sample_index$Year))
global_parts <- vector("list", length(years))
for (g in seq_along(years)) {
  yr <- years[g]
  rows <- global_sample_index[Year == yr]
  f <- file.path(era_dir, sprintf("era5_land_monthly_%d.nc", yr))
  climate <- rast(f); nm <- names(climate)
  t2 <- grep("t2m", nm, ignore.case=TRUE)
  sw1 <- grep("swvl1", nm, ignore.case=TRUE)
  sw2 <- grep("swvl2", nm, ignore.case=TRUE)
  if (length(t2)!=12L || length(sw1)!=12L || length(sw2)!=12L)
    stop("Unexpected ERA5-Land layer layout: ", basename(f))
  tavg <- climate[[t2]] - 273.15
  vswc <- (climate[[sw1]] + climate[[sw2]])/2
  names(tavg) <- sprintf("tavg_%02d",1:12)
  names(vswc) <- sprintf("vswc_%02d",1:12)
  env <- c(tavg,vswc,matmap)
  names(env)[25:26] <- c("MAT","MAP")
  xy <- as.matrix(rows[, .(x,y)])
  v <- as.data.table(terra::extract(env, xy, method="simple"))
  if ("ID" %in% names(v)) v[, ID := NULL]
  rows[, mean_ERA5_Tair_C := vapply(seq_len(.N),function(i) v[[sprintf("tavg_%02d",month[i])]][i],numeric(1))]
  rows[, mean_ERA5_VSWC := vapply(seq_len(.N),function(i) v[[sprintf("vswc_%02d",month[i])]][i],numeric(1))]
  z <- rows[, .(Year,month,EcoType,area_mha,mean_ERA5_Tair_C,mean_ERA5_VSWC)]
  # Extracted climatologies are in the same order as the production raster stack.
  vmap <- v[["MAP"]]; vmat <- v[["MAT"]]
  z[, MAP := vmap]; z[, MAT := vmat]
  # Match Stage 1's classification predictor; arid rows were excluded above.
  z[, aridity_index := fifelse(MAT > -9, MAP/(MAT+10), NA_real_)]
  global_parts[[g]] <- z[is.finite(aridity_index)]
  rm(climate,tavg,vswc,env,xy,v,z); invisible(gc())
  message(sprintf("Extracted global predictor sample %d/%d (%d)",g,length(years),yr))
}
global <- rbindlist(global_parts)
rm(global_parts, global_sample_index); invisible(gc())
fwrite(global,file.path(out,paste0("source_probability_global_predictor_sample_",support_tag,".csv")))
fwrite(neon,file.path(out,"source_probability_neon_predictors_for_support.csv"))

# Univariate overlap: report NEON 1-99% range, NEON min/max, and area-weighted
# global fractions outside those ranges for each ecosystem and predictor.
support_rows <- list(); k <- 1L
for (eco in sort(unique(neon$EcoType))) {
  nt <- neon[EcoType == eco]; gt <- global[EcoType == eco]
  for (v in numeric_predictors) {
    q <- quantile(nt[[v]], c(.01,.25,.5,.75,.99), na.rm=TRUE, names=FALSE)
    mn <- min(nt[[v]],na.rm=TRUE); mx <- max(nt[[v]],na.rm=TRUE)
    good <- is.finite(gt[[v]]) & is.finite(gt$area_mha)
    w <- gt$area_mha[good]; x <- gt[[v]][good]
    gq <- wquantile(x,w,c(.01,.99))
    support_rows[[k]] <- data.table(EcoType=eco,predictor=v,n_neon=nrow(nt),
      neon_min=mn,neon_q01=q[1],neon_q25=q[2],neon_median=q[3],neon_q75=q[4],neon_q99=q[5],neon_max=mx,
      global_area_weighted_q01=gq[1],global_area_weighted_q99=gq[2],
      global_area_fraction_outside_neon_minmax=sum(w[x<mn|x>mx])/sum(w),
      global_area_fraction_outside_neon_q01_q99=sum(w[x<q[1]|x>q[5]])/sum(w),
      neon_fraction_outside_global_q01_q99=mean(nt[[v]]<gq[1] | nt[[v]]>gq[2],na.rm=TRUE))
    k <- k+1L
  }
}
support <- rbindlist(support_rows)

# Joint support in five numeric dimensions. Scale by NEON IQR within ecosystem;
# set near-zero IQR to 1 to avoid unstable scaling. The cutoff is the 95th
# percentile distance of held-out-site NEON records to the remaining sites.
joint_rows <- list(); k <- 1L
for (eco in sort(unique(neon$EcoType))) {
  nt <- neon[EcoType == eco]
  gt <- global[EcoType == eco]
  if (uniqueN(nt$SITE_ID) < 2 || nrow(gt) == 0) next
  centers <- vapply(numeric_predictors, function(v) median(nt[[v]],na.rm=TRUE), numeric(1))
  scales <- vapply(numeric_predictors, function(v) IQR(nt[[v]],na.rm=TRUE), numeric(1))
  scales[!is.finite(scales) | scales < 1e-8] <- 1
  nmat <- as.matrix(nt[, ..numeric_predictors]); storage.mode(nmat) <- "double"
  gmat <- as.matrix(gt[, ..numeric_predictors]); storage.mode(gmat) <- "double"
  nscaled <- sweep(sweep(nmat,2,centers,"-"),2,scales,"/")
  gscaled <- sweep(sweep(gmat,2,centers,"-"),2,scales,"/")
  validn <- complete.cases(nscaled); validg <- complete.cases(gscaled)
  nt2 <- nt[validn]; nscaled <- nscaled[validn,,drop=FALSE]
  gt2 <- gt[validg]; gscaled <- gscaled[validg,,drop=FALSE]
  heldout_dist <- numeric()
  for (site in unique(nt2$SITE_ID)) {
    ix <- nt2$SITE_ID == site; tr <- nscaled[!ix,,drop=FALSE]
    if (!any(ix) || nrow(tr)==0) next
    heldout_dist <- c(heldout_dist, FNN::get.knnx(tr,nscaled[ix,,drop=FALSE],k=1)$nn.dist[,1])
  }
  cutoff <- as.numeric(quantile(heldout_dist,.95,na.rm=TRUE))
  gd <- FNN::get.knnx(nscaled,gscaled,k=1)$nn.dist[,1]
  joint_rows[[k]] <- data.table(EcoType=eco,n_sites=uniqueN(nt2$SITE_ID),n_neon=nrow(nt2),n_global_sample=nrow(gt2),
    heldout_neon_distance_q95=cutoff,
    global_area_fraction_beyond_loso95=sum(gt2$area_mha[gd>cutoff])/sum(gt2$area_mha),
    global_area_fraction_within_loso95=sum(gt2$area_mha[gd<=cutoff])/sum(gt2$area_mha),
    global_nearest_distance_median=median(gd),global_nearest_distance_q95=quantile(gd,.95))
  k <- k+1L
}
joint <- rbindlist(joint_rows)

writeLines(c("# Global vs NEON source-probability predictor support", "",
  sprintf("Global predictor samples are drawn uniformly within ecosystem from years %s and months %s (area weights are used for global fractions) and use the same ERA5-Land monthly and MAT/MAP climatology inputs as the spatial projection.", paste(support_years,collapse=","), paste(support_months,collapse=",")),
  "NEON prevalence is not used as a global prior. Univariate support reports global area outside the observed NEON range and central 98% range.",
  "Joint support uses five numeric predictors scaled by within-ecosystem NEON IQR. The support cutoff is the 95th percentile of NEON leave-one-site-out nearest-neighbor distances; this is a diagnostic threshold, not a universal extrapolation boundary.",
  "Global-area fractions summarize the sampled cells and are approximate Monte Carlo estimates.","",
  capture.output(print(support)),"",capture.output(print(joint))),
  file.path(out,sprintf("source_probability_global_neon_support_%s.md",support_tag)))
fwrite(support,file.path(out,sprintf("source_probability_global_neon_univariate_support_%s.csv",support_tag)))
fwrite(joint,file.path(out,sprintf("source_probability_global_neon_joint_support_%s.csv",support_tag)))
print(support)
print(joint)
