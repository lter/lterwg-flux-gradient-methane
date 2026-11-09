# Importance-weighted source-probability evaluation, restricted to predictor
# conditions with NEON leave-one-site-out support. This estimates transfer to
# the supported global covariate distribution only; it does not set global
# source prevalence or correct unsupported conditions.

library(data.table)
library(FNN)
library(mgcv)

root <- Sys.getenv("LOCALDIR_CH4",unset="/Volumes/MaloneLab/Research/FluxGradient/METHANE")
rfdir <- Sys.getenv("MONTHLY_RF_DIR",unset="/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF")
out <- file.path(rfdir,"OUTPUT")
tag <- "years_2000-2010-2020-2025_months_1-2-3-4-5-6-7-8-9-10-11-12"
num <- c("mean_ERA5_Tair_C","mean_ERA5_VSWC","MAP","MAT","aridity_index")
keys <- c("SITE_ID","EcoType","Year","month")

neon <- fread(file.path(out,"source_probability_neon_predictors_for_support.csv"))
global <- fread(file.path(out,paste0("source_probability_global_predictor_sample_",tag,".csv")))
pred <- fread(file.path(out,"source_probability_loso_predictions.csv"))
support <- fread(file.path(out,paste0("source_probability_global_neon_joint_support_",tag,".csv")))
neon[, SITE_ID:=as.character(SITE_ID)]; pred[, SITE_ID:=as.character(SITE_ID)]
neon[, EcoType:=as.character(EcoType)]; pred[, EcoType:=as.character(EcoType)]; global[,EcoType:=as.character(EcoType)]
setkeyv(neon,keys); setkeyv(pred,keys)
val <- merge(neon,pred[,.(SITE_ID,EcoType,Year,month,p_balanced)],by=keys,all=FALSE)
setorder(val,SITE_ID,Year,month)

metrics <- list(); weights_out <- list(); mi <- 1L; wi <- 1L
ecosystems <- strsplit(Sys.getenv("ECO_TYPES",unset="Forest,Grassland"),",",fixed=TRUE)[[1]]
site_limit_value <- Sys.getenv("SITE_LIMIT",unset="")
site_limit <- if (nzchar(site_limit_value)) as.integer(site_limit_value) else Inf
for (eco in ecosystems) {
  nall <- neon[EcoType==eco]
  gall <- global[EcoType==eco]
  if (!nrow(gall) || uniqueN(nall$SITE_ID)<3) next
  center <- vapply(num,function(v) median(nall[[v]],na.rm=TRUE),numeric(1))
  scale <- vapply(num,function(v) IQR(nall[[v]],na.rm=TRUE),numeric(1))
  scale[!is.finite(scale)|scale<1e-8] <- 1
  nmat <- sweep(sweep(as.matrix(nall[,..num]),2,center,"-"),2,scale,"/")
  gmat <- sweep(sweep(as.matrix(gall[,..num]),2,center,"-"),2,scale,"/")
  nvalid <- complete.cases(nmat); gvalid <- complete.cases(gmat)
  nall <- nall[nvalid]; nmat <- nmat[nvalid,,drop=FALSE]
  gall <- gall[gvalid]; gmat <- gmat[gvalid,,drop=FALSE]
  cutoff <- support[EcoType==eco,heldout_neon_distance_q95]
  test_all <- val[EcoType==eco]
  site_metrics <- list(); si <- 1L
  sites_to_eval <- sort(unique(test_all$SITE_ID))
  if (is.finite(site_limit)) sites_to_eval <- head(sites_to_eval,site_limit)
  for (site in sites_to_eval) {
    te <- test_all[SITE_ID==site]
    tr <- nall[SITE_ID!=site]
    trmat <- nmat[nall$SITE_ID!=site,,drop=FALSE]
    # Restrict the target distribution using only the training sites in this
    # fold, so the held-out site's covariates cannot define its own support.
    gdist_fold <- FNN::get.knnx(trmat,gmat,k=1)$nn.dist[,1]
    gkeep <- gdist_fold<=cutoff
    gtarget <- gall[gkeep]
    gscaled <- gmat[gkeep,,drop=FALSE]
    if (nrow(gtarget)<100) next
    tem <- sweep(sweep(as.matrix(te[,..num]),2,center,"-"),2,scale,"/")
    ok <- complete.cases(tem)
    te <- te[ok]; tem <- tem[ok,,drop=FALSE]
    if (!nrow(te) || !nrow(tr) || uniqueN(te$EcoType)!=1) next
    # Determine whether each held-out NEON predictor vector lies inside the
    # pooled leave-one-site-out NEON support envelope.
    dtest <- FNN::get.knnx(trmat,tem,k=1)$nn.dist[,1]
    keep <- dtest<=cutoff
    te <- te[keep]
    tem <- tem[keep,,drop=FALSE]
    if (!nrow(te)) next

    # Domain classifier estimates global/NEON covariate density ratio. It is
    # refit without the held-out site. Domain mass is balanced; global rows
    # use area weights and NEON training rows give equal total weight per site.
    ntr <- tr
    ntr[, row_weight := 0.5 / uniqueN(ntr$SITE_ID) / .N, by=SITE_ID]
    gwt <- 0.5 * gtarget$area_mha / sum(gtarget$area_mha)
    gd <- as.data.table(gscaled); setnames(gd,names(gd),make.names(num))
    gd[,domain:=1L]; gd[,weight:=gwt]
    nd <- as.data.table(sweep(sweep(as.matrix(ntr[,..num]),2,center,"-"),2,scale,"/"))
    setnames(nd,names(nd),make.names(num))
    nd[,domain:=0L]; nd[,weight:=ntr$row_weight]
    dtrain <- rbindlist(list(gd,nd),use.names=TRUE)
    # Rebuild a plain model frame in a fixed predictor order.
    form <- as.formula(paste("domain ~",paste(sprintf("s(%s,k=5)",make.names(num)),collapse=" + ")))
    fit <- tryCatch(mgcv::gam(form,data=dtrain,weights=weight,family=quasibinomial(),method="REML"),error=function(e) NULL)
    if (is.null(fit)) next
    xte <- as.data.table(tem); setnames(xte,names(xte),make.names(num))
    pg <- predict(fit,newdata=xte,type="response")
    raww <- pmin(pg/pmax(1-pg,1e-6),10)
    w <- raww/mean(raww)
    metrics[[mi]] <- data.table(EcoType=eco,heldout_site=site,n=nrow(te),
      n_supported_sites=uniqueN(test_all$SITE_ID),support_fraction=sum(keep)/nrow(test_all[SITE_ID==site]),
      weight_cap_fraction=mean(pg/pmax(1-pg,1e-6)>10),weight_ess=sum(w)^2/sum(w^2),
      unweighted_source_fraction=mean(te$observed_source),
      unweighted_mean_predicted_probability=mean(te$p_balanced),
      unweighted_brier=mean((te$p_balanced-te$observed_source)^2),
      unweighted_source_at_0_5=mean(te$p_balanced>=.5),
      unweighted_source_at_0_8=mean(te$p_balanced>=.8),
      observed_source_fraction=weighted.mean(te$observed_source,w),
      mean_predicted_probability=weighted.mean(te$p_balanced,w),
      brier=weighted.mean((te$p_balanced-te$observed_source)^2,w),
      predicted_source_fraction_at_0_5=weighted.mean(te$p_balanced>=.5,w),
      predicted_source_fraction_at_0_8=weighted.mean(te$p_balanced>=.8,w))
    weights_out[[wi]] <- data.table(SITE_ID=te$SITE_ID,EcoType=eco,Year=te$Year,month=te$month,
      observed_source=te$observed_source,p_balanced=te$p_balanced,nn_distance=dtest[keep],
      density_ratio_raw=pg/pmax(1-pg,1e-6),weight_capped=raww,weight_normalized=w)
    mi<-mi+1L;wi<-wi+1L
    message(sprintf("Weighted support evaluation %s, held-out site %s",eco,site))
  }
}

fold <- rbindlist(metrics,fill=TRUE)
setDT(fold)
overall <- fold[,.(n_site_folds=.N,n_month_rows=sum(n),n_sites=uniqueN(heldout_site),
  mean_site_unweighted_observed=mean(unweighted_source_fraction),
  mean_site_unweighted_prediction=mean(unweighted_mean_predicted_probability),
  mean_site_unweighted_brier=mean(unweighted_brier),
  mean_site_unweighted_source_at_0_5=mean(unweighted_source_at_0_5),
  mean_site_unweighted_source_at_0_8=mean(unweighted_source_at_0_8),
  mean_site_weighted_observed=mean(observed_source_fraction),
  mean_site_weighted_prediction=mean(mean_predicted_probability),
  mean_site_weighted_brier=mean(brier),
  mean_site_weighted_source_at_0_5=mean(predicted_source_fraction_at_0_5),
  mean_site_weighted_source_at_0_8=mean(predicted_source_fraction_at_0_8),
  mean_support_fraction=mean(support_fraction),mean_weight_ess=mean(weight_ess),
  mean_weight_cap_fraction=mean(weight_cap_fraction)),by=EcoType]
fwrite(fold,file.path(out,"source_probability_supported_weighted_loso_by_site.csv"))
fwrite(overall,file.path(out,"source_probability_supported_weighted_loso_summary.csv"))
fwrite(rbindlist(weights_out,fill=TRUE),file.path(out,"source_probability_supported_weighted_loso_weights.csv"))
print(overall)
