#!/usr/bin/env Rscript
# Nested leave-one-site-out threshold selection for monthly source classification.
# Outer site labels are never used to select the threshold or fit its model.
# Candidate cutoffs optimize equal-site macro balanced accuracy on held-out months.
# Current inputs imply a forest cutoff near 0.85; grassland has no held-out sink
# months and shrubland has only three sites, so neither supports threshold choice.
# Rerun after recalculating NEON gapfilled fluxes and regenerating the site summary.
# This audit does not change the production model bundle or its default cutoff.

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(ranger)
})
root <- Sys.getenv("LOCALDIR_CH4", unset="/Volumes/MaloneLab/Research/FluxGradient/METHANE")
rfdir <- Sys.getenv("MONTHLY_RF_DIR", unset=file.path(root,"Upscaling_Monthly_RF"))
out <- file.path(rfdir,"OUTPUT")
attrs <- fread(file.path(root,"OUTPUT/NEON_scale_annual_budget_summary.csv")) |>
  mutate(SITE_ID=as.character(SITE_ID)) |>
  filter(!is.na(EcoType), !str_detect(EcoType,regex("wetland|inundat|flood|marsh|swamp|bog|fen|lake|rice|crop|agri",ignore_case=TRUE))) |>
  distinct(SITE_ID,EcoType,MAP,MAT)
x <- fread(file.path(root,"OUTPUT/NEON_ERA5_gapfilled_30min.csv.gz"),
  select=c("SITE_ID","EcoType","Year","month","ERA5_Tair_C","ERA5_VSWC","gapfilled_CH4_mgC_30min")) |>
  mutate(SITE_ID=as.character(SITE_ID),Year=as.integer(Year),month=as.integer(month),
    across(c(ERA5_Tair_C,ERA5_VSWC,gapfilled_CH4_mgC_30min),as.numeric)) |>
  inner_join(attrs,by=c("SITE_ID","EcoType")) |>
  filter(is.finite(Year),is.finite(month),is.finite(ERA5_Tair_C),is.finite(ERA5_VSWC))
m <- x |> reframe(.by=c(SITE_ID,EcoType,Year,month),
  budget=sum(gapfilled_CH4_mgC_30min,na.rm=TRUE),
  mean_ERA5_Tair_C=mean(ERA5_Tair_C,na.rm=TRUE), mean_ERA5_VSWC=mean(ERA5_VSWC,na.rm=TRUE)) |>
  left_join(attrs,by=c("SITE_ID","EcoType")) |>
  mutate(y=as.integer(budget>0),aridity=if_else(MAT>-9,MAP/(MAT+10),NA_real_),
    is_arid=as.integer(is.finite(aridity)&aridity<15),
    EcoType=as.character(EcoType)) |>
  filter(is.finite(aridity),is.finite(mean_ERA5_Tair_C),is.finite(mean_ERA5_VSWC),
         is.finite(MAP),is.finite(MAT),is_arid==0)
predictors <- c("EcoType","mean_ERA5_Tair_C","mean_ERA5_VSWC","MAP","MAT","aridity")
thresholds <- seq(.1,.9,.01)
fit_predict <- function(tr,te,seed) {
  neg <- which(tr$y==0); pos <- which(tr$y==1)
  if (min(length(neg),length(pos))<2) return(NULL)
  set.seed(seed); bal <- tr[c(neg,sample(pos,length(neg))),]
  bal$yf <- factor(bal$y,levels=c(0,1))
  fit <- ranger(reformulate(predictors,"yf"),data=bal,num.trees=500,probability=TRUE,
    min.node.size=20,max.depth=8,replace=FALSE,sample.fraction=.7,
    importance="none",num.threads=4,seed=seed)
  raw <- fit$predictions[,"1"]
  ord <- order(raw); iso <- isoreg(raw[ord],as.numeric(bal$y[ord]))
  f <- approxfun(raw[ord],iso$yf,rule=2,ties="ordered")
  p <- pmin(pmax(f(predict(fit,data=te)$predictions[,"1"]),0),1)
  pi <- mean(tr$y); p <- plogis(qlogis(pmin(pmax(p,1e-5),1-1e-5))+qlogis(pi)-qlogis(.5))
  rm(fit,bal); invisible(gc(FALSE)); p
}
macro_ba <- function(d,t) {
  z <- d[,.(sens=if(any(y==1)) mean(p[y==1]>=t) else NA_real_,
            spec=if(any(y==0)) mean(p[y==0]<t) else NA_real_),by=SITE_ID]
  sens <- mean(z$sens,na.rm=TRUE); spec <- mean(z$spec,na.rm=TRUE)
  if(!is.finite(sens)||!is.finite(spec)) return(NA_real_)
  (sens+spec)/2
}
outer_rows <- list(); chosen <- list(); ri <- 1L; ci <- 1L
for (eco in c("Forest","Grassland","Shrubland")) {
  dat <- as.data.table(m)[EcoType==eco]
  sites <- sort(unique(dat$SITE_ID)); min_sites <- if (eco=="Shrubland") 5 else 3
  if(length(sites)<min_sites) { message("Skip threshold selection for ",eco,": only ",length(sites)," sites"); next }
  for (oi in seq_along(sites)) {
    outer <- sites[oi]; training <- dat[SITE_ID!=outer]; test <- dat[SITE_ID==outer]
    inner_preds <- list(); ii <- 1L
    inner_sites <- sort(unique(training$SITE_ID))
    for (inner in inner_sites) {
      itrain <- training[SITE_ID!=inner]; ival <- training[SITE_ID==inner]
      p <- fit_predict(itrain,ival,20261021L + oi*1000L + match(inner,inner_sites))
      if(!is.null(p)) { inner_preds[[ii]] <- data.table(SITE_ID=inner,Year=ival$Year,month=ival$month,
        y=ival$y,budget=ival$budget,p=p); ii <- ii+1L }
    }
    inner <- rbindlist(inner_preds,fill=TRUE)
    if(!nrow(inner)) next
    scores <- vapply(thresholds,function(t) macro_ba(inner,t),numeric(1))
    best <- if(any(is.finite(scores))) max(scores,na.rm=TRUE) else NA_real_
    if(is.finite(best)) {
      candidates <- thresholds[scores>=best-0.002]
      tstar <- candidates[which.min(abs(candidates-.5))]
    } else tstar <- .5 # cutoff is not identifiable without both source and sink years
    ptest <- fit_predict(training,test,20262021L+oi)
    if(is.null(ptest)) next
    outer_rows[[ri]] <- data.table(EcoType=eco,outer_site=outer,Year=test$Year,month=test$month,
      observed=test$y,budget=test$budget,p=ptest,selected_threshold=tstar); ri<-ri+1L
    chosen[[ci]] <- data.table(EcoType=eco,outer_site=outer,selected_threshold=tstar,
      inner_macro_balanced_accuracy=best,n_inner_sites=uniqueN(inner$SITE_ID),
      inner_source_months=sum(inner$y==1),inner_sink_months=sum(inner$y==0)); ci<-ci+1L
    message(sprintf("%s outer %d/%d %s: selected %.2f",eco,oi,length(sites),outer,tstar))
  }
}
pred <- rbindlist(outer_rows,fill=TRUE); th <- rbindlist(chosen,fill=TRUE)
fwrite(pred,file.path(out,"source_probability_nested_threshold_loso_predictions.csv"))
fwrite(th,file.path(out,"source_probability_nested_threshold_loso_selected.csv"))
calc_metric <- function(z,label,fixed=NULL) {
  eco <- unique(z$EcoType)
  cut <- if(is.null(fixed)) z$selected_threshold else rep(fixed,nrow(z))
  pr <- z$p>=cut; y <- z$observed
  sens <- if(any(y==1)) mean(pr[y==1]) else NA_real_
  spec <- if(any(y==0)) mean(!pr[y==0]) else NA_real_
  data.table(EcoType=eco,cutoff=label,n_months=nrow(z),n_sites=uniqueN(z$outer_site),
    observed_source_month_fraction=mean(y),predicted_source_month_fraction=mean(pr),
    sensitivity=sens,specificity=spec,
    balanced_accuracy=if(is.finite(sens)&&is.finite(spec))(sens+spec)/2 else NA_real_,
    accuracy=mean(pr==y),brier=mean((z$p-y)^2))
}
res <- rbindlist(lapply(unique(pred$EcoType),function(eco) {
  z <- pred[EcoType==eco]
  rbindlist(list(calc_metric(copy(z),"0.5",.5),calc_metric(copy(z),"0.8",.8),
                 calc_metric(copy(z),"nested selected")))
}))
fwrite(res,file.path(out,"source_probability_nested_threshold_loso_metrics.csv"))
print(res); print(th[,.(median=median(selected_threshold),q25=quantile(selected_threshold,.25),
  q75=quantile(selected_threshold,.75),min=min(selected_threshold),max=max(selected_threshold),
  source_months=sum(inner_source_months),sink_months=sum(inner_sink_months)),by=EcoType])
