library(data.table)
library(dplyr)
library(stringr)
library(ranger)
library(ggplot2)

root <- Sys.getenv("LOCALDIR_CH4", unset="/Volumes/MaloneLab/Research/FluxGradient/METHANE")
rfdir <- Sys.getenv("MONTHLY_RF_DIR", unset="/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF")
out <- file.path(rfdir,"OUTPUT"); figdir <- file.path(rfdir,"FIGURES")
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
  mean_ERA5_Tair_C=mean(ERA5_Tair_C,na.rm=TRUE),
  mean_ERA5_VSWC=mean(ERA5_VSWC,na.rm=TRUE)) |>
  left_join(attrs,by=c("SITE_ID","EcoType")) |>
  mutate(y=as.integer(budget>0),aridity=if_else(MAT>-9,MAP/(MAT+10),NA_real_),
    is_arid=as.integer(is.finite(aridity)&aridity<15),
    EcoType=factor(if_else(is_arid==1,"Arid",as.character(EcoType)),
      levels=c("Forest","Grassland","Shrubland","Arid"))) |>
  filter(is.finite(aridity),is.finite(mean_ERA5_Tair_C),is.finite(mean_ERA5_VSWC),is.finite(MAP),is.finite(MAT),is_arid==0)
predictors <- c("EcoType","mean_ERA5_Tair_C","mean_ERA5_VSWC","MAP","MAT","aridity")
sites <- sort(unique(m$SITE_ID))
iso_fit <- function(p,y) { o<-order(p); list(x=p[o],y=isoreg(p[o],as.numeric(y[o]))$yf) }
iso_apply <- function(f,p) pmin(pmax(approx(f$x,f$y,xout=p,method="linear",rule=2,ties="ordered")$y,0),1)
logit <- function(p) qlogis(pmin(pmax(p,1e-5),1-1e-5))
preds <- vector("list",length(sites))
for (i in seq_along(sites)) {
  s<-sites[i]; tr<-m[m$SITE_ID!=s,]; te<-m[m$SITE_ID==s,]
  neg<-which(tr$y==0); pos<-which(tr$y==1)
  if (min(length(neg),length(pos))<2) next
  set.seed(20261021+i); bal<-tr[c(neg,sample(pos,length(neg))),]
  bal$yf<-factor(bal$y,levels=c(0,1))
  fit<-ranger(reformulate(predictors,"yf"),data=bal,num.trees=500,probability=TRUE,
    min.node.size=20,max.depth=8,replace=FALSE,sample.fraction=.7,
    importance="none",num.threads=4,seed=20261021+i)
  iso<-iso_fit(fit$predictions[,"1"],bal$y)
  p<-iso_apply(iso,predict(fit,data=te)$predictions[,"1"])
  pi<-mean(tr$y)
  pprior<-plogis(logit(p)+qlogis(pi)-qlogis(.5))
  preds[[i]]<-data.frame(SITE_ID=s,EcoType=as.character(te$EcoType),Year=te$Year,month=te$month,
    observed_source=te$y,p_balanced=p,p_neon_prior=pprior,training_source_fraction=pi)
  message(sprintf("LOSO %d/%d %s; %d months",i,length(sites),s,nrow(te)))
  rm(fit,bal); invisible(gc(verbose=FALSE))
}
pdat<-bind_rows(preds)
auc<-function(y,p){n1<-sum(y==1);n0<-sum(y==0);if(!n1||!n0)return(NA_real_);r<-rank(p,ties.method="average");(sum(r[y==1])-n1*(n1+1)/2)/(n1*n0)}
score<-function(dat,col,label){
  y<-dat$observed_source;p<-dat[[col]];lp<-logit(p)
  slope<-tryCatch(coef(glm(y~lp,family=binomial()))[2],error=function(e)NA_real_)
  intercept<-tryCatch(coef(glm(y~offset(lp),family=binomial()))[1],error=function(e)NA_real_)
  data.frame(model=label,n=nrow(dat),sites=uniqueN(dat$SITE_ID),observed_source_fraction=mean(y),
    mean_probability=mean(p),calibration_intercept=intercept,calibration_slope=slope,
    auc=auc(y,p),brier=mean((p-y)^2),predicted_source_fraction_at_0_5=mean(p>=.5))
}
metrics<-bind_rows(score(pdat,"p_balanced","Balanced-prior isotonic"),
  score(pdat,"p_neon_prior","NEON-prevalence sensitivity"))
by_site<-pdat |> group_by(SITE_ID,EcoType) |> summarise(n=n(),observed=mean(observed_source),
  predicted_balanced=mean(p_balanced),predicted_neon_prior=mean(p_neon_prior),.groups="drop")
pdat$bin<-dplyr::ntile(pdat$p_balanced,10)
cal<-pdat |> group_by(bin) |> summarise(n=n(),predicted_balanced=mean(p_balanced),
  predicted_neon_prior=mean(p_neon_prior),observed=mean(observed_source),.groups="drop")
write.csv(pdat,file.path(out,"source_probability_loso_predictions.csv"),row.names=FALSE)
write.csv(metrics,file.path(out,"source_probability_loso_metrics.csv"),row.names=FALSE)
write.csv(by_site,file.path(out,"source_probability_loso_predictions_by_site.csv"),row.names=FALSE)
write.csv(cal,file.path(out,"source_probability_loso_calibration.csv"),row.names=FALSE)
plotdat<-bind_rows(
  cal |> transmute(x=predicted_balanced,y=observed,model="Balanced-prior isotonic"),
  cal |> transmute(x=predicted_neon_prior,y=observed,model="NEON-prior sensitivity"))
g<-ggplot(plotdat,aes(x,y,color=model))+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey50")+
  geom_line()+geom_point(size=2)+coord_equal(xlim=c(0,1),ylim=c(0,1))+theme_bw(base_size=12)+
  labs(title="Source probability on held-out NEON sites",
    subtitle="RF and isotonic calibration refit without each held-out site",
    x="Mean predicted P(source)",y="Observed source-month fraction",color=NULL,
    caption="Prior correction uses NEON prevalence as a sensitivity, not a global prior.")+
  theme(legend.position="bottom")
ggsave(file.path(figdir,"FIGURE_source_probability_loso_calibration.png"),g,width=8,height=5.2,dpi=250,bg="white")
writeLines(c("# Source-probability leave-one-site-out audit","",
  "Each site's full monthly labels were held out from both RF fitting and isotonic calibration. The NEON-prior sensitivity shifts balanced-prior odds to the training-fold source frequency; it is not a global prior estimate because NEON sites are not a probability sample of global land area.","",
  capture.output(print(metrics)),
  "","Monthly source labels inherit limitations of the existing NEON gapfilled flux. This tests transfer to unseen NEON sites, not global prevalence."),
  file.path(out,"source_probability_loso_audit.md"))
print(metrics)
cat("\nEqual-site means\n")
print(by_site |> summarise(n_sites=n(),observed=mean(observed),balanced=mean(predicted_balanced),neon_prior=mean(predicted_neon_prior)))

