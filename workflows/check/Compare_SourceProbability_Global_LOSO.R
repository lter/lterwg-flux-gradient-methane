library(data.table)
library(ggplot2)
rfdir <- Sys.getenv("MONTHLY_RF_DIR",unset="/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly_RF")
out <- file.path(rfdir,"OUTPUT"); figdir <- file.path(rfdir,"FIGURES")
spatial_file <- file.path(out,"monthly_cell_predictions_2000_2025.rds")
loso_file <- file.path(out,"source_probability_loso_predictions.csv")
p <- fread(loso_file)
siteyear <- p[,.(n_months=.N,observed_month_source_fraction=mean(observed_source),
  mean_predicted_P=mean(p_balanced),mean_predicted_P_neon_prior=mean(p_neon_prior)),
  by=.(SITE_ID,EcoType,Year)]
siteyear <- siteyear[n_months>=10]
val <- siteyear[,.(n_siteyears=.N,observed_source_year_fraction=mean(observed_month_source_fraction>=.5),
  loso_predicted_source_year_fraction=mean(mean_predicted_P>=.5),
  mean_loso_P=mean(mean_predicted_P),
  mean_observed_month_source_fraction=mean(observed_month_source_fraction)),by=EcoType]
message("Reading 1 GB spatial prediction file for annual area summaries...")
x <- readRDS(spatial_file); setDT(x)
annual_cells <- x[,.(area_mha=first(area_mha),mean_P=mean(source_prob_A,na.rm=TRUE)),
  by=.(Year,cell,EcoType,is_arid)]
area <- annual_cells[,.(annual_area_mha=sum(area_mha),
  source_area_share=sum(area_mha[mean_P>=.5])/sum(area_mha),
  area_weighted_mean_P=sum(area_mha*mean_P)/sum(area_mha)),by=.(Year,EcoType,is_arid)]
global <- area[,.(mean_annual_area_mha=mean(annual_area_mha),
  mean_annual_source_area_share=mean(source_area_share),
  mean_area_weighted_P=mean(area_weighted_mean_P)),by=.(EcoType,is_arid)]
global <- global[is_arid==0]
comparison <- merge(global,val,by="EcoType",all.x=TRUE)
write.csv(comparison,file.path(out,"source_probability_global_vs_loso_by_ecotype.csv"),row.names=FALSE)
print(comparison)
plotdata <- rbindlist(list(
  comparison[,.(EcoType,estimate=mean_annual_source_area_share,label="Global area classified source")],
  comparison[,.(EcoType,estimate=loso_predicted_source_year_fraction,label="LOSO predicted source site-years")],
  comparison[,.(EcoType,estimate=observed_source_year_fraction,label="LOSO observed source site-years")]),use.names=TRUE)
plotdata[,label:=factor(label,levels=c("Global area classified source","LOSO predicted source site-years","LOSO observed source site-years"))]
g<-ggplot(plotdata,aes(EcoType,estimate,color=label,group=label))+geom_point(size=3,position=position_dodge(width=.35))+
  geom_line(position=position_dodge(width=.35),linewidth=.5)+scale_y_continuous(limits=c(0,1),labels=scales::percent)+
  theme_bw(base_size=12)+labs(title="Global source-area classification vs held-out NEON sites",
    subtitle="A site-year is called source when mean monthly source probability/frequency is at least 0.5",
    x=NULL,y="Fraction classified source",color=NULL,
    caption="LOSO statistics use complete or near-complete site-years (≥10 months); shrubland has only 3 sites.")+
  theme(legend.position="bottom")
ggsave(file.path(figdir,"FIGURE_source_area_vs_loso_sites.png"),g,width=9,height=5.2,dpi=250,bg="white")

