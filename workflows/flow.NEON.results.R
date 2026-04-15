# Figures for parms:
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)
library(sf)
library(AOI)

# Import data and build the files needed: ####
localdir <- '/Volumes/MaloneLab/Research/FluxGradient/FluxData'

# Load parameter information:
load(file = fs::path(localdir,paste0("NEON_PARMS_DIEL_Q10.Rdata")))
# Load soil and root data
load(file=paste(localdir, "/", 'Soildata_YearMon.Rdata', sep=""))
# Load data obtained from Avni:
avni <- read.csv( paste("/Volumes/MaloneLab/Research/FluxGradient/Avni's Data/neondata_forSM.csv" , sep="" )) %>% rename( SITE_ID =site_code )

avni2 <- read.csv( paste("/Volumes/MaloneLab/Research/FluxGradient/Avni's Data/neondata_forsteve.csv" , sep="" )) %>% rename( SITE_ID =site_code ) %>% dplyr::select(
  SITE_ID, MAP, MAT, aridity.index,shallow_root, deep_root)

# Load site information:
soil.data.biogeo <- read.csv(file=paste(localdir, "/", 'Soil_Biogeochem_RootBiomass.csv', sep=""))
# Load canopy information: 
canopy <- read.csv(file.path(paste(localdir, "canopy_commbined.csv", sep="/"))) %>% distinct

canopy.Ht <- canopy %>% rename(SITE_ID = Site) %>% 
  reframe( .by= SITE_ID,
           canopyHeight_m = mean(canopyHeight_m))

metadata <- read.csv('/Volumes/MaloneLab/Research/FluxGradient/Ameriflux_NEON field-sites.csv') %>% 
  mutate(EcoType = case_when( Vegetation.Abbreviation..IGBP. == 'ENF' ~ 'Forest',
                              Vegetation.Abbreviation..IGBP. == 'DBF' ~ 'Forest', 
                              Vegetation.Abbreviation..IGBP. == 'MF' ~ 'Forest',
                              Vegetation.Abbreviation..IGBP. == 'EBF' ~ 'Forest',
                              Vegetation.Abbreviation..IGBP. == 'SAV' ~ 'Forest',
                              Vegetation.Abbreviation..IGBP. == 'WET' ~ 'Wetland',
                              
                              Vegetation.Abbreviation..IGBP. == 'GRA' ~ 'Grassland',
                              
                              Vegetation.Abbreviation..IGBP. == 'CVM' ~ 'Cropland',
                              Vegetation.Abbreviation..IGBP. == 'CRO' ~ 'Cropland',
                              Vegetation.Abbreviation..IGBP. == 'OSH' ~ 'Shrubland')) %>% 
  rename( site = Site_Id.NEON,
          SITE_ID = Site_Id.NEON) %>% 
  left_join(canopy.Ht , by = 'SITE_ID') %>% 
  left_join(soil.data.biogeo %>% rename(SITE_ID = siteID)  , by = 'SITE_ID')%>% 
  left_join(avni , by = 'SITE_ID')%>% 
  left_join(avni2 , by = 'SITE_ID')

# Create a simple feature
site.att.sf <- st_as_sf(x = metadata,                         
                        coords = c("Longitude..degrees.", "Latitude..degrees."), crs = 4326) 

#format data:
ENSEMBLE_DIELS_canopy <-  ENSEMBLE_DIELS %>% rename( SITE_ID = site) %>% left_join( metadata, by = "SITE_ID")

ENSEMBLE_Q10_eq4_canopy <-   ENSEMBLE_Q10_eq4 %>% 
  mutate( Date = as.Date(paste(idx, "-01", sep="")),
          month = format(Date, "%m") %>% as.numeric,
          Season = case_when(
            month %in% c(12, 1, 2) ~ "Winter",
            month %in% c(3, 4, 5) ~ "Spring",
            month %in% c(6, 7, 8) ~ "Summer",
            TRUE ~ "Autumn")) %>% left_join( metadata, by = "SITE_ID")

# This file has monthly data in it

parms <- ENSEMBLE_Q10_eq4_canopy %>% 
  dplyr::select(Q10.mean, SITE_ID, Date, Rref.mean, EcoType, Season, sulfurTot, dryMass, acidity, FRC.kg.C.m.2, SOC.kg.C.m.2, beta_roots, beta_soc, MAP, MAT, aridity.index,shallow_root, deep_root) %>% 
  left_join( ENSEMBLE_DIELS %>% 
               mutate(Date = as.Date(paste(YearMon, "-01", sep="")),
                      SITE_ID = site) %>% 
               reframe( .by = c(SITE_ID, Date), 
                        DIEL = sum(DIEL)), by=c('SITE_ID', 'Date')) %>% 
  left_join( Site_SoilData %>% 
               mutate(Date = as.Date(paste(YearMon, "-01", sep="")), SITE_ID = Site),by = c('SITE_ID', 'Date')) %>% 
  mutate( Source = case_when( DIEL > 0 ~ 1,
                              DIEL <= 0 ~ 0) %>% as.factor)

parms.site.season <- parms %>% reframe( .by = c(SITE_ID, Season),
                                        VSWCMean  = mean( VSWCMean, na.rm=TRUE), 
                                        Q10.mean = mean(Q10.mean, na.rm=TRUE) ,
                                        Rref.mean = mean(Rref.mean, na.rm=TRUE), 
                                        VSWCMin = mean(VSWCMin, na.rm=TRUE),
                                        VSWCVar= mean(VSWCVar, na.rm=TRUE),
                                        sulfurTot= mean(sulfurTot, na.rm=TRUE),
                                        dryMass= mean( dryMass, na.rm=TRUE),
                                        acidity= mean( acidity, na.rm=TRUE),
                                        DIEL = mean(DIEL, na.rm=T),
                                        FRC.kg.C.m.2 = mean(FRC.kg.C.m.2, na.rm=T), 
                                        SOC.kg.C.m.2 = mean(SOC.kg.C.m.2, na.rm=T),
                                        beta_roots= mean(beta_roots, na.rm=T),  
                                        beta_soc= mean(beta_soc, na.rm=T),
                                        MAT= mean(MAT, na.rm=T), 
                                        MAP= mean(MAP, na.rm=T))  %>% 
  mutate( Source = case_when( DIEL > 0 ~ "Source",
                              DIEL <= 0 ~ "Sink") %>% as.factor)

parms.site <- parms.site.season %>% mutate( Source = case_when( DIEL > 0 ~ 1,
                                                                DIEL <= 0 ~ 0) %>% as.numeric) %>% mutate(across(where(is.numeric), ~ ifelse(is.infinite(.), 0, .))) %>% 
  reframe( .by = SITE_ID, 
           Source = sum(Source ) %>% as.factor,
           VSWCMean  = mean( VSWCMean, na.rm=TRUE), 
           Q10.mean = mean(Q10.mean, na.rm=TRUE) ,
           Rref.mean = mean(Rref.mean, na.rm=TRUE), 
           VSWCMin = mean(VSWCMin, na.rm=TRUE),
           VSWCVar= mean(VSWCVar, na.rm=TRUE),
           sulfurTot= mean(sulfurTot, na.rm=TRUE),
           dryMass= mean( dryMass, na.rm=TRUE),
           acidity= mean( acidity, na.rm=TRUE),
           DIEL = mean(DIEL, na.rm=T),
           FRC.kg.C.m.2 = mean(FRC.kg.C.m.2, na.rm=T), 
           SOC.kg.C.m.2 = mean(SOC.kg.C.m.2, na.rm=T),
           beta_roots= mean(beta_roots, na.rm=T),  
           beta_soc= mean(beta_soc, na.rm=T),
           MAT= mean(MAT, na.rm=T), 
           MAP= mean(MAP, na.rm=T),
           Source_Sink = case_when( DIEL > 0 ~ "Source",
                                    DIEL <= 0 ~ "Sink") %>% as.factor)  %>% 
  left_join( metadata %>% dplyr::select(EcoType, SITE_ID), by = "SITE_ID")

parms.site.season.sf <-site.att.sf %>% left_join( parms.site.season , by="SITE_ID")

parms.site.sf <-site.att.sf %>% left_join( parms.site , by="SITE_ID")



# Figures: ####

# Seasonal methane pattern ####


NorthAmerica_AOI <- aoi.usa <- aoi_get(country = c('North America') )

map.total <- ggplot( ) + geom_sf(data= NorthAmerica_AOI, fill = "black")+
  geom_sf( data= parms.site.season.sf %>% arrange(desc(DIEL*1000)) ,
           aes( size = abs(DIEL*1000),
                col = DIEL*1000), alpha = 0.75) + theme_bw() +
  scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", 
                          direction=1,
                          limit = c(-0.015, 0.015),
                          oob = scales::squish) + 
  facet_wrap( ~ Season) +                
  scale_size_continuous(range = c(1, 5)) +
  theme(strip.background = element_rect(fill = "transparent")) + 
  labs(col=expression(paste("FCH"[4], "( mg m"^-2, "day"^-1, ")")),
       size =expression(paste("|FCH"[4], "| ( mg m"^-2, "day"^-1, ")")))

  density.total <-  parms.site.season.sf  %>% 
  ggplot(aes(DIEL*1000, fill = Season)) +
  geom_density( position = "stack", aes(col=Season), alpha = 0.5) + theme_bw()+
    scale_color_paletteer_d("ggthemes::Seattle_Grays") +
    scale_fill_paletteer_d("ggthemes::Seattle_Grays") +
    geom_vline( xintercept=0, linetype = "dashed")+ xlim( -0.05, 0.05) + 
    ylab("Density") + 
    xlab( expression(paste("Mean Seasonal FCH"[4], "( mg m"^-2, "day"^-1, ")")))
      
  total.season.final.plot <- ggarrange( map.total , 
                                        density.total, 
                                        labels=c("A", "B"), ncol=1,
                                        heights= c( 2, 1))
  
ggsave("FIGURES/Total_methane_season_plot.png", 
       plot =  total.season.final.plot, width = 7.6, 
       height = 7.3, units = "in")
  
# Q10 and base respiration by season and ecosystem type: ####
  
  source_colors <- c("Sink" = "red3", "Source" = "blue4")
  
 Q10.source.density.plots <-  parms.site.season  %>%   filter(!is.na(Source)) %>% 
    ggplot(aes( x = Q10.mean, y = after_stat(scaled), col= Source)) +
    geom_density(alpha=0.5)  +
    theme_bw() + ylab("Density") + xlab("Q10") + scale_color_manual(values = source_colors)
 
 
 Q10.source.box.plots <-  parms.site.season  %>%   filter(!is.na(Source)) %>% 
    ggplot(aes( x = Source, y = Q10.mean, col=Source)) +
    geom_boxplot(alpha=0.5)  + ylim(1,2) +
    theme_bw()  + stat_summary(
      fun = median, 
      geom = "text", 
      col = source_colors ,     # Color of the text
      vjust = -10, aes(label = round(..y.., digits = 2)) ) +
    scale_color_manual(values = source_colors)+ ylab("Q10") + 
    xlab("") +  theme(legend.position = "none")
  
 Q10.source.box.seasons.plots <- parms.site.season  %>%   filter(!is.na(Source)) %>% 
    ggplot(aes( x = Source, y = Q10.mean, col = Source)) +
    geom_boxplot(alpha=0.5)  + ylim(1,2) +
    theme_bw()  + stat_summary(
      fun = median, 
      geom = "text",     # Color of the text
      vjust = -13,
      size=3,
      aes(label = round(..y.., digits = 2)) )+ facet_wrap( ~ Season, ncol=4) + scale_color_manual(values = source_colors)+
   ylab("Q10") + xlab("") + theme(strip.background = element_rect(fill = "transparent", color = "black"),
     legend.position = "none")
  
  
DIEL.source.box.plots <- parms.site.season  %>%   filter(!is.na(Source)) %>% 
    ggplot(aes( x = Source, y = DIEL*1000, col = Source)) +
    geom_boxplot(alpha=0.5)  + ylim(-0.05,0.05) +
    theme_bw() +scale_color_manual(values = source_colors)+ theme(strip.background = element_rect(fill = "transparent", color = "black"), legend.position = "none") + ylab( expression(paste("Mean FCH"[4], " ( mg m"^-2, "day"^-1, ")")))
  
  
DIEL.source.box.seasons.plots <- parms.site.season  %>%   filter(!is.na(Source)) %>% 
    ggplot(aes( x = Source, y = DIEL*1000, col = Source)) +
    geom_boxplot(alpha=0.5)  + ylim(-0.05,0.05) +
    theme_bw() + facet_wrap( ~ Season, ncol=4)+scale_color_manual(values = source_colors)+ theme(strip.background = element_rect(fill = "transparent", color = "black"),
                                                                                                 legend.position = "none") + ylab( expression(paste("Mean FCH"[4], " ( mg m"^-2, "day"^-1, ")")))
  
DIEL.plots <-ggarrange(DIEL.source.box.plots,
DIEL.source.box.seasons.plots, ncol=1, labels= c("A", "B") )
  

Q10.plots <- ggarrange(
  ggarrange(Q10.source.box.plots,Q10.source.density.plots, ncol=2, labels= c("A", "B"), common.legend = TRUE),
ggarrange (Q10.source.box.seasons.plots,  labels= c("C") ), ncol=1)
  
  
ggsave("FIGURES/Q10_plots.png", plot =   Q10.plots, width =6, height =5, units = "in")

ggsave("FIGURES/DIEL_plot.png", plot =   DIEL.plots , width = 5, height = 6, units = "in")

  # Total Methane Figures Annual: ####
  library(ggpmisc)

map.DIEL.total <- ggplot( ) + geom_sf(data= NorthAmerica_AOI, fill = "black")+
  geom_sf( data= parms.site.sf %>% arrange(desc(DIEL )) ,
           aes( size = abs(DIEL*1000),
                col =DIEL*1000), alpha = 0.75) + theme_bw() +
  scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", 
                          direction=1,
                          limit = c(-0.015, 0.015),
                          oob = scales::squish,
                          breaks=pretty_breaks(n = 3)) +                
  scale_size_continuous(range = c(3, 5)) +
  theme(strip.background = element_rect(fill = "transparent")) + 
  labs(col=expression(paste("FCH"[4], "( mg m"^-2, "day"^-1, ")")),
       size =expression(paste("|FCH"[4], "| ( mg m"^-2, "day"^-1, ")")))

 MAT.MAP.plot <-  parms.site %>% filter(!is.na(Source)) %>% 
    ggplot(aes( x = MAT, y = MAP)) +
    geom_point(alpha=0.5, aes(col = DIEL*1000, size = DIEL*1000)) + 
    theme_bw()  +scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", 
                                         direction=1,
                                         limit = c(-0.015, 0.015),
                                         oob = scales::squish)+
   scale_size_continuous(range = c(1, 5))+
    theme(panel.background = element_rect(fill = "black"), # Plotting area
          panel.grid.major = element_line(color = "grey20"), # Optional: adjust grid visibility
      panel.grid.minor = element_line(color = "grey30") ) 
  
 DIEL.MAT.plot <-  parms.site  %>% 
   ggplot(aes( x = MAT, y = DIEL*1000)) +
   geom_point(alpha=0.75, aes(col = DIEL*1000), size = 3) + 
   theme_bw() + geom_smooth( method="loess", col='white') +
   stat_cor(method = "spearman", col="white",
            label.y.npc = "center") +
   facet_wrap( ~ Source_Sink + EcoType, ncol=3, scale="free") +
   scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", 
                           direction=1,
                           limit = c(-0.04, 0.04))  + 
   ylab( expression(paste("FCH"[4], " ( mg m"^-2, "day"^-1, ")"))) +
   theme(panel.background = element_rect(fill = "black"), # Plotting area
         panel.grid.major = element_line(color = "grey20"), 
         panel.grid.minor = element_line(color = "grey30"),
         legend.position = "none" ,
         strip.background = element_rect(fill = "transparent", color = "black"))
 
 
  annual.ch4.plots <- ggarrange( ggarrange(   map.DIEL.total,
               MAT.MAP.plot, labels=c("A", "B") , 
               common.legend = T, widths=c(2,1.25)), 
             DIEL.MAT.plot, labels=c(" ", "C"), ncol =1)

  
  ggsave("FIGURES/Total_Annual_CH4_plot.png", plot =     
           annual.ch4.plots , width = 10, height = 11, units = "in")
  
# MDS: ####
  
continious.variables <- c( 'VSWCMean' ,  'Q10.mean' , 'Rref.mean' , 'VSWCVar' , 'sulfurTot' , 'dryMass' , 'acidity','FRC.kg.C.m.2', 'SOC.kg.C.m.2', 'beta_roots' ,'beta_soc')
  
library(MASS)
library(ggrepel)
library(corrplot)
parms.mds <- parms.site %>% drop_na()
parms.mds %>% names
M = cor(  parms.mds[, 3: 17]  %>% drop_na())
corrplot(M, method = 'color', order = 'alphabet') # colorful number

library(vegan)

d <- vegdist(  parms.mds %>% dplyr::select( continious.variables))
fit <- metaMDS(d, k=2)

parms.mds$MDS1 <- fit$points[,1]
parms.mds$MDS2<- fit$points[,2]

cluster_res <- hclust(d, method = "ward.D2")

groups <- cutree(cluster_res, k = 4)
parms.mds$cluster <- groups # Add clusters to the df

library(paletteer)

plot.MDS <-   parms.mds %>%  ggplot(aes( x=MDS1, y = MDS2, col=cluster %>% as.factor)) +
  stat_ellipse()+
  geom_text_repel(aes(label = SITE_ID),
                   box.padding   = 0, 
                   point.padding = 0,
                  size =3) +
  geom_point( alpha = 0.25 ) +  theme_bw() + 
  scale_color_paletteer_d("ggsci::nrc_npg") +theme(legend.position = "none") +
  xlab( 'Dimension 1')+ ylab( 'Dimension 2')


# Additional plots:

plot.diel <-  parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =DIEL*1000, col=cluster %>% as.factor)) + geom_boxplot( )+ theme_bw() + geom_hline(yintercept = 0, col = "black", linetype="dashed") + 
  ylab(expression(paste(" Mean FCH"[4], "( mg m"^-2, "day"^-1, ")"))) + xlab("Cluster") + 
  scale_color_paletteer_d("ggsci::nrc_npg")+ theme(legend.position = "none") + ylim(-0.01, 0.01)

plot.source <-  parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =Source %>% as.numeric-1, col=cluster %>% as.factor)) + geom_violin( ) + theme_bw() +  stat_summary(fun = mean, geom = "crossbar", width = 0.2) + 
  scale_color_paletteer_d("ggsci::nrc_npg") + theme(legend.position = "none") + ylab( "Number of Month as a Source")+ xlab("Cluster")

plot.Q10 <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =Q10.mean , col=cluster %>% as.factor)) + geom_boxplot( ) + 
  theme_bw() + xlab("Cluster") + ylab('Q10') +
  stat_summary(fun = "median", 
               geom = "text", 
               aes(label = round(after_stat(y), 2)), 
               vjust = 1.5, # Moves text above the point
               size = 3) +  theme(legend.position = "none") +  scale_color_paletteer_d("ggsci::nrc_npg")

plot.sulfur <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =sulfurTot, col=cluster %>% as.factor)) + geom_boxplot( ) + theme_bw() + xlab("Cluster") + ylab('Soil Sulfur') +  
  scale_color_paletteer_d("ggsci::nrc_npg")+  theme(legend.position = "none") 


plot.rootmass <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =dryMass , col=cluster %>% as.factor)) + geom_boxplot( ) +  scale_color_paletteer_d("ggsci::nrc_npg")+  theme(legend.position = "none")  + xlab("Cluster") + 
  ylab('Root Dry Mass') + theme_bw() +  theme(legend.position = "none")

plot.RootMass <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =FRC.kg.C.m.2 , col=cluster %>% as.factor)) + geom_boxplot( )+  theme(legend.position = "none") +  scale_color_paletteer_d("ggsci::nrc_npg") +
  xlab("Cluster") +  ylab(expression(paste('Root Mass (kg C m'^-2,')'))) + theme_bw() +  theme(legend.position = "none")

plot.soc <-  parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =SOC.kg.C.m.2, col=cluster %>% as.factor)) + geom_boxplot( ) +  theme(legend.position = "none") +  scale_color_paletteer_d("ggsci::nrc_npg") +
  xlab("Cluster") +  ylab(expression(paste('SOC (kg C m'^-2,')'))) + theme_bw() +  theme(legend.position = "none")

plot.betaRoots <-   parms.mds%>% ggplot(aes( x=cluster %>% as.factor, y =beta_roots, col=cluster %>% as.factor)) +  
  theme(legend.position = "none") +  scale_color_paletteer_d("ggsci::nrc_npg") +
  xlab("Cluster") +  ylab(expression(paste(beta ,' Roots'))) + geom_boxplot( ) +  theme(legend.position = "none") + theme_bw() +  theme(legend.position = "none")

plot.betaSOC <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =beta_soc, col=cluster %>% as.factor)) +  
  theme(legend.position = "none") +  scale_color_paletteer_d("ggsci::nrc_npg") +
  xlab("Cluster") +  ylab(expression(paste(beta ,' SOC'))) + geom_boxplot( ) +  theme(legend.position = "none") + theme_bw() +  theme(legend.position = "none")

plot.soilAcidity <- parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =acidity, col=cluster %>% as.factor)) + geom_boxplot( ) +  scale_color_paletteer_d("ggsci::nrc_npg") + theme_bw()+  theme(legend.position = "none") + ylab("Soil Acidity") + xlab('Cluster')

plot.VSWCVAR <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =VSWCVar, col=cluster %>% as.factor)) + geom_boxplot( ) +  scale_color_paletteer_d("ggsci::nrc_npg") +
  ylab("SVWC Variability") + xlab('Cluster') + theme_bw()+  theme(legend.position = "none")

plot.VSWCMean <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =VSWCMean, col=cluster %>% as.factor)) + geom_boxplot( ) +  scale_color_paletteer_d("ggsci::nrc_npg") +
  ylab("Mean SVWC") + xlab('Cluster') + theme_bw()+  theme(legend.position = "none") 

# Maps: ####
site.att.sf.cluster <- site.att.sf %>% left_join(by='SITE_ID',   parms.mds )

map.cluster <- ggplot( ) + geom_sf(data= NorthAmerica_AOI, fill = "white")+
  geom_sf( data= site.att.sf.cluster ,
           aes( col = cluster %>% as.factor), alpha = 0.5, size = 4) + theme_bw() + scale_color_paletteer_d("ggsci::nrc_npg",na.translate = FALSE)  + theme(strip.background = element_rect(fill = "transparent"))+ labs(col="Cluster")

# Panels ####

# Patterns in FCH4


panel.methane1 <- ggarrange( map.cluster, plot.MDS ,
                             labels=c("A", "B"), nrow = 2, ncol =1, common.legend = TRUE)

panel.methane2 <- ggarrange( plot.diel, plot.source, plot.Q10 ,
                             labels=c("C", "D", "E"), nrow = 3, ncol =1)

final.MDS.plot <- ggarrange( panel.methane1, panel.methane2,
           widths = c(2,1), common.legend = TRUE)

ggsave("FIGURES/MDS_PLOT.png", plot =final.MDS.plot , width = 8, height = 8.7, units = "in")

# Patterns in Soil Conditions
panel.root_soil <- ggarrange(plot.RootMass,
                             plot.betaRoots,
                             plot.soc,
                             plot.sulfur ,
                             plot.soilAcidity ,
                             plot.VSWCMean ,
                             plot.VSWCVAR,
                             plot.betaSOC,
                             labels = c( "A", "B", "C", "D", "E", "F", "G", "H"))

ggsave("FIGURES/Roots_Soils_PLOT.png", plot =panel.root_soil , width = 8, height = 8, units = "in")

# RESULTS: ####

total.daily <- parms.site.season$DIEL*1000 
total.daily%>% summary

total.daily.annual <- parms.site$DIEL*1000 
total.daily.annual%>% summary


nrow( parms.site %>%  filter( Source_Sink == "Source"))
nrow( parms.site %>%  filter( Source_Sink == "Sink"))

parms.site$SITE_ID[ parms.site$Source == 4] %>% unique %>% length
parms.site$SITE_ID[ parms.site$Source == 0] %>% unique %>% length

parms.site.source <- parms.sf %>% filter( DIEL > 0)
# Measure the significance of difference between sinks and sources for Q10.
kruskal.test(Q10.mean ~ Source_Sink , data =parms.mds) 

kruskal.test(Q10.mean ~ Season , data =parms.site.season) 

kruskal.test(Q10.mean ~ Source, data = parms.site.season %>% filter(Season == "Winter")) 

# Explore differences in the clusters:

parms.mds.cluster1 <- parms.mds %>% filter(cluster == 1)
mean(parms.mds.cluster1$DIEL *1000)
sd(parms.mds.cluster1$DIEL *1000)

mean(parms.mds.cluster1$Q10.mean)
sd(parms.mds.cluster1$Q10.mean)

mean(parms.mds.cluster1$sulfurTot)
sd(parms.mds.cluster1$sulfurTot)

mean(parms.mds.cluster1$Q10.mean)
sd(parms.mds.cluster1$Q10.mean)

mean(parms.mds.cluster1$acidity)
sd(parms.mds.cluster1$acidity)

mean(parms.mds.cluster1$FRC.kg.C.m.2)
sd(parms.mds.cluster1$FRC.kg.C.m.2)


parms.mds.cluster2 <- parms.mds %>% filter(cluster == 2)
mean(parms.mds.cluster2$Q10.mean)
sd(parms.mds.cluster2$Q10.mean)

mean(parms.mds.cluster2$DIEL *1000)
sd(parms.mds.cluster2$DIEL*1000)

mean(parms.mds.cluster2$sulfurTot)
sd(parms.mds.cluster2$sulfurTot)

mean(parms.mds.cluster2$acidity)
sd(parms.mds.cluster2$acidity)

mean(parms.mds.cluster2$SOC.kg.C.m.2)
sd(parms.mds.cluster2$SOC.kg.C.m.2)

mean(parms.mds.cluster2$VSWCVar)
sd(parms.mds.cluster2$VSWCVar)

parms.mds.cluster3 <- parms.mds %>% filter(cluster == 3)

mean(parms.mds.cluster3$DIEL*1000)
sd(parms.mds.cluster3$DIEL*1000)

mean(parms.mds.cluster3$sulfurTot)
sd(parms.mds.cluster3$sulfurTot)

mean(parms.mds.cluster3$acidity)
sd(parms.mds.cluster3$acidity)

mean(parms.mds.cluster3$SOC.kg.C.m.2)
sd(parms.mds.cluster3$SOC.kg.C.m.2)

mean(parms.mds.cluster3$VSWCMean)
sd(parms.mds.cluster3$VSWCMean)

parms.mds.cluster4 <- parms.mds %>% filter(cluster == 4)

mean(parms.mds.cluster4$DIEL*1000)
sd(parms.mds.cluster4$DIEL*1000)

mean(parms.mds.cluster4$sulfurTot)
sd(parms.mds.cluster4$sulfurTot)

mean(parms.mds.cluster4$VSWCMean)
sd(parms.mds.cluster4$VSWCMean)

mean(parms.mds.cluster4$acidity)
sd(parms.mds.cluster4$acidity)

mean(parms.mds.cluster4$VSWCVar)
sd(parms.mds.cluster4$VSWCVar)

mean(parms.mds.cluster4$beta_soc)
sd(parms.mds.cluster4$beta_soc)

parms.mds %>% reframe( .by=cluster, 
                       DIEL = median(DIEL)*1000 %>% round(4))

                       