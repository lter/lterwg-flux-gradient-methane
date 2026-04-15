# Results 
library(tidyverse)
library( ggplot2)

load(file = fs::path(localdir,paste0("NEON_PARMS_DIEL_Q10.Rdata")))
load(file=paste(localdir, "/", 'Soildata_YearMon.Rdata', sep=""))
load( file = fs::path(localdir,paste0("NEON_Climate_Summary.Rdata")))

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
                              Vegetation.Abbreviation..IGBP. == 'OSH' ~ 'Shrubland')) %>% rename( Site = Site_Id.NEON)

 

ENSEMBLE_DIELS_Site <- ENSEMBLE_DIELS %>% 
  reframe( .by=c(site, gas, season), 
           total.FG = sum(DIEL),
           max.FG = max(DIEL),
           min.FG = min(DIEL),
           Night.DIFF = max.FG-max.FG,
           total.count = sum(count)) %>% rename(Site = site) %>%  
  left_join(metadata, by="Site") %>% 
  left_join(Climate_Summary_Season, by=c('Site', 'season') ) %>% 
  left_join(Site_SoilData_Season %>% rename(season = Season), by=c('Site', 'season')) %>% 
  mutate( source = case_when(total.FG > 0 ~ 1,
                             .default = 0) %>% as.factor)

# Join these files:

ENSEMBLE_Q10_eq5_soil <- ENSEMBLE_Q10_eq5 %>% rename( YearMon = idx, Site = SITE_ID) %>% 
  left_join(Site_SoilData, by = c( 'Site', 'YearMon')) %>% 
  left_join(Climate_Summary_YearMon,by=c('Site', 'YearMon'))
ENSEMBLE_Q10_eq4_soil <-ENSEMBLE_Q10_eq4  %>% rename( YearMon = idx, Site = SITE_ID) %>% 
  left_join(Site_SoilData, by = c( 'Site', 'YearMon'))%>% 
  left_join(Climate_Summary_YearMon,by=c('Site', 'YearMon'))

# DIEL PLOTS: ####

ENSEMBLE_DIELS_Site %>%  ggplot() + 
  geom_point(aes(x = total.FG , y = Site, col=VSWCMean )) +
  theme_bw() + geom_vline(xintercept = 0) 

ENSEMBLE_DIELS_Site %>%  ggplot(aes(x = VSWCMean , y = Tair_C_mean, col= source)) + 
  geom_point() + facet_wrap(~Vegetation.Abbreviation..IGBP.)
  theme_bw() + geom_smooth()
  ENSEMBLE_DIELS_Site %>% names()
  
  ggpairs(ENSEMBLE_DIELS_Site %>% filter("EcoType" != "WET" ), 
          columns=c(25,28, 30:34), 
          ggplot2::aes(color=source))


  ENSEMBLE_DIELS_Site  %>%  ggplot() + 
  geom_boxplot(aes(x = total.FG , y = Site )) +
  theme_bw() 



ENSEMBLE_Q10_eq4  %>%  ggplot() + 
  geom_boxplot(aes(x = Q10.mean, y = SITE_ID )) +
  theme_bw() 

ENSEMBLE_Q10_eq4  %>%  ggplot() + 
  geom_boxplot(aes(x = Rref.mean, y = SITE_ID )) +
  theme_bw() 

ENSEMBLE_Q10_eq4  %>%  ggplot() + 
  geom_point(aes(x = Rref.mean, y = Q10.mean )) +
  theme_bw() 


ENSEMBLE_Q10_eq4_soil  %>%  ggplot() + 
  geom_point(aes(x = Rref.mean, y = VSWCMax )) +
  theme_bw()


ENSEMBLE_Q10_eq4_soil  %>%  ggplot() + 
  geom_point(aes(x = VSWCMax, y = Site )) +
  theme_bw()

library(mgcv)

ENSEMBLE_DIELS_Site %>% names()
source.rf <- gam( source ~  season + EcoType + Climate.Class.Abbreviation..Koeppen. +
                    s(VSWCMean) + s(Tair_C_mean) + te(VSWCMax, Tair_C_mean)+ s(VSWCMax), te(Tair_C_mean, VSWCMean),
                           data =ENSEMBLE_DIELS_Site %>% filter("EcoType" != "WET" ),
                  method = "REML", random = ~ (1 | Site ),
                  family = binomial())
source.rf %>% summary

plot(source.rf)

 
  ggplot( data=ENSEMBLE_DIELS_Site, 
          aes( x=Tair_C_mean, y=total.FG)) +
  geom_point() +
  geom_smooth( method = "gam", formula = y  ~ s(x)) 
