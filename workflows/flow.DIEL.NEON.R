# The Diel analysis is currently set up by season for the ENSEMBLE data:

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)

localdir <- '/Volumes/MaloneLab/Research/FluxGradient/FluxData'

load( fs::path(localdir,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata")) )

DirRepo.CH4 <- "/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-methane"

source( fs::path(DirRepo.CH4, paste('/functions/calc_diel_ch4.R')))
source( fs::path(DirRepo.CH4, paste('/functions/calc_Q10.R')))

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

# Calculate diels by Season: ####
ENSEMBLE_DIELS <- data.frame()
ENSEMBLE_Q10_eq4 <- data.frame()
ENSEMBLE_Q10_eq5 <- data.frame()

site.list <- SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE %>% names 

options <- c('FG_ENSEMBLE_RSHP.CO2',
             'FG_ENSEMBLE_RSHP.H2O', 
             'FG_ENSEMBLE_RSHP.rf' )

for( site in site.list){
  print(site)
  
  for (opt in options){
    print(paste("Running option =", opt, sep=" "))
    
  data <- SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]] %>% 
    mutate(hour = format(time.rounded,'%H'),
           YearMon = format(time.rounded, "%Y-%m")) %>% distinct

  data$FG_mean <- data[, opt]
  message( paste("Running DIEL functions- CO2 for ", site))
  # Calculate Diurnal Patterns by Year-month:
  
  for( i in data$YearMon %>% unique){
    
    data.sub <- data %>% filter( YearMon == i)
    
    try(DIEL.CH4 <- DIEL.season.ch4( dataframe = data.sub, 
                                     flux = opt, 
                                     Gas = "CH4") %>% 
          mutate(gas= "CH4", 
                 YearMon = i,
                 site = site,
                 RSHP = opt), silent =TRUE)
    
    try(DIEL.CH4 %>%  ggplot() + geom_point(aes(x=Hour, y = DIEL, col=season)) + facet_wrap(~season) , silent =TRUE)
    
    try( ENSEMBLE_DIELS <- rbind( ENSEMBLE_DIELS, DIEL.CH4 ), silent =TRUE)
    
    try( rm( DIEL.CH4), silent=T)
  }
  
  # Q10:
 
  try(NEON_TRC_PARMS_04 <- TRC_PARMS_04(data.frame = data,
                                       iterations = 5000,
                                       priors.trc = brms::prior("normal(2.0, 0.3)", nlpar = "Q10", lb = 1.0, ub = 5) +
                                         brms::prior("normal(0.5, 0.3)", nlpar = "Rref", lb = 0.001, ub = 5),
                                       idx.colname = 'YearMon',
                                       NEE.colname = "FG_mean",
                                       TA.colname = 'Tair_C',
                                       Tref = 1)  %>%  mutate( SITE_ID = site, RSHP = opt), silent=T)
  
  try(ENSEMBLE_Q10_eq4 <- rbind(ENSEMBLE_Q10_eq4 ,NEON_TRC_PARMS_04
                                           ) , silent=T)
  }
  }

ENSEMBLE_DIELS_Site <- ENSEMBLE_DIELS %>% 
  reframe( .by=c(site, gas, season, RSHP), 
           total.FG = sum(DIEL, na.rm=T),
           max.FG = max(DIEL, na.rm=T),
           min.FG = min(DIEL, na.rm=T),
           total.count = sum(count))

save(ENSEMBLE_DIELS ,
     ENSEMBLE_Q10_eq4,
     ENSEMBLE_Q10_eq5 ,
     ENSEMBLE_DIELS_Site,
     file = fs::path(localdir,paste0("NEON_PARMS_DIEL_Q10.Rdata")))
