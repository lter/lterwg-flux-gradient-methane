# Create a climate summary by YearMon for NEON Sites:
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)

localdir <- '/Volumes/MaloneLab/Research/FluxGradient/METHANE'

load( fs::path(localdir,paste0("/SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata")) )

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
Climate_Summary_YearMon <- data.frame()
Climate_Summary_Season <- data.frame()
site.list <- SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage %>% names 

for( site in site.list){
  print(site)
  
  data.yearmon <- SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage[[site]] %>% mutate(
    hour = format(time.rounded,'%H'),
    YearMon = format(time.rounded, "%Y-%m")) %>% distinct %>% reframe( 
      .by=YearMon,
      Tair_C_mean = mean(Tair_C, na.rm=T ),
      Tair_C_min = min(Tair_C, na.rm=T ),
      Tair_C_max = max(Tair_C, na.rm=T ),
      Tair_C_var = var(Tair_C, na.rm=T )) %>% mutate(Site = site)
  
  try( Climate_Summary_YearMon <- rbind( Climate_Summary_YearMon, data.yearmon ), silent =TRUE)

  # Seasonal Summary:
  
  data.season <- SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage[[site]] %>% mutate(
    hour = format(time.rounded,'%H'),
    YearMon = format(time.rounded, "%Y-%m")) %>% distinct %>% reframe( 
      .by=season,
      Tair_C_mean = mean(Tair_C, na.rm=T ),
      Tair_C_min = min(Tair_C, na.rm=T ),
      Tair_C_max = max(Tair_C, na.rm=T ),
      Tair_C_var = var(Tair_C, na.rm=T )) %>% mutate(Site = site)
  
  try( Climate_Summary_Season <- rbind( Climate_Summary_Season, data.season ), silent =TRUE)
  
}

save(Climate_Summary_YearMon, 
     Climate_Summary_Season, 
     file = fs::path(localdir,paste0("NEON_Climate_Summary.Rdata")))

