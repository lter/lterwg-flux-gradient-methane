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

for( site in site.list){
  print(site)
  
  data <- SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]] %>% mutate(
                                             hour = format(time.rounded,'%H'),
                                             YearMon = format(time.rounded, "%Y-%m")) %>% distinct
  
  message( paste("Running DIEL functions- CO2 for ", site))
  # Calculate Diurnal Patterns by Year-month:
  
  try(DIEL.CH4 <- DIEL.season.ch4( dataframe = data, 
                                        flux = 'FG_ENSEMBLE', 
                                         Gas = "CH4") %>% mutate(gas= "CH4", 
                                                                 site = site), silent =TRUE)

  try(DIEL.CH4 %>%  ggplot() + geom_point(aes(x=Hour, y = DIEL, col=season)) + facet_wrap(~season) , silent =TRUE)
  
  try( ENSEMBLE_DIELS <- rbind( ENSEMBLE_DIELS, DIEL.CH4 ), silent =TRUE)
  
  try( rm( DIEL.CH4), silent=T)
  
  # Q10:
 
  try(NEON_TRC_PARMS_04 <- TRC_PARMS_04(data.frame = data,
                                       iterations = 5000,
                                       priors.trc = brms::prior("normal(2.0, 0.3)", nlpar = "Q10", lb = 1.0, ub = 5) +
                                         brms::prior("normal(0.5, 0.3)", nlpar = "Rref", lb = 0.001, ub = 5),
                                       idx.colname = 'YearMon',
                                       NEE.colname = 'FG_ENSEMBLE',
                                       TA.colname = 'Tair_C',
                                       Tref = 1)  %>%  mutate( SITE_ID = site), silent=T)
  
  try(ENSEMBLE_Q10_eq4 <- rbind(ENSEMBLE_Q10_eq4 ,NEON_TRC_PARMS_04
                                           ) , silent=T)

  try(NEON_TRC_PARMS_05 <- TRC_PARMS_05(data.frame = data,
                                       iterations = 5000,
                                       priors.trc = brms::prior("normal(0.2 , 1)", nlpar = "a", lb = 0.1, ub = 1) +
                                         brms::prior("normal(0.5, 0.03)", nlpar = "b", lb = 0.001, ub = 0.9),
                                       idx.colname = 'YearMon',
                                       NEE.colname = 'FG_ENSEMBLE',
                                       TA.colname = 'Tair_C') %>%  mutate( SITE_ID = site), silent=T)
  
  
  try(ENSEMBLE_Q10_eq5 <- rbind(   ENSEMBLE_Q10_eq5 ,
                               NEON_TRC_PARMS_05), silent=T)
  
}

save(ENSEMBLE_DIELS ,
     ENSEMBLE_Q10_eq4,
     ENSEMBLE_Q10_eq5 , file = fs::path(localdir,paste0("NEON_PARMS_DIEL_Q10.Rdata")))

ENSEMBLE_DIELS_Site <- ENSEMBLE_DIELS %>% 
  reframe( .by=c(site, gas, season), 
           total.FG = sum(FG),
           max.FG = max(FG),
           min.FG = min(FG),
           total.EC = sum(EC),
           max.EC = max(EC),
           min.EC = min(EC),
           Day.DIFF = min.EC-min.FG,
           Night.DIFF = max.EC-max.FG,
           total.diff = sum(DIFF.DIEL),
           percent.diff = total.diff/ total.EC *100,
           total.count = sum(count),
           Day.over.est.count = case_when( Day.DIFF > 0 ~ 1),
           Day.under.est.count = case_when( Day.DIFF < 0 ~ 1),
           Night.over.est.count = case_when( Night.DIFF > 0 ~ 1),
           Night.under.est.count = case_when( Night.DIFF < 0 ~ 1))

# DIEL PLOTS: ####

plot.diel.gas.season.fg <- ENSEMBLE_DIELS  %>%  ggplot()+ 
  geom_point(aes(x = Hour , y = DIEL, col=season)) +
  scale_colour_discrete_qualitative(palette = "Harmonic") +
  theme_bw() + facet_wrap(~ gas, , scales = "free_y")+ 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "transparent", linewidth = 0.5),
        legend.title = element_blank()) + ylab(expression(paste( "GF (g m"^2, ")")))


