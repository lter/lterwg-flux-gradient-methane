# Ensemble NEON:
library(tidyverse)

# Import the RSHP information:
load(fs::path(localdir,paste0("NEON_RSHP_CCC.Rdata")) )
# Import the filtered data:
load(fs::path(localdir,paste0("SITE_DATA_FILTERED_CH4.RDATA")))

# Now 
SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE <- list()


for( site in site.list){
  print( site)
  
  VAL_CGF_RSHP_sum_site <- NEON_RSHP_EVAL %>% filter( Site == site)
  VAL_RSHP_sum_site <- NEON_RSHP_EVAL %>% filter( Site == site)

  
  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]] <- SITE_DATA_FILTERED[[site]] %>% 
    full_join(VAL_RSHP_sum_site, by=c('gas', 'Approach', 'dLevelsAminusB' )) %>% 
    filter(gas == "CH4") %>% 
    mutate(month = format(timeEndA.local,'%m') %>% as.numeric,
           season = case_when(
             month %in% c(12, 1, 2) ~ "Winter",
             month %in% c(3, 4, 5) ~ "Spring",
             month %in% c(6, 7, 8) ~ "Summer",
             TRUE ~ "Autumn"),
           hour = format(timeEndA.local,'%H'),
           count= case_when( is.na(FG_mean) == FALSE ~ 1, TRUE ~ 0),
           FG_mean_RSHP.CO2 = case_when(RSHP.CO2 == '1' ~ FG_mean, TRUE ~ 0),
           FG_mean_RSHP.H2O = case_when(RSHP.H2O == '1' ~ FG_mean, TRUE ~ 0), 
           FG_mean_RSHP.rf = case_when(RSHP.rf == '1' ~FG_mean, TRUE ~ 0)) %>% 
    distinct  %>% 
    mutate( time.rounded = timeEndA.local %>% round_date( unit = "30 minutes") ) %>%
    reframe( .by= c(time.rounded), 
             FG_ENSEMBLE_RSHP.CO2= mean(FG_mean_RSHP.CO2, na.rm=T),
             FG_ENSEMBLE_RSHP.H2O= mean(FG_mean_RSHP.H2O, na.rm=T),
             FG_ENSEMBLE_RSHP.rf= mean(FG_mean_RSHP.rf, na.rm=T),
             count = sum(count),
             PAR = mean(PAR, na.rm=T),
             Tair_C = mean(Tair_C, na.rm=T)) %>% mutate(month = format(time.rounded,'%m') %>% as.numeric,
                                                               season = case_when(
                                                                 month %in% c(12, 1, 2) ~ "Winter",
                                                                 month %in% c(3, 4, 5) ~ "Spring",
                                                                 month %in% c(6, 7, 8) ~ "Summer",
                                                                 TRUE ~ "Autumn"),
                                                               hour = format(time.rounded,'%H')) 
  
  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]] %>% ggplot( ) + 
    geom_point( aes( x= FG_ENSEMBLE_RSHP.CO2, y=FG_ENSEMBLE_RSHP.H2O))
  
  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]] %>% ggplot( ) + 
    geom_point( aes( x= FG_ENSEMBLE_RSHP.CO2, y=FG_ENSEMBLE_RSHP.rf))
  
  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]] %>% ggplot( ) + 
    geom_point( aes( x= FG_ENSEMBLE_RSHP.H2O, y=FG_ENSEMBLE_RSHP.rf))
  
}

SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE # filtered for RSHP
# ENSEMBLE:


save( SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE ,file=fs::path(localdir,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata")))
googledrive::drive_upload(media = fs::path(localdir,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata")), overwrite = T, path = drive_url)
