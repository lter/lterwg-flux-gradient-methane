# Ensemble NEON:
library(tidyverse)


# Import the RSHP information:
load(fs::path(localdir,paste0("NEON_RSHP_CCC.Rdata")) )
# Import the filtered data:
load(fs::path(localdir,paste0("SITE_DATA_FILTERED_CH4.RDATA")))


# Separate the Approach and level: 
NEON_RSHP_EVAL <- NEON_RSHP_EVAL %>% separate(
    col = var, 
    into = c("Approach", "dLevelsAminusB"),
    sep = "-")



# Now 
SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE <- list()


for( site in site.list){
  print( site)
  
  VAL_CGF_RSHP_sum_site <- NEON_RSHP_EVAL %>% filter( Site == site)
  VAL_RSHP_sum_site <- NEON_RSHP_EVAL %>% filter( Site == site)

  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]] <- SITE_DATA_FILTERED[[site]] %>% 
    full_join(VAL_RSHP_sum_site, by=c('gas', 'Approach', 'dLevelsAminusB' )) %>% filter(RSHP.rf == 1 ) %>% filter(gas == "CH4") %>% 
    mutate(month = format(timeEndA.local,'%m') %>% as.numeric,
           season = case_when(
             month %in% c(12, 1, 2) ~ "Winter",
             month %in% c(3, 4, 5) ~ "Spring",
             month %in% c(6, 7, 8) ~ "Summer",
             TRUE ~ "Autumn"),
           hour = format(timeEndA.local,'%H'),
           count= case_when( is.na(FG_mean) == FALSE ~ 1,
                             TRUE ~ 0)) %>% distinct  %>% 
    mutate( time.rounded = timeEndA.local %>% round_date( unit = "30 minutes") ) %>% 
    filter( RSHP.rf == 1) %>%
    reframe( .by= c(time.rounded), 
             FG_ENSEMBLE = mean(FG_mean, na.rm=T),
             count = sum(count),
             PAR = mean(PAR, na.rm=T),
             Tair_C = mean(Tair_C, na.rm=T)) %>% mutate(month = format(time.rounded,'%m') %>% as.numeric,
                                                               season = case_when(
                                                                 month %in% c(12, 1, 2) ~ "Winter",
                                                                 month %in% c(3, 4, 5) ~ "Spring",
                                                                 month %in% c(6, 7, 8) ~ "Summer",
                                                                 TRUE ~ "Autumn"),
                                                               hour = format(time.rounded,'%H')) 
  

  
}

SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE # filtered for RSHP
# ENSEMBLE:

fileSave <- fs::path(localdir,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata"))
save( SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE ,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)
