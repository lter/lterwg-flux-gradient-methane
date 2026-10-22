# Ensemble NEON:
library(tidyverse)

# Import the RSHP information:
load(fs::path(localdir.ch4 ,paste0("NEON_RSHP_CCC.Rdata")) )

# Import the filtered data:
load(fs::path(localdir.ch4, paste0("/NEON_GradientFlux_Data_Filter/SITE_DATA_FILTERED_CH4.Rdata")))

# Now 
SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE <- list()


for( site in site.list){
  print( site)
  
  VAL_RSHP_sum_site <- NEON_RSHP_EVAL %>% filter( Site == site)

  
  site_flux <- SITE_DATA_FILTERED[[site]]
  if (!"gas" %in% names(site_flux) || anyNA(site_flux$gas)) {
    stop("Missing gas labels for ", site, "; regenerate the filtered inputs.")
  }
  # Existing saved files can contain CO2/H2O as well as CH4.
  site_flux <- site_flux %>% filter(gas == "CH4")
  if (nrow(site_flux) == 0L) stop("No CH4 observations for ", site)

  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]] <- site_flux %>%
    full_join(VAL_RSHP_sum_site, by=c('Approach', 'dLevelsAminusB' )) %>% 
    mutate(month = format(timeEndA.local,'%m') %>% as.numeric,
           season = case_when(
             month %in% c(12, 1, 2) ~ "Winter",
             month %in% c(3, 4, 5) ~ "Spring",
             month %in% c(6, 7, 8) ~ "Summer",
             TRUE ~ "Autumn"),
           hour = format(timeEndA.local,'%H'),
           count= case_when( is.na(FG_mean) == FALSE ~ 1, TRUE ~ 0),
           FG_mean_RSHP = case_when(RSHP.CO2 == '1' |RSHP.H2O == '1' ~ FG_mean,
                                    TRUE ~ NA)) %>% 
    distinct()  %>% 
    mutate( time.rounded = timeEndA.local %>% round_date( unit = "30 minutes") ) %>%
    reframe( .by= c(time.rounded), 
             FG_ENSEMBLE_RSHP= mean(FG_mean_RSHP, na.rm=T),
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
save( SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE ,file=fs::path(localdir.ch4 ,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata")))

googledrive::drive_upload(media = fs::path(localdir.ch4, paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata")), overwrite = T, path = drive_url)
