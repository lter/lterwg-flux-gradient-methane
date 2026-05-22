load("/Volumes/MaloneLab/Research/FluxGradient/NEON_Storage_Flux/NEON_storage.RData")
load(fs::path(localdir.ch4 ,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata")))

site.list <- names(SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE)


SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage <- list()

for(site in site.list){
  
  flux.site <- SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE[[site]]
  
 
  flux.storage.site <- NEON_storage [[site]] %>% mutate( time.rounded = time_halfhour_storage_end) %>% filter( gas == "CH4") %>% select( time.rounded, storage_flux_filled)
    
  
  SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage[[ site]] <-   flux.site %>% 
    full_join( flux.storage.site, by = 'time.rounded') %>% 
    mutate( flux_total =  case_when( !is.na(FG_ENSEMBLE_RSHP) ~ FG_ENSEMBLE_RSHP + storage_flux_filled,  TRUE ~ NA))

  print( paste0("Done with ", site))
  
  
}

save(   SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_storage,
        file=fs::path(localdir.ch4 ,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata")))
