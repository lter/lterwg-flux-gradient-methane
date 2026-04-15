# RSHP: ####
library( ggplot2)
source(fs::path(DirRepo.eval, 'functions/calc.One2One.CCC_testing.R'))

load(fs::path(localdir,paste0("SITEval_DATA_FILTERED_CH4.Rdata")) )
canopy.info <- read.csv( file.path(paste(localdir, "Val_canopy.csv", sep="/")))%>% rename(Site = SITE_ID)

VAL_RSHP <- CCC_SamplingHeightPairs(DATA = SITEval_DATA_FILTERED) 

file1 <- VAL_RSHP %>% 
  mutate(Combination = paste(var1, var2, sep ="-")) %>%
  rename( var= var1)%>% select(Site, gas, var, CCC )

VAL_RSHP_EVAL <- rbind( file1, file2) %>% 
  reframe( .by = c( Site, gas, var),
           CCC.GF.max = max( CCC, na.rm=T), 
           CCC.GF.min = min( CCC, na.rm=T), 
           CCC.GF.mean = mean(CCC, na.rm=T), 
           CCC.GF.median = median( CCC, na.rm=T),
           CCC.GF.range = CCC.GF.max - CCC.GF.min) %>% 
  separate(col = var, 
             into = c("Approach", "dLevelsAminusB"), 
             sep = "-") %>% 
  left_join(canopy.info, by=c('Site', 'dLevelsAminusB')) %>% filter( Canopy_L1 != "WW", Canopy_L1 != "WA" )

VAL_RSHP_EVAL$Canopy_L1 %>% unique

# Save the file:

save( VAL_RSHP,VAL_RSHP_EVAL,  file=fs::path(localdir,paste0("VAL_RSHP_CCC.Rdata")))
googledrive::drive_upload(media = fs::path(localdir,paste0("VAL_RSHP_CCC.Rdata")), overwrite = T, path = drive_url)

# Use the RSHP model to Filter THE DATA: #####

load(file= paste(localdir, 'SITE_RSHP_MODEL.Rdata', sep="") ) # The model is in here:

library(randomForest)
VAL_RSHP_EVAL$RSHP.rf <- predict(rf.good.ccc.ec, VAL_RSHP_EVAL)
VAL_RSHP_EVAL$RSHP.rf %>% summary

# Create Datasets with only the data in it from RSHP: 
site.list <- SITEval_DATA_FILTERED %>% names()

SITEval_DATA_FILTERED_RSHP <- list()

for( site in site.list){
  
  VAL_RSHP_sum_site <- VAL_RSHP_EVAL %>% filter( Site == site) %>% na.omit
 
  SITEval_DATA_FILTERED_RSHP[[site]] <- SITEval_DATA_FILTERED[[site]] %>% 
    full_join(VAL_RSHP_sum_site, by=c('gas', 'Approach', 'dLevelsAminusB' )) %>% 
    filter(RSHP.rf ==1) %>%  
    mutate(month = format(timeEndA.local,'%m') %>% as.numeric,
           Season = case_when(
             month %in% c(12, 1, 2) ~ "Winter",
             month %in% c(3, 4, 5) ~ "Spring",
             month %in% c(6, 7, 8) ~ "Summer",
             TRUE ~ "Autumn"),
           hour = format(timeEndA.local,'%H'))
  
}

SITEval_DATA_FILTERED_RSHP

save( SITEval_DATA_FILTERED_RSHP ,file=fs::path(localdir,paste0("SITEval_DATA_FILTERED_RSHP_ENSEMBLE.Rdata")))


# Compare FG to EC: ####
source(fs::path(DirRepo.eval, 'functions/calc.One2One.CCC_testing.R'))

SITEval_DATA_FILTERED <-  SITEval_DATA_FILTERED 

# Add canopy information:
SITEval_DATA_FILTEREDc <- list()
for( site in site.list){
  print(site)
  SITEval_DATA_FILTEREDc[[site]] <- SITEval_DATA_FILTERED[[site]] %>% mutate( Site = site) %>% left_join(canopy.info, by=c('Site', 'dLevelsAminusB')) %>% filter( Canopy_L1 != "WW", Canopy_L1 != "WA" ) 
}

# Data now only contains AW and AA:

CCC_RSHP <- ccc.val(sites.tibble =  SITEval_DATA_FILTERED_RSHP) 
CCC_VAL <- ccc.val(sites.tibble =  SITEval_DATA_FILTEREDc)


CCC_RSHP  %>%  ggplot( aes( x= CCC, y =Site, col=Approach)) + geom_point() + theme_bw() + facet_wrap(~ Season +gas, nrow=3) + xlim(-1, 1)

CCC_VAL %>%  ggplot( aes( x= CCC, y =Site, col=Approach)) + geom_point() + theme_bw() + facet_wrap(~ Season +gas, nrow=3) + xlim(-1, 1)

CCC_RSHP$CCC %>% summary

# Save the files: 
save( CCC_RSHP,CCC_VAL, file= fs::path(localdir,paste0("CCC_CH4.Rdata")))

# Add RSHP information to the fluxes:


SITEval_DATA_FILTERED_RSHPc <- list()

for( site in site.list){
  
  VAL_RSHP_sum_site <- VAL_RSHP_EVAL %>% filter( Site == site) %>% na.omit
  
  SITEval_DATA_FILTERED_RSHPc[[site]] <- SITEval_DATA_FILTEREDc[[site]] %>% 
    full_join(CCC_VAL, by=c('gas', 'Approach', 'dLevelsAminusB' )) %>% 
    filter(CCC > 0) %>%  
    mutate(month = format(timeEndA.local,'%m') %>% as.numeric,
           Season = case_when(
             month %in% c(12, 1, 2) ~ "Winter",
             month %in% c(3, 4, 5) ~ "Spring",
             month %in% c(6, 7, 8) ~ "Summer",
             TRUE ~ "Autumn"),
           hour = format(timeEndA.local,'%H'))
  
}

#Save files: 
  save( SITEval_DATA_FILTERED_RSHPc, SITEval_DATA_FILTERED_RSHP,
        file= fs::path(localdir,paste0("SITEval_DATA_FILTERED_RSHP.Rdata")))

# Harmonize data and test result:  ####

# Explore 2 options: CCC > 0 and RSHP
harmonize_val <- function( tibble){
  
  site.list <- names(tibble )
  harmonized <- list()
  for( i in site.list){
    print(i)
    
    df <- tibble[[i]]
    
    harmonized[[i]] <-   df %>% 
      mutate( time.rounded = round_date(timeEndA.local, unit = "30 minutes") ) %>% 
      reframe(.by= c(time.rounded, gas),
                                        FG_mean = median(FG_mean, na.rm=T),
                                        EC_mean = mean(EC_mean, na.rm=T),
                                        Tair_C = mean(Tair_C, na.rm=T),
                                        PAR = mean(PAR, na.rm=T))  %>% 
      mutate( Month = format(time.rounded , '%m') %>% as.numeric) %>% mutate(
        Season = case_when(
          Month %in% c(12, 1, 2) ~ "Winter",
          Month %in% 3:5 ~ "Spring",
          Month %in% 6:8 ~ "Summer",
          Month %in% 9:11 ~ "Autumn")) 
    
  }
  return(harmonized )
  
}

# RSHP model
SITEval_DATA_FILTERED_RSHP_H <- harmonize_val(tibble = SITEval_DATA_FILTERED_RSHP)

# CCC > 0
SITEval_DATA_FILTERED_RSHPc_H <- harmonize_val(tibble = SITEval_DATA_FILTERED_RSHPc)

#Save files: 
save( SITEval_DATA_FILTERED_RSHP_H, SITEval_DATA_FILTERED_RSHPc_H,
      file= fs::path(localdir,paste0("SITEval_DATA_FILTERED_RSHP_EnSEMBLE.Rdata")))

# CCC for Harmonized data:

harmonized_CCC <- function( tibble){
  
  CCCval_Harmonized <- data.frame()
  
  for( i in names( tibble)){
    print(i)
    
    seasons <- tibble[[i]] $Season %>% unique
    
    for(s in seasons){
      
      data <- tibble[[i]] %>% filter( Season == s) 
      GAS <- data$gas %>% unique
      
      for( g in GAS ){
        
        df.summary <-linear_terms2('FG_mean', 'EC_mean' , df = tibble[[i]] %>% filter(Season == s, gas == g), site=i) %>% mutate( Season = s, gas = g)
        
        CCCval_Harmonized  <- rbind(CCCval_Harmonized, df.summary)
      }
    }}
  
  return(CCCval_Harmonized )
}


CCC_RSHP_H <- harmonized_CCC(tibble=SITEval_DATA_FILTERED_RSHP_H )
CCC_RSHPc_H <- harmonized_CCC(tibble= SITEval_DATA_FILTERED_RSHPc_H )

# Data VIZ:
CCC_RSHP_H %>% ggplot() + geom_point( aes( x= CCC, y = Site, col=Season)) + facet_wrap( ~gas)

CCC_RSHPc_H %>% ggplot() + geom_point( aes( x= CCC, y = Site, col=Season)) + facet_wrap( ~gas)



