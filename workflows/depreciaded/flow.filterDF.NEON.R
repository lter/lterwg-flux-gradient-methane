# Builds the Filtered dataframe: CH4 only ####
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(gtools)
library(ggplot2)
library(ggpmisc)

source(fs::path(DirRepo.eval,'functions/calc.diel.R' ))
metadata <- read.csv('/Volumes/MaloneLab/Research/FluxGradient/Ameriflux_NEON field-sites.csv') 

# Need to bring in new canopy file with all the canopy_L1 options
canopy.info <- read.csv(file.path(paste(localdir, "canopy_commbined.csv", sep="/"))) %>% distinct

# Data Prep:
approach.df <- data.frame( Approach=c("MBR", "AE", "WP"))

canopy <- canopy.info %>% cross_join(approach.df) %>% mutate(Approach = factor(Approach, levels = c("MBR", "AE", "WP") ),
                                                             RelativeDistB = MeasurementHeight_m_B - CanopyHeight, 
                                                             RelativeDistA = MeasurementHeight_m_A - CanopyHeight, 
                                                             MeasurementDist = MeasurementHeight_m_A - MeasurementHeight_m_A)

site.list <- metadata$Site_Id.NEON %>% unique

SITE_DATA_FILTERED <- list() # Save all the data here:

SITES_MBR_9min_FILTER <- list()
SITES_AE_9min_FILTER <- list()
SITES_WP_9min_FILTER <- list()

# Filtered Data into one object: ####
for( site in site.list){
  
  print(paste("Working on " , site))
  
  message( paste("Importing the data for ", site))
  
  localdir.site <- paste(localdir,"/", site, sep = "")
  load(paste(localdir.site, "/", site, "_FILTER.Rdata", sep=""))
  
  canopy.sub <- canopy %>% filter( Site == site) %>% select(Site, Canopy_L1, dLevelsAminusB, Approach)
  

  MBR_9min_FILTER_canopy <- canopy.sub %>% filter( Approach == "MBR") %>% 
    full_join( MBR_9min_FILTER , by = c('dLevelsAminusB'), relationship = "many-to-many")  
  %>% select( timeEndA.local,FG_mean , Approach, Canopy_L1, gas, 
            TowerPosition_A, TowerPosition_B,  EC_mean, cross_grad_flag, dLevelsAminusB,
            RH, PAR, Tair_K) %>% mutate( Tair_C = Tair_K -273.15)
  
  WP_9min_FILTER_canopy <- canopy.sub %>% filter( Approach == "WP") %>% 
    full_join( WP_9min_FILTER , by = c('dLevelsAminusB'), relationship = "many-to-many" )  %>% 
    select( timeEndA.local,FG_mean ,Approach, Canopy_L1, gas, 
            TowerPosition_A, TowerPosition_B,  EC_mean, cross_grad_flag, dLevelsAminusB,
            RH, PAR, Tair_K) %>% mutate( Tair_C = Tair_K -273.15)
  
  AE_9min_FILTER_canopy <- canopy.sub %>% 
    filter( Approach == "AE") %>% 
    full_join( AE_9min_FILTER , by = c('dLevelsAminusB'), 
               relationship = "many-to-many")  %>% 
    select( timeEndA.local,FG_mean, Approach, Canopy_L1, gas, 
            TowerPosition_A, TowerPosition_B,  EC_mean, cross_grad_flag, dLevelsAminusB,
            RH, PAR, Tair_K) %>% mutate( Tair_C = Tair_K -273.15)
  
  Data <- rbind( MBR_9min_FILTER_canopy,  AE_9min_FILTER_canopy, WP_9min_FILTER_canopy) %>% as.data.frame() 
  
  SITE_DATA_FILTERED[[site]] <- Data
  
  SITES_MBR_9min_FILTER[[site]] <- MBR_9min_FILTER_canopy
  SITES_AE_9min_FILTER[[site]] <- AE_9min_FILTER_canopy
  SITES_WP_9min_FILTER[[site]] <- WP_9min_FILTER_canopy
}

print("Create Site List by approach and reports")

# REPORTS: ####
filter.report.CH4 <- data.frame()
filter.report.stability.CH4 <- data.frame()
#  filter report
for( site in site.list){
  
  print( site)
  
  # Load the files:
  localdir.site <- paste(localdir,"/", site, sep = "")
  
  files.CH4 <- paste(site, "_9min.report.CH4.csv", sep = "")
  
  report.CH4 <- read.csv( paste( localdir.site,"/", files.CH4, sep="" ))
  
  filter.report.CH4 <- rbind(filter.report.CH4,  report.CH4 )
  
  message("Done with ", site)
}
# Report with Stability
for( site in site.list){
  
  print( site)
  
  # Load the files:
  localdir.site <- paste(localdir,"/", site, sep = "")
  
  files.CH4 <- paste(site, "_9min.report.stability.CH4.csv", sep = "")
  report.CH4 <- read.csv( paste( localdir.site,"/", files.CH4, sep="" ))
  filter.report.stability.CH4 <- rbind(filter.report.stability.CH4,  report.CH4 )
  
  message("Done with ", site)
}