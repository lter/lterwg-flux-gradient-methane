# Filter the data and calculate any needed flags:
library(lutz)
library(sf)

site.list <- metadata$Site_Id.NEON %>% unique

source(fs::path(DirRepo.eval,'functions/calc.filter_FG.R' ))
source(fs::path(DirRepo.eval,'functions/calc.SITELIST_FORMATTING.R'))

# Creates the files for the filter report and Filters the data #####
for( site in site.list){
  
  print( site)
  
  # Load the files:
  localdir.site <- paste(localdir,"/", site, sep = "")
  
  files <- paste(site, "_Evaluation.Rdata", sep = "")
  
  load(paste(localdir.site, "/", files, sep=""))
  
  # Change the time to local:
  
  # Get NEON sites from the server and find the time zones: https://cran.r-project.org/web/packages/lutz/readme/README.html
  sites.location <- metadata %>%  st_as_sf(coords = c("Longitude..degrees.", "Latitude..degrees."),
                                           crs = "+proj=longlat +datum=WGS84")
  
  sites.location$TZ <- tz_lookup(sites.location, method = "accurate")
  sites.location.sub <- sites.location %>%  select( "Site_Id.NEON" , "TZ")
  
  site.tz <- sites.location.sub$TZ[which( sites.location.sub$Site_Id.NEON == site)]
 
  MBR_9min.df.final$timeEndA.local <- MBR_9min.df.final$timeEndA %>%  as.POSIXlt( tz = site.tz)
  AE_9min.df.final$timeEndA.local <- AE_9min.df.final$timeEnd_A %>%  as.POSIXlt( tz = site.tz)
  WP_9min.df.final$timeEndA.local <- WP_9min.df.final$timeEnd_A %>%  as.POSIXlt( tz = site.tz)
  
  MBR_9min.df.final$Month.local <- MBR_9min.df.final$timeEndA.local %>% format("%m")
  MBR_9min.df.final$Month.local <- MBR_9min.df.final$timeEndA.local %>% format("%m")
  MBR_9min.df.final$Month.local <- MBR_9min.df.final$timeEndA.local %>% format("%m")
  
  MBR_9min.df.final$time.local <- MBR_9min.df.final$timeEndA.local %>% format("%H:%M")
  AE_9min.df.final$time.local <- AE_9min.df.final$timeEndA.local %>% format("%H:%M")
  WP_9min.df.final$time.local <- WP_9min.df.final$timeEndA.local %>% format("%H:%M")
  
  MBR_9min.df.final$hour.local <- MBR_9min.df.final$timeEndA.local %>% format("%H")
  AE_9min.df.final$hour.local <- AE_9min.df.final$timeEndA.local %>% format("%H")
  WP_9min.df.final$hour.local <- WP_9min.df.final$timeEndA.local %>% format("%H")
  
  # Run the filter functions... Report:
  # CH4
  WP_9min.report.CH4 <-  filter_report(df = WP_9min.df.final %>% filter(gas == "CH4"),
                                       dConcSNR.min = 3,
                                       approach = "WP") %>% mutate(approach = "WP",
                                                                   site =site)
  
  AE_9min.report.CH4  <-  filter_report( df = AE_9min.df.final %>% filter(gas == "CH4"),
                                         dConcSNR.min = 3,
                                         approach = "AE") %>% mutate(approach = "AE",
                                                                     site =site)
  
  MBR_9min.report.CH4  <-  filter_report( df = MBR_9min.df.final %>% filter(gas == "CH4"),
                                          dConcSNR.min = 3,
                                          approach = "MBR") %>% mutate(approach = "MBR",
                                                                       site =site)
  
  
  WP_9min.report.stability.CH4 <-  filter_report_stability(df = WP_9min.df.final %>% filter(gas == "CH4"),
                                                           dConcSNR.min = 3,
                                                           approach = "WP") %>% mutate(approach = "WP",
                                                                                       site =site)
  
  AE_9min.report.stability.CH4 <-  filter_report_stability( df = AE_9min.df.final%>% filter(gas == "CH4"),
                                                            dConcSNR.min = 3,
                                                            approach = "AE") %>% mutate(approach = "AE",
                                                                                        site =site)
  
  MBR_9min.report.stability.CH4 <-  filter_report_stability( df = MBR_9min.df.final%>% filter(gas == "CH4"),
                                                             dConcSNR.min = 3,
                                                             approach = "MBR") %>% mutate(approach = "MBR",
                                                                                          site =site)
  
  SITE_9min.report.CH4 <- rbind( WP_9min.report.CH4, AE_9min.report.CH4, MBR_9min.report.CH4)
  SITE_9min.report.stability.CH4 <- rbind( WP_9min.report.stability.CH4, AE_9min.report.stability.CH4, MBR_9min.report.stability.CH4)
  
  # CO2
  WP_9min.report.CO2 <-  filter_report(df = WP_9min.df.final %>% filter(gas == "CO2"),
                                       dConcSNR.min = 3,
                                       approach = "WP") %>% mutate(approach = "WP",
                                                                   site =site)
  
  AE_9min.report.CO2  <-  filter_report( df = AE_9min.df.final %>% filter(gas == "CO2"),
                                         dConcSNR.min = 3,
                                         approach = "AE") %>% mutate(approach = "AE",
                                                                     site =site)
  
  MBR_9min.report.CO2  <-  filter_report( df = MBR_9min.df.final %>% filter(gas == "CO2"),
                                          dConcSNR.min = 3,
                                          approach = "MBR") %>% mutate(approach = "MBR",
                                                                       site =site)
  
  
  WP_9min.report.stability.CO2 <-  filter_report_stability(df = WP_9min.df.final %>% filter(gas == "CO2"),
                                                           dConcSNR.min = 3,
                                                           approach = "WP") %>% mutate(approach = "WP",
                                                                                       site =site)
  
  AE_9min.report.stability.CO2 <-  filter_report_stability( df = AE_9min.df.final%>% filter(gas == "CO2"),
                                                            dConcSNR.min = 3,
                                                            approach = "AE") %>% mutate(approach = "AE",
                                                                                        site =site)
  
  MBR_9min.report.stability.CO2 <-  filter_report_stability( df = MBR_9min.df.final%>% filter(gas == "CO2"),
                                                             dConcSNR.min = 3,
                                                             approach = "MBR") %>% mutate(approach = "MBR",
                                                                                          site =site)
  
  SITE_9min.report.CO2 <- rbind( WP_9min.report.CO2, AE_9min.report.CO2, MBR_9min.report.CO2)
  SITE_9min.report.stability.CO2 <- rbind( WP_9min.report.stability.CO2, AE_9min.report.stability.CO2, MBR_9min.report.stability.CO2)
  
  # H2O
  
  WP_9min.report.H2O <-  filter_report(df = WP_9min.df.final %>% filter(gas == "H2O"),
                                       dConcSNR.min = 3,
                                       approach = "WP") %>% mutate(approach = "WP",
                                                                   site =site)
  
  AE_9min.report.H2O  <-  filter_report( df = AE_9min.df.final %>% filter(gas == "H2O"),
                                         dConcSNR.min = 3,
                                         approach = "AE") %>% mutate(approach = "AE",
                                                                     site =site)
  
  MBR_9min.report.H2O  <-  filter_report( df = MBR_9min.df.final %>% filter(gas == "H2O"),
                                          dConcSNR.min = 3,
                                          approach = "MBR") %>% mutate(approach = "MBR",
                                                                       site =site)
  
  
  WP_9min.report.stability.H2O <-  filter_report_stability(df = WP_9min.df.final%>% filter(gas == "H2O"),
                                                           dConcSNR.min = 3,
                                                           approach = "WP") %>% mutate(approach = "WP",
                                                                                       site =site)
  
  AE_9min.report.stability.H2O <-  filter_report_stability( df = AE_9min.df.final%>% filter(gas == "H2O"),
                                                            dConcSNR.min = 3,
                                                            approach = "AE") %>% mutate(approach = "AE",
                                                                                        site =site)
  
  MBR_9min.report.stability.H2O <-  filter_report_stability( df = MBR_9min.df.final%>% filter(gas == "H2O"),
                                                             dConcSNR.min = 3,
                                                             approach = "MBR") %>% mutate(approach = "MBR",
                                                                                          site =site)
  
  SITE_9min.report.H2O <- rbind( WP_9min.report.H2O, AE_9min.report.H2O, MBR_9min.report.H2O)
  SITE_9min.report.stability.H2O <- rbind( WP_9min.report.stability.H2O, AE_9min.report.stability.H2O, MBR_9min.report.stability.H2O)
  

  # Run the filter functions... FILTER data and adjust the GF sign:

  MBR_9min_FILTER <- filter_fluxes( df = MBR_9min.df.final,
                                    dConcSNR.min = 3,
                                    approach = "MBR") %>% 
    mutate( EC_mean = case_when(gas == "CO2" ~  FC_turb_interp,
                                gas == "H2O" ~ FH2O_interp,
                                gas == "CH4"~ NA))
  
  AE_9min_FILTER <- filter_fluxes( df = AE_9min.df.final,
                                   dConcSNR.min = 3,
                                   approach = "AE")%>% 
    mutate( EC_mean = case_when(gas == "CO2" ~   FC_turb_interp,
                                gas == "H2O" ~ FH2O_interp,
                                gas == "CH4"~ NA))
  
  WP_9min_FILTER <- filter_fluxes ( df = WP_9min.df.final,
                                    dConcSNR.min = 3,
                                    approach = "WP")%>% 
    mutate( EC_mean = case_when(gas == "CO2" ~   FC_turb_interp,
                                gas == "H2O" ~ FH2O_interp,
                                gas == "CH4"~ NA))
  
  
  # Output the files
  localdir.site <- paste(localdir,"/", site, sep = "")
  
  write.csv( SITE_9min.report.stability.CO2,  paste(localdir.site, "/", site,"_9min.report.stability.CO2.csv", sep=""))
  write.csv( SITE_9min.report.CO2,  paste(localdir.site, "/", site,"_9min.report.CO2.csv", sep=""))
  
  write.csv( SITE_9min.report.stability.CH4,  paste(localdir.site, "/", site,"_9min.report.stability.CH4.csv", sep=""))
  write.csv( SITE_9min.report.CH4,  paste(localdir.site, "/", site,"_9min.report.CH4.csv", sep=""))
  
  write.csv( SITE_9min.report.stability.H2O,  paste(localdir.site, "/", site,"_9min.report.stability.H2O.csv", sep=""))
  write.csv( SITE_9min.report.H2O,  paste(localdir.site, "/", site,"_9min.report.H2O.csv", sep=""))
  
  
  save( AE_9min_FILTER,
        WP_9min_FILTER,
        MBR_9min_FILTER, file = paste(localdir.site, "/", site, "_FILTER.Rdata", sep=""))
  
  
  # Upload files to the google
  
  site_folder <-data_folder$id[data_folder$name==site]
  
  fileSave <- paste(localdir.site, "/", site, "_FILTER.Rdata", sep="")
  googledrive::drive_upload(media = fileSave, overwrite = T, path = site_folder)
  
  fileSave <- paste(localdir.site, "/", site,"_9min.report.CO2.csv", sep="")
  googledrive::drive_upload(media = fileSave, overwrite = T, path =site_folder)
  
  fileSave <- paste(localdir.site, "/", site,"_9min.report.stability.CO2.csv", sep="")
  googledrive::drive_upload(media = fileSave, overwrite = T, path =site_folder)
  
  fileSave <- paste(localdir.site, "/", site,"_9min.report.H2O.csv", sep="")
  googledrive::drive_upload(media = fileSave, overwrite = T, path =site_folder)
  
  fileSave <- paste(localdir.site, "/", site,"_9min.report.stability.H2O.csv", sep="")
  googledrive::drive_upload(media = fileSave, overwrite = T, path =site_folder)
  
  fileSave <- paste(localdir.site, "/", site,"_9min.report.CH4.csv", sep="")
  googledrive::drive_upload(media = fileSave, overwrite = T, path =site_folder)
  
  fileSave <- paste(localdir.site, "/", site,"_9min.report.stability.CH4.csv", sep="")
  googledrive::drive_upload(media = fileSave, overwrite = T, path =site_folder)
  
  message( paste("Done with filtering", site))
  
  rm( AE_9min_FILTER,
      WP_9min_FILTER,
      MBR_9min_FILTER,
      SITE_9min.report.CO2, 
      WP_9min.report.CO2, 
      AE_9min.report.CO2, 
      MBR_9min.report.CO2,
      SITE_9min.report.stability.CO2, 
      WP_9min.report.stability.CO2, 
      AE_9min.report.stability.CO2, 
      MBR_9min.report.stability.CO2,
      SITE_9min.report.H2O, 
      WP_9min.report.H2O, 
      AE_9min.report.H2O, 
      MBR_9min.report.H2O,
      SITE_9min.report.stability.H2O, 
      WP_9min.report.stability.H2O, 
      AE_9min.report.stability.H2O, 
      MBR_9min.report.stability.H2O,
      SITE_9min.report.H2O,
      sample.diel, sample.month,
      SITE_9min.report.CH4, 
      WP_9min.report.CH4, 
      AE_9min.report.CH4, 
      MBR_9min.report.CH4,
      SITE_9min.report.stability.CH4, 
      WP_9min.report.stability.CH4, 
      AE_9min.report.stability.CH4, 
      MBR_9min.report.stability.CH4 )
}

# Builds the Filtered dataframe: CH4 only ####
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(gtools)
library(ggplot2)
library(ggpmisc)


site.list <- metadata$Site_Id.NEON %>% unique

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
    full_join( MBR_9min_FILTER , by = c('dLevelsAminusB'), relationship = "many-to-many")  %>% 
    select( timeEndA.local,FG_mean , Approach, Canopy_L1, gas, 
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

# RSHP: ####
source(fs::path(DirRepo.ch4,'functions/calc_CCC_SHP.R' ))
load( fs::path(localdir,paste0("SITE_DATA_FILTERED_CH4.Rdata")))

NEON_RSHP <- CCC_SamplingHeightPairs(DATA = SITE_DATA_FILTERED) 

NEON_RSHP %>% filter(gas == "CH4") %>%  ggplot( aes( x= CCC, y = Site)) + geom_boxplot() + theme_bw()  + xlim(-1, 1) 

# select the good pairs
file1 <-NEON_RSHP %>% 
  mutate(Combination = paste(var1, var2, sep ="-")) %>%
  rename( var= var1)%>% select(Site, gas, var, CCC )
  
  
 file2 <- NEON_RSHP %>% 
          mutate(Combination = paste(var1, var2, sep ="-")) %>% 
          rename( var= var2) %>% select(Site, gas, var, CCC )

 NEON_RSHP_EVAL <- rbind( file1, file2) %>% 
  reframe( .by = c( Site, gas, var),
           CCC.GF.max = max( CCC, na.rm=T), 
           CCC.GF.min = min( CCC, na.rm=T), 
           CCC.GF.mean = mean(CCC, na.rm=T), 
           CCC.GF.median = median( CCC, na.rm=T),
           CCC.GF.range = mean( CCC, na.rm=T)) %>% 
  mutate(Good.CCC.GF = case_when(
    CCC.GF.max >= 0.95 ~ 1,.default = 0) %>% as.factor) 

load(file= paste(localdir, 'SITE_RSHP_MODEL.Rdata', sep="") ) # The model is in here:

NEON_RSHP_EVAL$RSHP.rf <- predict(rf.good.ccc.ec, NEON_RSHP_EVAL)
NEON_RSHP_EVAL$RSHP.rf %>% summary
NEON_RSHP_EVAL$Good.CCC.GF %>% summary

library(caret)

confusionMatrix(NEON_RSHP_EVAL$RSHP.rf, NEON_RSHP_EVAL$Good.CCC.GF)

NEON_RSHP_EVAL %>% filter(RSHP.rf == 1) %>% ggplot() + geom_boxplot( aes(x=CCC.GF.max, y = Site)) + 
  facet_wrap(~gas) + theme_bw() + xlim(-1,1)

NEON_RSHP_EVAL %>% filter(Good.CCC.GF == 1) %>% ggplot() + geom_boxplot( aes(x=CCC.GF.max, y = Site)) + 
  facet_wrap(~gas) + theme_bw() + xlim(-1,1)

fileSave <- fs::path(localdir,paste0("NEON_RSHP_CCC.Rdata"))
save( NEON_RSHP, NEON_RSHP_EVAL, file=fileSave)
