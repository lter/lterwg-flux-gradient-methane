#flow.filter.validation
# -------------- Change this stuff -------------
# Add all sites here:
metadata <- read.csv('/Volumes/MaloneLab/Research/FluxGradient/metadata_validation.csv') # has a list of all the sites

# Add local directory for downloaded data here:
DirRepo <-  "/Users/sm3466/Library/CloudStorage/Dropbox-YSE/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient"

setwd(DirRepo.ch4)
datadir <- '/Volumes/MaloneLab/Research/FluxGradient/FluxData' # MaloneLab Server

DnldFromGoogleDrive <- FALSE # Enter TRUE to grab files listed in dnld_files from Google Drive. Enter FALSE if you have the most up-to-date versions locally in localdir
email <- 'sparklelmalone@gmail.com'
googledrive::drive_auth(email = TRUE) 
drive_url <- googledrive::as_id("https://drive.google.com/drive/folders/1Q99CT77DnqMl2mrUtuikcY47BFpckKw3") # The Data 
data_folder <- googledrive::drive_ls(path = drive_url)

# Filter the data and calculate any needed flags:
library(lutz)
library(sf)

source(fs::path(DirRepo.eval ,'functions/calc.filter_FG.R' ))
source(fs::path(DirRepo.eval ,'functions/calc.SITELIST_FORMATTING.R'))

site.list <- c('SE-Sto', 'SE-Svb', 'US-Uaf')
# Creates the files for the filter report and Filters the data #####
for( site in site.list){
  
  print( site)
  
  # Load the files:
  localdir.site <- paste(datadir,"/", site, sep = "")
  
  files <- paste(site, "_Evaluation.Rdata", sep = "")
  
  load(paste(localdir.site, "/", files, sep=""))
  
  # Change the time to local:
  
  # Get NEON sites from the server and find the time zones: https://cran.r-project.org/web/packages/lutz/readme/README.html
  sites.location <- metadata %>%  st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
                                           crs = "+proj=longlat +datum=WGS84")
  
  sites.location$TZ <- tz_lookup(sites.location, method = "accurate")
  sites.location.sub <- sites.location %>%  select( "SITE_ID" , "TZ")
  
  site.tz <- sites.location.sub$TZ[which( sites.location.sub$SITE_ID == site)]
  
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
                                gas == "CH4"~ FCH4_turb_interp))
  
  AE_9min_FILTER <- filter_fluxes( df = AE_9min.df.final,
                                   dConcSNR.min = 3,
                                   approach = "AE")%>% 
    mutate( EC_mean = case_when(gas == "CO2" ~   FC_turb_interp,
                                gas == "H2O" ~ FH2O_interp,
                                gas == "CH4"~ FCH4_turb_interp))
  
  WP_9min_FILTER <- filter_fluxes ( df = WP_9min.df.final,
                                    dConcSNR.min = 3,
                                    approach = "WP")%>% 
    mutate( EC_mean = case_when(gas == "CO2" ~   FC_turb_interp,
                                gas == "H2O" ~ FH2O_interp,
                                gas == "CH4"~ FCH4_turb_interp))
  
  
  # Output the files
  localdir.site <- paste(datadir ,"/", site, sep = "")
  
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

source(fs::path(DirRepo.eval,'functions/calc.diel.R' ))

SITEval_DATA_FILTERED <- list() # Save all the data here:
SITESval_MBR_9min_FILTER <- list()
SITESval_AE_9min_FILTER <- list()
SITESval_WP_9min_FILTER <- list()

# Filtered Data into one object: ####
for( site in site.list){
  
  print(paste("Working on " , site))
  
  message( paste("Importing the data for ", site))
  
  localdir.site <- paste(datadir ,"/", site, sep = "")
  load(paste(localdir.site, "/", site, "_FILTER.Rdata", sep=""))
 
  MBR_9min_FILTER_canopy <-  MBR_9min_FILTER %>% mutate(Approach = "MBR") %>% 
    select( timeEndA.local,FG_mean , Approach, gas, 
            TowerPosition_A, TowerPosition_B,  EC_mean, cross_grad_flag, dLevelsAminusB,
            RH, PAR, Tair_K ) %>% mutate( Tair_C = Tair_K -273.15)
  
  WP_9min_FILTER_canopy <- WP_9min_FILTER %>% mutate( Approach = "WP") %>%
    select( timeEndA.local,FG_mean ,Approach, gas, 
            TowerPosition_A, TowerPosition_B,  EC_mean, cross_grad_flag, dLevelsAminusB,
            RH, PAR, Tair_K) %>% mutate( Tair_C = Tair_K -273.15)
  
  AE_9min_FILTER_canopy <- AE_9min_FILTER %>% mutate( Approach = "AE") %>%
    select( timeEndA.local,FG_mean ,Approach, gas, 
            TowerPosition_A, TowerPosition_B,  EC_mean, cross_grad_flag, dLevelsAminusB,
            RH, PAR, Tair_K) %>% mutate( Tair_C = Tair_K -273.15)
  
  Data <- rbind( MBR_9min_FILTER_canopy,  AE_9min_FILTER_canopy, WP_9min_FILTER_canopy) %>% as.data.frame() %>% 
    mutate(month = format(timeEndA.local,'%m') %>% as.numeric,
           Season = case_when(
             month %in% c(12, 1, 2) ~ "Winter",
             month %in% c(3, 4, 5) ~ "Spring",
             month %in% c(6, 7, 8) ~ "Summer",
             TRUE ~ "Autumn"),
           hour = format(timeEndA.local,'%H'))
  
  SITEval_DATA_FILTERED[[site]] <- Data
  
  SITESval_MBR_9min_FILTER[[site]] <- MBR_9min_FILTER_canopy
  SITESval_AE_9min_FILTER[[site]] <- AE_9min_FILTER_canopy
  SITESval_WP_9min_FILTER[[site]] <- WP_9min_FILTER_canopy
}

print("Create Site List by approach and reports")

# REPORTS: ####
filter.report.CH4 <- data.frame()
filter.report.stability.CH4 <- data.frame()
#  filter report
for( site in site.list){
  
  print( site)
  
  # Load the files:
  localdir.site <- paste(datadir ,"/", site, sep = "")
  
  files.CH4 <- paste(site, "_9min.report.CH4.csv", sep = "")
  
  report.CH4 <- read.csv( paste( localdir.site,"/", files.CH4, sep="" ))
  
  filter.report.CH4 <- rbind(filter.report.CH4,  report.CH4 )
  
  message("Done with ", site)
}
# Report with Stability
for( site in site.list){
  
  print( site)
  
  # Load the files:
  localdir.site <- paste(datadir ,"/", site, sep = "")
  
  files.CH4 <- paste(site, "_9min.report.stability.CH4.csv", sep = "")
  report.CH4 <- read.csv( paste( localdir.site,"/", files.CH4, sep="" ))
  filter.report.stability.CH4 <- rbind(filter.report.stability.CH4,  report.CH4 )
  
  message("Done with ", site)
}

# Save Data: 


