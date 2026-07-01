# Send Files to google:

library(ggpubr)
library(ggplot2)
library(colorspace)
library(sf)
library(AOI)
library(tidyverse)

# Import data and build the files needed: ####
localdir <- '/Volumes/MaloneLab/Research/FluxGradient/FluxData'

# Load parameter information:
load(file = fs::path(localdir,paste0("NEON_PARMS_DIEL_Q10.Rdata")))

# Load soil and root data
load(file=paste(localdir, "/", 'Soildata_YearMon.Rdata', sep=""))
# Load data obtained from Avni:
avni <- read.csv( paste("/Volumes/MaloneLab/Research/FluxGradient/Avni's Data/neondata_forSM.csv" , sep="" )) %>% rename( SITE_ID =site_code )

avni2 <- read.csv( paste("/Volumes/MaloneLab/Research/FluxGradient/Avni's Data/neondata_forsteve.csv" , sep="" )) %>% rename( SITE_ID =site_code ) %>% dplyr::select(
  SITE_ID, MAP, MAT, aridity.index,shallow_root, deep_root)

# Load site information:
soil.data.biogeo <- read.csv(file=paste(localdir, "/", 'Soil_Biogeochem_RootBiomass.csv', sep=""))
# Load canopy information: 
canopy <- read.csv(file.path(paste(localdir, "canopy_commbined.csv", sep="/"))) %>% distinct

canopy.Ht <- canopy %>% rename(SITE_ID = Site) %>% 
  reframe( .by= SITE_ID,
           canopyHeight_m = mean(canopyHeight_m))

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
                              Vegetation.Abbreviation..IGBP. == 'OSH' ~ 'Shrubland')) %>% 
  rename( site = Site_Id.NEON,
          SITE_ID = Site_Id.NEON) %>% 
  left_join(canopy.Ht , by = 'SITE_ID') %>% 
  left_join(soil.data.biogeo %>% rename(SITE_ID = siteID)  , by = 'SITE_ID')%>% 
  left_join(avni , by = 'SITE_ID')%>% 
  left_join(avni2 , by = 'SITE_ID')

# Create a simple feature
site.att.sf <- st_as_sf(x = metadata,                         
                        coords = c("Longitude..degrees.", "Latitude..degrees."), crs = 4326) 

#format data:
ENSEMBLE_DIELS_canopy <-  ENSEMBLE_DIELS %>% rename( SITE_ID = site) %>% left_join( metadata, by = "SITE_ID")

ENSEMBLE_Q10_eq4_canopy <-   ENSEMBLE_Q10_eq4 %>% 
  mutate( Date = as.Date(paste(idx, "-01", sep="")),
          month = format(Date, "%m") %>% as.numeric,
          Season = case_when(
            month %in% c(12, 1, 2) ~ "Winter",
            month %in% c(3, 4, 5) ~ "Spring",
            month %in% c(6, 7, 8) ~ "Summer",
            TRUE ~ "Autumn")) %>% left_join( metadata, by = "SITE_ID")

# This file has monthly data in it

parms <- ENSEMBLE_Q10_eq4_canopy %>% 
  dplyr::select(Q10.mean, SITE_ID, Date, Rref.mean, EcoType, Season, sulfurTot, dryMass, acidity, FRC.kg.C.m.2, SOC.kg.C.m.2, beta_roots, beta_soc, MAP, MAT, aridity.index,shallow_root, deep_root) %>% 
  left_join( ENSEMBLE_DIELS %>% 
               mutate(Date = as.Date(paste(YearMon, "-01", sep="")),
                      SITE_ID = site) %>% 
               reframe( .by = c(SITE_ID, Date), 
                        DIEL = sum(DIEL)), by=c('SITE_ID', 'Date')) %>% 
  left_join( Site_SoilData %>% 
               mutate(Date = as.Date(paste(YearMon, "-01", sep="")), SITE_ID = Site),by = c('SITE_ID', 'Date')) %>% 
  mutate( Source = case_when( DIEL > 0 ~ 1,
                              DIEL <= 0 ~ 0) %>% as.factor)

parms.site.season <- parms %>% reframe( .by = c(SITE_ID, Season),
                                        VSWCMean  = mean( VSWCMean, na.rm=TRUE), 
                                        Q10.mean = mean(Q10.mean, na.rm=TRUE) ,
                                        Rref.mean = mean(Rref.mean, na.rm=TRUE), 
                                        VSWCMin = mean(VSWCMin, na.rm=TRUE),
                                        VSWCVar= mean(VSWCVar, na.rm=TRUE),
                                        sulfurTot= mean(sulfurTot, na.rm=TRUE),
                                        dryMass= mean( dryMass, na.rm=TRUE),
                                        acidity= mean( acidity, na.rm=TRUE),
                                        DIEL = mean(DIEL, na.rm=T),
                                        FRC.kg.C.m.2 = mean(FRC.kg.C.m.2, na.rm=T), 
                                        SOC.kg.C.m.2 = mean(SOC.kg.C.m.2, na.rm=T),
                                        beta_roots= mean(beta_roots, na.rm=T),  
                                        beta_soc= mean(beta_soc, na.rm=T),
                                        MAT= mean(MAT, na.rm=T), 
                                        MAP= mean(MAP, na.rm=T))  %>% 
  mutate( Source = case_when( DIEL > 0 ~ "Source",
                              DIEL <= 0 ~ "Sink") %>% as.factor)

parms.site <- parms.site.season %>% mutate( Source = case_when( DIEL > 0 ~ 1,
                                                                DIEL <= 0 ~ 0) %>% as.numeric) %>% mutate(across(where(is.numeric), ~ ifelse(is.infinite(.), 0, .))) %>% 
  reframe( .by = SITE_ID, 
           Source = sum(Source ) %>% as.factor,
           VSWCMean  = mean( VSWCMean, na.rm=TRUE), 
           Q10.mean = mean(Q10.mean, na.rm=TRUE) ,
           Rref.mean = mean(Rref.mean, na.rm=TRUE), 
           VSWCMin = mean(VSWCMin, na.rm=TRUE),
           VSWCVar= mean(VSWCVar, na.rm=TRUE),
           sulfurTot= mean(sulfurTot, na.rm=TRUE),
           dryMass= mean( dryMass, na.rm=TRUE),
           acidity= mean( acidity, na.rm=TRUE),
           DIEL = mean(DIEL, na.rm=T),
           FRC.kg.C.m.2 = mean(FRC.kg.C.m.2, na.rm=T), 
           SOC.kg.C.m.2 = mean(SOC.kg.C.m.2, na.rm=T),
           beta_roots= mean(beta_roots, na.rm=T),  
           beta_soc= mean(beta_soc, na.rm=T),
           MAT= mean(MAT, na.rm=T), 
           MAP= mean(MAP, na.rm=T),
           Source_Sink = case_when( DIEL > 0 ~ "Source",
                                    DIEL <= 0 ~ "Sink") %>% as.factor)  %>% 
  left_join( metadata %>% dplyr::select(EcoType, SITE_ID), by = "SITE_ID")

parms.site.season.sf <-site.att.sf %>% left_join( parms.site.season , by="SITE_ID")

parms.site.sf <-site.att.sf %>% left_join( parms.site , by="SITE_ID")

parms.mds <- parms.site %>% drop_na()

library(vegan)

d <- vegdist(  parms.mds %>% dplyr::select( continious.variables))
fit <- metaMDS(d, k=2)

parms.mds$MDS1 <- fit$points[,1]
parms.mds$MDS2<- fit$points[,2]

cluster_res <- hclust(d, method = "ward.D2")

groups <- cutree(cluster_res, k = 4)
parms.mds$cluster <- groups # Add clusters to the df

site.att.sf.cluster <- site.att.sf %>% left_join(by='SITE_ID',   parms.mds )

DnldFromGoogleDrive <- FALSE # Enter TRUE to grab files listed in dnld_files from Google Drive. Enter FALSE if you have the most up-to-date versions locally in localdir

email <- 'sparklelmalone@gmail.com'
googledrive::drive_auth(email = TRUE) 
drive_url <- googledrive::as_id("https://drive.google.com/drive/folders/1Q99CT77DnqMl2mrUtuikcY47BFpckKw3") # The Data 
data_folder <- googledrive::drive_ls(path = drive_url)


fileSave <- fs::path(localdir,paste0("PARMS_CH4_FINAL.Rdata"))
save( parms.mds, parms.site.season.sf,
      parms.site.season,
      parms.site.sf,
      parms.site,
      file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)


load( fs::path(localdir,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata")) )
fileSave <- fs::path(localdir,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE.Rdata")) 
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)
