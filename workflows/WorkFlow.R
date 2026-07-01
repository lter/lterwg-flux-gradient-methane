# flow.evaluation.batch.ch4

# Runs the analysis for NEON files:

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(sf)

# -------------- Change this stuff -------------
#DirRepo <- 'C:/Users/csturtevant/Documents/Git/lterwg-flux-gradient' # Relative or absolute path to lterwg-flux-gradient git repo on your local machine. Make sure you've pulled the latest from main!
DirRepo.ch4 <-"/Users/sm3466/Dropbox-YSE/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-methane"

# Laptop
DirRepo.ch4 <-"/Users/sm3466/Library/CloudStorage/Dropbox-YSE/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-methane"
DirRepo.eval <-"/Users/sm3466/Dropbox-YSE/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-eval"

setwd(DirRepo.ch4)
localdir <- '/Volumes/MaloneLab/Research/FluxGradient'
localdir.ch4 <- '/Volumes/MaloneLab/Research/FluxGradient/METHANE'
DnldFromGoogleDrive <- FALSE # Enter TRUE to grab files listed in dnld_files from Google Drive. Enter FALSE if you have the most up-to-date versions locally in localdir

email <- 'sparklelmalone@gmail.com'
googledrive::drive_auth(email = TRUE) 
drive_url <- googledrive::as_id("https://drive.google.com/drive/folders/1Q99CT77DnqMl2mrUtuikcY47BFpckKw3") # The Data 
data_folder <- googledrive::drive_ls(path = drive_url)

metadata <- read.csv('/Volumes/MaloneLab/Research/FluxGradient/Ameriflux_NEON field-sites.csv') # has a list of all the sites

site.list <- metadata$Site_Id.NEON %>% unique
# -------------------------------------------------------
# Application of Filter Functions & Compiles Data frames into one list: : #### 
message('Running Filter...')
source(fs::path(DirRepo.ch4,'workflows/flow.filter.NEON.R'))

fileSave <- fs::path(localdir.ch4,paste0("/NEON_GradientFlux_Data_Filter/SITE_DATA_FILTERED_CH4.Rdata"))
save( SITE_DATA_FILTERED,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir.ch4 ,paste0("/NEON_GradientFlux_Data_Filter/SITES_MBR_9min_FILTER_CH4.Rdata"))
save( SITES_MBR_9min_FILTER,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir.ch4 ,paste0("/NEON_GradientFlux_Data_Filter/SITES_AE_9min_FILTER_CH4.Rdata"))
save( SITES_AE_9min_FILTER ,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir.ch4 ,paste0("/NEON_GradientFlux_Data_Filter/SITES_WP_9min_FILTER_CH4.Rdata"))
save( SITES_WP_9min_FILTER,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

# RSHP:
source(fs::path(DirRepo.ch4,'workflows/flow.RSHP.NEON_ID.R'))

# Ensemble GF:
source(fs::path(DirRepo.ch4,'workflows/flow.ENSEMBLE.NEON.R'))

# Combine ensembled gradient flux with storage flux:
source(fs::path(DirRepo.ch4,'workflows/flow.TotalFlux.R'))

# load(fs::path(localdir.ch4 ,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata")))

# Summarize Data Availability here!!!

# DIELS (NEON): ####
source(fs::path(DirRepo.ch4,'workflows/flow.DIEL.NEON.R'))

# see which DIELS are better CO2 or H2O

#Explains why sites show specific patterns and how those patterns scale temporally-
source(fs::path(DirRepo.ch4,'workflows/flow.30min.analysis.R'))  # makes files needed below
#source(fs::path(DirRepo.ch4,'workflows/NEON.30min.Gapfill.R')) 
source(fs::path(DirRepo.ch4,'workflows/NEON.ERA5.HalfHourlyGapfill.R')) # Gap-fills the data


source(fs::path(DirRepo.ch4,'workflows/flow.plots.R')) # produces figures based NEON Fluxes


source(fs::path(DirRepo.ch4,'workflows/NEON.FLUXNET.CH4FluxComparison.R')) # Ensure this is using the ERA 5 data for the comparisons!
## OUTPUT/CH4_flux_medians_by_source_and_behavior.csv 


# Attempts at upscale for a global budget comparison:
source(fs::path(DirRepo.ch4,'workflows/ERA-Upscaling.R')) # Not usig this version!

source(fs::path(DirRepo.ch4,'workflows/ERA-SpatialProbability.R')) 
source(fs::path(DirRepo.ch4,'workflows/Download-ERA5Land-Monthly.R')) 
source(fs::path(DirRepo.ch4,'workflows/Download-Process-MODIS-WAD2M.R')) 

# Random Forest Models 
source("/Users/sm3466/Library/CloudStorage/Dropbox-YSE/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-methane/workflows/ERA-SpatialUpscaling-Monthly-RF.R")
source(fs::path(DirRepo.ch4, '/workflows/ERA-SpatialUpscaling-Figures-RF.R'))

# Supplemental 
source(fs::path(DirRepo.ch4,'workflows/Flux_Justification_Supplement.R'))
source(fs::path(DirRepo.ch4,'/workflows/Supplemental_GF_Uncertainty.R'))
source(fs::path(DirRepo.ch4,'/workflows/Supplemental_Budget_Uncertainty.R'))

# Methods
source(fs::path(DirRepo.ch4,'workflows/flow.map.R'))

# Depreciated Files:
# Exploratory:
# GLM Models
source(fs::path(DirRepo.ch4,'workflows/NEON.DriveScale.Analysis.R')) # produces needed files

source(fs::path(DirRepo.ch4,'workflows/ERA-SpatialUpscaling-Monthly.R')) 
source(fs::path(DirRepo.ch4,'workflows/NEON.StrongSink.DriverComparison.R'))
source(fs::path(DirRepo.ch4,'workflows/NEON.DriveScale.Analysis.R')) # fix figure colors!
source(fs::path(DirRepo.ch4,'workflows/NEON.DIEL.Analysis2.R'))
source(fs::path(DirRepo.ch4,'workflows/NEON.MonthlySinkBehavior.Analysis.R'))
source(fs::path(DirRepo.ch4,'workflows/NEON.Supplementary.ModelDriverPlots.R'))
source(fs::path(DirRepo.ch4,'workflows/NEON.30min.Gapfill.R')) # I am not sure I need this file!
source(fs::path(DirRepo.ch4,'workflows/NEON.ConsistencyMagnitude.Analysis.R')) # fix figure colors!
source(fs::path(DirRepo.ch4,'workflows/NEON.TotalFlux.AnnualBudget.R'))
source(fs::path(DirRepo.ch4,'workflows/NEON.DIEL.SynthesisFigure.R'))
source(fs::path(DirRepo.ch4,'workflows/flow.climate.NEON.R'))
source(fs::path(DirRepo.ch4,'workflows/flow.soilmoisture.NEON.R'))
source(fs::path(DirRepo.ch4,'workflows/flow.parms_results.R'))


