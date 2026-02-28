# flow.evaluation.batch.ch4

rm(list=ls())

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(sf)

# -------------- Change this stuff -------------
#DirRepo <- 'C:/Users/csturtevant/Documents/Git/lterwg-flux-gradient' # Relative or absolute path to lterwg-flux-gradient git repo on your local machine. Make sure you've pulled the latest from main!
DirRepo.ch4 <-"/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-methane"

DirRepo.eval <-"/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-eval"

setwd(DirRepo.ch4)
#localdir <- 'C:/Users/csturtevant/OneDrive - Battelle Ecology/FluxGradient/filterTesting' # We'll deposit output files here prior to uploading to Google Drive

localdir <- '/Volumes/MaloneLab/Research/FluxGradient/FluxData'
DnldFromGoogleDrive <- FALSE # Enter TRUE to grab files listed in dnld_files from Google Drive. Enter FALSE if you have the most up-to-date versions locally in localdir

email <- 'sparklelmalone@gmail.com'
googledrive::drive_auth(email = TRUE) 
drive_url <- googledrive::as_id("https://drive.google.com/drive/folders/1Q99CT77DnqMl2mrUtuikcY47BFpckKw3") # The Data 
data_folder <- googledrive::drive_ls(path = drive_url)

metadata <- read.csv('/Volumes/MaloneLab/Research/FluxGradient/Ameriflux_NEON field-sites.csv') # has a list of all the sites

# -------------------------------------------------------
# Application of Filter Functions & Compiles Data frames into one list: : #### 
message('Running Filter...')
source(fs::path(DirRepo.ch4,'workflows/flow.filter.NEON.R'))

fileSave <- fs::path(localdir,paste0("SITE_DATA_FILTERED_CH4.RDS"))
saveRDS( SITE_DATA_FILTERED,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir,paste0("SITES_MBR_9min_FILTER_CH4.Rdata"))
save( SITES_MBR_9min_FILTER,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir,paste0("SITES_AE_9min_FILTER_CH4.Rdata"))
save( SITES_AE_9min_FILTER ,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir,paste0("SITES_WP_9min_FILTER_CH4.Rdata"))
save( SITES_WP_9min_FILTER,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

message('Running Filter...')
source(fs::path(DirRepo.eval,'workflows/flow.filter.validation.R'))


fileSave <- fs::path(localdir,paste0("SITEval_DATA_FILTERED_CH4.Rdata"))
save( SITEval_DATA_FILTERED,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir,paste0("SITESval_MBR_9min_FILTER_CH4.Rdata"))
save( SITESval_MBR_9min_FILTER,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir,paste0("SITESval_AE_9min_FILTER_CH4.Rdata"))
save( SITESval_AE_9min_FILTER ,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir,paste0("SITESval_WP_9min_FILTER_CH4.Rdata"))
save( SITESval_WP_9min_FILTER,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

# Counter Gradient (VAL): ####

# Gradient Flux Validation (VAL): ####

# Reliable sampling Height pairs (NEON): ####

# DIELS (NEON): ####

# Q10 FIT (NEON): ####

# Visualizations and Analysis (NEON): ####
