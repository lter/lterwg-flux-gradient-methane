DirRepo.ch4 <-"/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-methane"

DirRepo.eval <-"/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-eval"

setwd(DirRepo.ch4)
localdir <- '/Volumes/MaloneLab/Research/FluxGradient'
localdir.ch4 <- '/Volumes/MaloneLab/Research/FluxGradient/Validation_Sites'

DnldFromGoogleDrive <- FALSE # Enter TRUE to grab files listed in dnld_files from Google Drive. Enter FALSE if you have the most up-to-date versions locally in localdir

email <- 'sparklelmalone@gmail.com'
googledrive::drive_auth(email = TRUE) 
drive_url <- googledrive::as_id("https://drive.google.com/drive/folders/1Q99CT77DnqMl2mrUtuikcY47BFpckKw3") # The Data 
data_folder <- googledrive::drive_ls(path = drive_url)


site.list <- c('SE-Sto', 'SE-Svb', 'US-Uaf')


message('Running Filter...')
source(fs::path(DirRepo.ch4,'workflows/Validation/flow.filter.validation.R'))

fileSave <- fs::path(localdir.ch4, paste0("SITEval_DATA_FILTERED_CH4.Rdata"))
save( SITEval_DATA_FILTERED,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)


fileSave <- fs::path(localdir.ch4 ,paste0("SITESval_MBR_9min_FILTER_CH4.Rdata"))
save( SITESval_MBR_9min_FILTER,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir.ch4, paste0("SITESval_AE_9min_FILTER_CH4.Rdata"))
save( SITESval_AE_9min_FILTER ,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

fileSave <- fs::path(localdir.ch4, paste0("SITESval_WP_9min_FILTER_CH4.Rdata"))
save( SITESval_WP_9min_FILTER,file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

source(fs::path(DirRepo.ch4,'workflows/flow.map.R')) # Makes Map but also create a canopy file needed below.

source(fs::path(DirRepo.ch4,'workflows/Validation/flow.RSHP_VAL.R'))

source(fs::path(DirRepo.ch4,'workflows/Validation/flow.DIEL.VAL.R'))

source(fs::path(DirRepo.ch4,'workflows/Validation/flow.Val.plots.R'))