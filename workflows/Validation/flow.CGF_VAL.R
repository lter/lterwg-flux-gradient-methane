# Put the counter gradient with the CH4 and explore use as a filter: 
library( ggplot2)
source(fs::path(DirRepo.eval, 'functions/calc.One2One.CCC_testing.R'))

fileSave <- fs::path(localdir,paste0("SITEval_DATA_FILTERED_CH4.Rdata"))
load( fileSave )

site.list <- names(SITEval_DATA_FILTERED)

SITEval_DATA_FILTERED_Final <- list()

# Extracts CGF for H2O and merges it by timestamp to the other gases! ####
for(site in site.list){
  
  print(site)
  
  SITEval_DATA_FILTERED_H2O <- SITEval_DATA_FILTERED[[site]] %>% filter( gas == "H2O") %>% select(timeEndA.local, cross_grad_flag) %>% rename(cross_grad_flag_H2O = cross_grad_flag )  %>% 
    mutate( Month = format(timeEndA.local, '%m') %>% as.numeric) %>% mutate(
      Season = case_when(
        Month %in% c(12, 1, 2) ~ "Winter",
        Month %in% 3:5 ~ "Spring",
        Month %in% 6:8 ~ "Summer",
        Month %in% 9:11 ~ "Autumn"))
  
  
  SITEval_DATA_FILTERED_Final[[site]] <- SITEval_DATA_FILTERED[[site]] %>% left_join( SITEval_DATA_FILTERED_H2O, by= 'timeEndA.local')
  
  rm(SITEval_DATA_FILTERED_H2O)

  print("complete")
  
}

# Filters by CGF ####
SITEval_DATA_FILTERED_Final_CGF <- list()
for(site in site.list){
  
  print(site)
  
  SITEval_DATA_FILTERED_H2O <- SITEval_DATA_FILTERED[[site]] %>% filter( gas == "H2O") %>% select(timeEndA.local, cross_grad_flag) %>% rename(cross_grad_flag_H2O = cross_grad_flag )
  
  
  SITEval_DATA_FILTERED_Final_CGF[[site]] <- SITEval_DATA_FILTERED[[site]] %>% left_join( SITEval_DATA_FILTERED_H2O, by= 'timeEndA.local') %>% filter(cross_grad_flag_H2O != 1) %>% 
    mutate( Month = format(timeEndA.local, '%m') %>% as.numeric) %>% mutate(
      Season = case_when(
        Month %in% c(12, 1, 2) ~ "Winter",
        Month %in% 3:5 ~ "Spring",
        Month %in% 6:8 ~ "Summer",
        Month %in% 9:11 ~ "Autumn"))
  
  rm(SITEval_DATA_FILTERED_H2O)
  
  print("complete")
  
}

fileSave <- fs::path(localdir,paste0("SITEval_DATA_FILTERED_Final_CH4.Rdata"))
load( fileSave )

fileSave <- fs::path(localdir,paste0("SITEval_DATA_FILTERED_Final_CH4.Rdata"))
save( SITEval_DATA_FILTERED_Final,
      SITEval_DATA_FILTERED_Final_CGF,
      file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)


