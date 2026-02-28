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

# RSHP ####
fileSave <- fs::path(localdir,paste0("SITEval_DATA_FILTERED_Final_CH4.Rdata"))
load( fileSave )

# This is for the CGF file: ####
VAL_CGF_RSHP <- CCC_SamplingHeightPairs(DATA = SITEval_DATA_FILTERED_Final_CGF) %>% mutate(Combination = paste(var1, var2, sep ="-"),
         CCC = as.numeric(CCC)) %>% 
  reframe( .by= c(Combination, Site, gas, var1), 
           CCC.1.max = max( CCC, na.rm=T), 
           CCC.1.min = min( CCC, na.rm=T)) %>% 
  rename( var= var1) %>% 
  rbind(VAL_CGF_RSHP  %>% 
          mutate(Combination = paste(var1, var2, sep ="-")) %>% 
          reframe( .by=c(Combination, Site, gas, var2), 
                   CCC.1.max = max(CCC, na.rm=T),
                   CCC.1.min = min(CCC, na.rm=T)) %>% 
          rename( var= var2) )  %>% 
  mutate(Good.CCC.GF = case_when(
    CCC.1.max >= 0.75 ~ 1,.default = 0) %>% as.factor) 

VAL_CGF_RSHP %>% filter(Good.CCC.GF == 1) %>% ggplot() + geom_boxplot( aes(x=CCC.1.max, y = Site)) + 
  facet_wrap(~gas) + theme_bw() + xlim(-1,1)


VAL_CGF_RSHP_sum <- VAL_RSHP_EVAL %>% select("Site","gas","var", 'Good.CCC.GF' ) %>% reframe( .by= c("Site","gas","var"), Good.CCC.GF = sum(Good.CCC.GF %>% as.numeric)) %>% separate(
  col = var, 
  into = c("Approach", "dLevelsAminusB"), 
  sep = "-")

# This is for the non CGF file: ##

RSHP <- CCC_SamplingHeightPairs(DATA = SITEval_DATA_FILTERED_Final)

VAL_RSHP <-RSHP %>% mutate(Combination = paste(var1, var2, sep ="-"),                                                   CCC = as.numeric(CCC)) %>% 
  reframe( .by= c(Combination, Site, gas, var1), 
           CCC.GF.max = max( CCC, na.rm=T), 
           CCC.GF.min = min( CCC, na.rm=T), 
           CCC.GF.mean = mean( CCC, na.rm=T), 
           CCC.GF.median = median( CCC, na.rm=T),
           CCC.GF.range = range( CCC, na.rm=T)) %>% 
  rename( var= var1) %>% 
  rbind(RSHP %>% 
          mutate(Combination = paste(var1, var2, sep ="-")) %>% 
          reframe( .by=c(Combination, Site, gas, var2), 
                   CCC.GF.max = max( CCC, na.rm=T), 
                   CCC.GF.min = min( CCC, na.rm=T), 
                   CCC.GF.mean = mean( CCC, na.rm=T), 
                   CCC.GF.median = median( CCC, na.rm=T),
                   CCC.GF.range = range( CCC, na.rm=T)) %>% 
          rename( var= var2) )  %>% 
  mutate(Good.CCC.GF = case_when(
    CCC.GF.max >= 0.75 ~ 1,.default = 0) %>% as.factor) 

# Save the file:
fileSave <- fs::path(localdir,paste0("CCC_CH4_RSHP.Rdata"))
save( VAL_RSHP, file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

# Use the RSHP model from eval:

load(file= paste(localdir, 'SITE_RSHP_MODEL.Rdata', sep="") ) # The model is in here:

#load(fs::path(localdir,paste0("CCC_CH4_RSHP.Rdata")))# The data on the RSHP for CH4 is in there
VAL_RSHP$RSHP.rf <- predict(rf.good.ccc.ec, VAL_RSHP)

VAL_RSHP$RSHP.rf %>% summary
VAL_RSHP$Good.CCC.GF %>% summary

caret::confusionMatrix(VAL_RSHP$RSHP.rf, VAL_RSHP$Good.CCC.GF)


# This now uses RSHP.rf!
VAL_RSHP_sum <- VAL_RSHP %>% select("Site","gas","var", 'Good.CCC.GF', 'RSHP.rf' ) %>% reframe( .by= c("Site","gas","var"), Good.CCC.GF = sum(RSHP.rf %>% as.numeric)) %>% separate(
  col = var, 
  into = c("Approach", "dLevelsAminusB"), 
  sep = "-")
site.list <- SITEval_DATA_FILTERED_Final_CGF %>% names()

# Now 
SITEval_DATA_FILTERED_Final_CGF_RSHP <- list()
SITEval_DATA_FILTERED_Final_RSHP <- list()

for( site in site.list){
  
  VAL_CGF_RSHP_sum_site <- VAL_CGF_RSHP_sum %>% filter( Site == site)
  VAL_RSHP_sum_site <- VAL_RSHP_sum %>% filter( Site == site)
  
  SITEval_DATA_FILTERED_Final_RSHP[[site]] <- SITEval_DATA_FILTERED_Final[[site]] %>% full_join(VAL_RSHP_sum_site, by=c('gas', 'Approach', 'dLevelsAminusB' )) %>% filter(Good.CCC.GF >1 )
  
  SITEval_DATA_FILTERED_Final_CGF_RSHP[[site]] <- SITEval_DATA_FILTERED_Final_CGF[[site]] %>% full_join(VAL_CGF_RSHP_sum_site, by=c('gas', 'Approach', 'dLevelsAminusB' )) %>% filter(Good.CCC.GF >1 )
  
  SITEval_DATA_FILTERED_Final_RSHP[[site]] <- SITEval_DATA_FILTERED_Final[[site]] %>% full_join(VAL_RSHP_sum_site, by=c('gas', 'Approach', 'dLevelsAminusB' )) %>% filter(Good.CCC.GF >1 )
  
}

SITEval_DATA_FILTERED_Final_CGF_RSHP # filtered for CGF and RSHP
SITEval_DATA_FILTERED_Final_RSHP # filtered for RSHP

# Compare FG to EC: ####
source(fs::path(DirRepo.eval, 'functions/calc.One2One.CCC_testing.R'))

CCC_RSHP <- ccc.val(sites.tibble =  SITEval_DATA_FILTERED_Final_RSHP)
CCC_CGF_RSHP <- ccc.val(sites.tibble =  SITEval_DATA_FILTERED_Final_CGF_RSHP) %>% rename(CCC_CGF = CCC, R2_CGF = R2) 


CCC_RSHP_FINAL <- CCC_RSHP %>% left_join( CCC_CGF_RSHP , by =c('Approach', 'dLevelsAminusB', 'gas',   'Site', 'Season') )

CCC_RSHP_FINAL %>% filter(gas == "CH4") %>%  ggplot( aes( x= CCC_CGF, y =Site, col=Approach)) + geom_point() + theme_bw() + facet_wrap(~ Season, nrow=1) + xlim(-1, 1)

CCC_RSHP_FINAL %>% filter(gas == "CH4") %>%  ggplot( aes( x= CCC_CGF, y =CCC, col=Approach)) + geom_point() + theme_bw()  + xlim(-1, 1) + ylim(-1, 1) +  geom_abline(intercept = 0, slope = 1,linetype = "dashed") 


# Save the files: 

fileSave <- fs::path(localdir,paste0("CCC_CH4.Rdata"))
save( CCC_RSHP_FINAL, file=fileSave)
googledrive::drive_upload(media = fileSave, overwrite = T, path = drive_url)

load(fileSave)

# Harmonize data and test result: 

harmonize_val <- function( tibble){
  
  site.list <- names(tibble )
  harmonized <- list()
  for( i in site.list){
    print(i)
    
    df <- tibble[[i]]
    
    df %>% names
    
    harmonized[[i]] <-   df %>% reframe(.by= c(timeEndA.local, gas),
                                      FG_mean = median(FG_mean, na.rm=T),
                                      EC_mean = mean(EC_mean, na.rm=T),
                                      Tair_C = mean(Tair_C, na.rm=T),
                                      PAR = mean(PAR, na.rm=T))  %>% 
      mutate( Month = format(timeEndA.local, '%m') %>% as.numeric) %>% mutate(
        Season = case_when(
          Month %in% c(12, 1, 2) ~ "Winter",
          Month %in% 3:5 ~ "Spring",
          Month %in% 6:8 ~ "Summer",
          Month %in% 9:11 ~ "Autumn"))
  
  }
  return(harmonized )

}

SITEval_DATA_FILTERED_Final_CGF_RSHP_H <- harmonize_val(tibble = SITEval_DATA_FILTERED_Final_CGF_RSHP)

SITEval_DATA_FILTERED_Final_RSHP_H <- harmonize_val(tibble = SITEval_DATA_FILTERED_Final_RSHP)

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


CCC_RSHP_FINAL_H <- harmonized_CCC(tibble=SITEval_DATA_FILTERED_Final_RSHP_H )
CCC_RSHP_CGF_FINAL_H <- harmonized_CCC(tibble=SITEval_DATA_FILTERED_Final_CGF_RSHP )

# Data VIZ:
CCC_RSHP_FINAL_H %>% ggplot() + geom_point( aes( x= CCC, y = Site, col=Season)) + facet_wrap( ~gas)

CCC_RSHP_CGF_FINAL_H %>% ggplot() + geom_point( aes( x= CCC, y = Site, col=Season)) + facet_wrap( ~gas)


# DEVLOP DIEL AND Q10 DATASETS

