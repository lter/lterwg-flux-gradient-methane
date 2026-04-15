# RSHP: ####
source(fs::path(DirRepo.ch4,'functions/calc_CCC_SHP.R' ))
load( fs::path(localdir,paste0("SITE_DATA_FILTERED_CH4.Rdata")))
load(file= paste(localdir, 'SITE_RSHP_MODEL.Rdata', sep="") )

NEON_RSHP <- CCC_SamplingHeightPairs(DATA = SITE_DATA_FILTERED) 
NEON_RSHP %>% names()

# select the good pairs
file1 <-NEON_RSHP %>% 
  mutate(Combination = paste(var1, var2, sep ="-")) %>%
  rename( var= var1)%>% select(Site, gas, var, CCC ) %>% separate(
    col = var, 
    into = c("Approach", "dLevelsAminusB"), 
    sep = "-")


file2 <- NEON_RSHP %>% 
  mutate(Combination = paste(var1, var2, sep ="-")) %>% 
  rename( var= var2) %>% select(Site, gas, var, CCC ) %>% separate(
    col = var, 
    into = c("Approach", "dLevelsAminusB"), 
    sep = "-")

NEON_RSHP_EVAL <- rbind( file1, file2) %>% 
  reframe( .by = c( Site, gas, Approach, dLevelsAminusB),
           CCC.GF.max = max( CCC, na.rm=T), 
           CCC.GF.min = min( CCC, na.rm=T), 
           CCC.GF.mean = mean(CCC, na.rm=T), 
           CCC.GF.median = median( CCC, na.rm=T),
           CCC.GF.range = CCC.GF.max - CCC.GF.min) %>% distinct


save(NEON_RSHP_EVAL, file=fs::path(localdir,paste0("NEON_RSHP_CCC.Rdata")))

load(file=fs::path(localdir,paste0("NEON_RSHP_CCC.Rdata")))
     
# Get a list of good levels at each site for CO2 and H2O:
load(fs::path(localdir,paste0("SITES_One2One_AA_AW.Rdata")))

## Make a list of good levels for CO2 and H20:

RSHP.CO2 <- SITES_One2One %>% 
  filter(Good.CCC == 1, gas == "CO2") %>% 
  select( Site, dLevelsAminusB, Approach) %>% mutate( RSHP.CO2 = 1)

RSHP.H2O <- SITES_One2One %>% 
  filter(Good.CCC == 1, gas == "H2O") %>% 
  select( Site, dLevelsAminusB, Approach) %>% mutate( RSHP.H2O = 1)


NEON_RSHP_EVAL <- NEON_RSHP_EVAL %>% 
  left_join( RSHP.CO2, by = c('Site', 'dLevelsAminusB', 'Approach')) %>% 
  left_join( RSHP.H2O, by = c('Site', 'dLevelsAminusB', 'Approach')) 


NEON_RSHP_EVAL <- NEON_RSHP_EVAL %>% mutate(RSHP.CO2 = replace_na(RSHP.CO2, 0) %>% as.factor,
                                            RSHP.H2O = replace_na(RSHP.H2O, 0) %>% as.factor)

library(caret)

NEON_RSHP_EVAL$RSHP.rf <- predict( rf.good.ccc.ec,NEON_RSHP_EVAL)
confusionMatrix(NEON_RSHP_EVAL$RSHP.rf, NEON_RSHP_EVAL$RSHP.CO2)
confusionMatrix(NEON_RSHP_EVAL$RSHP.rf, NEON_RSHP_EVAL$RSHP.H2O)
confusionMatrix(NEON_RSHP_EVAL$RSHP.CO2, NEON_RSHP_EVAL$RSHP.H2O)

NEON_RSHP_EVAL$RSHP.rf %>% summary
NEON_RSHP_EVAL$RSHP.CO2 %>% summary
NEON_RSHP_EVAL$RSHP.H2O %>% summary

save( NEON_RSHP, NEON_RSHP_EVAL, 
      file=fs::path(localdir,paste0("NEON_RSHP_CCC.Rdata")))