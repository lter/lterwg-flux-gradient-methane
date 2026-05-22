# Determine the RSHP based on CCC >= 0.5 for CO2 and H2O:

source(fs::path(DirRepo.ch4,'functions/calc_CCC_SHP.R' ))

# Get a list of good levels at each site for CO2 and H2O:
load(fs::path(localdir,paste0("/NEON_GradientFlux_Data_Filter/SITES_One2One_AA_AW.Rdata")))

## Make a list of good levels for CO2 and H20:

RSHP.CO2 <- SITES_One2One %>% 
  filter(Good.CCC == 1, gas == "CO2") %>% 
  select( Site, dLevelsAminusB, Approach) %>% mutate( RSHP.CO2 = 1)

RSHP.H2O <- SITES_One2One %>% 
  filter(Good.CCC == 1, gas == "H2O") %>% 
  select( Site, dLevelsAminusB, Approach) %>% mutate( RSHP.H2O = 1)


NEON_RSHP_EVAL <- RSHP.CO2 %>% 
  left_join( RSHP.H2O, by = c('Site', 'dLevelsAminusB', 'Approach')) 


NEON_RSHP_EVAL <- NEON_RSHP_EVAL %>% mutate(RSHP.CO2 = replace_na(RSHP.CO2, 0) %>% as.factor,
                                            RSHP.H2O = replace_na(RSHP.H2O, 0) %>% as.factor)

library(caret)


confusionMatrix(NEON_RSHP_EVAL$RSHP.CO2, NEON_RSHP_EVAL$RSHP.H2O)

NEON_RSHP_EVAL$RSHP.CO2 %>% summary
NEON_RSHP_EVAL$RSHP.H2O %>% summary

save( NEON_RSHP_EVAL, 
      file=fs::path(localdir.ch4 ,paste0("NEON_RSHP_CCC.Rdata")))