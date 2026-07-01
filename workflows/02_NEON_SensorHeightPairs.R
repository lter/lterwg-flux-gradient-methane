# Determine the RSHP based on CCC >= 0.5 for CO2 and H2O:

source(fs::path(DirRepo.ch4,'functions/calc_CCC_SHP.R' ))

# Get a list of good levels at each site for CO2 and H2O:
load(fs::path(localdir,paste0("/NEON_GradientFlux_Data_Filter/SITES_One2One_AA_AW.Rdata")))


# Determine the maximum CCC at a site 

SITES_One2One <- SITES_One2One %>% group_by(Site, gas) %>% mutate(Site.max.CCC= max(CCC)) %>% ungroup() %>% as.data.frame()

SITES_One2One$Site.max.CCC %>% summary
SITES_One2One$CCC %>% summary
## Make a list of good levels for CO2 and H20:

all.sites <- SITES_One2One$Site %>% unique

RSHP.CO2 <- SITES_One2One %>% 
  filter(Good.CCC == 1, gas == "CO2") %>% 
  select( Site, dLevelsAminusB, Approach) %>% mutate( RSHP.CO2 = 1)

missing.CO2 <- setdiff(all.sites, RSHP.CO2$Site %>% unique) 

RSHP.CO2.add <- SITES_One2One %>%
  filter( Site %in% missing.CO2, gas == "CO2", CCC == Site.max.CCC ) %>%
  select( Site, dLevelsAminusB, Approach) %>% mutate( RSHP.CO2 = 1)


RSHP.H2O <- SITES_One2One %>% 
  filter(Good.CCC == 1, gas == "H2O") %>% 
  select( Site, dLevelsAminusB, Approach) %>% mutate( RSHP.H2O = 1)

missing.H2O  <- setdiff(all.sites, RSHP.H2O$Site %>% unique) 

RSHP.H2O.add <- SITES_One2One %>% 
  filter(Site %in% missing.H2O , gas == "H2O", CCC == Site.max.CCC ) %>% 
  select( Site, dLevelsAminusB, Approach) %>% mutate( RSHP.H2O = 1)


RSHP.CO2.all <- rbind( RSHP.CO2, RSHP.CO2.add)
RSHP.H2O.all <- rbind( RSHP.H2O, RSHP.H2O.add)

NEON_RSHP_EVAL <- RSHP.CO2.all %>%
  full_join( RSHP.H2O.all, by = c('Site', 'dLevelsAminusB', 'Approach'))


NEON_RSHP_EVAL <- NEON_RSHP_EVAL %>% mutate(RSHP.CO2 = replace_na(RSHP.CO2, 0) %>% as.factor,
                                            RSHP.H2O = replace_na(RSHP.H2O, 0) %>% as.factor)

NEON_RSHP_EVAL$Site %>% unique # All sites should be present
library(caret)


confusionMatrix(NEON_RSHP_EVAL$RSHP.CO2, NEON_RSHP_EVAL$RSHP.H2O)

NEON_RSHP_EVAL$RSHP.CO2 %>% summary
NEON_RSHP_EVAL$RSHP.H2O %>% summary




save( NEON_RSHP_EVAL, 
      file=fs::path(localdir.ch4 ,paste0("NEON_RSHP_CCC.Rdata")))

