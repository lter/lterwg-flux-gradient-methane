# Determine the RSHP based on CCC >= 0.5 for CO2 and H2O, and the CH4 pairs (= CO2 RSHPs):
# (calc_CCC_SHP.R not needed: CCC comes from the evaluation run's SITES_One2One)

# Get a list of good levels at each site for CO2 and H2O:
# v3b evaluation run (lterwg-flux-gradient-eval, run.tag "_AA_AW_2021_2026_v3b"): CO2/H2O CCC against EC per pair
if (!exists("one2one.file")) one2one.file <- fs::path(localdir, "NEON_GradientFlux_Data_Filter/SITES_One2One_AA_AW_2021_2026_v3b.RData")
load(one2one.file)


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

# CH4 pairs (no CH4 EC at NEON): the pairs that are reliable for CO2 (CCC >= 0.5 against CO2 EC).
# Chosen by testing the rules on CO2 and H2O, where EC can score them (essd-2026-413, Supplement S2):
#  - CO2-RSHP pairs used for H2O keep the H2O benchmark (median site CCC 0.67 vs 0.69 with H2O's own RSHPs);
#  - pairs reliable only for H2O ("either gas"), or extra pairs screened against a gradient-flux reference
#    (even at CCC_ref >= 0.8, or only to fill gaps), lower agreement with EC at many sites.
# The max-CCC CO2 pair is NOT used for CH4: sites with no CO2 RSHP get the no-EC (S2) fallback in step 03
# (AA MBR reference + AE/WP pairs with CCC >= 0.5 against it, rescaled to it).
RSHP.CH4 <- RSHP.CO2 %>% distinct(Site, dLevelsAminusB, Approach) %>% mutate(RSHP.CH4 = 1)
CH4.fallback.sites <- setdiff(all.sites, RSHP.CH4$Site)
message("CH4 pairs from CO2 RSHPs at ", length(unique(RSHP.CH4$Site)), " sites; S2 fallback at: ",
        paste(CH4.fallback.sites, collapse = ", "))

NEON_RSHP_EVAL <- RSHP.CO2.all %>%
  full_join( RSHP.H2O.all, by = c('Site', 'dLevelsAminusB', 'Approach')) %>%
  full_join( RSHP.CH4, by = c('Site', 'dLevelsAminusB', 'Approach'))


NEON_RSHP_EVAL <- NEON_RSHP_EVAL %>% mutate(RSHP.CO2 = replace_na(RSHP.CO2, 0) %>% as.factor,
                                            RSHP.H2O = replace_na(RSHP.H2O, 0) %>% as.factor,
                                            RSHP.CH4 = replace_na(RSHP.CH4, 0) %>% as.factor)

NEON_RSHP_EVAL$Site %>% unique # All sites should be present
if (requireNamespace("caret", quietly = TRUE)) print(caret::confusionMatrix(NEON_RSHP_EVAL$RSHP.CO2, NEON_RSHP_EVAL$RSHP.H2O))

NEON_RSHP_EVAL$RSHP.CO2 %>% summary
NEON_RSHP_EVAL$RSHP.H2O %>% summary

save( NEON_RSHP_EVAL, CH4.fallback.sites,
      file=fs::path(localdir.ch4 ,paste0("NEON_RSHP_CCC.Rdata")))

