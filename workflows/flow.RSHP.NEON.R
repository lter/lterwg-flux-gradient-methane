
# RSHP: ####
library(randomForest)
source(fs::path(DirRepo.eval,'functions/calc_validation.R' ))

load( fs::path(localdir,paste0("SITE_DATA_FILTERED_CH4.Rdata")))

NEON_RSHP <- CCC_SamplingHeightPairs(DATA = SITE_DATA_FILTERED) 

NEON_RSHP %>% filter(gas == "CH4") %>%  ggplot( aes( x= CCC, y = Site)) + geom_boxplot() + theme_bw()  + xlim(-1, 1) 

# select the good pairs
file1 <-NEON_RSHP %>% 
  mutate(Combination = paste(var1, var2, sep ="-")) %>%
  rename( var= var1)%>% select(Site, gas, var, CCC )


file2 <- NEON_RSHP %>% 
  mutate(Combination = paste(var1, var2, sep ="-")) %>% 
  rename( var= var2) %>% select(Site, gas, var, CCC )

NEON_RSHP_EVAL <- rbind( file1, file2) %>% 
  reframe( .by = c( Site, gas, var),
           CCC.GF.max = max( CCC, na.rm=T), 
           CCC.GF.min = min( CCC, na.rm=T), 
           CCC.GF.mean = mean(CCC, na.rm=T), 
           CCC.GF.median = median( CCC, na.rm=T),
           CCC.GF.range = mean( CCC, na.rm=T)) %>% 
  mutate(Good.CCC.GF = case_when(
    CCC.GF.max >= 0.95 ~ 1,.default = 0) %>% as.factor) 

load(file= paste(localdir, 'SITE_RSHP_MODEL.Rdata', sep="") ) # The model is in here:

NEON_RSHP_EVAL$RSHP.rf <- predict(rf.good.ccc.ec, NEON_RSHP_EVAL)
NEON_RSHP_EVAL$RSHP.rf %>% summary
NEON_RSHP_EVAL$Good.CCC.GF %>% summary

library(caret)

confusionMatrix(NEON_RSHP_EVAL$RSHP.rf, NEON_RSHP_EVAL$Good.CCC.GF)

NEON_RSHP_EVAL %>% filter(RSHP.rf == 1) %>% ggplot() + geom_boxplot( aes(x=CCC.GF.max, y = Site)) + 
  facet_wrap(~gas) + theme_bw() + xlim(-1,1)

NEON_RSHP_EVAL %>% filter(Good.CCC.GF == 1) %>% ggplot() + geom_boxplot( aes(x=CCC.GF.max, y = Site)) + 
  facet_wrap(~gas) + theme_bw() + xlim(-1,1)

save( NEON_RSHP, NEON_RSHP_EVAL, 
      file=fs::path(localdir,paste0("NEON_RSHP_CCC.Rdata")))
