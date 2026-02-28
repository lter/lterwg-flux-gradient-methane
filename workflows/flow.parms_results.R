# Results 
library(tidyverse)
library( ggplot2)

load(file = fs::path(localdir,paste0("NEON_PARMS_DIEL_Q10.Rdata")))

ENSEMBLE_DIELS %>% 
ENSEMBLE_Q10_eq4 
ENSEMBLE_Q10_eq5 

ENSEMBLE_DIELS_Site <- ENSEMBLE_DIELS %>% 
  reframe( .by=c(site, gas, season), 
           total.FG = sum(DIEL),
           max.FG = max(DIEL),
           min.FG = min(DIEL),
           Night.DIFF = max.FG-max.FG,
           total.count = sum(count))

# DIEL PLOTS: ####

ENSEMBLE_DIELS_Site  %>%  ggplot() + 
  geom_point(aes(x = total.FG , y = site )) +
  theme_bw() + geom_vline(xintercept = 0) + facet_wrap( ~season)


ENSEMBLE_DIELS_Site  %>%  ggplot() + 
  geom_boxplot(aes(x = total.FG , y = site )) +
  theme_bw() + geom_vline(xintercept = 0, linetype='dashed',
                          col='red') 

ENSEMBLE_Q10_eq4$idx

ENSEMBLE_Q10_eq4  %>%  ggplot() + 
  geom_boxplot(aes(x = Q10.mean, y = SITE_ID )) +
  theme_bw() 

ENSEMBLE_Q10_eq4  %>%  ggplot() + 
  geom_boxplot(aes(x = Rref.mean, y = SITE_ID )) +
  theme_bw() 

ENSEMBLE_Q10_eq4  %>%  ggplot() + 
  geom_point(aes(x = Rref.mean, y = Q10.mean )) +
  theme_bw() 

