# Diel for Validation Data:

# The Diel analysis is currently set up by season for the ENSEMBLE data:

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)

load( fs::path(localdir,paste0("/Validation_Sites/SITEval_DATA_FILTERED_RSHP_EnSEMBLE.Rdata")) )

source( fs::path(DirRepo.ch4, paste('/functions/calc_diel_ch4.R')))
source( fs::path(DirRepo.ch4, paste('/functions/calc_Q10.R')))

canopy.info <- read.csv( file.path(paste(localdir, "Val_canopy.csv", sep="/")))%>% rename(Site = SITE_ID) %>% 
  mutate(EcoType = case_when( IGBP == 'ENF' ~ 'Forest',
                              IGBP == 'DBF' ~ 'Forest', 
                              IGBP == 'MF' ~ 'Forest',
                              IGBP == 'EBF' ~ 'Forest',
                              IGBP == 'SAV' ~ 'Forest',
                              IGBP == 'WET' ~ 'Wetland',
                              
                              IGBP == 'GRA' ~ 'Grassland',
                              
                              IGBP == 'CVM' ~ 'Cropland',
                              IGBP == 'CRO' ~ 'Cropland',
                              IGBP == 'OSH' ~ 'Shrubland')) %>% 
  filter(Canopy_L1 == "AA"| Canopy_L1 == "AW")


# Calculate diels by Season: ####
# RSHP CCC 0:

# Calculate diels by Season: ####
ENSEMBLE_DIELSc <- data.frame()
ENSEMBLE_Q10_eq4c <- data.frame()

site.list <- SITEval_DATA_FILTERED_RSHPc_H %>% names 

options <- c('FG_mean', 'EC_mean' )

for( site in site.list){
  print(site)
  
  for (opt in options){
    print(paste("Running option =", opt, sep=" "))
    
    data <- SITEval_DATA_FILTERED_RSHPc_H[[site]] %>% 
      mutate(hour = format(time.rounded,'%H'),
             YearMon = format(time.rounded, "%Y-%m"),
             season = Season) %>% 
      distinct %>% filter(gas == "CH4")
    
    data$FG_mean <- data[[opt]]
    message( paste("Running DIEL functions- CH4 for ", site))
    # Calculate Diurnal Patterns by Year-month:
    
    for( i in data$YearMon %>% unique){
      
      data.sub <- data %>% filter( YearMon == i)
      
      try(DIEL.CH4 <- DIEL.season.ch4( dataframe = data.sub, 
                                       flux = opt, 
                                       Gas = "CH4") %>% 
            mutate(gas= "CH4", 
                   YearMon = i,
                   site = site,
                   RSHP = opt), silent =TRUE)
      
      try(DIEL.CH4 %>%  ggplot() + geom_point(aes(x=Hour, y = DIEL, col=season)) + facet_wrap(~season) , silent =TRUE)
      
      try( ENSEMBLE_DIELSc <- rbind( ENSEMBLE_DIELSc, DIEL.CH4 ), silent =TRUE)
      
      try( rm( DIEL.CH4), silent=T)
    }
    
    # Q10:
    
    try(NEON_TRC_PARMS_04 <- TRC_PARMS_04(data.frame = data,
                                          iterations = 5000,
                                          priors.trc = brms::prior("normal(2.0, 0.3)", nlpar = "Q10", lb = 1.0, ub = 5) +
                                            brms::prior("normal(0.5, 0.3)", nlpar = "Rref", lb = 0.001, ub = 5),
                                          idx.colname = 'YearMon',
                                          NEE.colname = "FG_mean",
                                          TA.colname = 'Tair_C',
                                          Tref = 1)  %>%  mutate( SITE_ID = site, RSHP = opt), silent=T)
    
    try( ENSEMBLE_Q10_eq4c <- rbind(ENSEMBLE_Q10_eq4c ,NEON_TRC_PARMS_04) , silent=T)
  }}

ENSEMBLE_DIELS_Sitec <- ENSEMBLE_DIELSc %>% 
  reframe( .by=c(site, gas, season, RSHP), 
           total.Flux = sum(DIEL),
           max.Flux = max(DIEL),
           min.Flux = min(DIEL),
           total.count = sum(count))

save(ENSEMBLE_DIELSc ,
     ENSEMBLE_Q10_eq4c,
     ENSEMBLE_DIELS_Sitec,
     file = fs::path(localdir,paste0("/Validation_Sites/Val_PARMS_DIEL_Q10_RSHPc.Rdata")))


# DIEL PLOTS: ####

load(file = fs::path(localdir,paste0("/Validation_Sites/Val_PARMS_DIEL_Q10_RSHPc.Rdata")))

ENSEMBLE_DIELSc  %>%  ggplot()+ 
  geom_line(aes(x = Hour , y = DIEL, col=season, linetype = RSHP), size=0.25) +
  scale_colour_discrete_qualitative(palette = "Harmonic") +
  theme_bw() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "transparent", linewidth = 0.5),
        legend.title = element_blank()) + ylab(expression(paste( "GF (g m"^-2, ")"))) + 
  facet_wrap(~site, scales = "free_y")

# Comparison Dataframes & validation

ENSEMBLE_DIELSc %>% summary

ENSEMBLE_DIELSc_FG <- 
  ENSEMBLE_DIELSc %>% filter( RSHP == "FG_mean") %>% 
  rename(FG_DIEL =DIEL, FG_DIEL.SE = DIEL.SE,
         FG_Peak.Hour =Peak.Hour , FG_Min.Hour = Min.Hour, FG_count = count )

ENSEMBLE_DIELSc_EC <- 
  ENSEMBLE_DIELSc %>% filter( RSHP == "EC_mean") %>% 
  rename(EC_DIEL =DIEL, EC_DIEL.SE = DIEL.SE,
         EC_Peak.Hour =Peak.Hour , EC_Min.Hour = Min.Hour, EC_count = count )


ENSEMBLE_DIELSc_wide <- ENSEMBLE_DIELSc_FG %>% left_join(ENSEMBLE_DIELSc_EC,
                                                            by=c("Hour",'YearMon', 'site', 'gas')) %>% 
  reframe(.by=c('YearMon', 'site', 'gas'),
          EC_DIEL= sum(EC_DIEL),
          FG_DIEL = sum(FG_DIEL)) %>% 
  mutate( Date = paste0(YearMon, "-01") %>% as.Date( format = "%Y-%m-%d") ,
          month = format(  Date,'%m') %>% as.numeric,
                 Season = case_when(
                   month %in% c(12, 1, 2) ~ "Winter",
                   month %in% c(3, 4, 5) ~ "Spring",
                   month %in% c(6, 7, 8) ~ "Summer",
                   TRUE ~ "Autumn"))


lm( ENSEMBLE_DIELSc_wide$EC_DIEL ~
      ENSEMBLE_DIELSc_wide$FG_DIEL) %>% summary

ENSEMBLE_DIELSc_wide %>% 
  ggplot(aes(x= EC_DIEL*1000, y=FG_DIEL*1000, col=Season)) +
  geom_point(alpha=0.5, size=3) +
  geom_smooth(method = "lm", se = FALSE, size=0.5, col="grey40",
              formula = 'y ~ x + 0') +
  stat_regline_equation(aes(label = ..eq.label..), 
                        col='black', size =3.5,
                        formula = 'y ~ x + 0', fontface=2,
                        label.x.npc = "left",
                        label.y.npc = "top") +
  stat_regline_equation(aes(label = ..rr.label..), 
                        col='black', size =3.5,
                        formula = 'y ~ x + 0', fontface=2,
                        label.x.npc = "center",
                        label.y.npc = "bottom") +
  theme_bw() + 
  geom_abline(slope=1, col="grey40", linetype="dashed", size =1) +
  scale_colour_discrete_qualitative(palette = "Harmonic") + facet_wrap(~site, scales = "free")

# Q10: ####

ENSEMBLE_Q10_eq4c_FG <- ENSEMBLE_Q10_eq4c %>% filter( RSHP == 'FG_mean') %>% rename(Q10.mean_FG =Q10.mean, Q10.se_FG= Q10.se, Q10.Bulk_ESS_FG =Q10.Bulk_ESS,
                                                                                  Q10.Tail_ESS_FG =Q10.Tail_ESS, Q10.Rhat_FG =Q10.Rhat, 
                                                                                  Rref.mean_FG = Rref.mean, Rref.se_FG = Rref.se, 
                                                                                  Rref.Bulk_ESS_FG =Rref.Bulk_ESS, Rref.Tail_ESS_FG = Rref.Tail_ESS,
                                                                                  Rref.Rhat_FG=Rref.Rhat, Tref_FG = Tref, samples_FG = samples)

ENSEMBLE_Q10_eq4c_EC <- ENSEMBLE_Q10_eq4c %>% filter( RSHP == 'EC_mean') %>% rename(Q10.mean_EC =Q10.mean, Q10.se_EC= Q10.se, Q10.Bulk_ESS_EC =Q10.Bulk_ESS,
                                                                                  Q10.Tail_ESS_EC =Q10.Tail_ESS, Q10.Rhat_EC =Q10.Rhat, 
                                                                                  Rref.mean_EC = Rref.mean, Rref.se_EC = Rref.se, 
                                                                                  Rref.Bulk_ESS_EC =Rref.Bulk_ESS, Rref.Tail_ESS_EC = Rref.Tail_ESS,
                                                                                  Rref.Rhat_EC=Rref.Rhat, Tref_EC = Tref, samples_EC = samples)

ENSEMBLE_Q10_eq4c_wide <- ENSEMBLE_Q10_eq4c_FG %>% 
  full_join(ENSEMBLE_Q10_eq4c_EC, by = c('idx', 'SITE_ID'))


ENSEMBLE_Q10_eq4c_wide %>% ggplot( aes( x= Q10.mean_EC, y=Q10.mean_FG, col=SITE_ID)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE, size=0.5, col="grey40",
                           formula = 'y ~ x + 0') +
  stat_regline_equation(aes(label = ..eq.label..), 
                        col='black', size =3.5,
                        formula = 'y ~ x + 0', fontface=2,
                        label.x.npc = "left",
                        label.y.npc = "top") +
  stat_regline_equation(aes(label = ..rr.label..), 
                        col='black', size =3.5, fontface=2,
                        label.x.npc = "center",
                        label.y.npc = "bottom") +
  geom_abline( slope=1, col="grey40", linetype="dashed", size =1) +
  theme_bw() 



ENSEMBLE_Q10_eq4c_wide %>% ggplot( aes( x= Rref.mean_EC, y=Rref.mean_FG, col=SITE_ID)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, size=0.5, col="grey40",formula = 'y ~ x + 0' ) +
  stat_regline_equation(aes(label = ..eq.label..), 
                        col='black', size =3.5,
                        formula = 'y ~ x + 0', fontface=2,
                        label.x.npc = "left",
                        label.y.npc = "top") +
  stat_regline_equation(aes(label = ..rr.label..), 
                        col='black', size =3.5, fontface=2,
                        label.x.npc = "center",
                        label.y.npc = "bottom") +
  geom_abline( slope=1, col="grey40", linetype="dashed", size =1) +
  theme_bw() 

# Save all the files needed for plots:
save(ENSEMBLE_DIELSc ,
     ENSEMBLE_Q10_eq4c,
     ENSEMBLE_DIELSc_wide ,
     ENSEMBLE_Q10_eq4c_wide,
     ENSEMBLE_DIELS_Sitec,
     file = fs::path(localdir,paste0("/Validation_Sites/Val_PARMS_DIEL_Q10_RSHPc_wide.Rdata")))
