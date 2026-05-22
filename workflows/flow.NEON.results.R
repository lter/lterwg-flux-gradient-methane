# Figures for parms:
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)
library(sf)
library(AOI)

# Import data and build the files needed: ####
localdir <- '/Volumes/MaloneLab/Research/FluxGradient/METHANE'

# Load parameter information:
load(file = fs::path(localdir,paste0("NEON_PARMS_DIEL_Q10.Rdata")))
# Load soil and root data
load(file=paste(localdir, "/", 'Soildata_YearMon.Rdata', sep=""))
# Load data obtained from Avni:
avni <- read.csv( paste("/Volumes/MaloneLab/Research/FluxGradient/Avni's Data/neondata_forSM.csv" , sep="" )) %>% rename( SITE_ID =site_code )

avni2 <- read.csv( paste("/Volumes/MaloneLab/Research/FluxGradient/Avni's Data/neondata_forsteve.csv" , sep="" )) %>% rename( SITE_ID =site_code ) %>% dplyr::select(
  SITE_ID, MAP, MAT, aridity.index,shallow_root, deep_root)

# Load site information:
soil.data.biogeo <- read.csv(file=paste(localdir, "/", 'Soil_Biogeochem_RootBiomass.csv', sep=""))
# Load canopy information: 
canopy <- read.csv("/Volumes/MaloneLab/Research/FluxGradient/canopy_commbined.csv") %>% distinct

canopy.Ht <- canopy %>% rename(SITE_ID = Site) %>% 
  reframe( .by= SITE_ID,
           canopyHeight_m = mean(canopyHeight_m))

metadata <- read.csv('/Volumes/MaloneLab/Research/FluxGradient/Ameriflux_NEON field-sites.csv') %>% 
  mutate(EcoType = case_when( Vegetation.Abbreviation..IGBP. == 'ENF' ~ 'Forest',
                              Vegetation.Abbreviation..IGBP. == 'DBF' ~ 'Forest', 
                              Vegetation.Abbreviation..IGBP. == 'MF' ~ 'Forest',
                              Vegetation.Abbreviation..IGBP. == 'EBF' ~ 'Forest',
                              Vegetation.Abbreviation..IGBP. == 'SAV' ~ 'Forest',
                              Vegetation.Abbreviation..IGBP. == 'WET' ~ 'Wetland',
                              
                              Vegetation.Abbreviation..IGBP. == 'GRA' ~ 'Grassland',
                              
                              Vegetation.Abbreviation..IGBP. == 'CVM' ~ 'Cropland',
                              Vegetation.Abbreviation..IGBP. == 'CRO' ~ 'Cropland',
                              Vegetation.Abbreviation..IGBP. == 'OSH' ~ 'Shrubland')) %>% 
  rename( site = Site_Id.NEON,
          SITE_ID = Site_Id.NEON) %>% 
  left_join(canopy.Ht , by = 'SITE_ID') %>% 
  left_join(soil.data.biogeo %>% rename(SITE_ID = siteID)  , by = 'SITE_ID')%>% 
  left_join(avni , by = 'SITE_ID')%>% 
  left_join(avni2 , by = 'SITE_ID')

# Create a simple feature
site.att.sf <- st_as_sf(x = metadata,                         
                        coords = c("Longitude..degrees.", "Latitude..degrees."), crs = 4326) 

#format data:
ENSEMBLE_DIELS_canopy <-  ENSEMBLE_DIELS %>% rename( SITE_ID = site) %>% left_join( metadata, by = "SITE_ID")

ENSEMBLE_Q10_eq4_canopy <-   ENSEMBLE_Q10_eq4 %>% 
  mutate( Date = as.Date(paste(idx, "-01", sep="")),
          month = format(Date, "%m") %>% as.numeric,
          Season = case_when(
            month %in% c(12, 1, 2) ~ "Winter",
            month %in% c(3, 4, 5) ~ "Spring",
            month %in% c(6, 7, 8) ~ "Summer",
            TRUE ~ "Autumn")) %>% left_join( metadata, by = "SITE_ID")

# This file has monthly data in it

parms <- ENSEMBLE_Q10_eq4_canopy %>% 
  dplyr::select(Q10.mean, SITE_ID, Date, Rref.mean, EcoType, Season, sulfurTot, dryMass, acidity, FRC.kg.C.m.2, SOC.kg.C.m.2, beta_roots, beta_soc, MAP, MAT, aridity.index,shallow_root, deep_root) %>% 
  left_join( ENSEMBLE_DIELS %>% 
               mutate(Date = as.Date(paste(YearMon, "-01", sep="")),
                      SITE_ID = site) %>% 
               reframe( .by = c(SITE_ID, Date), 
                        DIEL = sum(DIEL)), by=c('SITE_ID', 'Date')) %>% 
  left_join( Site_SoilData %>% 
               mutate(Date = as.Date(paste(YearMon, "-01", sep="")), SITE_ID = Site),by = c('SITE_ID', 'Date')) %>% 
  mutate( Source = case_when( DIEL > 0 ~ 1,
                              DIEL <= 0 ~ 0) %>% as.factor)

parms.site.season <- parms %>% reframe( .by = c(SITE_ID, Season),
                                        VSWCMean  = mean( VSWCMean, na.rm=TRUE), 
                                        Q10.mean = mean(Q10.mean, na.rm=TRUE) ,
                                        Rref.mean = mean(Rref.mean, na.rm=TRUE), 
                                        VSWCMin = mean(VSWCMin, na.rm=TRUE),
                                        VSWCVar= mean(VSWCVar, na.rm=TRUE),
                                        sulfurTot= mean(sulfurTot, na.rm=TRUE),
                                        dryMass= mean( dryMass, na.rm=TRUE),
                                        acidity= mean( acidity, na.rm=TRUE),
                                        DIEL = mean(DIEL, na.rm=T),
                                        FRC.kg.C.m.2 = mean(FRC.kg.C.m.2, na.rm=T), 
                                        SOC.kg.C.m.2 = mean(SOC.kg.C.m.2, na.rm=T),
                                        beta_roots= mean(beta_roots, na.rm=T),  
                                        beta_soc= mean(beta_soc, na.rm=T),
                                        MAT= mean(MAT, na.rm=T), 
                                        MAP= mean(MAP, na.rm=T))  %>% 
  mutate( Source = case_when( DIEL > 0 ~ "Source",
                              DIEL <= 0 ~ "Sink") %>% as.factor)

parms.site <- parms.site.season %>% mutate( Source = case_when( DIEL > 0 ~ 1,
                                                                DIEL <= 0 ~ 0) %>% as.numeric) %>% mutate(across(where(is.numeric), ~ ifelse(is.infinite(.), 0, .))) %>% 
  reframe( .by = SITE_ID, 
           Source = sum(Source ) %>% as.factor,
           VSWCMean  = mean( VSWCMean, na.rm=TRUE), 
           Q10.mean = mean(Q10.mean, na.rm=TRUE) ,
           Rref.mean = mean(Rref.mean, na.rm=TRUE), 
           VSWCMin = mean(VSWCMin, na.rm=TRUE),
           VSWCVar= mean(VSWCVar, na.rm=TRUE),
           sulfurTot= mean(sulfurTot, na.rm=TRUE),
           dryMass= mean( dryMass, na.rm=TRUE),
           acidity= mean( acidity, na.rm=TRUE),
           DIEL = mean(DIEL, na.rm=T),
           FRC.kg.C.m.2 = mean(FRC.kg.C.m.2, na.rm=T), 
           SOC.kg.C.m.2 = mean(SOC.kg.C.m.2, na.rm=T),
           beta_roots= mean(beta_roots, na.rm=T),  
           beta_soc= mean(beta_soc, na.rm=T),
           MAT= mean(MAT, na.rm=T), 
           MAP= mean(MAP, na.rm=T),
           Source_Sink = case_when( DIEL > 0 ~ "Source",
                                    DIEL <= 0 ~ "Sink") %>% as.factor)  %>% 
  left_join( metadata %>% dplyr::select(EcoType, SITE_ID), by = "SITE_ID")

parms.site.season.sf <-site.att.sf %>% left_join( parms.site.season , by="SITE_ID")

parms.site.sf <-site.att.sf %>% left_join( parms.site , by="SITE_ID")

# Drivers of CH4, stored as DIEL in parms: ####

parms.diel <- parms %>%
  mutate(DIEL_mgC_m2_day = DIEL * 1000) %>%
  filter(is.finite(DIEL_mgC_m2_day))

# Soil-moisture variables vary monthly within sites, so evaluate them as anomalies.
monthly.driver.variables <- c('VSWCMean', 'VSWCMin', 'VSWCVar') %>%
  intersect(names(parms.diel))

DIEL.monthly.driver.data <- parms.diel %>%
  dplyr::select(SITE_ID, Date, Season, DIEL_mgC_m2_day, dplyr::all_of(monthly.driver.variables)) %>%
  group_by(SITE_ID) %>%
  mutate(
    DIEL_site_mean = mean(DIEL_mgC_m2_day, na.rm = TRUE),
    DIEL_anomaly = DIEL_mgC_m2_day - DIEL_site_mean,
    across(dplyr::all_of(monthly.driver.variables), ~ .x - mean(.x, na.rm = TRUE), .names = "{.col}_anomaly")
  ) %>%
  ungroup()

DIEL.monthly.driver.summary <- purrr::map_dfr(monthly.driver.variables, function(driver) {
  driver.df <- DIEL.monthly.driver.data %>%
    dplyr::select(SITE_ID, Season, DIEL_anomaly, driver_anomaly = dplyr::all_of(paste0(driver, "_anomaly"))) %>%
    drop_na() %>%
    filter(is.finite(DIEL_anomaly), is.finite(driver_anomaly))
  
  if (nrow(driver.df) < 10 || length(unique(driver.df$driver_anomaly)) < 4) {
    return(tibble(driver = driver, scale = "within_site_monthly", n = nrow(driver.df),
                  spearman_rho = NA_real_, p.value = NA_real_, deviance_explained = NA_real_,
                  delta_AIC = NA_real_))
  }
  
  fit <- mgcv::gam(
    DIEL_anomaly ~ s(driver_anomaly, k = 5) + Season,
    data = driver.df,
    method = "REML"
  )
  null_fit <- mgcv::gam(DIEL_anomaly ~ Season, data = driver.df, method = "REML")
  test <- suppressWarnings(cor.test(driver.df$driver_anomaly, driver.df$DIEL_anomaly, method = "spearman"))
  
  tibble(
    driver = driver,
    scale = "within_site_monthly",
    n = nrow(driver.df),
    spearman_rho = unname(test$estimate),
    p.value = test$p.value,
    deviance_explained = summary(fit)$dev.expl,
    delta_AIC = AIC(null_fit) - AIC(fit)
  )
})

plot.DIEL.monthly.drivers <- DIEL.monthly.driver.data %>%
  dplyr::select(SITE_ID, Season, DIEL_anomaly, dplyr::all_of(paste0(monthly.driver.variables, "_anomaly"))) %>%
  pivot_longer(
    cols = dplyr::all_of(paste0(monthly.driver.variables, "_anomaly")),
    names_to = "driver",
    values_to = "driver_anomaly"
  ) %>%
  filter(is.finite(DIEL_anomaly), is.finite(driver_anomaly)) %>%
  mutate(driver = str_remove(driver, "_anomaly$")) %>%
  left_join(DIEL.monthly.driver.summary, by = "driver") %>%
  mutate(driver_label = paste0(driver, "\nwithin-site rho = ", round(spearman_rho, 2))) %>%
  ggplot(aes(x = driver_anomaly, y = DIEL_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_point(aes(color = Season), alpha = 0.4, size = 1.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), color = "black", linewidth = 0.8, se = TRUE) +
  facet_wrap(~driver_label, scales = "free_x", ncol = 3) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
  labs(
    x = "Within-site driver anomaly",
    y = expression(paste("Within-site CH"[4], " anomaly (mg C ", m^-2, " ", day^-1, ")")),
    color = "Season"
  ) +
  theme(strip.background = element_rect(fill = "transparent", color = "black"))

ggsave("FIGURES/DIEL_monthly_driver_anomalies.png", plot = plot.DIEL.monthly.drivers, width = 9, height = 4.6, units = "in")

# Site-level variables explain across-site differences, so evaluate them on site means.
site.driver.variables <- c(
  'sulfurTot', 'dryMass', 'acidity',
  'FRC.kg.C.m.2', 'SOC.kg.C.m.2',
  'beta_roots', 'beta_soc',
  'MAP', 'MAT'
) %>% intersect(names(parms.site))

DIEL.site.driver.data <- parms.site %>%
  dplyr::select(SITE_ID, EcoType, DIEL, dplyr::all_of(site.driver.variables)) %>%
  mutate(DIEL_mgC_m2_day = DIEL * 1000) %>%
  filter(is.finite(DIEL_mgC_m2_day))

DIEL.site.driver.summary <- purrr::map_dfr(site.driver.variables, function(driver) {
  driver.df <- DIEL.site.driver.data %>%
    dplyr::select(SITE_ID, EcoType, DIEL_mgC_m2_day, driver_value = dplyr::all_of(driver)) %>%
    drop_na() %>%
    filter(is.finite(DIEL_mgC_m2_day), is.finite(driver_value))
  
  if (nrow(driver.df) < 6 || length(unique(driver.df$driver_value)) < 4) {
    return(tibble(driver = driver, scale = "across_site", n = nrow(driver.df),
                  spearman_rho = NA_real_, p.value = NA_real_, deviance_explained = NA_real_,
                  delta_AIC = NA_real_))
  }
  
  fit <- mgcv::gam(DIEL_mgC_m2_day ~ s(driver_value, k = min(5, length(unique(driver.df$driver_value)) - 1)),
                   data = driver.df, method = "REML")
  null_fit <- lm(DIEL_mgC_m2_day ~ 1, data = driver.df)
  test <- suppressWarnings(cor.test(driver.df$driver_value, driver.df$DIEL_mgC_m2_day, method = "spearman"))
  
  tibble(
    driver = driver,
    scale = "across_site",
    n = nrow(driver.df),
    spearman_rho = unname(test$estimate),
    p.value = test$p.value,
    deviance_explained = summary(fit)$dev.expl,
    delta_AIC = AIC(null_fit) - AIC(fit)
  )
})

DIEL.driver.summary <- bind_rows(DIEL.monthly.driver.summary, DIEL.site.driver.summary) %>%
  arrange(desc(delta_AIC), desc(abs(spearman_rho)))

print(DIEL.driver.summary)
write.csv(DIEL.driver.summary, file = "OUTPUT/DIEL_driver_summary.csv", row.names = FALSE)

important.site.DIEL.drivers <- DIEL.site.driver.summary %>%
  filter(is.finite(delta_AIC)) %>%
  slice_max(order_by = delta_AIC, n = 6, with_ties = FALSE) %>%
  pull(driver)

plot.DIEL.site.drivers <- DIEL.site.driver.data %>%
  dplyr::select(SITE_ID, EcoType, DIEL_mgC_m2_day, dplyr::all_of(important.site.DIEL.drivers)) %>%
  pivot_longer(cols = dplyr::all_of(important.site.DIEL.drivers), names_to = "driver", values_to = "driver_value") %>%
  filter(is.finite(DIEL_mgC_m2_day), is.finite(driver_value)) %>%
  left_join(DIEL.site.driver.summary, by = "driver") %>%
  mutate(driver_label = paste0(driver, "\nsite rho = ", round(spearman_rho, 2))) %>%
  ggplot(aes(x = driver_value, y = DIEL_mgC_m2_day)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_point(aes(color = EcoType), alpha = 0.75, size = 2.2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), color = "black", linewidth = 0.8, se = TRUE) +
  facet_wrap(~driver_label, scales = "free_x", ncol = 3) +
  theme_bw() +
  labs(
    x = "Site-level driver value",
    y = expression(paste("Mean CH"[4], " as DIEL (mg C ", m^-2, " ", day^-1, ")")),
    color = "EcoType"
  ) +
  theme(strip.background = element_rect(fill = "transparent", color = "black"))

ggsave("FIGURES/DIEL_site_driver_relationships.png", plot = plot.DIEL.site.drivers, width = 9, height = 6.5, units = "in")

plot.DIEL.driver.evidence <- DIEL.driver.summary %>%
  mutate(driver = factor(driver, levels = rev(driver))) %>%
  ggplot(aes(x = driver, y = delta_AIC, fill = scale)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(x = "", y = "AIC improvement over null model", fill = "Analysis scale")

ggsave("FIGURES/DIEL_driver_model_evidence.png", plot = plot.DIEL.driver.evidence, width = 7, height = 5.5, units = "in")

# Site consistency as CH4 source, sink, or fluctuating: ####

site.consistency.threshold <- 0.75

DIEL.site.consistency <- parms.diel %>%
  dplyr::select(SITE_ID, Date, Season, DIEL_mgC_m2_day) %>%
  arrange(SITE_ID, Date) %>%
  group_by(SITE_ID) %>%
  mutate(
    source_month = DIEL_mgC_m2_day > 0,
    sign_change = source_month != lag(source_month)
  ) %>%
  reframe(
    n_months = sum(!is.na(DIEL_mgC_m2_day)),
    prop_source_months = mean(source_month, na.rm = TRUE),
    prop_sink_months = 1 - prop_source_months,
    mean_DIEL_mgC_m2_day = mean(DIEL_mgC_m2_day, na.rm = TRUE),
    median_DIEL_mgC_m2_day = median(DIEL_mgC_m2_day, na.rm = TRUE),
    sd_DIEL_mgC_m2_day = sd(DIEL_mgC_m2_day, na.rm = TRUE),
    sign_changes = sum(sign_change, na.rm = TRUE)
  ) %>%
  mutate(
    CH4_behavior = case_when(
      prop_source_months >= site.consistency.threshold ~ "Consistent source",
      prop_source_months <= 1 - site.consistency.threshold ~ "Consistent sink",
      TRUE ~ "Fluctuating"
    ),
    CH4_behavior = factor(CH4_behavior, levels = c("Consistent sink", "Fluctuating", "Consistent source"))
  ) %>%
  left_join(
    parms.site %>%
      dplyr::select(SITE_ID, EcoType, dplyr::all_of(site.driver.variables)),
    by = "SITE_ID"
  )

write.csv(DIEL.site.consistency, file = "OUTPUT/DIEL_site_consistency.csv", row.names = FALSE)

DIEL.site.behavior.tests <- purrr::map_dfr(site.driver.variables, function(driver) {
  test.df <- DIEL.site.consistency %>%
    dplyr::select(CH4_behavior, value = dplyr::all_of(driver)) %>%
    drop_na() %>%
    filter(is.finite(value))
  
  if (nrow(test.df) < 6 || n_distinct(test.df$CH4_behavior) < 2 || n_distinct(test.df$value) < 4) {
    return(tibble(driver = driver, n = nrow(test.df), p.value = NA_real_))
  }
  
  test <- kruskal.test(value ~ CH4_behavior, data = test.df)
  
  tibble(
    driver = driver,
    n = nrow(test.df),
    p.value = test$p.value
  )
}) %>%
  arrange(p.value)

DIEL.site.behavior.medians <- DIEL.site.consistency %>%
  dplyr::select(CH4_behavior, dplyr::all_of(site.driver.variables)) %>%
  pivot_longer(cols = dplyr::all_of(site.driver.variables), names_to = "driver", values_to = "value") %>%
  filter(is.finite(value)) %>%
  reframe(
    .by = c(driver, CH4_behavior),
    n = n(),
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE)
  )

DIEL.site.behavior.summary <- DIEL.site.behavior.medians %>%
  left_join(DIEL.site.behavior.tests, by = "driver") %>%
  arrange(p.value, driver, CH4_behavior)

print(DIEL.site.consistency %>% count(CH4_behavior))
print(DIEL.site.behavior.summary)

write.csv(DIEL.site.behavior.summary, file = "OUTPUT/DIEL_site_behavior_driver_summary.csv", row.names = FALSE)

important.site.behavior.drivers <- DIEL.site.behavior.tests %>%
  filter(!is.na(p.value)) %>%
  slice_min(order_by = p.value, n = 6, with_ties = FALSE) %>%
  pull(driver)

plot.DIEL.site.behavior.map <- ggplot() +
  geom_sf(data = site.att.sf, fill = "grey90", color = "white") +
  geom_sf(
    data = site.att.sf %>% left_join(DIEL.site.consistency, by = "SITE_ID"),
    aes(color = CH4_behavior, size = abs(mean_DIEL_mgC_m2_day)),
    alpha = 0.8
  ) +
  theme_bw() +
  scale_color_manual(values = c("Consistent sink" = "red3", "Fluctuating" = "grey35", "Consistent source" = "blue4"), na.translate = FALSE) +
  scale_size_continuous(range = c(2, 6)) +
  labs(color = "CH4 behavior", size = expression(paste("|mean CH"[4], "|")))

ggsave("FIGURES/DIEL_site_behavior_map.png", plot = plot.DIEL.site.behavior.map, width = 7, height = 5, units = "in")

plot.DIEL.site.behavior.phase <- DIEL.site.consistency %>%
  ggplot(aes(x = prop_source_months, y = mean_DIEL_mgC_m2_day)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_vline(xintercept = c(1 - site.consistency.threshold, site.consistency.threshold), linetype = "dotted", color = "grey45") +
  geom_point(aes(color = CH4_behavior, size = n_months), alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = SITE_ID), size = 2.6, max.overlaps = 30) +
  theme_bw() +
  scale_color_manual(values = c("Consistent sink" = "red3", "Fluctuating" = "grey35", "Consistent source" = "blue4"), na.translate = FALSE) +
  labs(
    x = "Fraction of months with CH4 source behavior",
    y = expression(paste("Mean CH"[4], " as DIEL (mg C ", m^-2, " ", day^-1, ")")),
    color = "CH4 behavior",
    size = "Months"
  )

ggsave("FIGURES/DIEL_site_behavior_phase_space.png", plot = plot.DIEL.site.behavior.phase, width = 8, height = 6, units = "in")

plot.DIEL.site.behavior.drivers <- DIEL.site.consistency %>%
  dplyr::select(SITE_ID, CH4_behavior, dplyr::all_of(important.site.behavior.drivers)) %>%
  pivot_longer(cols = dplyr::all_of(important.site.behavior.drivers), names_to = "driver", values_to = "value") %>%
  filter(is.finite(value)) %>%
  left_join(DIEL.site.behavior.tests, by = "driver") %>%
  mutate(driver_label = paste0(driver, "\nKruskal p = ", signif(p.value, 2))) %>%
  ggplot(aes(x = CH4_behavior, y = value, color = CH4_behavior)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.65, size = 1.8) +
  facet_wrap(~driver_label, scales = "free_y", ncol = 3) +
  theme_bw() +
  scale_color_manual(values = c("Consistent sink" = "red3", "Fluctuating" = "grey35", "Consistent source" = "blue4"), na.translate = FALSE) +
  labs(x = "", y = "Site-level driver value", color = "CH4 behavior") +
  theme(
    strip.background = element_rect(fill = "transparent", color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave("FIGURES/DIEL_site_behavior_driver_boxplots.png", plot = plot.DIEL.site.behavior.drivers, width = 10, height = 6.5, units = "in")



# Figures: ####

# Seasonal methane pattern ####


NorthAmerica_AOI <- aoi.usa <- aoi_get(country = c('North America') )


map.total <- ggplot( ) + geom_sf(data= NorthAmerica_AOI, fill = "black")+
  geom_sf( data= parms.site.season.sf %>% arrange(desc(DIEL*1000)) ,
           aes( size = abs(DIEL*1000),
                col = DIEL*1000), alpha = 0.75) + theme_bw() +
  scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", 
                          direction=1,
                          limit = c(-5, 5),
                          oob = scales::squish) + 
  facet_wrap( ~ Season) +                
  scale_size_continuous(range = c(1, 5)) +
  theme(strip.background = element_rect(fill = "transparent")) + 
  labs(col=expression(paste("FCH"[4], "( mg C m"^-2, "day"^-1, ")")),
       size =expression(paste("|FCH"[4], "| ( mg C m"^-2, "day"^-1, ")")))

  density.total <-  parms.site.season.sf  %>% 
  ggplot(aes(DIEL*1000, fill = Season)) +
  geom_density( position = "stack", aes(col=Season), alpha = 0.5) + theme_bw()+
    scale_color_paletteer_d("ggthemes::Seattle_Grays") +
    scale_fill_paletteer_d("ggthemes::Seattle_Grays") +
    geom_vline( xintercept=0, linetype = "dashed")+ xlim( -20, 20) + 
    ylab("Density") + 
    xlab( expression(paste("Mean Seasonal FCH"[4], "( mg C m"^-2, "day"^-1, ")")))
      
  total.season.final.plot <- ggarrange( map.total , 
                                        density.total, 
                                        labels=c("A", "B"), ncol=1,
                                        heights= c( 2, 1))
  
ggsave("FIGURES/Total_methane_season_plot.png", 
       plot =  total.season.final.plot, width = 7.6, 
       height = 7.3, units = "in")
  
# Q10 and base respiration by season and ecosystem type: ####
  
  source_colors <- c("Sink" = "red3", "Source" = "blue4")
  
 Q10.source.density.plots <-  parms.site.season  %>%   filter(!is.na(Source)) %>% 
    ggplot(aes( x = Q10.mean, y = after_stat(scaled), col= Source)) +
    geom_density(alpha=0.5)  +
    theme_bw() + ylab("Density") + xlab("Q10") + scale_color_manual(values = source_colors)
 
 
 Q10.source.box.plots <-  parms.site.season  %>%   filter(!is.na(Source)) %>% 
    ggplot(aes( x = Source, y = Q10.mean, col=Source)) +
    geom_boxplot(alpha=0.5)  + ylim(1,2) +
    theme_bw()  + stat_summary(
      fun = median, 
      geom = "text", 
      col = source_colors ,     # Color of the text
      vjust = -10, aes(label = round(..y.., digits = 2)) ) +
    scale_color_manual(values = source_colors)+ ylab("Q10") + 
    xlab("") +  theme(legend.position = "none")
  
 Q10.source.box.seasons.plots <- parms.site.season  %>%   filter(!is.na(Source)) %>% 
    ggplot(aes( x = Source, y = Q10.mean, col = Source)) +
    geom_boxplot(alpha=0.5)  + ylim(1,2) +
    theme_bw()  + stat_summary(
      fun = median, 
      geom = "text",     # Color of the text
      vjust = -12,
      size=3,
      aes(label = round(..y.., digits = 2)) )+ facet_wrap( ~ Season, ncol=4) + scale_color_manual(values = source_colors)+
   ylab("Q10") + xlab("") + theme(strip.background = element_rect(fill = "transparent", color = "black"),
     legend.position = "none")
  
  
DIEL.source.box.plots <- parms.site.season  %>%   filter(!is.na(Source)) %>% 
    ggplot(aes( x = Source, y = DIEL*1000, col = Source)) +
    geom_boxplot(alpha=0.5)  + ylim(-20,20) +
    theme_bw() +scale_color_manual(values = source_colors)+ 
  theme(strip.background = element_rect(fill = "transparent", color = "black"), legend.position = "none") +
  ylab( expression(paste("FCH"[4], " ( mg C m"^-2, "day"^-1, ")")))
  
  
DIEL.source.box.seasons.plots <- parms.site.season  %>%   filter(!is.na(Source)) %>% 
    ggplot(aes( x = Source, y = DIEL*1000, col = Source)) +
    geom_boxplot(alpha=0.5)  + ylim(-20,20) +
    theme_bw() + facet_wrap( ~ Season, ncol=4)+scale_color_manual(values = source_colors)+
  theme(strip.background = element_rect(fill = "transparent", color = "black"),
        legend.position = "none") + ylab( expression(paste("FCH"[4], " ( mg C m"^-2, "day"^-1, ")")))
  
DIEL.plots <-ggarrange(DIEL.source.box.plots,
DIEL.source.box.seasons.plots, ncol=1, labels= c("A", "B") )
  

Q10.plots <- ggarrange(
  ggarrange(Q10.source.box.plots,Q10.source.density.plots, ncol=2, labels= c("A", "B"), common.legend = TRUE),
ggarrange (Q10.source.box.seasons.plots,  labels= c("C") ), ncol=1)
  
  
ggsave("FIGURES/Q10_plots.png", plot =   Q10.plots, width =6, height =5, units = "in")

ggsave("FIGURES/DIEL_plot.png", plot =   DIEL.plots , width = 5, height = 6, units = "in")

  # Total Methane Figures Annual: ####
  library(ggpmisc)
library(scales)

map.DIEL.total <- ggplot( ) + geom_sf(data= NorthAmerica_AOI, fill = "black")+
  geom_sf( data= parms.site.sf %>% arrange(desc(DIEL )) ,
           aes( size = abs(DIEL*1000),
                col =DIEL*1000), alpha = 0.75) + theme_bw() +
  scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", 
                          direction=1,
                          limit = c(-10, 10),
                          oob = scales::squish,
                          breaks=pretty_breaks(n = 3)) +                
  scale_size_continuous(range = c(3, 5)) +
  theme(strip.background = element_rect(fill = "transparent")) + 
  labs(col=expression(paste("FCH"[4], "( mg C m"^-2, "day"^-1, ")")),
       size =expression(paste("|FCH"[4], "| ( mg C m"^-2, "day"^-1, ")")))

 MAT.MAP.plot <-  parms.site %>% filter(!is.na(Source)) %>% 
    ggplot(aes( x = MAT, y = MAP)) +
    geom_point(alpha=0.5, aes(col = DIEL*1000, size = DIEL*1000)) + 
    theme_bw()  +scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", 
                                         direction=1,
                                         limit = c(-10, 10),
                                         oob = scales::squish)+
   scale_size_continuous(range = c(1, 5))+
    theme(panel.background = element_rect(fill = "black"), # Plotting area
          panel.grid.major = element_line(color = "grey20"), # Optional: adjust grid visibility
      panel.grid.minor = element_line(color = "grey30") ) 
  
 DIEL.MAT.plot <-  parms.site  %>% 
   ggplot(aes( x = MAT, y = DIEL*1000)) +
   geom_point(alpha=0.75, aes(col = DIEL*1000), size = 3) + 
   theme_bw() + geom_smooth( method="loess", col='white') +
   stat_cor(method = "spearman", col="white",
            label.y.npc = "center") +
   facet_wrap( ~ Source_Sink + EcoType, ncol=3, scale="free") +
   scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", 
                           direction=1,
                           limit = c(-10, 10))  + 
   ylab( expression(paste("FCH"[4], " ( mg C m"^-2, "day"^-1, ")"))) +
   theme(panel.background = element_rect(fill = "black"), # Plotting area
         panel.grid.major = element_line(color = "grey20"), 
         panel.grid.minor = element_line(color = "grey30"),
         legend.position = "none" ,
         strip.background = element_rect(fill = "transparent", color = "black"))
 
 
  annual.ch4.plots <- ggarrange( ggarrange(   map.DIEL.total,
               MAT.MAP.plot, labels=c("A", "B") , 
               common.legend = T, widths=c(2,1.25)), 
             DIEL.MAT.plot, labels=c(" ", "C"), ncol =1)

  
  ggsave("FIGURES/Total_Annual_CH4_plot.png", plot =     
           annual.ch4.plots , width = 10, height = 11, units = "in")
  
# MDS: ####
  
continious.variables <- c( 'VSWCMean' ,  'Q10.mean' , 'Rref.mean' , 'VSWCVar' , 'sulfurTot' , 'dryMass' , 'acidity','FRC.kg.C.m.2', 'SOC.kg.C.m.2', 'beta_roots' ,'beta_soc')
  
library(MASS)
library(ggrepel)
library(corrplot)
parms.mds <- parms.site %>% drop_na()
parms.mds %>% names
M = cor(  parms.mds[, 3: 17]  %>% drop_na())
corrplot(M, method = 'color', order = 'alphabet') # colorful number

library(vegan)

d <- vegdist(  parms.mds %>% dplyr::select( continious.variables))
fit <- metaMDS(d, k=2)

parms.mds$MDS1 <- fit$points[,1]
parms.mds$MDS2<- fit$points[,2]

cluster_res <- hclust(d, method = "ward.D2")

groups <- cutree(cluster_res, k = 4)
parms.mds$cluster <- groups # Add clusters to the df

library(paletteer)

plot.MDS <-   parms.mds %>%  ggplot(aes( x=MDS1, y = MDS2, col=cluster %>% as.factor)) +
  stat_ellipse()+
  geom_text_repel(aes(label = SITE_ID),
                   box.padding   = 0, 
                   point.padding = 0,
                  size =3) +
  geom_point( alpha = 0.25 ) +  theme_bw() + 
  scale_color_paletteer_d("ggsci::nrc_npg") +theme(legend.position = "none") +
  xlab( 'Dimension 1')+ ylab( 'Dimension 2')


# Additional plots:

plot.diel <-  parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =DIEL*1000, col=cluster %>% as.factor)) + geom_boxplot( )+ theme_bw() + geom_hline(yintercept = 0, col = "black", linetype="dashed") + 
  ylab(expression(paste(" Mean FCH"[4], "( mg C m"^-2, "day"^-1, ")"))) + xlab("Cluster") + 
  scale_color_paletteer_d("ggsci::nrc_npg")+ theme(legend.position = "none") + ylim(-10, 10)

plot.source <-  parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =Source %>% as.numeric-1, col=cluster %>% as.factor)) + geom_violin( ) + theme_bw() +  stat_summary(fun = mean, geom = "crossbar", width = 0.2) + 
  scale_color_paletteer_d("ggsci::nrc_npg") + theme(legend.position = "none") + ylab( "Number of Month as a Source")+ xlab("Cluster")

plot.Q10 <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =Q10.mean , col=cluster %>% as.factor)) + geom_boxplot( ) + 
  theme_bw() + xlab("Cluster") + ylab('Q10') +
  stat_summary(fun = "median", 
               geom = "text", 
               aes(label = round(after_stat(y), 2)), 
               vjust = 1.5, # Moves text above the point
               size = 3) +  theme(legend.position = "none") +  scale_color_paletteer_d("ggsci::nrc_npg")

plot.sulfur <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =sulfurTot, col=cluster %>% as.factor)) + geom_boxplot( ) + theme_bw() + xlab("Cluster") + ylab('Soil Sulfur') +  
  scale_color_paletteer_d("ggsci::nrc_npg")+  theme(legend.position = "none") 


plot.rootmass <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =dryMass , col=cluster %>% as.factor)) + geom_boxplot( ) +  scale_color_paletteer_d("ggsci::nrc_npg")+  theme(legend.position = "none")  + xlab("Cluster") + 
  ylab('Root Dry Mass') + theme_bw() +  theme(legend.position = "none")

plot.RootMass <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =FRC.kg.C.m.2 , col=cluster %>% as.factor)) + geom_boxplot( )+  theme(legend.position = "none") +  scale_color_paletteer_d("ggsci::nrc_npg") +
  xlab("Cluster") +  ylab(expression(paste('Root Mass (kg C m'^-2,')'))) + theme_bw() +  theme(legend.position = "none")

plot.soc <-  parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =SOC.kg.C.m.2, col=cluster %>% as.factor)) + geom_boxplot( ) +  theme(legend.position = "none") +  scale_color_paletteer_d("ggsci::nrc_npg") +
  xlab("Cluster") +  ylab(expression(paste('SOC (kg C m'^-2,')'))) + theme_bw() +  theme(legend.position = "none")

plot.betaRoots <-   parms.mds%>% ggplot(aes( x=cluster %>% as.factor, y =beta_roots, col=cluster %>% as.factor)) +  
  theme(legend.position = "none") +  scale_color_paletteer_d("ggsci::nrc_npg") +
  xlab("Cluster") +  ylab(expression(paste(beta ,' Roots'))) + geom_boxplot( ) +  theme(legend.position = "none") + theme_bw() +  theme(legend.position = "none")

plot.betaSOC <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =beta_soc, col=cluster %>% as.factor)) +  
  theme(legend.position = "none") +  scale_color_paletteer_d("ggsci::nrc_npg") +
  xlab("Cluster") +  ylab(expression(paste(beta ,' SOC'))) + geom_boxplot( ) +  theme(legend.position = "none") + theme_bw() +  theme(legend.position = "none")

plot.soilAcidity <- parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =acidity, col=cluster %>% as.factor)) + geom_boxplot( ) +  scale_color_paletteer_d("ggsci::nrc_npg") + theme_bw()+  theme(legend.position = "none") + ylab("Soil Acidity") + xlab('Cluster')

plot.VSWCVAR <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =VSWCVar, col=cluster %>% as.factor)) + geom_boxplot( ) +  scale_color_paletteer_d("ggsci::nrc_npg") +
  ylab("SVWC Variability") + xlab('Cluster') + theme_bw()+  theme(legend.position = "none")

plot.VSWCMean <-   parms.mds %>% ggplot(aes( x=cluster %>% as.factor, y =VSWCMean, col=cluster %>% as.factor)) + geom_boxplot( ) +  scale_color_paletteer_d("ggsci::nrc_npg") +
  ylab("Mean SVWC") + xlab('Cluster') + theme_bw()+  theme(legend.position = "none") 

# Maps: ####
site.att.sf.cluster <- site.att.sf %>% left_join(by='SITE_ID',   parms.mds )

map.cluster <- ggplot( ) + geom_sf(data= NorthAmerica_AOI, fill = "white")+
  geom_sf( data= site.att.sf.cluster ,
           aes( col = cluster %>% as.factor), alpha = 0.5, size = 4) + theme_bw() + scale_color_paletteer_d("ggsci::nrc_npg",na.translate = FALSE)  + theme(strip.background = element_rect(fill = "transparent"))+ labs(col="Cluster")

# Panels ####

# Patterns in FCH4


panel.methane1 <- ggarrange( map.cluster, plot.MDS ,
                             labels=c("A", "B"), nrow = 2, ncol =1, common.legend = TRUE)

panel.methane2 <- ggarrange( plot.diel, plot.source, plot.Q10 ,
                             labels=c("C", "D", "E"), nrow = 3, ncol =1)

final.MDS.plot <- ggarrange( panel.methane1, panel.methane2,
           widths = c(2,1), common.legend = TRUE)

ggsave("FIGURES/MDS_PLOT.png", plot =final.MDS.plot , width = 8, height = 8.7, units = "in")

# Patterns in Soil Conditions
panel.root_soil <- ggarrange(plot.RootMass,
                             plot.betaRoots,
                             plot.soc,
                             plot.sulfur ,
                             plot.soilAcidity ,
                             plot.VSWCMean ,
                             plot.VSWCVAR,
                             plot.betaSOC,
                             labels = c( "A", "B", "C", "D", "E", "F", "G", "H"))

ggsave("FIGURES/Roots_Soils_PLOT.png", plot =panel.root_soil , width = 8, height = 8, units = "in")

# RESULTS: ####

total.daily <- parms.site.season$DIEL*1000 
total.daily%>% summary

total.daily.annual <- parms.site$DIEL*1000 
total.daily.annual%>% summary


nrow( parms.site %>%  filter( Source_Sink == "Source"))
nrow( parms.site %>%  filter( Source_Sink == "Sink"))

parms.site$SITE_ID[ parms.site$Source == 4] %>% unique %>% length
parms.site$SITE_ID[ parms.site$Source == 0] %>% unique %>% length


# Measure the significance of difference between sinks and sources for Q10.
kruskal.test(Q10.mean ~ Source_Sink , data =parms.mds) 

kruskal.test(Q10.mean ~ Season , data =parms.site.season) 

kruskal.test(Q10.mean ~ Source, data = parms.site.season %>% filter(Season == "Winter")) 

# Explore differences in the clusters:

parms.mds.cluster1 <- parms.mds %>% filter(cluster == 1)
mean(parms.mds.cluster1$DIEL *1000)
sd(parms.mds.cluster1$DIEL *1000)

mean(parms.mds.cluster1$Q10.mean)
sd(parms.mds.cluster1$Q10.mean)

mean(parms.mds.cluster1$sulfurTot)
sd(parms.mds.cluster1$sulfurTot)

mean(parms.mds.cluster1$Q10.mean)
sd(parms.mds.cluster1$Q10.mean)

mean(parms.mds.cluster1$acidity)
sd(parms.mds.cluster1$acidity)

mean(parms.mds.cluster1$FRC.kg.C.m.2)
sd(parms.mds.cluster1$FRC.kg.C.m.2)


parms.mds.cluster2 <- parms.mds %>% filter(cluster == 2)
mean(parms.mds.cluster2$Q10.mean)
sd(parms.mds.cluster2$Q10.mean)

mean(parms.mds.cluster2$DIEL *1000)
sd(parms.mds.cluster2$DIEL*1000)

mean(parms.mds.cluster2$sulfurTot)
sd(parms.mds.cluster2$sulfurTot)

mean(parms.mds.cluster2$acidity)
sd(parms.mds.cluster2$acidity)

mean(parms.mds.cluster2$SOC.kg.C.m.2)
sd(parms.mds.cluster2$SOC.kg.C.m.2)

mean(parms.mds.cluster2$VSWCVar)
sd(parms.mds.cluster2$VSWCVar)

parms.mds.cluster3 <- parms.mds %>% filter(cluster == 3)

mean(parms.mds.cluster3$DIEL*1000)
sd(parms.mds.cluster3$DIEL*1000)

mean(parms.mds.cluster3$sulfurTot)
sd(parms.mds.cluster3$sulfurTot)

mean(parms.mds.cluster3$acidity)
sd(parms.mds.cluster3$acidity)

mean(parms.mds.cluster3$SOC.kg.C.m.2)
sd(parms.mds.cluster3$SOC.kg.C.m.2)

mean(parms.mds.cluster3$VSWCMean)
sd(parms.mds.cluster3$VSWCMean)

parms.mds.cluster4 <- parms.mds %>% filter(cluster == 4)

mean(parms.mds.cluster4$DIEL*1000)
sd(parms.mds.cluster4$DIEL*1000)

mean(parms.mds.cluster4$sulfurTot)
sd(parms.mds.cluster4$sulfurTot)

mean(parms.mds.cluster4$VSWCMean)
sd(parms.mds.cluster4$VSWCMean)

mean(parms.mds.cluster4$acidity)
sd(parms.mds.cluster4$acidity)

mean(parms.mds.cluster4$VSWCVar)
sd(parms.mds.cluster4$VSWCVar)

mean(parms.mds.cluster4$beta_soc)
sd(parms.mds.cluster4$beta_soc)

summary.4 <- parms.mds.cluster4$DIEL*1000
summary.4 %>% summary
parms.mds %>% reframe( .by=cluster, 
                       DIEL = median(DIEL)*1000 %>% round(4),
                       DIEL.SD = sd(DIEL)*1000 %>% round(4))

parms.mds %>% names
parms.mds %>% ggplot( aes( x=VSWCMin, y=Q10.mean, col = cluster %>% as.factor)) + geom_smooth(method = 'loess') +
  facet_wrap( ~cluster)
randomForest::randomForest(cluster ~  Q10.mean + acidity + sulfurTot + DIEL+ MAT+ MAP ,
                           data=parms.mds %>% mutate(cluster = as.factor(cluster)))

randomForest::randomForest(cluster ~  Q10.mean + acidity + sulfurTot + DIEL+ MAT+ MAP +EcoType +beta_roots+ beta_soc,
  data=parms.mds %>% mutate(cluster = as.factor(cluster))) 

write.csv( parms.mds, file = "parms_mds.csv")

# Copilot:

# Cluster Analysis of Q10, Rref, and DIEL

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Simulating some data
# In practice, you would want to load your actual dataset here
set.seed(123)
clusters <- rep(1:5, each = 20)
Q10 <- rnorm(100, mean=2, sd=0.5)
Rref <- rnorm(100, mean=5, sd=1.5)
DIEL <- rnorm(100, mean=10, sd=2)

# Creating a data frame
data <- parms.mds

# 1. Q10 by cluster
q10_plot <- ggplot(data, aes(x=factor(cluster), y=Q10.mean)) + 
  geom_boxplot() + 
  labs(title='Q10 by Cluster', x='Cluster', y='Q10')
print(q10_plot)

# 2. Rref by cluster
rref_plot <- ggplot(data, aes(x=factor(cluster), y=Rref.mean)) + 
  geom_boxplot() + 
  labs(title='Rref by Cluster', x='Cluster', y='Rref')
print(rref_plot)

# 3. DIEL by cluster
diel_plot <- ggplot(data, aes(x=factor(cluster), y=DIEL*1000)) + 
  geom_boxplot() + 
  labs(title='DIEL by Cluster', x='Cluster', y='DIEL')
print(diel_plot) + ylim(-20, 20)

# 4. Scatter plot of Q10 vs DIEL
q10_diel_plot <- ggplot(data, aes(x=Q10.mean, y=DIEL*1000, col=factor(cluster))) + 
  geom_point(aes(color=factor(cluster))) + 
  labs(title='Q10 vs DIEL', x='Q10', y='DIEL') + ylim( -20,20) + geom_smooth(method="lm")

print(q10_diel_plot) 

# 5. Scatter plot of Rref vs DIEL
rref_diel_plot <- ggplot(data, aes(x=Rref.mean, y=DIEL)) + 
  geom_point(aes(color=factor(cluster))) + 
  labs(title='Rref vs DIEL', x='Rref', y='DIEL')
print(rref_diel_plot)

# 6. Relationship between Q10/Rref and DIEL patterns
combined_plot <- ggplot(data, aes(x=Q10.mean/Rref.mean, y=DIEL*1000, color=factor(cluster))) + 
  geom_point(aes(color=factor(cluster))) + 
  labs(title='Relationship between Q10/Rref and DIEL', x='Q10/Rref', y='DIEL')+ ylim( -20,20) + geom_smooth(method="lm")
print(combined_plot)

                       
