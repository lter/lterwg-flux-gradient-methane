
library(tidyverse)
library(sf)
library(AOI)
library(ggplot2)
library(colorspace)

load(  fs::path(localdir,paste0("SITEval_DATA_FILTERED_CH4.Rdata")))

# See what is in the attribute table for the datasets:

# Location and IGBP are needed
meta.val <- read.csv( fs::path('/Volumes/MaloneLab/Research/FluxGradient/',paste0("/metadata_validation.csv"))) %>% rename(Tower_Measurement_HT = Measurement_HT)

val.sites.shp <- meta.val %>% st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
                                       crs = 4326) 

Study_AOI <- aoi.usa <- aoi_get(country = c('North America', 'Europe'))
aoi.usa <- aoi_get(country = c('PR', 'USA'))


ggplot() + geom_sf(data = Study_AOI) + 
  geom_sf(data = val.sites.shp) + 
  theme_bw() +theme(legend.position="top")


# Information from the attribute tables:
adjust.attr <- function(attr){
  attr2 <- attr %>% t() %>% as.data.frame()
  colnames(attr2) <- attr2[1,]
  attr3 <- attr2[-c(1:2), c('DistZaxsLvlMeasTow', 'DistZaxsTow', 'ElevRefeTow',
                            'LatTow' , 'LonTow', 'LvlMeasTow', 'TowerPosition', 'Site')] %>% 
    rename( Mes_Ht = DistZaxsLvlMeasTow,
            Elevation = ElevRefeTow,
            Lat = LatTow,
            Long =LonTow )
  
  
  return(attr3)
  
}

measLevel <- function( attr){
  sub1 <- attr %>% select(Mes_Ht, TowerPosition)
  sub2 <- sub1 %>% cross_join( sub1 %>% rename(Mes_Ht_B = Mes_Ht,TowerPosition_B =TowerPosition )) %>% 
    mutate( Mes_Dist = Mes_Ht %>% as.numeric - Mes_Ht_B%>% as.numeric,
            dLevelsAminusB =  paste(TowerPosition, TowerPosition_B, sep="_")) %>% filter( Mes_Dist >0) %>% select ( Mes_Dist , dLevelsAminusB, Mes_Ht, Mes_Ht_B)
}

usuaf.att <- read.csv(fs::path(localdir,paste0("/US-Uaf/usuaf_attr.csv"))) %>% adjust.attr () %>% measLevel %>% mutate( SITE_ID = "US-Uaf")
sesvb.att <- read.csv(fs::path(localdir,paste0("/SE-Svb/sesvb_attr.csv")))%>% adjust.attr () %>% measLevel %>% mutate( SITE_ID = "SE-Svb")
sesto.att <- read.csv(fs::path(localdir,paste0("/SE-Sto/sesto_attr.csv")))%>% adjust.attr () %>% measLevel %>%mutate( SITE_ID = "SE-Sto")


attributes_val <- rbind( usuaf.att, sesvb.att, sesto.att) %>% left_join(val.sites.shp , by="SITE_ID")%>% mutate( canopy_level_A = case_when(Mes_Ht > CHM ~ "A",
                                                      Mes_Ht <= CHM ~ "W"),
                           canopy_level_B = case_when(Mes_Ht_B > CHM ~ "A",
                                                      Mes_Ht_B <= CHM ~ "W"),
                           Canopy_L1 = paste( canopy_level_A , canopy_level_B, sep=""))


# Save Canopy Information:
canopy.info <- attributes_val %>% select( SITE_ID,  IGBP, dLevelsAminusB, Mes_Dist, CHM, Tower_Measurement_HT, Canopy_L1)
  
write.csv(canopy.info, file.path(paste(localdir, "Val_canopy.csv", sep="/")))
# View site attributes and create a site visualization:

library(tidyverse)
library(sf)
library(AOI)

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
                              Vegetation.Abbreviation..IGBP. == 'OSH' ~ 'Shrubland')) %>% rename( Site = Site_Id.NEON)

site.att.sf <- st_as_sf(x = metadata,                         
                        coords = c("Longitude..degrees.", "Latitude..degrees."),
                        crs = 4326) 

site.att.sf %>% filter(EcoType == "Cropland")

summary.igbp <-  metadata %>% reframe( .by = c(Vegetation.Abbreviation..IGBP., EcoType), 
                                       towers = length(Vegetation.Abbreviation..IGBP.))  %>% 
  reframe( .by = EcoType, towers = sum(towers))


canopy <- read.csv(file.path(paste(localdir, "canopy_commbined.csv", sep="/"))) %>% distinct

canopy.Ht <- canopy %>% reframe( .by= Site,
                                 canopyHeight_m = mean(canopyHeight_m))

canopy.Ht$canopyHeight_m %>% range

plot.tower.counts <- ggplot(data=site.att.sf, aes(x=EcoType , col=EcoType, fill=EcoType)) +
  geom_bar(stat="count", width=0.7) + theme_bw() + ylab('NEON Towers (Counts)') +
  xlab('') + scale_y_continuous(breaks = scales::breaks_pretty(n = 5)) +
  scale_color_discrete_sequential(palette = "Hawaii") +
  scale_fill_discrete_sequential(palette = "Hawaii") +
  theme(legend.position = "none",text = element_text(size = 20))+ theme(
    panel.background = element_rect(fill = "black", colour = "black"), # Black panel background
    plot.background = element_rect(fill = "black", colour = "black"),   # Black overall plot background
    panel.grid.major = element_line(colour = "grey30"),                 # Darker grey major grid lines
    panel.grid.minor = element_line(colour = "grey10"),                 # Even darker grey minor grid lines
    axis.text = element_text(colour = "white"),                         # White axis text
    axis.title = element_text(colour = "white")                         # White axis titles
  )

site.att.sf.ht <- site.att.sf %>% left_join(canopy.Ht, by = 'Site')

World_AOI <- aoi.usa <- aoi_get(country = c('North America', 'South America', "Asia", "Africa", "Austrailia", "Europe") )

library(colorspace)

map <- ggplot() + geom_sf(data = World_AOI, col='white', fill="black") + 
  geom_sf(data = val.sites.shp %>% filter(SITE_ID != 'FI-Hyy'), fill='transparent', alpha=0.75, col="red", size=3.5) + 
  geom_sf(data = site.att.sf.ht, size=2.5, alpha =0.85, aes(col = EcoType)) +
  theme_bw()  + coord_sf(xlim = c(-160, 30), ylim = c(20, 75))+
  scale_color_discrete_sequential(palette = "Hawaii", name="") + theme(text = element_text(size = 20))

# See palattes:
hcl_palettes(type = "sequential")

library(ggpubr)

final.plot <- ggarrange(map , plot.tower.counts, ncol=1, common.legend = TRUE)

ggsave("FIGURES/Map_plot.png", plot = final.plot, width = 7.6, height = 7.3, units = "in")

