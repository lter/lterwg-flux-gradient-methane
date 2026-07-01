
library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(colorspace)

load( fs::path(localdir.ch4 ,paste0("SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE_TotalFlux.Rdata")) )

# See what is in the attribute table for the datasets:

# Location and IGBP are needed
meta.val <- read.csv( fs::path('/Volumes/MaloneLab/Research/FluxGradient/',paste0("/metadata_validation.csv"))) %>% rename(Tower_Measurement_HT = Measurement_HT)

val.sites.shp <- meta.val %>% st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
                                       crs = 4326) 

Study_AOI <- rnaturalearth::ne_countries(
  continent = c("North America"),
  scale = 50,
  returnclass = "sf")

aoi.usa <- rnaturalearth::ne_countries(
  country = c("Puerto Rico", "United States of America"),
  scale = 50,
  returnclass = "sf")


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

usuaf.att <- read.csv(fs::path(localdir,paste0("/Validation_Sites/US-Uaf/usuaf_attr.csv"))) %>% adjust.attr () %>% measLevel %>% mutate( SITE_ID = "US-Uaf")
sesvb.att <- read.csv(fs::path(localdir,paste0("/Validation_Sites/SE-Svb/sesvb_attr.csv")))%>% adjust.attr () %>% measLevel %>% mutate( SITE_ID = "SE-Svb")
sesto.att <- read.csv(fs::path(localdir,paste0("/Validation_Sites/SE-Sto/sesto_attr.csv")))%>% adjust.attr () %>% measLevel %>%mutate( SITE_ID = "SE-Sto")


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
library(rnaturalearth)

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

World_AOI <- rnaturalearth::ne_countries(
  continent = c("North America", "South America", "Asia", "Africa", "Europe", "Oceania"),
  scale = 50,
  returnclass = "sf"
)

library(colorspace)

map <- ggplot() + geom_sf(data = World_AOI, col='white', fill="black") +
  #geom_sf(data = val.sites.shp %>% filter(SITE_ID != 'FI-Hyy'), fill='transparent', alpha=0.75, col="red", size=3.5) +
  geom_sf(data = site.att.sf.ht, size=2.5, alpha =0.85, aes(col = EcoType)) +
  theme_bw()  + coord_sf(xlim = c(-160, -40), ylim = c(10, 75))+
  scale_color_discrete_sequential(palette = "Hawaii", name="") +
  labs(tag = "A") +
  theme(text = element_text(size = 20),
        plot.tag = element_text(color = "black", face = "bold", size = 14))
map
# See palattes:
hcl_palettes(type = "sequential")

library(ggpubr)

# ── Panel C: MAP vs MAT climate space with Aridity Index ──────────────────────
# De Martonne Aridity Index: AI = MAP / (MAT + 10)
# AI < 20  → Arid   (aligns with Koeppen B-class sites)
# AI 20-30 → Semi-arid
# AI ≥ 30  → Humid
#
# Boundary curve in MAP/MAT space: MAP = 20 * (MAT + 10)

climate.data <- metadata %>%
  rename(
    MAP_mm  = Mean.Average.Precipitation..mm.,
    MAT_C   = Mean.Average.Tempurature..degrees.C.,
    Koeppen = Climate.Class.Abbreviation..Koeppen.
  ) %>%
  mutate(
    MAP_mm  = as.numeric(MAP_mm),
    MAT_C   = as.numeric(MAT_C),
    # de Martonne AI is undefined when MAT <= -10 (denominator <= 0); set to NA
    AI      = if_else(MAT_C + 10 > 0, MAP_mm / (MAT_C + 10), NA_real_),
    Aridity = case_when(
      !is.na(AI) & AI <= 15 ~ "Arid",
      AI <  30              ~ "Semi-arid",
      TRUE                  ~ "Humid"
    ),
    Aridity = factor(Aridity, levels = c("Arid", "Semi-arid", "Humid"))
  ) %>%
  filter(!is.na(MAP_mm), !is.na(MAT_C))

# Shade the arid region (MAP < 20*(MAT+10)) across the observed MAT range
mat_seq    <- seq(min(climate.data$MAT_C, na.rm = TRUE) - 2,
                  max(climate.data$MAT_C, na.rm = TRUE) + 2, length.out = 200)
arid_bound <- data.frame(MAT_C = mat_seq, MAP_bound = 15 * (mat_seq + 10))

plot.climate <- ggplot(climate.data, aes(x = MAT_C, y = MAP_mm)) +
  # Arid shading below the AI = 20 boundary
  geom_ribbon(data  = arid_bound,
              aes(x = MAT_C, ymin = 0, ymax = pmax(0, MAP_bound)),
              inherit.aes = FALSE,
              fill = "#c8a96e", alpha = 0.18) +
  # AI = 20 boundary line
  geom_line(data  = arid_bound,
            aes(x = MAT_C, y = MAP_bound),
            inherit.aes = FALSE,
            linetype = "dashed", color = "white", linewidth = 0.7) +
  annotate("text", x = max(mat_seq) - 1, y = 15 * (max(mat_seq) + 10) + 80,
           label = "AI = 15 (arid boundary)", hjust = 1,
           size = 3.2, color = "white") +
  # All sites: filled circle colored by EcoType
  geom_point(aes(color = EcoType), shape = 16, size = 3, alpha = 0.85) +
  # Arid sites (AI < 20): white outline ring drawn on top to highlight them
  geom_point(data = filter(climate.data, Aridity == "Arid"),
             aes(color = EcoType), shape = 21, size = 4,
             stroke = 1.4, fill = NA, color = "white") +
  scale_color_discrete_sequential(palette = "Hawaii", name = "") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "black", colour = "black"),
    plot.background  = element_rect(fill = "black", colour = "black"),
    panel.grid.major = element_line(colour = "grey30"),
    panel.grid.minor = element_line(colour = "grey10"),
    axis.text        = element_text(colour = "white"),
    axis.title       = element_text(colour = "white"),
    legend.background = element_rect(fill = "black"),
    legend.text       = element_text(colour = "white"),
    legend.title      = element_text(colour = "white"),
    text             = element_text(size = 16)
  ) +
  labs(
    x = "Mean Annual Temperature (°C)",
    y = "Mean Annual Precipitation (mm)"
  )

# ── Assemble 3-panel figure ───────────────────────────────────────────────────
final.plot <- ggarrange(
  map,
  ggarrange(plot.tower.counts, plot.climate,
            ncol = 2, common.legend = FALSE,
            labels = c("B", "C"), font.label = list(color = "white", size = 14)),
  nrow          = 2,
  common.legend = TRUE,
  legend        = "bottom"
)

ggsave("FIGURES/Map_plot.png", plot = final.plot, width = 12, height = 10, units = "in")
