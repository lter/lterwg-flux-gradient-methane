
# This flow extract NEON soil chemical and phisical properties for all towersites.
# Load the packages
library(neonUtilities)
library(tidyverse)

# Import the metdata file for tower sites to make a sitelist:
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


site.list <- metadata$Site %>% unique

# Soil data is currently summarized across horizons. 

for ( site in site.list){
  # Download the data
  # Setting check.size = FALSE prevents a user prompt for large downloads
  print(paste('Downloading Data for:', site))
 
  soil.chem <- loadByProduct(
    dpID = 'DP1.10047.001', # Megapit information,
    site = site,
    #startdate = start_date,
    #enddate = end_date,
    package = "basic", # You can choose "basic" or "expanded"
    timeIndex = "30", # 30-minute averaging interval
    check.size = FALSE,
    nCores=7)
  
  root.data <- loadByProduct(
    dpID = 'DP1.10067.001',# Root biomass and chemistry peroiodic
    site = site,
    #startdate = start_date,
    #enddate = end_date,
    package = "basic", # You can choose "basic" or "expanded"
    timeIndex = "30", # 30-minute averaging interval
    check.size = FALSE,
    nCores=7)
  
  
  
  localdir.site <- paste(localdir,"/", site, sep = "")
  
  files <- paste(site, "_soilChem_Root.Rdata", sep = "")
  
  save( soil.chem,
        root.data, 
        file=paste(localdir.site, "/", files, sep=""))
  
  
}

soil.data <- data.frame()
for ( site in site.list){
  print(paste('Building soil chemical and rootbiomass file for:', site))
# Build the site file that contains all the data needed:
localdir.site <- paste(localdir,"/", site, sep = "")
files <- paste(site, "_soilChem_Root.Rdata", sep = "")
load(file=paste(localdir.site, "/", files, sep=""))
  

bd.vars <- c( 'collectDate', 'horizonName', 'siteID', 
              'bulkDensCenterDepth', 'bulkDensTopDepth',
              'bulkDensBottomDepth' ,'bulkDensThirdBar', 'bulkDensOvenDry')
  

chem.vars <- c( 'siteID', 'collectDate', 'horizonName', 'sulfurTot',
          'biogeoTopDepth','biogeoBottomDepth', 'carbonTot','nitrogenTot',
          'ctonRatio', 'estimatedOC', 'acidity')

particle.vars <- c('siteID', 'collectDate', 'horizonName', 
                   'sandTotal', 'siltTotal', 'clayTotal')

root.biomass.vars <- c( 'siteID', 'collectDate', "sizeCategory",
                        "rootStatus","dryMass","mycorrhizaeVisible" )  

soil.chem.biogeochem <- soil.chem$spc_biogeochem[chem.vars] %>% 
  reframe( .by = c('siteID'), 
           sulfurTot = mean(sulfurTot, na.rm=T ),
           biogeoTopDepth = mean(biogeoTopDepth, na.rm=T ),
           biogeoBottomDepth = mean( biogeoBottomDepth, na.rm=T ),
           carbonTot = mean(sulfurTot, na.rm=T ),
           nitrogenTot = mean(nitrogenTot, na.rm=T ),
           ctonRatio = mean( ctonRatio, na.rm=T ),
           estimatedOC = mean( estimatedOC, na.rm=T ),
           acidity= mean(acidity, na.rm=T ))

soil.chem.particlesize <- soil.chem$spc_particlesize [particle.vars] %>% 
  reframe( .by = c('siteID'), 
           sandTotal= mean(sandTotal, na.rm=T ),
           siltTotal= mean(siltTotal, na.rm=T ),
           clayTotal= mean(clayTotal, na.rm=T )) 

root.biomass <- root.data$bbc_rootmass[root.biomass.vars  ] %>% 
  reframe( .by = c('siteID','sizeCategory'), 
           dryMass= mean(dryMass, na.rm=T )) %>%  
  reframe( .by = c('siteID'), 
           dryMass= sum(dryMass, na.rm=T ))

if( 'spc_bulkdensity' %in% names(soil.chem)){
  
  soil.chem.bulkdensity <- soil.chem$spc_bulkdensity[bd.vars] %>% 
    reframe( .by = c('siteID'),
             bulkDensOvenDry = mean(bulkDensOvenDry, na.rm=T ))
  
  soil.root.total <- soil.chem.biogeochem %>% 
    full_join(soil.chem.bulkdensity, by = "siteID" ) %>%
    full_join( soil.chem.particlesize , by = "siteID"  ) %>% 
    full_join(root.biomass , by = "siteID" )} else{
  
      soil.root.total <- soil.chem.biogeochem %>% mutate(bulkDensOvenDry = NA) %>% 
    full_join( soil.chem.particlesize , by = "siteID"  ) %>% 
    full_join(root.biomass , by = "siteID" ) 
}

 soil.data <- gtools::smartbind( soil.data , soil.root.total )

}

# Save the file created:
write.csv( soil.data,
      file=paste(localdir, "/", 'Soil_Biogeochem_RootBiomass.csv', sep=""))
  
