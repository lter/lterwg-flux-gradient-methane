# Get soil moisture data for NEON tower sites:

# Install the packages if not already installed
# install.packages("neonUtilities")
# install.packages("tidyverse")

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
start_date <- "2020-01"
end_date <- "2024-12"
data_product_id <- "DP1.00094.001" 

# Download Soil moisture data for all sites

for ( site in site.list){
  # Download the data
  # Setting check.size = FALSE prevents a user prompt for large downloads
  print(paste('Downloading Data for:', site))
  soil.data <- loadByProduct(
    dpID = data_product_id,
    site = site,
    startdate = start_date,
    enddate = end_date,
    package = "basic", # You can choose "basic" or "expanded"
    timeIndex = "30", # 30-minute averaging interval
    check.size = FALSE,
    nCores=7)
  
  localdir.site <- paste(localdir,"/", site, sep = "")
  
  files <- paste(site, "_soildata.Rdata", sep = "")
  
  save( soil.data, file=paste(localdir.site, "/", files, sep=""))
  
  
}

# EOF