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

# Download data: ####

start_date <- "2020-01"
end_date <- "2024-12"
data_product_id <- "DP1.00094.001" 

# Download Soil moisture data for all sites
# there is no data for PUUM (31)
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

# Build soilmoisture Dataframe:(YearMon and Season) ####

Site_SoilData <- data.frame()
Site_SoilData_Season <- data.frame()

for ( site in site.list[c(1:30, 32:47)]){
  print(site)
  localdir.site <- paste(localdir,"/", site, sep = "")
  
  files <- paste(site, "_soildata.Rdata", sep = "")
  
  load(file=paste(localdir.site, "/", files, sep=""))
  
  soil.data.filter <- soil.data$SWS_30_minute %>% filter(verticalPosition == '501') %>% mutate(
    TIMESTAMP = as.POSIXct( startDateTime),
    YearMon = format(TIMESTAMP,'%Y-%m')) %>% reframe(.by= c(YearMon, siteID, domainID),
                                                     VSWCMean = mean(VSWCMean, na.rm=T),
                                                     VSWCMin = min(VSWCMinimum, na.rm=T),
                                                     VSWCMax = max(VSWCMaximum, na.rm=T),
                                                     VSWCVar = mean(VSWCVariance, na.rm=T) ) %>% 
    rename( Site = siteID,
            domain = domainID)
                                                                                                         
  Site_SoilData <- rbind( Site_SoilData, soil.data.filter )
  
  soil.data.filter.seaon <- soil.data$SWS_30_minute %>% filter(verticalPosition == '501') %>% mutate(
    TIMESTAMP = as.POSIXct( startDateTime),
    Month = format(TIMESTAMP,'%m')%>% as.numeric,
    Season = case_when(
      Month %in% c(12, 1, 2) ~ "Winter",
      Month %in% 3:5 ~ "Spring",
      Month %in% 6:8 ~ "Summer",
      Month %in% 9:11 ~ "Autumn")) %>% 
    reframe(.by= c(Season, siteID, domainID),
            VSWCMean = mean(VSWCMean, na.rm=T),
            VSWCMin = min(VSWCMinimum, na.rm=T),
            VSWCMax = max(VSWCMaximum, na.rm=T),
            VSWCVar = mean(VSWCVariance, na.rm=T) ) %>%                                          
    rename( Site = siteID,
            domain = domainID)
  Site_SoilData_Season <- rbind( Site_SoilData_Season, soil.data.filter.seaon )
  
}



save( Site_SoilData,
      Site_SoilData_Season,
      file=paste(localdir.ch4, "/", 'Soildata_YearMon.Rdata', sep=""))
