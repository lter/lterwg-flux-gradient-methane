# Diel Analysis Functions for CH4:

DIEL.season.ch4 <- function( dataframe, flux, Gas){
  
  dataframe <- dataframe %>% as.data.frame %>% rename( timeEndA.local = time.rounded)
  
  if(Gas == "CH4"){
    dataframe$flux <- ((dataframe[, flux]* 16.042)/1000000)*1800 
  }
  
  
  if(Gas == "CO2"){
    dataframe$flux <- ((dataframe[, flux]* 44.01)/1000000)*1800 
  }
  
  if(Gas == "H2O"){
    dataframe$flux <- ((dataframe[, flux]* 16)/1000000)*1800 
  }
  
  # Make sure season is in the dataframe:
  dataframe.gs <- dataframe %>% 
    mutate(  Month = timeEndA.local %>% format("%m") %>% as.numeric,
             YearMon = timeEndA.local %>% format("%Y-%m"),
             Year = timeEndA.local %>% format("%Y"),
             Hour = timeEndA.local %>% format("%H") %>% as.numeric) %>% 
    filter(!is.na( flux) == TRUE)
  
  season <- unique(dataframe.gs$season)
  
  message(" Ready to summarize by season")
  Final.data.all <- data.frame()
  
  for( i in season){
    print(i)
    try({
      
      # Remove outliers:
      subset <- dataframe.gs %>% filter(season == i, 
                                        flux > -50, flux < 50)
      
      count <- subset$flux %>% na.omit %>% length
      
      if(count > 48){
        model <- loess( flux ~ Hour , data = subset )
        model %>% plot
        Diel.df <- data.frame(Hour = seq(0, 23), season = i)
        pred <- predict(model, newdata = Diel.df, se=TRUE)
        
        new.data  <- Diel.df %>% mutate(DIEL = pred$fit, 
                                        DIEL.SE =pred$fit* qt(0.95 / 2 + 0.5, 
                                                              pred$df)) %>% mutate(data="all")
        
        
        Peak <- new.data$DIEL %>% max(na.rm=T)
        MIN <- new.data$DIEL %>% min(na.rm=T)
        
        new.data$Peak.Hour <-  new.data$Hour[new.data$DIEL == Peak] %>% as.numeric %>% mean(na.rm=T)
        new.data$Min.Hour <-  new.data$Hour[new.data$DIEL == MIN] %>% as.numeric %>% mean(na.rm=T)
        new.data$count <- count
        Final.data.all <- rbind(Final.data.all, new.data) 
        rm(new.data)
        } },silent=T)}
  
 
  message(" Done fitting loess for all data")
  
  Final.data <- Final.data.all
  
  if( exists('Final.data' )) { 
    return( Final.data)
  }
  
}
