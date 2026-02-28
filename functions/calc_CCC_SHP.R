source(fs::path(DirRepo.eval,"./functions/calc.linear.terms.R"))

CCC_SamplingHeightPairs <- function( DATA){
  
  site.list <- DATA %>% names()
  RSHP <- data.frame()
  
  for( site in site.list) {
    
    print(site)
    
    df <-DATA[[site]]  %>% 
      mutate(month = format(timeEndA.local,'%m') %>% as.numeric,
             season = case_when(
               month %in% c(12, 1, 2) ~ "Winter",
               month %in% c(3, 4, 5) ~ "Spring",
               month %in% c(6, 7, 8) ~ "Summer",
               TRUE ~ "Autumn"),
             hour = format(timeEndA.local,'%H'),
             count= case_when( is.na(FG_mean) == FALSE ~ 1,
                               TRUE ~ 0)) %>% distinct %>% 
      mutate( time.rounded = timeEndA.local %>% round_date( unit = "30 minutes") )  %>%
      reframe( .by= c( dLevelsAminusB, gas, Approach, time.rounded), 
               EC_mean = mean(EC_mean, na.rm=T), 
               FG_mean = mean(FG_mean, na.rm=T)) %>% mutate( levels = paste(Approach, dLevelsAminusB, sep="-"))
    
    
    # CO2 : ####
    df_CO2 <- df %>% filter( gas == "CO2") 
    
    if(df_CO2$levels %>% unique %>% length > 1){
      
      df_level_wide <- df_CO2 %>% select(levels, time.rounded, FG_mean ) %>% 
        pivot_wider(names_from = levels, values_from = c(FG_mean)) %>% as.data.frame() 
      
      canopy_level.list.co2 <- df_CO2$levels %>% unique
      
      combinations.canopy_level.list.co2 <- combn(canopy_level.list.co2,2)
      
      for(i in 1:ncol(combinations.canopy_level.list.co2)){
        print( paste( "Working on combination", combinations.canopy_level.list.co2[1, i], combinations.canopy_level.list.co2[2, i], sep = " "))
        
        new.data <- linear_terms(  df = df_level_wide,
                                   x_col=combinations.canopy_level.list.co2[1, i],
                                   y_col= combinations.canopy_level.list.co2[2, i],
                                   site= site,
                                   var1= combinations.canopy_level.list.co2[1, i],
                                   var2=combinations.canopy_level.list.co2[2, i])%>% 
          mutate( gas = "CO2")
        
        RSHP <- rbind( RSHP, new.data)
      }
      
      
    }
    
    # H2O : ####
    df_H2O <- df %>% filter( gas == "H2O")
    
    if(df_H2O$levels %>% unique %>% length > 1){
      
      df_level_wide <- df_H2O %>% select(levels, time.rounded, FG_mean ) %>% 
        pivot_wider(names_from = levels, values_from = c(FG_mean)) %>% as.data.frame() 
      
      canopy_level.list.h2o <- df_H2O$levels %>% unique
      
      combinations.canopy_level.list.h2o <- combn(canopy_level.list.h2o,2)
      
      for(i in 1:ncol(combinations.canopy_level.list.h2o)){
        print( paste( "Working on combination", combinations.canopy_level.list.h2o[1, i], combinations.canopy_level.list.h2o[2, i], sep = " "))
        
        new.data <- linear_terms(  df = df_level_wide,
                                   x_col= combinations.canopy_level.list.h2o[1, i],
                                   y_col= combinations.canopy_level.list.h2o[2, i],
                                   site = site,
                                   var1= combinations.canopy_level.list.h2o[1, i],
                                   var2=combinations.canopy_level.list.h2o[2, i])%>% 
          mutate( gas = "H2O")
        
        RSHP <- rbind( RSHP, new.data)
      }
      
    }
    
    # CH4 : ####
    df_CH4 <- df %>% filter( gas == "CH4")
    
    if(df_CH4$levels %>% unique %>% length > 1){
      
      df_level_wide <- df_CH4 %>% select(levels, time.rounded, FG_mean ) %>% 
        pivot_wider(names_from = levels, values_from = c(FG_mean)) %>% as.data.frame() 
      
      canopy_level.list.CH4 <- df_CH4$levels %>% unique
      
      combinations.canopy_level.list.CH4 <- combn(canopy_level.list.CH4,2)
      
      for(i in 1:ncol(combinations.canopy_level.list.CH4)){
        print( paste( "Working on combination", combinations.canopy_level.list.CH4[1, i], combinations.canopy_level.list.CH4[2, i], sep = " "))
        
        new.data <- linear_terms(  df = df_level_wide,
                                   x_col= combinations.canopy_level.list.CH4[1, i],
                                   y_col= combinations.canopy_level.list.CH4[2, i],
                                   site = site,
                                   var1= combinations.canopy_level.list.CH4[1, i],
                                   var2=combinations.canopy_level.list.CH4[2, i])%>% 
          mutate( gas = "CH4")
        
        RSHP <- rbind( RSHP, new.data)
      }
      
    }
    
  }
  
  return( RSHP )
}


