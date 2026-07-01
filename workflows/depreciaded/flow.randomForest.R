# Look at differences in Q10 and Rbase for sites that are sinks and sources:

Site_SoilData %>% mutate(Date = as.Date(paste(YearMon, "-01", sep="")),
                         SITE_ID = Site)


library(randomForest)

rf.model <- randomForest( Source ~  VSWCMean + 
                            Q10.mean + Rref.mean + VSWCVar + 
                            sulfurTot + dryMass + acidity + 
                            FRC.kg.C.m.2 + SOC.kg.C.m.2 +
                            beta_roots + beta_soc, data = parms.site.season %>% na.omit)

varImpPlot(rf.model)


# Sensitivity analysis 
continious.variables <- c( 'VSWCMean' ,  'Q10.mean' , 'Rref.mean' , 
                           'VSWCVar' , 'sulfurTot' , 'dryMass' , 
                           'acidity','FRC.kg.C.m.2', 'SOC.kg.C.m.2', 
                           'beta_roots' ,'beta_soc')

min.conditions <- parms %>% dplyr::select(continious.variables ) %>% na.omit %>%  
  summarise(across(where(is.numeric), ~ min(.x, na.rm = TRUE)))

max.conditions <- parms %>% dplyr::select(continious.variables )  %>% na.omit  %>%  
  summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))

median.conditions <- parms %>% dplyr::select(continious.variables ) %>% na.omit %>%  
  summarise(across(where(is.numeric), ~ (min(.x, na.rm = TRUE) - max(.x, na.rm = TRUE))/2))

conditions <- rbind( min.conditions %>% mutate(level = "min"), max.conditions%>% mutate(level = "max"),
                     median.conditions%>% mutate(level = "median"))

condition.ranges <- data.frame()

for( i in continious.variables){
  print(i)
  
  min <- min.conditions[i] %>% as.numeric
  max <- max.conditions[i]%>% as.numeric
  range <- (max - min)/20
  seq.cond <- seq(min , max, range)
  target.df <- data.frame()
  target.df <- data.frame( target = as.numeric(seq.cond))
  target.df[, i] <- target.df$target
  
  
  conditions.sub <- conditions %>% dplyr::select(-c(i))
  
  conditions.sub.join <- cross_join( conditions.sub, target.df %>% dplyr::select(i)) %>% mutate( target = i)
  condition.ranges <- rbind(  conditions.sub.join ,  condition.ranges )
  
  
}


condition.ranges$predictions <- predict(rf.model, condition.ranges, type="prob" )[,2]

condition.ranges <- condition.ranges %>% distinct()

condition.ranges %>% filter (target == "VSWCMean") %>% ggplot( aes( x= VSWCMean, y = predictions, col=level)) + geom_point() + geom_smooth()

condition.ranges %>% filter (target == "dryMass") %>% ggplot( aes( x= dryMass, y = predictions, col=level)) + geom_point() + geom_smooth()

condition.ranges %>% filter (target == "sulfurTot") %>% ggplot( aes( x= sulfurTot, y = predictions, col=level)) + geom_point() + geom_smooth()

condition.ranges %>% filter (target == "Q10.mean") %>% ggplot( aes( x= Q10.mean, y = predictions, col=level)) + geom_point() + geom_smooth()

condition.ranges %>% filter (target == "acidity") %>% ggplot( aes( x= acidity, y = predictions, col=level)) + geom_point() + geom_smooth()

condition.ranges %>% filter (target == "Rref.mean") %>% ggplot( aes( x= Rref.mean, y = predictions, col=level)) + geom_point() + geom_smooth()

condition.ranges %>% filter (target == "FRC.kg.C.m.2") %>% ggplot( aes( x= FRC.kg.C.m.2, y = predictions, col=level)) + geom_point() + geom_smooth()

condition.ranges %>% filter (target == "SOC.kg.C.m.2") %>% ggplot( aes( x= SOC.kg.C.m.2, y = predictions, col=level)) + geom_point() + geom_smooth()

condition.ranges %>% filter (target == "beta_roots") %>% ggplot( aes( x= beta_roots, y = predictions, col=level)) + geom_point() + geom_smooth()

condition.ranges %>% filter (target == "beta_soc") %>% ggplot( aes( x= beta_soc, y = predictions, col=level)) + geom_point() + geom_smooth()

# Multidemensional Scaling 