# see how well the RSHP scheme form eval works here:
library(tidyverse)
library(randomForest)

load(file= paste(localdir, 'SITE_RSHP_MODEL.Rdata', sep="") ) # The model is in here:

load(fs::path(localdir,paste0("CCC_CH4_RSHP.Rdata")))# The data on the RSHP for CH4 is in there
VAL_RSHP$RSHP.rf <- predict(rf.good.ccc.ec, VAL_RSHP)

VAL_RSHP$RSHP.rf %>% summary
VAL_RSHP$Good.CCC.GF %>% summary

library(caret)

confusionMatrix(VAL_RSHP$RSHP.rf, VAL_RSHP$Good.CCC.GF)
