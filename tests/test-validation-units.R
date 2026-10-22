# Run from repository root: Rscript tests/test-validation-units.R
suppressPackageStartupMessages(library(tidyverse))
source("functions/validation.units.R")
source("functions/ch4.units.R")
eq <- function(x,y) stopifnot(isTRUE(all.equal(x,y,tolerance=1e-12)))
DirRepo.ch4 <- getwd()
localdir.ch4 <- tempfile(); dir.create(localdir.ch4)
dir.create(file.path(localdir.ch4,"Storage_Flux/SinglePoint"),recursive=TRUE)
t <- as.POSIXct("2024-01-01",tz="UTC") + 0:1*1800
SITEval_DATA_FILTERED_RSHPc_H <- list()
Validation_storage_singlepoint <- list()
for (site in c("SE-Sto","SE-Svb","US-Uaf")) {
  k <- if (site=="SE-Sto") 1 else 1000
  d <- data.frame(time.rounded=t,gas="CH4",FG_mean=c(1,-1)*k,EC_mean=c(2,-2)*k)
  SITEval_DATA_FILTERED_RSHPc_H[[site]] <- d
  Validation_storage_singlepoint[[site]] <- data.frame(time.rounded=t,storage_flux_1pt=c(500,NA))
  normalized <- validation_flux_to_umol(d,site)
  eq(normalized$FG_mean,c(1,-1))
  stopifnot(identical(validation_flux_to_umol(normalized,site),normalized))
  stopifnot(inherits(try(validation_flux_to_umol(transform(d,FG_total=FG_mean),site),silent=TRUE),"try-error"))
}
save(SITEval_DATA_FILTERED_RSHPc_H,file=file.path(localdir.ch4,"SITEval_DATA_FILTERED_RSHP_EnSEMBLE.Rdata"))
save(Validation_storage_singlepoint,file=file.path(localdir.ch4,"Storage_Flux/SinglePoint/Validation_storage_singlepoint.RData"))
source("workflows/Validation/02b_VAL_TotalFlux.R")
for (d in SITEval_DATA_FILTERED_RSHPc_H_total) {
  eq(d$FG_total,c(1.5,-1));eq(d$EC_mean,c(2,-2))
  stopifnot(identical(d$storage_added,c(TRUE,FALSE)),all(d$CH4_flux_unit=="umol m-2 s-1"))
  eq(ch4_rate_to_mgC(d$FG_total,"umol"),c(32.4297,-21.6198))
}
unlink(localdir.ch4,recursive=TRUE)
message("Validation unit normalization and storage-addition checks passed.")
