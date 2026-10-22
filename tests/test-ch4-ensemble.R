# Run from repository root: Rscript tests/test-ch4-ensemble.R
suppressPackageStartupMessages(library(tidyverse))
script <- parse("workflows/03_NEON_EnsembleGapfill.R")
# Run the actual calculation, with temporary inputs; skip only the Drive upload.
script <- script[!vapply(script, function(x) grepl("googledrive::drive_upload", paste(deparse(x), collapse=" ")), logical(1))]
run <- function(gases, rates) {
  e <- new.env(parent=globalenv())
  e$localdir.ch4 <- tempfile(); dir.create(e$localdir.ch4)
  on.exit(unlink(e$localdir.ch4, recursive=TRUE))
  dir.create(file.path(e$localdir.ch4, "NEON_GradientFlux_Data_Filter"))
  e$site.list <- "TEST"
  NEON_RSHP_EVAL <- data.frame(Site="TEST", Approach="AE", dLevelsAminusB="2_1", RSHP.CO2="1", RSHP.H2O="0")
  SITE_DATA_FILTERED <- list(TEST=data.frame(gas=gases, FG_mean=rates,
    Approach="AE", dLevelsAminusB="2_1", timeEndA.local=as.POSIXct("2024-01-01", tz="UTC"), PAR=1, Tair_C=10))
  save(NEON_RSHP_EVAL, file=file.path(e$localdir.ch4,"NEON_RSHP_CCC.Rdata"))
  save(SITE_DATA_FILTERED, file=file.path(e$localdir.ch4,"NEON_GradientFlux_Data_Filter/SITE_DATA_FILTERED_CH4.Rdata"))
  eval(script, envir=e)
  e$SITE_DATA_FILTERED_Final_RSHP_ENSEMBLE$TEST$FG_ENSEMBLE_RSHP
}
stopifnot(identical(run(c("CH4","CO2","H2O"), c(-2,100,10000)), -2),
          identical(run("CH4", -2), -2),
          inherits(try(run("CO2", 100), silent=TRUE), "try-error"),
          inherits(try(run(NA_character_, 1), silent=TRUE), "try-error"))
message("CH4 ensemble gas-selection checks passed.")
