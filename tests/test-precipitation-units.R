# Check the actual rolling-expression used in the workflow.
expressions <- parse("workflows/06_NEON_ERA5Gapfill.R")
find_roll <- function(x) {
  if (!is.call(x)) return(NULL)
  if (identical(x[[1]],quote(slider::slide_dbl)) &&
      grepl("ERA5_precip_mm",paste(deparse(x),collapse=" "))) return(x)
  for (child in as.list(x)[-1]) { if (missing(child)) next; result <- find_roll(child); if (!is.null(result)) return(result) }
  NULL
}
roll <- NULL
for (x in expressions) { roll <- find_roll(x); if (!is.null(roll)) break }
stopifnot(!is.null(roll))
ERA5_precip_mm <- rep(1,336)
x <- eval(roll)
stopifnot(tail(x,1)==168, x[48]==24, x[1]==0.5)
ERA5_precip_mm <- rep(0,336)
stopifnot(all(eval(roll)==0))
ERA5_precip_mm <- c(NA_real_, rep(1,335))
stopifnot(tail(eval(roll),1)==167.5)
message("Precipitation accumulation checks passed.")
