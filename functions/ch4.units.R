# CH4 molar rates -> mg carbon per interval. Choose the input prefix explicitly.
ch4_rate_to_mgC <- function(x, input_unit, seconds = 1800) {
  if (length(input_unit) != 1L || is.na(input_unit) ||
      !input_unit %in% c("umol", "nmol")) stop("input_unit must be umol or nmol")
  if (length(seconds) != 1L || !is.finite(seconds) || seconds <= 0) {
    stop("seconds must be a positive finite interval")
  }
  mgC_per_unit <- if (input_unit == "umol") 12.011 / 1000 else 12.011 / 1e6
  x * seconds * mgC_per_unit
}
