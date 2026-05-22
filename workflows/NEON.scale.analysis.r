# Compatibility wrapper.
# Native-scale flux products now live in flow.30min.analysis.R so that mean
# 30-minute flux, lookup-filled daily flux, and ERA5 annual budgets are created
# in one methods-consistent place.

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  "workflows"
}

source(file.path(script_dir, "flow.30min.analysis.R"))
