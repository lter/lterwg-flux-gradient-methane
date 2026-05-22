# Compatibility wrapper for the earlier filename used in notes.
# The canonical source/sink and driver workflow is NEON.DriverScale.Analysis.R.

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  "workflows"
}

source(file.path(script_dir, "NEON.DriverScale.Analysis.R"))
