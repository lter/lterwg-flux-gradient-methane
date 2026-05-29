# Download ERA5-Land monthly gridded climate inputs for monthly CH4 upscaling.
#
# Output directory:
#   /Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly/DATA/era5_land_monthly
#
# Required CDS setup:
#   1. Create/login to a Copernicus CDS account.
#   2. Accept the ERA5-Land monthly means licence in the CDS web portal.
#   3. Install the Python package `cdsapi`.
#   4. Put credentials in ~/.cdsapirc or set CDSAPI_URL and CDSAPI_KEY.

library(jsonlite)

spatial_dir <- Sys.getenv(
  "MONTHLY_UPSCALING_DIR",
  unset = "/Volumes/MaloneLab/Research/FluxGradient/METHANE/Upscaling_Monthly"
)

era5_dir <- file.path(spatial_dir, "DATA", "era5_land_monthly")
dir.create(era5_dir, recursive = TRUE, showWarnings = FALSE)

years <- 2000:2025
months <- sprintf("%02d", 1:12)

cdsapirc <- path.expand("~/.cdsapirc")
cds_config_present <- file.exists(cdsapirc) && file.info(cdsapirc)$size > 0
cds_env_present <- nzchar(Sys.getenv("CDSAPI_URL")) && nzchar(Sys.getenv("CDSAPI_KEY"))

if (!cds_config_present && !cds_env_present) {
  stop(
    "No CDS credentials found. ~/.cdsapirc exists but is empty, or CDSAPI_URL/CDSAPI_KEY are not set. ",
    "Add CDS credentials and accept the ERA5-Land monthly means licence, then rerun this script."
  )
}

python_bin <- Sys.getenv("CDSAPI_PYTHON", unset = "")
if (!nzchar(python_bin)) {
  candidate_python <- "/private/tmp/era5_cdsapi_venv/bin/python"
  python_bin <- if (file.exists(candidate_python)) candidate_python else "python3"
}

has_cdsapi <- system2(
  python_bin,
  c("-c", shQuote("import cdsapi")),
  stdout = TRUE,
  stderr = TRUE
)

if (!identical(attr(has_cdsapi, "status"), NULL)) {
  stop(
    "Python package `cdsapi` is not installed. Install it with `python3 -m pip install --user cdsapi`, ",
    "then rerun this script."
  )
}

request_script <- file.path(era5_dir, "retrieve_era5_land_monthly.py")

python_code <- paste(
  "import cdsapi",
  "from pathlib import Path",
  "",
  sprintf("out_dir = Path(%s)", toJSON(era5_dir, auto_unbox = TRUE)),
  sprintf("years = %s", toJSON(as.character(years), auto_unbox = TRUE)),
  sprintf("months = %s", toJSON(months, auto_unbox = TRUE)),
  "",
  "client = cdsapi.Client()",
  "",
  "for year in years:",
  "    target = out_dir / f'era5_land_monthly_{year}.nc'",
  "    if target.exists() and target.stat().st_size > 0:",
  "        print(f'Skipping existing {target}')",
  "        continue",
  "    request = {",
  "        'product_type': 'monthly_averaged_reanalysis',",
  "        'variable': [",
  "            '2m_temperature',",
  "            'volumetric_soil_water_layer_1',",
  "            'volumetric_soil_water_layer_2',",
  "            'total_precipitation',",
  "        ],",
  "        'year': year,",
  "        'month': months,",
  "        'time': ['00:00'],",
  "        'data_format': 'netcdf',",
  "        'download_format': 'unarchived',",
  "        'area': [90, -180, -90, 180],",
  "        'grid': [0.5, 0.5],",
  "    }",
  "    print(f'Downloading {target}')",
  "    client.retrieve('reanalysis-era5-land-monthly-means', request, str(target))",
  sep = "\n"
)

writeLines(python_code, request_script)

message("Wrote CDS request script to ", request_script)
message("Downloading ERA5-Land monthly files to ", era5_dir)

status <- system2(python_bin, request_script)
if (!identical(status, 0L)) {
  stop("ERA5-Land download failed with status ", status)
}

downloaded <- list.files(era5_dir, pattern = "^era5_land_monthly_[0-9]{4}\\.nc$", full.names = TRUE)
message("Downloaded/found ", length(downloaded), " annual ERA5-Land NetCDF files.")
