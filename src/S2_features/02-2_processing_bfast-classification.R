#!/usr/bin/env Rscript

############################################################
# Script:   02-2_processing_bfast-classification.R
# Author:   [Your Name]
# Project:  [Your Project Name]
#
# Purpose
# -------
# Detect per-pixel structural breaks in vegetation dynamics using
# BFAST on a **kNDVI time series** derived from a Sentinel-2
# gdalcubes NetCDF time cube.
#
# Conceptual workflow
# -------------------
# 1) Open an existing Sentinel-2 monthly NetCDF cube with gdalcubes.
# 2) Use `reduce_time()` with a **custom R reducer**:
#      - for each pixel, build a kNDVI time series from bands B08/B04,
#      - run `bfastmonitor()` on that series,
#      - return:
#          * `change_date`     – breakpoint time (if any),
#          * `change_magnitude` – magnitude of change at breakpoint.
# 3) Write the resulting per-pixel change information to a new
#    NetCDF file for further analysis and mapping.
#
# INPUT
# -----
#   - Monthly Sentinel-2 cube:
#       file.path(root_folder, "data", "burgwald_2018_2022_all.nc")
#
# OUTPUT
# ------
#   - BFAST results (per-pixel):
#       file.path(root_folder, "data", "burgwald_bfast_results.nc")
#
# Assumptions / requirements
# --------------------------
#   - `01_setup-burgwald.R` has been sourced so that:
#       * `root_folder <- here::here()` is defined.
#   - The NetCDF cube contains at least bands:
#       * "B08" (NIR)
#       * "B04" (RED)
#   - The cube is **monthly** from 2018–2022:
#       * time index used as `ts(..., start = c(2018, 1), frequency = 12)`
#   - Reflectances are scaled to 0–10000 (or 0–1; the /10000 is harmless).
#
# Technical notes
# ---------------
#   - `gdalcubes::reduce_time()` executes the custom reducer on the
#     *server side* / cube engine; only the resulting summary bands
#     are written out.
#   - The `FUN` passed to `reduce_time()` must be **self-contained**:
#       * all packages needed inside must be loaded within `FUN`
#         (here: `library(bfast)`).
#   - Heavy work (kNDVI computation + BFAST) happens inside the
#     gdalcubes processing chain; R only orchestrates the call.
############################################################

source(here::here("_core", "00-setup-burgwald.R"))

# Set parallelisation for gdalcubes (adjust cores to your machine)
gdalcubes_options(parallel = 12)

# Path to your existing Sentinel-2 monthly cube
cube_file  <- file.path(root_folder, "data", "burgwald_2018_2022_all.nc")

# Output NetCDF with BFAST results
bfast_file <- file.path(root_folder, "data", "burgwald_bfast_results.nc")

system.time({
  ncdf_cube(cube_file) |>
    # reduce_time with a custom R reducer:
    # For each pixel, we:
    #   1) compute a kNDVI time series from bands B08 and B04
    #   2) run bfastmonitor() on that time series
    #   3) return change_date and change_magnitude
    reduce_time(
      names = c("change_date", "change_magnitude"),
      FUN = function(x) {
        # x is a matrix (bands x time) for ONE pixel.
        # We must work INSIDE this function only.
        
        # ---- 1) Compute kNDVI time series --------------------
        # Here we use a kNDVI variant as in your MOF example:
        #   kNDVI = tanh(((NIR - RED) / (NIR + RED))^2)
        # Assumes Sentinel-2 reflectances are scaled to 0–10000
        # – if already 0–1, the /10000 step is harmless but unnecessary.
        nir <- x["B08", ] / 10000
        red <- x["B04", ] / 10000
        
        kndvi <- tanh(((nir - red) / (nir + red))^2)
        
        # If everything is NA → no information, return NA
        if (all(is.na(kndvi))) {
          return(c(NA, NA))
        }
        
        # ---- 2) Build time series object for bfastmonitor ----
        # Time dimension:
        #   - monthly cube 2018–2022 → start c(2018, 1), frequency = 12
        kndvi_ts <- ts(kndvi, start = c(2018, 1), frequency = 12)
        
        # ---- 3) Run bfastmonitor() and handle errors ---------
        # IMPORTANT:
        # - Packages must be loaded INSIDE the reducer,
        #   because gdalcubes runs this in a fresh R process.
        library(bfast)
        
        out <- tryCatch({
          # start = when monitoring begins (tail period)
          # history = "all" = use full pre-monitoring history
          # level = 0.01 = significance threshold
          res <- bfastmonitor(
            kndvi_ts,
            start   = c(2020, 1),   # monitoring phase start (adapt if needed)
            history = "all",
            level   = 0.01
          )
          
          # bfastmonitor() returns:
          #   res$breakpoint  - time of first structural break
          #   res$magnitude   - magnitude of change at breakpoint
          c(res$breakpoint, res$magnitude)
        },
        error = function(e) {
          # On any problem: return NA values
          c(NA, NA)
        })
        
        out
      }
    ) |>
    # Write the per-pixel change results to a NetCDF file
    write_ncdf(bfast_file, overwrite = TRUE)
})

cat("BFAST results written to:\n  ", bfast_file, "\n")


bfast_star <- gdalcubes::ncdf_cube(bfast_file) |>
  stars::st_as_stars()

tmap_mode("view")
tm_shape(bfast_star["change_date"]) +
  tm_raster() +
  tm_layout(title = "BFAST change date (kNDVI)")

tmap_mode("plot")
tm_shape(bfast_star["change_magnitude"]) +
  tm_raster() +
  tm_layout(title = "BFAST change magnitude (kNDVI)")
