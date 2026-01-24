#!/usr/bin/env Rscript

############################################################
# Script: 02-4_processing_bfast-classification.R
#
# Purpose:
# Detect per-pixel structural breaks in vegetation dynamics using
# BFAST on a kNDVI time series derived from a Sentinel-2 NetCDF cube.
#
# Contract:
# - Input cube path and output result path come from metadata/outputs.tsv
#   compiled in setup as `paths`.
# - No hard-coded output file paths.
#
# Context (what this is, technically)
# ----------------------------------
# This script runs a *per-pixel* time-series break detection on a pre-built
# Sentinel-2 data cube stored as NetCDF (gdalcubes cube). For each pixel:
#   - extract the temporal sequence of reflectances (here: B08, B04),
#   - derive kNDVI as a scalar vegetation signal per timestamp,
#   - run bfast::bfastmonitor() to detect a structural break,
#   - return two numeric outputs per pixel:
#       (1) change_date      : breakpoint position in the ts object
#       (2) change_magnitude : estimated magnitude at the breakpoint
# The result is written back as a NetCDF raster cube with two variables
# ("change_date", "change_magnitude") aligned to the input cube grid.
#
# Key dependencies / system requirements
# --------------------------------------
# R packages:
#   - gdalcubes  : reads NetCDF cube and executes reduce_time() in parallel
#   - bfast      : bfastmonitor() for breakpoint + magnitude
# Optional for quick visualization:
#   - stars      : reads NetCDF result as stars object
#   - tmap       : quick raster viewing / plotting
#
# External / system side:
#   - gdalcubes depends on GDAL (and often netCDF/HDF backends) installed on the system.
#   - You need a working GDAL/netCDF stack such that gdalcubes can open the cube_file
#     and write NetCDF output.
#
# Important implementation note (gdalcubes worker processes)
# ---------------------------------------------------------
# gdalcubes executes the reducer FUN in worker processes.
# Therefore:
#   - any package used inside FUN must be loaded inside FUN
#     (hence: library(bfast) inside the function).
#   - only base objects shipped into FUN are available; avoid relying on
#     global environment state.
############################################################

source(here::here("src", "_core", "01-setup-burgwald.R"))

# ---- Preconditions ---------------------------------------------------------
# These checks stop early if the required packages are missing.
# This is a hard requirement because reduce_time() and bfastmonitor() are central.
if (!requireNamespace("gdalcubes", quietly = TRUE)) stop("Missing package: gdalcubes")
if (!requireNamespace("bfast", quietly = TRUE)) stop("Missing package: bfast")

# Optional (only used for visualization at the end)
# If these packages are absent, the script still produces the NetCDF output.
has_stars <- requireNamespace("stars", quietly = TRUE)
has_tmap  <- requireNamespace("tmap",  quietly = TRUE)

# gdalcubes parallelism (adjust to your machine)
# parallel = number of worker processes used by gdalcubes (not threads in one process).
gdalcubes::gdalcubes_options(parallel = 12)

# ---- Canonical I/O (outputs.tsv) ------------------------------------------
# cube_file:
#   NetCDF cube containing the time series per pixel.
#   This cube is assumed to contain at least bands "B08" and "B04".
# bfast_file:
#   Output NetCDF where each pixel stores the reducer output variables:
#     - change_date
#     - change_magnitude
cube_file <- paths[["s2_gdalcubes_cube_2018_2022_all"]]
if (is.null(cube_file)) stop("Missing outputs.tsv entry for key: s2_gdalcubes_cube_2018_2022_all")
if (!file.exists(cube_file)) stop("Input cube does not exist: ", cube_file)

bfast_file <- paths[["s2_bfast_kndvi_results"]]
if (is.null(bfast_file)) stop("Missing outputs.tsv entry for key: s2_bfast_kndvi_results")

# Ensure output directory exists (only parent dir creation; no extra structure logic)
#fs::dir_create(dirname(bfast_file), recurse = TRUE)

# ---- BFAST reducer notes ---------------------------------------------------
# gdalcubes executes the reducer FUN in an R worker process.
# Therefore, any packages used inside FUN must be loaded inside FUN.

system.time({
  # ncdf_cube(cube_file):
  #   Opens the NetCDF cube as a gdalcubes cube object (spatio-temporal).
  gdalcubes::ncdf_cube(cube_file) |>
    # reduce_time():
    #   Applies FUN along the time dimension for each pixel.
    #   FUN receives x: a matrix of dimension (bands x time) for one pixel.
    #   It must return a numeric vector of fixed length == length(names).
    gdalcubes::reduce_time(
      names = c("change_date", "change_magnitude"),
      FUN = function(x) {
        # x: matrix (bands x time) for one pixel
        
        # Load inside worker process
        # Required because the reducer runs in a separate process that does not
        # inherit the caller's attached packages.
        suppressPackageStartupMessages(library(bfast))
        
        # Guard: required bands
        # If the cube does not provide B08 and B04 for this pixel, return NA outputs.
        if (!all(c("B08", "B04") %in% rownames(x))) return(c(NA_real_, NA_real_))
        
        # Convert digital numbers to reflectance [0..1] by dividing by 10000.
        # Assumption: cube stores scaled integer reflectances (S2 common convention).
        nir <- x["B08", ] / 10000
        red <- x["B04", ] / 10000
        
        # kNDVI definition used here:
        #   kNDVI = tanh( NDVI^2 )
        # where NDVI = (nir - red) / (nir + red)
        # This stabilizes extremes and yields a bounded vegetation signal.
        kndvi <- tanh(((nir - red) / (nir + red))^2)
        
        # If no usable data -> NA
        # This checks the derived signal, not the raw bands.
        if (all(is.na(kndvi))) return(c(NA_real_, NA_real_))
        
        # Build an R ts object for BFAST.
        # frequency=12 implies monthly series.
        # start=c(2018,1) anchors index and aligns breakpoint reporting to a timeline.
        #
        # NOTE: This assumes the time axis in the cube corresponds to a monthly
        # cadence aligned with this ts definition. If the cube is not monthly or
        # has missing months, interpretation of breakpoint index becomes ambiguous.
        kndvi_ts <- ts(kndvi, start = c(2018, 1), frequency = 12)
        
        # Run bfastmonitor:
        # - start=c(2020,1): monitoring period begins Jan 2020
        # - history="all"  : use all available pre-monitoring history
        # - level=0.01     : significance level; smaller -> fewer detections (stricter)
        #
        # Output extracted:
        # - res$breakpoint : index/time location of detected break (ts units)
        # - res$magnitude  : estimated magnitude of change at break
        out <- tryCatch({
          res <- bfastmonitor(
            kndvi_ts,
            start   = c(2020, 1),
            history = "all",
            level   = 0.01
          )
          c(res$breakpoint, res$magnitude)
        }, error = function(e) {
          # Pixel-level failures (e.g., too few observations, singularities, etc.)
          # are handled by returning NA outputs for that pixel.
          c(NA_real_, NA_real_)
        })
        
        out
      }
    ) |>
    # write_ncdf():
    #   Writes the resulting cube (2 variables) to a NetCDF file.
    # overwrite=TRUE replaces an existing result file.
    gdalcubes::write_ncdf(bfast_file, overwrite = TRUE)
})

cat("BFAST results written to:\n  ", bfast_file, "\n")

# ---- Optional quick look (only if packages exist) --------------------------
# This block is non-essential. It:
#   - reads the NetCDF result into a stars object,
#   - displays change_date in interactive view mode,
#   - plots change_magnitude in static plot mode.
if (has_stars && has_tmap) {
  bfast_star <- gdalcubes::ncdf_cube(bfast_file) |>
    stars::st_as_stars()
  
  tmap::tmap_mode("view")
  print(
    tmap::tm_shape(bfast_star["change_date"]) +
      tmap::tm_raster() +
      tmap::tm_layout(title = "BFAST change date (kNDVI)")
  )
  
  tmap::tmap_mode("plot")
  print(
    tmap::tm_shape(bfast_star["change_magnitude"]) +
      tmap::tm_raster() +
      tmap::tm_layout(title = "BFAST change magnitude (kNDVI)")
  )
}
