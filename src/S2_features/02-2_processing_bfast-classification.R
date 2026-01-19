#!/usr/bin/env Rscript

############################################################
# Script: 02-2_processing_bfast-classification.R
#
# Purpose:
# Detect per-pixel structural breaks in vegetation dynamics using
# BFAST on a kNDVI time series derived from a Sentinel-2 NetCDF cube.
#
# Contract:
# - Input cube path and output result path come from metadata/outputs.tsv
#   compiled in setup as `paths`.
# - No hard-coded output file paths.
############################################################

source(here::here("src", "_core", "01-setup-burgwald.R"))

# ---- Preconditions ---------------------------------------------------------

if (!requireNamespace("gdalcubes", quietly = TRUE)) stop("Missing package: gdalcubes")
if (!requireNamespace("bfast", quietly = TRUE)) stop("Missing package: bfast")

# Optional (only used for visualization at the end)
has_stars <- requireNamespace("stars", quietly = TRUE)
has_tmap  <- requireNamespace("tmap",  quietly = TRUE)

# gdalcubes parallelism (adjust to your machine)
gdalcubes::gdalcubes_options(parallel = 12)

# ---- Canonical I/O (outputs.tsv) ------------------------------------------

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
  gdalcubes::ncdf_cube(cube_file) |>
    gdalcubes::reduce_time(
      names = c("change_date", "change_magnitude"),
      FUN = function(x) {
        # x: matrix (bands x time) for one pixel
        
        # Load inside worker process
        suppressPackageStartupMessages(library(bfast))
        
        # Guard: required bands
        if (!all(c("B08", "B04") %in% rownames(x))) return(c(NA_real_, NA_real_))
        
        nir <- x["B08", ] / 10000
        red <- x["B04", ] / 10000
        
        kndvi <- tanh(((nir - red) / (nir + red))^2)
        
        # If no usable data -> NA
        if (all(is.na(kndvi))) return(c(NA_real_, NA_real_))
        
        kndvi_ts <- ts(kndvi, start = c(2018, 1), frequency = 12)
        
        out <- tryCatch({
          res <- bfastmonitor(
            kndvi_ts,
            start   = c(2020, 1),
            history = "all",
            level   = 0.01
          )
          c(res$breakpoint, res$magnitude)
        }, error = function(e) {
          c(NA_real_, NA_real_)
        })
        
        out
      }
    ) |>
    gdalcubes::write_ncdf(bfast_file, overwrite = TRUE)
})

cat("BFAST results written to:\n  ", bfast_file, "\n")

# ---- Optional quick look (only if packages exist) --------------------------

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
