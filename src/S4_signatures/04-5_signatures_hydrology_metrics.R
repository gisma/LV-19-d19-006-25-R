#!/usr/bin/env Rscript
############################################################
# Script:  04-3_signatures_hydro_metrics.R
# Project: Burgwald
#
# Purpose:
#   Segment-based hydrological signatures (S4_signatures).
#   Pure zonal aggregation of precomputed hydrology rasters.
#
#   - NO thresholds
#   - NO classification
#   - NO decision logic
#
# Inputs (productive, from outputs.tsv via paths[]):
#   - layer0_segments            (S3_structure, gpkg)
#   - watershed_id_10m           (S2_features, tif)
#   - strahler_order_10m         (S2_features, tif)
#   - flowacc_10m                (S2_features, tif)
#   - dist_to_stream_10m         (S2_features, tif) [optional]
#
# Output (productive, from outputs.tsv via paths[]):
#   - layer0_attr_hydro_metrics  (S4_signatures, gpkg)
#
# Output columns (segment-level):
#   segment_id
#   watershed_major
#   watershed_frac
#   strahler_max
#   strahler_mean
#   strahler_frac_stream
#   flowacc_mean
#   flowacc_p90
#   flowacc_max
#   dist_stream_mean   (if available)
#   dist_stream_min    (if available)
############################################################

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(terra)
  library(dplyr)
  library(tibble)
  library(exactextractr)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))
source(here::here("src", "r-libs", "metrics-fun.R"))

## -------------------------------------------------------------------
## 1) Productive input (from outputs.tsv)
## -------------------------------------------------------------------

seg_file        <- paths[["layer0_segments"]]
watershed_file <- paths[["watershed_id_10m"]]
strahler_file  <- paths[["strahler_order_10m"]]
flowacc_file   <- paths[["flowacc_10m"]]

stopifnot(file.exists(seg_file))
stopifnot(file.exists(watershed_file))
stopifnot(file.exists(strahler_file))
stopifnot(file.exists(flowacc_file))

has_dist_stream <- "dist_to_stream_10m" %in% names(paths) &&
  file.exists(paths[["dist_to_stream_10m"]])

if (has_dist_stream) {
  dist_stream_file <- paths[["dist_to_stream_10m"]]
  stopifnot(file.exists(dist_stream_file))
}

## -------------------------------------------------------------------
## 2) Productive outputs (from outputs.tsv)
## -------------------------------------------------------------------

out_file <- paths[["layer0_attr_hydro_metrics"]]

## -------------------------------------------------------------------
## 2b) tmp output folder (NON-productive)
## -------------------------------------------------------------------

tmp_root <- here::here("data", "tmp", "layer0_segments", "tmp_hydro_metrics")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Read inputs
## -------------------------------------------------------------------

segments_sf <- sf::read_sf(seg_file)
stopifnot("segment_id" %in% names(segments_sf))

watershed_r <- terra::rast(watershed_file)
strahler_r  <- terra::rast(strahler_file)
flowacc_r   <- terra::rast(flowacc_file)

if (has_dist_stream) {
  dist_stream_r <- terra::rast(dist_stream_file)
}

## -------------------------------------------------------------------
## 4) Harmonise CRS
## -------------------------------------------------------------------

ref_crs <- terra::crs(watershed_r, proj = TRUE)
stopifnot(!is.na(ref_crs), ref_crs != "")

if (!identical(sf::st_crs(segments_sf)$wkt, ref_crs)) {
  segments_sf <- sf::st_transform(segments_sf, ref_crs)
}

if (any(!sf::st_is_valid(segments_sf))) {
  segments_sf <- sf::st_make_valid(segments_sf)
}

## -------------------------------------------------------------------
## 5) Segment-wise zonal extraction
## -------------------------------------------------------------------

message("Extracting hydrological metrics per segment ...")

# --- Watershed: dominant class + fraction ----------------------------

watershed_major <- exactextractr::exact_extract(
  watershed_r,
  segments_sf,
  function(values, coverage_fraction) {
    values <- values[!is.na(values)]
    if (length(values) == 0) return(NA_integer_)
    as.integer(names(sort(table(values), decreasing = TRUE)[1]))
  }
)

watershed_frac <- exactextractr::exact_extract(
  watershed_r,
  segments_sf,
  function(values, coverage_fraction) {
    values <- values[!is.na(values)]
    if (length(values) == 0) return(NA_real_)
    tab <- table(values)
    max(tab) / sum(tab)
  }
)

# --- Strahler --------------------------------------------------------

strahler_max  <- exactextractr::exact_extract(strahler_r, segments_sf, "max")
strahler_mean <- exactextractr::exact_extract(strahler_r, segments_sf, "mean")

strahler_frac_stream <- exactextractr::exact_extract(
  strahler_r,
  segments_sf,
  function(values, coverage_fraction) {
    values <- values[!is.na(values)]
    if (length(values) == 0) return(NA_real_)
    mean(values > 0)
  }
)

# --- Flow accumulation ----------------------------------------------

flowacc_mean <- exactextractr::exact_extract(flowacc_r, segments_sf, "mean")
flowacc_p90  <- exactextractr::exact_extract(flowacc_r, segments_sf, "quantile", probs = 0.9)
flowacc_max  <- exactextractr::exact_extract(flowacc_r, segments_sf, "max")

# --- Distance to stream (optional) ----------------------------------

if (has_dist_stream) {
  dist_stream_mean <- exactextractr::exact_extract(dist_stream_r, segments_sf, "mean")
  dist_stream_min  <- exactextractr::exact_extract(dist_stream_r, segments_sf, "min")
} else {
  dist_stream_mean <- rep(NA_real_, nrow(segments_sf))
  dist_stream_min  <- rep(NA_real_, nrow(segments_sf))
}

## -------------------------------------------------------------------
## 6) Assemble table + geometry
## -------------------------------------------------------------------

hydro_df <- tibble::tibble(
  segment_id             = segments_sf$segment_id,
  watershed_major        = as.integer(watershed_major),
  watershed_frac         = as.numeric(watershed_frac),
  strahler_max            = as.numeric(strahler_max),
  strahler_mean           = as.numeric(strahler_mean),
  strahler_frac_stream    = as.numeric(strahler_frac_stream),
  flowacc_mean             = as.numeric(flowacc_mean),
  flowacc_p90              = as.numeric(flowacc_p90),
  flowacc_max              = as.numeric(flowacc_max),
  dist_stream_mean         = as.numeric(dist_stream_mean),
  dist_stream_min          = as.numeric(dist_stream_min)
)

out_sf <- segments_sf %>%
  dplyr::select(segment_id, geometry) %>%
  dplyr::left_join(hydro_df, by = "segment_id")

## -------------------------------------------------------------------
## 7) Write productive output
## -------------------------------------------------------------------

out_dir <- dirname(out_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sf::st_write(out_sf, out_file, delete_dsn = TRUE, quiet = TRUE)
