#!/usr/bin/env Rscript
############################################################
# Script:  04-4_signatures_biostructure_metrics.R
# Project: Burgwald
#
# Purpose:
#   Segment-based bio-structure signatures (S4_signatures) from
#   precomputed raster products (assumed generated in S2).
#
#   - NO clustering
#   - NO candidates
#   - NO decision logic
#
# Inputs (productive, from outputs.tsv via paths[]):
#   - layer0_segments                 (S3_structure, gpkg)
#   - chm_mean_10m                    (S2_features, tif)
#   - chm_p95_10m                     (S2_features, tif)
#   - chm_sd_10m                      (S2_features, tif)
#   - canopy_fraction_10m             (S2_features, tif)
#
# Output (productive, from outputs.tsv via paths[]):
#   - layer0_attr_biostructure_metrics (S4_signatures, gpkg)
#
# Output columns (segment-level):
#   segment_id
#   chm_mean_mean
#   chm_p95_mean
#   chm_sd_mean
#   canopy_fraction_mean
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
chm_mean_file   <- paths[["chm_mean_10m"]]
chm_p95_file    <- paths[["chm_p95_10m"]]
chm_sd_file     <- paths[["chm_sd_10m"]]
can_frac_file   <- paths[["canopy_fraction_10m"]]

stopifnot(file.exists(seg_file))
stopifnot(file.exists(chm_mean_file))
stopifnot(file.exists(chm_p95_file))
stopifnot(file.exists(chm_sd_file))
stopifnot(file.exists(can_frac_file))

## -------------------------------------------------------------------
## 2) Productive outputs (from outputs.tsv)
## -------------------------------------------------------------------

out_file <- paths[["layer0_attr_biostructure_metrics"]]

## -------------------------------------------------------------------
## 2b) tmp output folder (NON-productive)
## -------------------------------------------------------------------

tmp_root <- here::here("data", "tmp", "layer0_segments", "tmp_biostructure_metrics")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Read inputs
## -------------------------------------------------------------------

segments_sf <- sf::read_sf(seg_file)
stopifnot("segment_id" %in% names(segments_sf))

chm_mean_r <- terra::rast(chm_mean_file)
chm_p95_r  <- terra::rast(chm_p95_file)
chm_sd_r   <- terra::rast(chm_sd_file)
can_frac_r <- terra::rast(can_frac_file)

## -------------------------------------------------------------------
## 4) Harmonise CRS
## -------------------------------------------------------------------

ref_crs <- terra::crs(chm_mean_r, proj = TRUE)
stopifnot(!is.na(ref_crs), ref_crs != "")

if (!identical(sf::st_crs(segments_sf)$wkt, ref_crs)) {
  segments_sf <- sf::st_transform(segments_sf, ref_crs)
}

if (any(!sf::st_is_valid(segments_sf))) {
  segments_sf <- sf::st_make_valid(segments_sf)
}

## -------------------------------------------------------------------
## 5) Segment-wise zonal extraction (pure signatures)
## -------------------------------------------------------------------

bio_df <- tibble::tibble(
  segment_id            = segments_sf$segment_id,
  chm_mean_mean         = as.numeric(exactextractr::exact_extract(chm_mean_r, segments_sf, "mean")),
  chm_p95_mean          = as.numeric(exactextractr::exact_extract(chm_p95_r,  segments_sf, "mean")),
  chm_sd_mean           = as.numeric(exactextractr::exact_extract(chm_sd_r,   segments_sf, "mean")),
  canopy_fraction_mean  = as.numeric(exactextractr::exact_extract(can_frac_r, segments_sf, "mean"))
)

out_sf <- segments_sf %>%
  dplyr::select(segment_id, geometry) %>%
  dplyr::left_join(bio_df, by = "segment_id")

## -------------------------------------------------------------------
## 6) Write productive output
## -------------------------------------------------------------------

out_dir <- dirname(out_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sf::st_write(out_sf, out_file, delete_dsn = TRUE, quiet = TRUE)
