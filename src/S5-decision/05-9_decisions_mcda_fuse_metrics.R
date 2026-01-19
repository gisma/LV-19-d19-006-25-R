#!/usr/bin/env Rscript
############################################################
# Script:  05-9_decisions_mcda_fuse_metrics.R
# Project: Burgwald
#
# Purpose:
#   S5 MCDA fusion: combine per-domain representativeness distances
#   into a single weighted score on segment level, then derive points
#   in selected segments.
#
# Inputs (productive):
#   - layer0_it_strata_segments
#   - layer0_physio_strata_segments
#   - layer0_hydro_strata_segments
#   - layer0_biostructure_strata_segments
#   - layer0_coverage_strata_segments
#
# Outputs (productive):
#   - layer0_mcda_segments_scored
#   - layer0_mcda_candidates_pts
############################################################

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(dplyr)
  library(tidyr)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))
source(here::here("src", "r-libs", "metrics-fun.R"))

## -------------------------------------------------------------------
## 0) Parameters (weights + selection)
## -------------------------------------------------------------------
w_it       <- 1.0
w_physio   <- 1.0
w_hydro    <- 1.0
w_bio      <- 1.0
w_coverage <- 1.0

top_n <- 50L
set.seed(42)

## -------------------------------------------------------------------
## 1) Productive inputs
## -------------------------------------------------------------------
it_seg_file     <- paths[["layer0_it_strata_segments"]]
physio_seg_file <- paths[["layer0_physio_strata_segments"]]
hydro_seg_file  <- paths[["layer0_hydro_strata_segments"]]
bio_seg_file    <- paths[["layer0_biostructure_strata_segments"]]
cov_seg_file    <- paths[["layer0_coverage_strata_segments"]]

stopifnot(file.exists(it_seg_file))
stopifnot(file.exists(physio_seg_file))
stopifnot(file.exists(hydro_seg_file))
stopifnot(file.exists(bio_seg_file))
stopifnot(file.exists(cov_seg_file))

## -------------------------------------------------------------------
## 2) Productive outputs
## -------------------------------------------------------------------
mcda_seg_out <- paths[["layer0_mcda_segments_scored"]]
mcda_pts_out <- paths[["layer0_mcda_candidates_pts"]]

## -------------------------------------------------------------------
## 2b) tmp output folder (NON-productive)
## -------------------------------------------------------------------
tmp_root <- here::here("data", "tmp", "layer0_decisions", "tmp_mcda_fuse")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Read inputs (geometry backbone: IT layer)
## -------------------------------------------------------------------
it_sf     <- sf::read_sf(it_seg_file)
physio_sf <- sf::read_sf(physio_seg_file)
hydro_sf  <- sf::read_sf(hydro_seg_file)
bio_sf    <- sf::read_sf(bio_seg_file)
cov_sf    <- sf::read_sf(cov_seg_file)

stopifnot("segment_id" %in% names(it_sf))

## expect distance columns (from scripts)
stopifnot("it_dist_to_center" %in% names(it_sf))
stopifnot("physio_dist_to_center" %in% names(physio_sf))
stopifnot("hydro_dist_to_center" %in% names(hydro_sf))
stopifnot("bio_dist_to_center" %in% names(bio_sf))
stopifnot("cov_dist_to_center" %in% names(cov_sf))

## -------------------------------------------------------------------
## 4) Join distances onto one segment table
## -------------------------------------------------------------------
base_sf <- it_sf %>%
  dplyr::select(segment_id, geometry, it_stratum, it_dist_to_center)

physio_df <- physio_sf %>% sf::st_drop_geometry() %>%
  dplyr::select(segment_id, physio_stratum, physio_dist_to_center)

hydro_df <- hydro_sf %>% sf::st_drop_geometry() %>%
  dplyr::select(segment_id, hydro_stratum, hydro_dist_to_center)

bio_df <- bio_sf %>% sf::st_drop_geometry() %>%
  dplyr::select(segment_id, bio_stratum, bio_dist_to_center)

cov_df <- cov_sf %>% sf::st_drop_geometry() %>%
  dplyr::select(segment_id, cov_stratum, cov_dist_to_center)

seg_sf <- base_sf %>%
  dplyr::left_join(physio_df, by = "segment_id") %>%
  dplyr::left_join(hydro_df,  by = "segment_id") %>%
  dplyr::left_join(bio_df,    by = "segment_id") %>%
  dplyr::left_join(cov_df,    by = "segment_id")

## -------------------------------------------------------------------
## 5) Normalize distances (0..1) + MCDA score
## -------------------------------------------------------------------
norm01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[2] == r[1]) return(rep(NA_real_, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

seg_sf <- seg_sf %>%
  dplyr::mutate(
    it_d_n       = norm01(it_dist_to_center),
    physio_d_n   = norm01(physio_dist_to_center),
    hydro_d_n    = norm01(hydro_dist_to_center),
    bio_d_n      = norm01(bio_dist_to_center),
    coverage_d_n = norm01(cov_dist_to_center),
    mcda_score   = (w_it       * it_d_n) +
      (w_physio   * physio_d_n) +
      (w_hydro    * hydro_d_n) +
      (w_bio      * bio_d_n) +
      (w_coverage * coverage_d_n)
  )

## -------------------------------------------------------------------
## 6) Select top-N segments + derive points in segments
## -------------------------------------------------------------------
sel_sf <- seg_sf %>%
  dplyr::filter(!is.na(mcda_score)) %>%
  dplyr::arrange(mcda_score) %>%
  dplyr::slice_head(n = top_n)

sel_pts <- sf::st_point_on_surface(sel_sf) %>%
  dplyr::mutate(space = "mcda") %>%
  dplyr::select(segment_id, space, mcda_score, geometry)

## -------------------------------------------------------------------
## 7) Write outputs
## -------------------------------------------------------------------
dir.create(dirname(mcda_seg_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(mcda_pts_out), recursive = TRUE, showWarnings = FALSE)

sf::st_write(seg_sf, mcda_seg_out, delete_dsn = TRUE, quiet = TRUE)
sf::st_write(sel_pts, mcda_pts_out, delete_dsn = TRUE, quiet = TRUE)
