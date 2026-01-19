#!/usr/bin/env Rscript
############################################################
# Script:  05-5_decisions_coverage_metrics.R
# Project: Burgwald
#
# Purpose:
#   S5 decisions (Coverage/LUCC): cluster segments in cover-fraction space
#   and select representative candidates (segments + points in segments).
#
# Input (productive):
#   - layer0_segments_attrstack_metrics   (S4_signatures, gpkg)
#
# Outputs (productive):
#   - layer0_coverage_strata_segments     (S5_decisions, gpkg)
#   - layer0_coverage_candidates_pts      (S5_decisions, gpkg)
############################################################

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stats)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))
source(here::here("src", "r-libs", "metrics-fun.R"))

## -------------------------------------------------------------------
## 0) Parameters
## -------------------------------------------------------------------
n_per_stratum <- 3L
k_range <- 2:12
set.seed(42)

## -------------------------------------------------------------------
## 1) Productive input
## -------------------------------------------------------------------
attr_file <- paths[["layer0_segments_attrstack_metrics"]]
stopifnot(file.exists(attr_file))

## -------------------------------------------------------------------
## 2) Productive outputs
## -------------------------------------------------------------------
cov_seg_out <- paths[["layer0_coverage_strata_segments"]]
cov_pts_out <- paths[["layer0_coverage_candidates_pts"]]

## -------------------------------------------------------------------
## 2b) tmp output folder (NON-productive)
## -------------------------------------------------------------------
tmp_root <- here::here("data", "tmp", "layer0_decisions", "tmp_coverage_decisions")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Read input
## -------------------------------------------------------------------
segs_sf <- sf::read_sf(attr_file)
stopifnot("segment_id" %in% names(segs_sf))

## -------------------------------------------------------------------
## 4) Coverage decision space (subset from FULL S4 attrstack)
## -------------------------------------------------------------------
cover_cols <- c("cov_forest", "cov_agri", "cov_grass", "cov_built", "cov_water")
stopifnot(all(cover_cols %in% names(segs_sf)))

cov_df <- segs_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(segment_id, dplyr::all_of(cover_cols)) %>%
  tidyr::drop_na()

k_cov <- choose_k_silhouette(cov_df[, cover_cols, drop = FALSE], k_range = k_range)$k
cov_x <- scale(as.matrix(cov_df[, cover_cols, drop = FALSE]))
cov_km <- kmeans(cov_x, centers = k_cov, nstart = 50)

cov_dist <- numeric(nrow(cov_x))
for (i in seq_len(nrow(cov_x))) {
  c_id <- cov_km$cluster[i]
  cov_dist[i] <- sqrt(sum((cov_x[i, ] - cov_km$centers[c_id, ])^2))
}

cov_assign <- tibble::tibble(
  segment_id          = cov_df$segment_id,
  cov_stratum         = as.integer(cov_km$cluster),
  cov_dist_to_center  = as.numeric(cov_dist),
  cov_k               = as.integer(k_cov)
)

cov_seg_sf <- segs_sf %>%
  dplyr::left_join(cov_assign, by = "segment_id")

cov_candidates_ids <- cov_assign %>%
  dplyr::group_by(cov_stratum) %>%
  dplyr::slice_min(cov_dist_to_center, n = n_per_stratum, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(segment_id, cov_stratum, cov_dist_to_center)

cov_candidates_seg <- cov_seg_sf %>%
  dplyr::inner_join(cov_candidates_ids, by = "segment_id") %>%
  dplyr::select(segment_id, cov_stratum, cov_dist_to_center, geometry)

cov_candidates_pts <- sf::st_point_on_surface(cov_candidates_seg) %>%
  dplyr::mutate(space = "coverage") %>%
  dplyr::select(segment_id, space, stratum = cov_stratum, dist_to_center = cov_dist_to_center, geometry)

## -------------------------------------------------------------------
## 5) Write outputs
## -------------------------------------------------------------------
dir.create(dirname(cov_seg_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cov_pts_out), recursive = TRUE, showWarnings = FALSE)

sf::st_write(cov_seg_sf, cov_seg_out, delete_dsn = TRUE, quiet = TRUE)
sf::st_write(cov_candidates_pts, cov_pts_out, delete_dsn = TRUE, quiet = TRUE)
