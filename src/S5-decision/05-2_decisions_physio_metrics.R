#!/usr/bin/env Rscript
############################################################
# Script:  05-2_decisions_physio_metrics.R
# Project: Burgwald
#
# Purpose:
#   S5 decisions (Physio): cluster segments in physio-space and select
#   representative candidates (segments + points in segments).
#
# Input (productive):
#   - layer0_segments_attrstack_metrics   (S4_signatures, gpkg)
#
# Outputs (productive):
#   - layer0_physio_strata_segments       (S5_decisions, gpkg)
#   - layer0_physio_candidates_pts        (S5_decisions, gpkg)
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
physio_seg_out <- paths[["layer0_physio_strata_segments"]]
physio_pts_out <- paths[["layer0_physio_candidates_pts"]]

## -------------------------------------------------------------------
## 2b) tmp output folder (NON-productive)
## -------------------------------------------------------------------
tmp_root <- here::here("data", "tmp", "layer0_decisions", "tmp_physio_decisions")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Read input
## -------------------------------------------------------------------
segs_sf <- sf::read_sf(attr_file)
stopifnot("segment_id" %in% names(segs_sf))

## -------------------------------------------------------------------
## 4) Physio decision space (subset from FULL S4 attrstack)
## -------------------------------------------------------------------
physio_cols <- c("elev_mean", "slope_mean_deg", "southness_mean")
stopifnot(all(physio_cols %in% names(segs_sf)))

physio_df <- segs_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(segment_id, dplyr::all_of(physio_cols)) %>%
  tidyr::drop_na()

k_physio <- choose_k_silhouette(physio_df[, physio_cols, drop = FALSE], k_range = k_range)$k
physio_x <- scale(as.matrix(physio_df[, physio_cols, drop = FALSE]))
physio_km <- kmeans(physio_x, centers = k_physio, nstart = 50)

physio_dist <- numeric(nrow(physio_x))
for (i in seq_len(nrow(physio_x))) {
  c_id <- physio_km$cluster[i]
  physio_dist[i] <- sqrt(sum((physio_x[i, ] - physio_km$centers[c_id, ])^2))
}

physio_assign <- tibble::tibble(
  segment_id            = physio_df$segment_id,
  physio_stratum        = as.integer(physio_km$cluster),
  physio_dist_to_center = as.numeric(physio_dist),
  physio_k              = as.integer(k_physio)
)

physio_seg_sf <- segs_sf %>%
  dplyr::left_join(physio_assign, by = "segment_id")

physio_candidates_ids <- physio_assign %>%
  dplyr::group_by(physio_stratum) %>%
  dplyr::slice_min(physio_dist_to_center, n = n_per_stratum, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(segment_id, physio_stratum, physio_dist_to_center)

physio_candidates_seg <- physio_seg_sf %>%
  dplyr::inner_join(physio_candidates_ids, by = "segment_id") %>%
  dplyr::select(segment_id, physio_stratum, physio_dist_to_center, geometry)

physio_candidates_pts <- sf::st_point_on_surface(physio_candidates_seg) %>%
  dplyr::mutate(space = "physio") %>%
  dplyr::select(segment_id, space, stratum = physio_stratum, dist_to_center = physio_dist_to_center, geometry)

## -------------------------------------------------------------------
## 5) Write outputs
## -------------------------------------------------------------------
dir.create(dirname(physio_seg_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(physio_pts_out), recursive = TRUE, showWarnings = FALSE)

sf::st_write(physio_seg_sf, physio_seg_out, delete_dsn = TRUE, quiet = TRUE)
sf::st_write(physio_candidates_pts, physio_pts_out, delete_dsn = TRUE, quiet = TRUE)
