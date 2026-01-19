#!/usr/bin/env Rscript
############################################################
# Script:  05-4_decisions_biostructure_metrics.R
# Project: Burgwald
#
# Purpose:
#   S5 decisions (Bio-Structure): cluster segments in bio-structure space
#   and select representative candidates (segments + points in segments).
#
# Input (productive):
#   - layer0_segments_attrstack_metrics   (S4_signatures, gpkg)
#
# Outputs (productive):
#   - layer0_biostructure_strata_segments (S5_decisions, gpkg)
#   - layer0_biostructure_candidates_pts  (S5_decisions, gpkg)
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
bio_seg_out <- paths[["layer0_biostructure_strata_segments"]]
bio_pts_out <- paths[["layer0_biostructure_candidates_pts"]]

## -------------------------------------------------------------------
## 2b) tmp output folder (NON-productive)
## -------------------------------------------------------------------
tmp_root <- here::here("data", "tmp", "layer0_decisions", "tmp_biostructure_decisions")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Read input
## -------------------------------------------------------------------
segs_sf <- sf::read_sf(attr_file)
stopifnot("segment_id" %in% names(segs_sf))

## -------------------------------------------------------------------
## 4) Bio-Structure decision space (subset from FULL S4 attrstack)
## -------------------------------------------------------------------
bio_cols <- c("chm_p95_mean", "canopy_fraction_mean", "chm_sd_mean")
stopifnot(all(bio_cols %in% names(segs_sf)))

bio_df <- segs_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(segment_id, dplyr::all_of(bio_cols)) %>%
  tidyr::drop_na()

k_bio <- choose_k_silhouette(bio_df[, bio_cols, drop = FALSE], k_range = k_range)$k
bio_x <- scale(as.matrix(bio_df[, bio_cols, drop = FALSE]))
bio_km <- kmeans(bio_x, centers = k_bio, nstart = 50)

bio_dist <- numeric(nrow(bio_x))
for (i in seq_len(nrow(bio_x))) {
  c_id <- bio_km$cluster[i]
  bio_dist[i] <- sqrt(sum((bio_x[i, ] - bio_km$centers[c_id, ])^2))
}

bio_assign <- tibble::tibble(
  segment_id          = bio_df$segment_id,
  bio_stratum         = as.integer(bio_km$cluster),
  bio_dist_to_center  = as.numeric(bio_dist),
  bio_k               = as.integer(k_bio)
)

bio_seg_sf <- segs_sf %>%
  dplyr::left_join(bio_assign, by = "segment_id")

bio_candidates_ids <- bio_assign %>%
  dplyr::group_by(bio_stratum) %>%
  dplyr::slice_min(bio_dist_to_center, n = n_per_stratum, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(segment_id, bio_stratum, bio_dist_to_center)

bio_candidates_seg <- bio_seg_sf %>%
  dplyr::inner_join(bio_candidates_ids, by = "segment_id") %>%
  dplyr::select(segment_id, bio_stratum, bio_dist_to_center, geometry)

bio_candidates_pts <- sf::st_point_on_surface(bio_candidates_seg) %>%
  dplyr::mutate(space = "biostructure") %>%
  dplyr::select(segment_id, space, stratum = bio_stratum, dist_to_center = bio_dist_to_center, geometry)

## -------------------------------------------------------------------
## 5) Write outputs
## -------------------------------------------------------------------
dir.create(dirname(bio_seg_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(bio_pts_out), recursive = TRUE, showWarnings = FALSE)

sf::st_write(bio_seg_sf, bio_seg_out, delete_dsn = TRUE, quiet = TRUE)
sf::st_write(bio_candidates_pts, bio_pts_out, delete_dsn = TRUE, quiet = TRUE)
