#!/usr/bin/env Rscript
############################################################
# Script:  05-3_decisions_hydro_metrics.R
# Project: Burgwald
#
# Purpose:
#   S5 decisions (Hydro): cluster segments in hydro-space and select
#   representative candidates (segments + points in segments).
#
# Input (productive):
#   - layer0_segments_attrstack_metrics   (S4_signatures, gpkg)
#
# Outputs (productive):
#   - layer0_hydro_strata_segments        (S5_decisions, gpkg)
#   - layer0_hydro_candidates_pts         (S5_decisions, gpkg)
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
hydro_seg_out <- paths[["layer0_hydro_strata_segments"]]
hydro_pts_out <- paths[["layer0_hydro_candidates_pts"]]

## -------------------------------------------------------------------
## 2b) tmp output folder (NON-productive)
## -------------------------------------------------------------------
tmp_root <- here::here("data", "tmp", "layer0_decisions", "tmp_hydro_decisions")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Read input
## -------------------------------------------------------------------
segs_sf <- sf::read_sf(attr_file)
stopifnot("segment_id" %in% names(segs_sf))

## -------------------------------------------------------------------
## 4) Hydro decision space (subset from FULL S4 attrstack)
## -------------------------------------------------------------------
hydro_cols <- c("strahler_max", "flowacc_p90", "dist_stream_mean")
stopifnot(all(hydro_cols %in% names(segs_sf)))

hydro_df <- segs_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(segment_id, dplyr::all_of(hydro_cols)) %>%
  tidyr::drop_na()

k_hydro <- choose_k_silhouette(hydro_df[, hydro_cols, drop = FALSE], k_range = k_range)$k
hydro_x <- scale(as.matrix(hydro_df[, hydro_cols, drop = FALSE]))
hydro_km <- kmeans(hydro_x, centers = k_hydro, nstart = 50)

hydro_dist <- numeric(nrow(hydro_x))
for (i in seq_len(nrow(hydro_x))) {
  c_id <- hydro_km$cluster[i]
  hydro_dist[i] <- sqrt(sum((hydro_x[i, ] - hydro_km$centers[c_id, ])^2))
}

hydro_assign <- tibble::tibble(
  segment_id           = hydro_df$segment_id,
  hydro_stratum        = as.integer(hydro_km$cluster),
  hydro_dist_to_center = as.numeric(hydro_dist),
  hydro_k              = as.integer(k_hydro)
)

hydro_seg_sf <- segs_sf %>%
  dplyr::left_join(hydro_assign, by = "segment_id")

hydro_candidates_ids <- hydro_assign %>%
  dplyr::group_by(hydro_stratum) %>%
  dplyr::slice_min(hydro_dist_to_center, n = n_per_stratum, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(segment_id, hydro_stratum, hydro_dist_to_center)

hydro_candidates_seg <- hydro_seg_sf %>%
  dplyr::inner_join(hydro_candidates_ids, by = "segment_id") %>%
  dplyr::select(segment_id, hydro_stratum, hydro_dist_to_center, geometry)

hydro_candidates_pts <- sf::st_point_on_surface(hydro_candidates_seg) %>%
  dplyr::mutate(space = "hydro") %>%
  dplyr::select(segment_id, space, stratum = hydro_stratum, dist_to_center = hydro_dist_to_center, geometry)

## -------------------------------------------------------------------
## 5) Write outputs
## -------------------------------------------------------------------
dir.create(dirname(hydro_seg_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(hydro_pts_out), recursive = TRUE, showWarnings = FALSE)

sf::st_write(hydro_seg_sf, hydro_seg_out, delete_dsn = TRUE, quiet = TRUE)
sf::st_write(hydro_candidates_pts, hydro_pts_out, delete_dsn = TRUE, quiet = TRUE)
