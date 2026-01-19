#!/usr/bin/env Rscript
############################################################
# Script:  05-1_decisions_it_metrics.R
# Project: Burgwald
#
# Purpose:
#   S5 decisions (IT): cluster segments in IT-space and select
#   representative candidates (segments + points in segments).
#
# Input (productive):
#   - layer0_segments_attrstack_metrics   (S4_signatures, gpkg)
#
# Outputs (productive):
#   - layer0_it_strata_segments           (S5_decisions, gpkg)
#   - layer0_it_candidates_pts            (S5_decisions, gpkg)
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
it_seg_out <- paths[["layer0_it_strata_segments"]]
it_pts_out <- paths[["layer0_it_candidates_pts"]]

## -------------------------------------------------------------------
## 2b) tmp output folder (NON-productive)
## -------------------------------------------------------------------
tmp_root <- here::here("data", "tmp", "layer0_decisions", "tmp_it_decisions")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Read input
## -------------------------------------------------------------------
segs_sf <- sf::read_sf(attr_file)
stopifnot("segment_id" %in% names(segs_sf))

## -------------------------------------------------------------------
## 4) IT decision space (subset from FULL S4 attrstack)
## -------------------------------------------------------------------
it_cols <- c("H_norm", "U")
if (!all(it_cols %in% names(segs_sf))) it_cols <- c("H", "U")
stopifnot(all(it_cols %in% names(segs_sf)))

it_df <- segs_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(segment_id, dplyr::all_of(it_cols)) %>%
  tidyr::drop_na()

k_it <- choose_k_silhouette(it_df[, it_cols, drop = FALSE], k_range = k_range)$k
it_x <- scale(as.matrix(it_df[, it_cols, drop = FALSE]))
it_km <- kmeans(it_x, centers = k_it, nstart = 50)

## distances to own cluster center
it_dist <- numeric(nrow(it_x))
for (i in seq_len(nrow(it_x))) {
  c_id <- it_km$cluster[i]
  it_dist[i] <- sqrt(sum((it_x[i, ] - it_km$centers[c_id, ])^2))
}

it_assign <- tibble::tibble(
  segment_id        = it_df$segment_id,
  it_stratum        = as.integer(it_km$cluster),
  it_dist_to_center = as.numeric(it_dist),
  it_k              = as.integer(k_it)
)

it_seg_sf <- segs_sf %>%
  dplyr::left_join(it_assign, by = "segment_id")

it_candidates_ids <- it_assign %>%
  dplyr::group_by(it_stratum) %>%
  dplyr::slice_min(it_dist_to_center, n = n_per_stratum, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(segment_id, it_stratum, it_dist_to_center)

it_candidates_seg <- it_seg_sf %>%
  dplyr::inner_join(it_candidates_ids, by = "segment_id") %>%
  dplyr::select(segment_id, it_stratum, it_dist_to_center, geometry)

it_candidates_pts <- sf::st_point_on_surface(it_candidates_seg) %>%
  dplyr::mutate(space = "it") %>%
  dplyr::select(segment_id, space, stratum = it_stratum, dist_to_center = it_dist_to_center, geometry)

## -------------------------------------------------------------------
## 5) Write outputs
## -------------------------------------------------------------------
dir.create(dirname(it_seg_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(it_pts_out), recursive = TRUE, showWarnings = FALSE)

sf::st_write(it_seg_sf, it_seg_out, delete_dsn = TRUE, quiet = TRUE)
sf::st_write(it_candidates_pts, it_pts_out, delete_dsn = TRUE, quiet = TRUE)
