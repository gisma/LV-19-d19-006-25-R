#!/usr/bin/env Rscript
############################################################
# Script:  05-0_decisions_candidate_selection_layer0.R
# Project: Burgwald
#
# Purpose:
#   Decision space (S5): derive candidate selections from the
#   full segment signature stack (S4_signatures).
#
#   This script produces:
#     1) IT-strata clustering (pattern-based)
#     2) Physio-strata clustering (process-based)
#     3) Representative segment candidates per stratum (centroid-nearest)
#
#   - NO changes to segmentation geometry
#   - Uses S4 attrstack as the single truth source
#
# Inputs (productive, from outputs.tsv via paths[]):
#   - layer0_segments_attrstack_metrics   (S4_signatures, gpkg)
#
# Outputs (productive, from outputs.tsv via paths[]):
#   - layer0_it_strata_segments           (S5_decisions, gpkg)
#   - layer0_physio_strata_segments       (S5_decisions, gpkg)
#   - layer0_it_candidates_pts            (S5_decisions, gpkg)
#   - layer0_physio_candidates_pts        (S5_decisions, gpkg)
#
# Notes:
#   - The concrete S5 output keys MUST exist in outputs.tsv.
#   - Cluster count (k) and candidates per stratum (n) are explicit
#     parameters below.
############################################################

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stats)
  library(units)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))
source(here::here("src", "r-libs", "metrics-fun.R"))

## -------------------------------------------------------------------
## 0) Parameters (explicit decision knobs)
## -------------------------------------------------------------------

k_it      <- 6L   # number of IT strata
k_physio  <- 6L   # number of Physio strata
n_per_stratum <- 3L  # candidates per stratum (centroid-nearest)

set.seed(42)

## -------------------------------------------------------------------
## 1) Productive input (from outputs.tsv)
## -------------------------------------------------------------------

attr_file <- paths[["layer0_segments_attrstack_metrics"]]
stopifnot(file.exists(attr_file))

## -------------------------------------------------------------------
## 2) Productive outputs (from outputs.tsv)
## -------------------------------------------------------------------

it_seg_out     <- paths[["layer0_it_strata_segments"]]
physio_seg_out <- paths[["layer0_physio_strata_segments"]]
it_pts_out     <- paths[["layer0_it_candidates_pts"]]
physio_pts_out <- paths[["layer0_physio_candidates_pts"]]

## -------------------------------------------------------------------
## 2b) tmp output folder (NON-productive)
## -------------------------------------------------------------------

tmp_root <- here::here("data", "tmp", "layer0_decisions", "tmp_candidate_selection")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Read inputs
## -------------------------------------------------------------------

segs_sf <- sf::read_sf(attr_file)
stopifnot("segment_id" %in% names(segs_sf))

## -------------------------------------------------------------------
## 4) IT decision space: choose columns + cluster
## -------------------------------------------------------------------
# Expected columns from S4 IT metrics:
#   H_norm, U  (or H, U if you kept non-normalized)
#
# We use: H_norm + U if available, else fall back to H + U

it_cols <- c("H_norm", "U")
if (!all(it_cols %in% names(segs_sf))) {
  it_cols <- c("H", "U")
}
stopifnot(all(it_cols %in% names(segs_sf)))

it_df <- segs_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(segment_id, dplyr::all_of(it_cols)) %>%
  tidyr::drop_na()

it_x <- scale(as.matrix(it_df[, it_cols, drop = FALSE]))
it_km <- kmeans(it_x, centers = k_it, nstart = 50)

it_assign <- tibble::tibble(
  segment_id = it_df$segment_id,
  it_stratum = as.integer(it_km$cluster)
)

segs_it_sf <- segs_sf %>%
  dplyr::left_join(it_assign, by = "segment_id")

## -------------------------------------------------------------------
## 5) Physio decision space: choose columns + cluster
## -------------------------------------------------------------------
# Expected columns from S4 physio metrics:
#   elev_mean, slope_mean_deg, southness_mean, forest_fraction
#
# forest_fraction comes from cover metrics; if not present, stop.

physio_cols <- c("elev_mean", "slope_mean_deg", "southness_mean", "forest_fraction")
stopifnot(all(physio_cols %in% names(segs_sf)))

physio_df <- segs_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(segment_id, dplyr::all_of(physio_cols)) %>%
  tidyr::drop_na()

physio_x <- scale(as.matrix(physio_df[, physio_cols, drop = FALSE]))
physio_km <- kmeans(physio_x, centers = k_physio, nstart = 50)

physio_assign <- tibble::tibble(
  segment_id      = physio_df$segment_id,
  physio_stratum  = as.integer(physio_km$cluster)
)

segs_physio_sf <- segs_sf %>%
  dplyr::left_join(physio_assign, by = "segment_id")

## -------------------------------------------------------------------
## 6) Representative candidates per stratum (centroid-nearest)
## -------------------------------------------------------------------

# helper: select n closest to cluster center in scaled space
select_nearest <- function(df, x_scaled, cluster, centers, n) {
  idx <- which(cluster == centers)
  if (length(idx) == 0) return(integer(0))
  # distance to that cluster center
  c_id <- centers
  mu <- as.numeric(attr(cluster, "centers")[c_id, , drop = TRUE]) # not used
}

# IT candidates
it_dist <- matrix(NA_real_, nrow = nrow(it_x), ncol = 1)
for (i in seq_len(nrow(it_x))) {
  c_id <- it_km$cluster[i]
  it_dist[i, 1] <- sqrt(sum((it_x[i, ] - it_km$centers[c_id, ])^2))
}
it_df2 <- it_df %>%
  dplyr::mutate(
    it_stratum = as.integer(it_km$cluster),
    dist_to_center = as.numeric(it_dist[, 1])
  )

it_candidates_ids <- it_df2 %>%
  dplyr::group_by(it_stratum) %>%
  dplyr::slice_min(dist_to_center, n = n_per_stratum, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(segment_id, it_stratum, dist_to_center)

it_candidates_sf <- segs_it_sf %>%
  dplyr::inner_join(it_candidates_ids, by = "segment_id") %>%
  dplyr::mutate(candidate_type = "it") %>%
  dplyr::select(segment_id, candidate_type, it_stratum, dist_to_center, geometry)

it_candidates_pts <- sf::st_centroid(it_candidates_sf)

# Physio candidates
physio_dist <- matrix(NA_real_, nrow = nrow(physio_x), ncol = 1)
for (i in seq_len(nrow(physio_x))) {
  c_id <- physio_km$cluster[i]
  physio_dist[i, 1] <- sqrt(sum((physio_x[i, ] - physio_km$centers[c_id, ])^2))
}
physio_df2 <- physio_df %>%
  dplyr::mutate(
    physio_stratum = as.integer(physio_km$cluster),
    dist_to_center = as.numeric(physio_dist[, 1])
  )

physio_candidates_ids <- physio_df2 %>%
  dplyr::group_by(physio_stratum) %>%
  dplyr::slice_min(dist_to_center, n = n_per_stratum, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(segment_id, physio_stratum, dist_to_center)

physio_candidates_sf <- segs_physio_sf %>%
  dplyr::inner_join(physio_candidates_ids, by = "segment_id") %>%
  dplyr::mutate(candidate_type = "physio") %>%
  dplyr::select(segment_id, candidate_type, physio_stratum, dist_to_center, geometry)

physio_candidates_pts <- sf::st_centroid(physio_candidates_sf)

## -------------------------------------------------------------------
## 7) Write productive outputs
## -------------------------------------------------------------------

dir.create(dirname(it_seg_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(physio_seg_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(it_pts_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(physio_pts_out), recursive = TRUE, showWarnings = FALSE)

sf::st_write(segs_it_sf,     it_seg_out,     delete_dsn = TRUE, quiet = TRUE)
sf::st_write(segs_physio_sf, physio_seg_out, delete_dsn = TRUE, quiet = TRUE)
sf::st_write(it_candidates_pts,     it_pts_out,     delete_dsn = TRUE, quiet = TRUE)
sf::st_write(physio_candidates_pts, physio_pts_out, delete_dsn = TRUE, quiet = TRUE)
