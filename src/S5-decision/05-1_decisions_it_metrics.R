#!/usr/bin/env Rscript
############################################################
# Script:  05-1_decisions_it_metrics.R
# Project: Burgwald
#
# Purpose
# -------
# IT-only decision preview via ordination:
#   - restrict the effective search domain to a polygon (gauge-station hull)
#   - extract IT columns (H, K, I, U) from the segment attrstack
#   - compute a simple ordination (PCA on z-scored variables)
#   - derive a compact "candidate set" for downstream decision stages
#
# Conceptual role
# --------------
# This is not a classifier and not a clustering workflow.
# The IT space is semantically defined by its axes:
#   H, K, I, U  = information-theoretic / structural descriptors
#
# PCA is used strictly as a coordinate system:
#   - z-scoring makes H/K/I/U comparable
#   - PC1/PC2 provide a low-dimensional map of IT variability
#
# Candidate logic (two modes)
# ---------------------------
# Mode A (pc1_bins):
#   - discretise PC1 into q_pc1 bins (quantiles or ntile fallback)
#   - per bin select n_per_stratum segments closest to PCA origin
#
# Mode B (top_n):
#   - ignore bins; select global Top-N segments closest to PCA origin
#
# Why "closest to origin"?
# ------------------------
# PCA is computed on z-scored variables: the origin corresponds to the
# multivariate mean ("typical IT profile"). Selecting near-origin segments
# yields representative, non-extreme exemplars.
#
# IMPORTANT CHANGE (requested)
# ----------------------------
# Candidate search is constrained to a polygon derived from gauge stations
# (the hull). The "within" predicate is applied BEFORE PCA and BEFORE
# candidate selection. This ensures candidates are only sought inside.
#
# Outputs
# -------
# - out_seg: all segments with IT fields + PCA coordinates attached
#            (outside hull or incomplete IT -> NA)
# - out_pts: candidate points (one point per selected segment)
############################################################

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(dplyr)
  library(tidyr)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))

# -------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------

# Candidate selection mode:
#   "pc1_bins" -> q_pc1 bins * n_per_stratum candidates (up to)
#   "top_n"    -> global N candidates (no bins)
candidate_mode <- "pc1_bins"

# Mode A: stratified along PC1
q_pc1 <- 9L
n_per_stratum <- 10L

# Mode B: global Top-N (used only if candidate_mode == "top_n")
top_n <- 100L

# Point placement should be done in projected CRS (planar geometry).
crs_proj <- 25832

# -------------------------------------------------------------------
# IO
# -------------------------------------------------------------------
in_file  <- paths[["layer0_segments_attrstack_metrics"]]
out_seg  <- paths[["layer0_it_strata_segments"]]
out_pts  <- paths[["layer0_it_candidates_pts"]]

dir.create(dirname(out_seg), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_pts), recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Domain polygon (gauge-station hull)
# -------------------------------------------------------------------
# Use the key you actually have in your project. If it differs, replace here.
hull_file <- paths[["gauges_hull"]]

stopifnot(file.exists(in_file))
stopifnot(file.exists(hull_file))

# -------------------------------------------------------------------
# Load segments (truth source) + hull
# -------------------------------------------------------------------
segs_all <- sf::read_sf(in_file)
stopifnot(inherits(segs_all, "sf"))
stopifnot("segment_id" %in% names(segs_all))

# Your pipeline uses "geom" as geometry column name; keep it.
stopifnot("geom" %in% names(segs_all))

hull <- sf::read_sf(hull_file)
stopifnot(inherits(hull, "sf"))
hull <- sf::st_make_valid(hull)
hull <- sf::st_union(hull)   # ensure single polygon/multipolygon

# Ensure same CRS before spatial predicate.
if (sf::st_crs(segs_all) != sf::st_crs(hull)) {
  hull <- sf::st_transform(hull, sf::st_crs(segs_all))
}

# -------------------------------------------------------------------
# 1) EARLY domain restriction: ONLY segments fully inside hull
# -------------------------------------------------------------------
# This is the requested behavior: PCA + candidates are derived only inside.
segs_in <- segs_all[sf::st_within(segs_all, hull, sparse = FALSE), ]
message("Segments within hull: ", nrow(segs_in), " / ", nrow(segs_all))
stopifnot(nrow(segs_in) > 0)

# -------------------------------------------------------------------
# 2) IT-only extraction (inside hull)
# -------------------------------------------------------------------
it_cols <- c("H", "K", "I", "U")
stopifnot(all(it_cols %in% names(segs_all)))

# Keep only the minimum IT state required for ordination.
# drop_na() defines which segments can be assigned PCA coordinates.
it_seg <- segs_in %>%
  dplyr::select(segment_id, dplyr::all_of(it_cols), geom) %>%
  dplyr::mutate(segment_id = as.integer(segment_id)) %>%
  tidyr::drop_na(dplyr::all_of(it_cols))

stopifnot(nrow(it_seg) > 0)

# -------------------------------------------------------------------
# 3) Ordination (PCA on z-scored H/K/I/U)
# -------------------------------------------------------------------
X  <- as.matrix(sf::st_drop_geometry(it_seg)[, it_cols, drop = FALSE])
Xz <- scale(X)

# prcomp() assumes centered/scaled input if you already did it.
pca <- prcomp(Xz, center = FALSE, scale. = FALSE)

scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
names(scores) <- c("it_pc1", "it_pc2")

it_seg$it_pc1 <- scores$it_pc1
it_seg$it_pc2 <- scores$it_pc2

# Variance explained (useful later for interpretation and reporting).
attr(it_seg, "it_pca_var") <- summary(pca)$importance[2, 1:2]

# Distance to PCA origin (selection score for "typical" segments).
it_seg$it_dist_to_center <- sqrt(it_seg$it_pc1^2 + it_seg$it_pc2^2)

# -------------------------------------------------------------------
# 4) Candidate selection (two explicit modes)
# -------------------------------------------------------------------
# Candidate selection operates ONLY on it_seg (inside hull + complete IT).
if (candidate_mode == "pc1_bins") {
  
  # Quantile binning along PC1 yields roughly balanced occupancy.
  # If breaks collapse (ties), ntile() provides a robust rank-based fallback.
  pc1_breaks <- as.numeric(quantile(
    it_seg$it_pc1,
    probs = seq(0, 1, length.out = q_pc1 + 1),
    na.rm = TRUE
  ))
  pc1_breaks <- unique(pc1_breaks)
  
  if (length(pc1_breaks) < 3L) {
    it_seg$it_stratum <- dplyr::ntile(it_seg$it_pc1, q_pc1)
  } else {
    it_seg$it_stratum <- cut(
      it_seg$it_pc1,
      breaks = pc1_breaks,
      include.lowest = TRUE,
      labels = FALSE
    )
  }
  it_seg$it_stratum <- as.integer(it_seg$it_stratum)
  
  # For each PC1 bin: choose n_per_stratum closest to PCA origin.
  cand_ids <- it_seg %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(it_stratum) %>%
    dplyr::slice_min(it_dist_to_center, n = n_per_stratum, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(segment_id)
  
} else if (candidate_mode == "top_n") {
  
  # Global Top-N closest to origin (no stratification).
  it_seg$it_stratum <- NA_integer_
  
  cand_ids <- it_seg %>%
    sf::st_drop_geometry() %>%
    dplyr::slice_min(it_dist_to_center, n = top_n, with_ties = FALSE) %>%
    dplyr::select(segment_id)
  
} else {
  stop("Unknown candidate_mode: ", candidate_mode)
}

# Candidate polygons: filter it_seg via IDs only (prevents .x/.y suffixes).
cand_poly <- it_seg %>%
  dplyr::semi_join(cand_ids, by = "segment_id")

# -------------------------------------------------------------------
# 5) Candidate points (geometry-safe point placement)
# -------------------------------------------------------------------
# st_point_on_surface() is correct in projected CRS; transform back afterwards.
cand_pts <- cand_poly %>%
  sf::st_make_valid() %>%
  sf::st_transform(crs_proj) %>%
  sf::st_point_on_surface() %>%
  sf::st_transform(sf::st_crs(cand_poly)) %>%
  dplyr::mutate(space = "it_pca") %>%
  dplyr::transmute(
    segment_id,
    space,
    stratum = it_stratum,
    dist_to_center = it_dist_to_center,
    it_pc1,
    it_pc2,
    geom
  )

# -------------------------------------------------------------------
# 6) Attach PCA/stratum fields back to ALL segments for output
# -------------------------------------------------------------------
# Outside hull -> no ordination intended -> NA.
# Inside hull but incomplete IT -> NA (because drop_na removed them).
it_assign <- it_seg %>%
  sf::st_drop_geometry() %>%
  dplyr::select(segment_id, dplyr::all_of(it_cols), it_pc1, it_pc2, it_stratum, it_dist_to_center)

it_seg_out <- segs_all %>%
  dplyr::select(segment_id, dplyr::all_of(it_cols), geom) %>%
  dplyr::left_join(it_assign, by = c("segment_id", it_cols))

# -------------------------------------------------------------------
# 7) Write outputs
# -------------------------------------------------------------------
sf::st_write(it_seg_out, out_seg, driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)
sf::st_write(cand_pts,   out_pts, driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)

message("Wrote segments:   ", out_seg, " (rows=", nrow(it_seg_out), ")")
message("Wrote candidates: ", out_pts, " (rows=", nrow(cand_pts), ")")
