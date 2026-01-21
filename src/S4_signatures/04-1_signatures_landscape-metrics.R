#!/usr/bin/env Rscript

############################################################
# Script:   03-2_analysis_landscape-metrics-burgwald.R
# Project:  Burgwald
#
# Purpose:
# --------
# Segment-domain landscape metrics on the *stabilized MeanShift segments*.
#
# This is the S3-consistent replacement of the old grid-based (2×2 km) approach:
#   - OLD (deprecated): compute H/U on arbitrary grid cells (didactic only).
#   - NEW (S3): compute H/U *per segmentation unit* (your optimal segments).
#
# Outputs are written ONLY via outputs.tsv (paths[...] from 01-setup-burgwald.R).
#
# Required outputs.tsv keys (add them if missing):
#   S3_structure  layer0_segments_it_metrics   metrics   gpkg
# Optional (only if you want raster debug fields; otherwise ignore/remove blocks):
#   S3_structure  layer0_it_entropy_raster     metrics   tif
#   S3_structure  layer0_it_relmutinf_raster   metrics   tif
#
# Conceptual notes (short, but explicit):
# --------------------------------------
# - H(x) (Shannon entropy) is a *compositional* measure (class mixture within a segment).
# - U (relative mutual information, "relmutinf") is a *configurational* measure:
#   it uses adjacency (co-occurrence) of classes *within a segment*.
#
# Circularity warning (important for the Reader logic):
# -----------------------------------------------------
# Segmentation scale ↔ process/landscape scale ↔ observable spectral scale is circular:
# you pick a segmentation scale partly from data-driven stability (ARI_prev),
# but the "right" segmentation scale is also a hypothesis about relevant processes.
# Therefore H/U on segments is not "ground truth" — it's a structured diagnostic
# that depends on (and feeds back into) earlier modelling choices (S0–S3 recursion).
############################################################

## -------------------------------------------------------------------
## 0) Libraries + setup
## -------------------------------------------------------------------

library(here)
library(sf)
library(terra)
library(dplyr)
library(exactextractr)

# NOTE: For the configurational metric U we need fast aggregation over many pairs.
# data.table is used only for counting/aggregation (no workflow logic).
library(data.table)

source(here::here("src", "_core", "01-setup-burgwald.R"))

## -------------------------------------------------------------------
## 1) Inputs (canonical)
## -------------------------------------------------------------------

# (A) Segments: stabilized MeanShift output (S3)
seg_file <- paths[["layer0_segments"]]

# (B) Categorical landcover / classification raster for H/U
# Prefer a 10 m classification (RF/MLC/kmeans) if available; else use CLC.
# You decide by changing this one line.
lc_file  <- paths[["classification_rf"]]  # <- preferred (10 m, categorical)
# lc_file <- paths[["aoi_clc"]]           # <- fallback (CLC 100 m)

stopifnot(file.exists(seg_file))
stopifnot(file.exists(lc_file))

## -------------------------------------------------------------------
## 2) Output (canonical)
## -------------------------------------------------------------------

out_gpkg <- paths[["layer0_attr_it_metrics"]]

## -------------------------------------------------------------------
## 3) Load data
## -------------------------------------------------------------------

segs <- sf::read_sf(seg_file)

# Ensure there is a stable segment id field.
# Your segmentation export uses "segment_id" (from the label raster value).
stopifnot("segment_id" %in% names(segs))

# Ensure integer ID (important for rasterization and joins).
segs$segment_id <- as.integer(segs$segment_id)

# Landcover / classification raster (categorical)
lc <- terra::rast(lc_file)

crs_lc <- terra::crs(lc, proj = TRUE)
if (!identical(sf::st_crs(segs)$wkt, crs_lc)) segs <- sf::st_transform(segs, crs_lc)

## -------------------------------------------------------------------
## 4) Metric 1: H(x) per segment (composition) using exactextractr
## -------------------------------------------------------------------
# exactextractr returns per-polygon class coverage fast and with partial-pixel weights.
# For entropy we only need class proportions p_k within each segment.

# Helper: Shannon entropy from a named numeric vector of class weights
shannon_entropy <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0) return(NA_real_)
  p <- w / sum(w)
  -sum(p * log(p))
}

# exact_extract(..., fun = "freq") returns a data.frame per polygon:
# value = class value; count = weighted cell count (coverage-weighted).
# We compute H from coverage-weighted counts.
freq_list <- exactextractr::exact_extract(lc, segs, fun = "freq")

# Convert list-of-data.frames to entropy per segment_id
H_df <- rbindlist(
  lapply(seq_along(freq_list), function(i) {
    df <- freq_list[[i]]
    # exactextractr::freq output columns: "value", "count"
    # (count is coverage-weighted cell count or area-weighted depending on raster)
    if (is.null(df) || nrow(df) == 0) {
      return(data.table(segment_id = segs$segment_id[i], H = NA_real_, K = 0L))
    }
    w <- df$count
    H <- shannon_entropy(w)
    K <- sum(is.finite(w) & w > 0)
    data.table(segment_id = segs$segment_id[i], H = H, K = as.integer(K))
  }),
  use.names = TRUE, fill = TRUE
)

# Optional normalized entropy (0..1), defined as H / log(K) for K>=2
H_df[, H_norm := fifelse(is.finite(H) & K >= 2L, H / log(K), NA_real_)]

## -------------------------------------------------------------------
## 5) Metric 2: U per segment (configuration) from adjacency within segments
## -------------------------------------------------------------------
# Goal: relative mutual information U = I / H, where:
#   I = sum_{i,j} p(i,j) log( p(i,j) / (p(i)*p(j)) )
# computed from within-segment neighbor pairs (4-neighborhood).
#
# Implementation strategy (single raster scan; no per-segment loops):
#   1) Rasterize segment_id onto the landcover raster grid (same extent/resolution).
#   2) Build adjacency pairs (right + down neighbors) where segment_id matches.
#   3) Count pairs per segment and class combination.
#   4) Compute MI and U per segment.
#
# Note: This configurational computation is pixel-adjacency based.
# It does NOT use fractional pixel overlaps (unlike H via exactextractr).
# That is a deliberate compromise: exact fractional adjacency is not defined
# cleanly without a sub-pixel model. For 10 m categorical maps the effect is
# typically small compared to segmentation-scale uncertainty.

# Rasterize segments to the landcover grid
seg_v  <- terra::vect(segs)
seg_r  <- terra::rasterize(seg_v, lc, field = "segment_id", touches = TRUE)
names(seg_r) <- "segment_id"

# Extract matrices (fast adjacency construction)
lc_m  <- terra::as.matrix(lc, wide = TRUE)
seg_m <- terra::as.matrix(seg_r, wide = TRUE)

nr <- nrow(lc_m)
nc <- ncol(lc_m)

# Right-neighbor pairs (cell -> cell east)
a_seg <- seg_m[, 1:(nc - 1), drop = FALSE]
b_seg <- seg_m[, 2:nc,       drop = FALSE]
a_lc  <- lc_m[, 1:(nc - 1),  drop = FALSE]
b_lc  <- lc_m[, 2:nc,        drop = FALSE]

keep_r <- !is.na(a_seg) & !is.na(b_seg) & (a_seg == b_seg) &
  !is.na(a_lc) & !is.na(b_lc)

# Down-neighbor pairs (cell -> cell south)
c_seg <- seg_m[1:(nr - 1), , drop = FALSE]
d_seg <- seg_m[2:nr,       , drop = FALSE]
c_lc  <- lc_m[1:(nr - 1),  , drop = FALSE]
d_lc  <- lc_m[2:nr,        , drop = FALSE]

keep_d <- !is.na(c_seg) & !is.na(d_seg) & (c_seg == d_seg) &
  !is.na(c_lc) & !is.na(d_lc)

# Build long tables of adjacency pairs (segment_id, class_i, class_j)
pairs_dt <- rbind(
  data.table(
    segment_id = as.integer(a_seg[keep_r]),
    i = as.integer(a_lc[keep_r]),
    j = as.integer(b_lc[keep_r])
  ),
  data.table(
    segment_id = as.integer(c_seg[keep_d]),
    i = as.integer(c_lc[keep_d]),
    j = as.integer(d_lc[keep_d])
  )
)

# Count co-occurrences per segment and class pair
pairs_cnt <- pairs_dt[, .N, by = .(segment_id, i, j)]
setnames(pairs_cnt, "N", "n_ij")

# MI per segment from p(i,j), p(i), p(j)
# p(i,j) from pair counts; marginals computed from joint.
pairs_cnt[, n_seg := sum(n_ij), by = segment_id]
pairs_cnt[, p_ij := n_ij / n_seg]

pairs_cnt[, p_i := sum(p_ij), by = .(segment_id, i)]
pairs_cnt[, p_j := sum(p_ij), by = .(segment_id, j)]

pairs_cnt[, mi_term := p_ij * log(p_ij / (p_i * p_j))]

U_df <- pairs_cnt[, .(I = sum(mi_term, na.rm = TRUE)), by = segment_id]
# Join H to define U = I / H
U_df <- merge(U_df, H_df[, .(segment_id, H)], by = "segment_id", all.x = TRUE)
U_df[, U := fifelse(is.finite(I) & is.finite(H) & H > 0, I / H, 0)]
U_df <- U_df[, .(segment_id, I, U)]

## -------------------------------------------------------------------
## 6) Attach metrics to segments + write canonical output
## -------------------------------------------------------------------

# Add segment area (metric CRS if available; otherwise compute in native CRS)
# If segments are lon/lat, st_area returns spherical area; still OK for relative size.
segs$area_m2 <- as.numeric(sf::st_area(segs$geometry))

# Join metrics
segs <- segs %>%
  left_join(as.data.frame(H_df), by = "segment_id") %>%
  left_join(as.data.frame(U_df), by = "segment_id")

# Write GPKG (canonical)
if (file.exists(out_gpkg)) {
  try(sf::st_delete_dsn(out_gpkg), silent = TRUE)
}
sf::st_write(segs, out_gpkg, delete_dsn = TRUE, quiet = TRUE)

message("Segment-based IT metrics written to: ", out_gpkg)

############################################################
# Result columns in output:
#   segment_id  : integer label id (stable)
#   area_m2     : segment area
#   H           : Shannon entropy (composition) per segment
#   K           : number of classes with p>0 in segment
#   H_norm      : normalized entropy (H/log(K))
#   I           : mutual information from within-segment adjacency pairs
#   U           : relative mutual information U = I/H (configuration)
############################################################
