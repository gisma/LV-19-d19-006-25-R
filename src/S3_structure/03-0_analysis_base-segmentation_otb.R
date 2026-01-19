#!/usr/bin/env Rscript

############################################################
# Script:  03-0_analysis_base-segmentation_otb.R
# Project: Burgwald
#
# Purpose:
# --------
# Base segmentation of Sentinel predictor stack using OTB (link2GI):
#   - Z-score standardisation (BandMathX)
#   - PCA feature extraction (DimensionalityReduction)
#   - Multi-scale MeanShift segmentation + ARI_prev stability metric
#   - Export final segmentation:
#        (1) run MeanShift in raster-label mode
#        (2) compute label-size distribution
#        (3) choose min-pixel threshold from distribution
#        (4) sieve (<min_px) and reassign to dominant neighbors
#        (5) polygonize the cleaned label raster and write GPKG
#
# Rules implemented:
# ------------------
# - Productive outputs come ONLY from outputs.tsv via `paths[...]`.
# - tmp stays non-productive (not in outputs.tsv).
# - OTB is linked via a fixed local searchLocation (no system-wide search).
############################################################

## -------------------------------------------------------------------
## 0) Setup + libraries
## -------------------------------------------------------------------

library(here)
library(link2GI)

library(terra)
library(sf)
library(dplyr)
library(tibble)
library(readr)

source(here::here("src", "_core", "01-setup-burgwald.R"))
source(here::here("src", "r-libs", "metrics-fun.R"))

## -------------------------------------------------------------------
## 1) Input definition (productive input via outputs.tsv)
## -------------------------------------------------------------------

pred_stack_file <- paths[["s2_pred_stack_2021"]]
stopifnot(file.exists(pred_stack_file))

## -------------------------------------------------------------------
## 2) Productive outputs (from outputs.tsv via `paths`)
## -------------------------------------------------------------------

zscore_file  <- paths[["layer0_zscore_stack"]]
pca_file     <- paths[["layer0_pca_features"]]
metrics_file <- paths[["layer0_meanshift_metrics"]]
seg_file     <- paths[["layer0_segments"]]

## -------------------------------------------------------------------
## 2b) tmp output folder (explicitly NOT productive / not in outputs.tsv)
## -------------------------------------------------------------------

tmp_root <- here::here("data", "processed", "layer0_segments", "tmp_meanshift")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Link OTB
## -------------------------------------------------------------------

otb <- link2GI::linkOTB(searchLocation = "~/apps/otb911")

## -------------------------------------------------------------------
## 4) Z-score standardisation (BandMathX) — Self-describing API
## -------------------------------------------------------------------
#
# GOAL
# ----
# Convert a multi-band predictor raster into a standardized feature stack
# where each band is transformed to:
#
#   z = (x - mean(x)) / sqrt(var(x))
#
# This removes scale differences between bands and makes them comparable
# for downstream PCA and segmentation.
#
# IMPORTANT TECHNICAL NOTE
# ------------------------
# BandMathX expressions contain semicolons and parentheses; therefore the
# expression must be shell-quoted. This is handled inside
# run_otb_bandmathx_zscore() via build_bandmathx_zscore_expr().
#

run_otb_bandmathx_zscore(
  otb      = otb,
  in_file  = pred_stack_file,
  out_file = zscore_file,
  ram      = 256,
  quiet    = FALSE
)
stopifnot(file.exists(zscore_file))

## -------------------------------------------------------------------
## 5) PCA feature extraction (DimensionalityReduction) — Self-describing API
## -------------------------------------------------------------------

run_otb_pca(
  otb      = otb,
  in_file  = zscore_file,
  out_file = pca_file,
  ncomp    = 6L,
  ram      = 256,
  quiet    = FALSE
)
stopifnot(file.exists(pca_file))

## -------------------------------------------------------------------
## 6) Multi-scale MeanShift + ARI_prev evaluation
## -------------------------------------------------------------------
# NOTE ON SCALE SELECTION
# -----------------------
# The parameters (spatialr, ranger, minsize) define operational segmentation scales,
# not physical process scales. They control neighborhood size, feature similarity
# tolerance, and internal region merging within the MeanShift algorithm.
#
# There is an inherent circularity between:
#   (1) landscape scale (what we assume to be meaningful spatial units),
#   (2) process scale (at which physical/ecological processes operate), and
#   (3) spectrally derivable scale (what can be extracted from RS features + PCA).
#
# Segmentation scale influences the observed structure, which is later used to
# justify the segmentation scale itself. Therefore parameters are treated as
# exploratory control parameters and evaluated via stability (ARI_prev),
# not interpreted as physically meaningful scales.
#
# Final minimum segment size is enforced separately via raster sieve post-processing.

scales <- tibble::tibble(
  scale_id = c("s20m", "s30m", "s40m", "s60m", "s80m", "s120m"),
  spatialr = c(2,     3,     4,     6,     8,     12),
  ranger   = c(0.06,  0.08,  0.10,  0.12,  0.14,  0.16),
  minsize  = c(30,    40,    50,    80,    100,   150)
)


metrics_list <- vector("list", nrow(scales))

for (i in seq_len(nrow(scales))) {
  
  s <- scales[i, ]
  message("Processing scale: ", s$scale_id)
  
  metrics_list[[i]] <- compute_ari_prev(
    otblink     = otb,
    feat_raster = pca_file,
    out_dir_tmp = tmp_root,
    spatialr    = s$spatialr,
    ranger      = s$ranger,
    minsize     = s$minsize,
    perturb     = list(dr = 0.02, ds = 1L, dm = 20L, K = 8L),
    sample_n    = 2e5,
    seed        = 1L,
    tilesize    = 256L,
    ram         = 256,
    quiet       = TRUE
  ) %>%
    dplyr::mutate(scale_id = s$scale_id)
}

metrics_df <- bind_rows(metrics_list)

# Normalize naming to the column used downstream (avoid silent mismatch)
# compute_ari_prev() returns 'ari_prev' (lowercase); downstream expects 'ARI_prev'
if ("ari_prev" %in% names(metrics_df) && !"ARI_prev" %in% names(metrics_df)) {
  metrics_df <- dplyr::rename(metrics_df, ARI_prev = ari_prev)
}

# write metrics CSV to productive path (outputs.tsv)
readr::write_csv(metrics_df, metrics_file)

## -------------------------------------------------------------------
## 7) Select best scale
## -------------------------------------------------------------------

best_id <- metrics_df %>%
  arrange(desc(ARI_prev)) %>%
  slice(1) %>%
  pull(scale_id)

message("Best scale selected: ", best_id)

best_params <- scales %>%
  filter(scale_id == best_id)

best_spatialr <- best_params$spatialr[[1]]
best_ranger   <- best_params$ranger[[1]]
best_minsize  <- best_params$minsize[[1]]

## -------------------------------------------------------------------
## 8) Final segmentation: raster labels -> stats -> sieve -> polygonize
## -------------------------------------------------------------------
#
# Rationale
# ---------
# OTB's 'minsize' influences the internal merging step but does NOT guarantee
# that the final segmentation contains no tiny regions. We therefore enforce
# a hard postcondition:
#
#   - compute the label size distribution
#   - choose a min pixel threshold from that distribution
#   - absorb labels smaller than min_px into the dominant neighbor label
#   - polygonize the cleaned label raster and write the productive GPKG
#

# (8.1) Produce label raster in tmp (non-productive)
label_file <- file.path(tmp_root, sprintf("layer0_labels_%s.tif", best_id))

run_otb_meanshift_labels(
  otb              = otb,
  in_raster        = pca_file,
  out_label_raster = label_file,
  spatialr         = best_spatialr,
  ranger           = best_ranger,
  minsize          = best_minsize,
  tilesize         = 256L,
  ram              = 256,
  quiet            = TRUE
)
stopifnot(file.exists(label_file))

# (8.2) Compute label size distribution and choose min_px from distribution
st <- label_size_stats(label_file)

message("Label-size quantiles (pixels per label):")
print(st$quantile)
message("Shares of small labels (by count):")
print(st$shares)

# Explicit, audit-friendly decision rule (see metrics-fun.R)
min_px <- choose_min_px(
  st,
  onepx_threshold     = 0.05,  # trigger cleanup if >=5% of labels are ~1 pixel
  min_px_if_triggered = 5L,
  min_px_else         = 0L
)
message("Chosen min_px for sieve/reassign: ", min_px)

# (8.3) Sieve / reassign if requested
label_file_clean <- label_file

if (min_px > 0) {
  lab_clean <- sieve_labels(label_file, min_px = min_px, max_iter = 50L)
  
  label_file_clean <- file.path(tmp_root, sprintf("layer0_labels_%s_sieved_minpx%d.tif", best_id, min_px))
  terra::writeRaster(lab_clean, label_file_clean, overwrite = TRUE)
  
  stopifnot(file.exists(label_file_clean))
  message("Sieved label raster written to (tmp): ", label_file_clean)
} else {
  message("No sieve applied (min_px == 0).")
}

# (8.4) Polygonize cleaned label raster and write productive GPKG
#
# We polygonize AFTER sieve so that the productive segmentation output
# reflects the hard postcondition (no tiny regions).
#
# dissolve=TRUE ensures one polygon per label.
lab_r <- terra::rast(label_file_clean)
pol   <- terra::as.polygons(lab_r, dissolve = TRUE, values = TRUE, na.rm = TRUE)

# The polygon attribute holding the label value will have the raster layer name.
# Make it explicit and stable.
lab_name <- names(pol)[1]
names(pol)[names(pol) == lab_name] <- "segment_id"

# Convert to sf and ensure segment_id is integer
seg_sf <- sf::st_as_sf(pol)
seg_sf$segment_id <- as.integer(seg_sf$segment_id)

# Write as GPKG to productive path (from outputs.tsv)
# Overwrite if exists.
if (file.exists(seg_file)) {
  try(sf::st_delete_dsn(seg_file), silent = TRUE)
}

sf::st_write(seg_sf, seg_file, delete_dsn = TRUE, quiet = TRUE)

message("Final segmentation written to (productive):")
message(seg_file)
stopifnot(file.exists(seg_file))

