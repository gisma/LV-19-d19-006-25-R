#!/usr/bin/env Rscript

#######################################################################
# Script:   03-0_analysis_base-segmentation_otb.R
# Project:  Burgwald
#
# Purpose
# -------
# Base segmentation of a Sentinel predictor stack (S2-derived features) into
# spatial units (S3 structure) using OTB via link2GI:
#
#   (1) Z-score standardisation (OTB BandMathX) on the multi-band predictor stack
#   (2) PCA feature extraction (OTB DimensionalityReduction) -> low-dim feature raster
#   (3) Multi-scale MeanShift segmentation stability screening via ARI_prev:
#         - run MeanShift label raster for a base parameter set
#         - perturb parameters locally (spatialr/ranger/minsize) and re-run
#         - compute ARI between base vs. perturbed labelings on a sampled pixel set
#   (4) Final segmentation export (productive GPKG):
#         - run MeanShift in raster-label mode for the best parameter set
#         - compute label-size distribution (pixels per segment)
#         - choose a minimum size cutoff from the distribution (jump/knee)
#         - mask segments smaller than cutoff (set to NA)  [postcondition]
#         - polygonize cleaned label raster and write GPKG
#
# Project rules
# -------------
# - Productive inputs/outputs are resolved ONLY via outputs.tsv -> `paths[...]`
# - tmp outputs are non-productive and explicitly kept outside outputs.tsv
# - OTB is linked via a fixed local searchLocation (no system-wide search)
#
# Notes on concepts + key references
# ---------------------------------
# - Z-score standardisation: brings bands to comparable scale before PCA/segmentation.
#   (General statistical standardisation; used here as preprocessing convention.)
# - PCA: dimensionality reduction and de-correlation of features.
#   Jolliffe, I. (2002). Principal Component Analysis. Springer.
# - MeanShift: nonparametric mode-seeking / clustering; OTB implements a large-scale
#   mean-shift segmentation with region merging.
#   Comaniciu, D. & Meer, P. (2002). Mean shift: A robust approach toward feature
#   space analysis. IEEE TPAMI.
# - ARI: adjusted Rand index for comparing two partitions/labelings.
#   Hubert, L. & Arabie, P. (1985). Comparing partitions. Journal of Classification.
# - Stability screening: parameter choice is treated as an operational control problem
#   (exploratory), evaluated via perturbation stability rather than interpreted as
#   physical process scale (scale circulari#!/usr/bin/env Rscript
#######################################################################

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
## 1) Productive input (from outputs.tsv)
## -------------------------------------------------------------------

pred_stack_file <- paths[["s2_pred_stack_2021"]]
stopifnot(file.exists(pred_stack_file))

## -------------------------------------------------------------------
## 2) Productive outputs (from outputs.tsv)
## -------------------------------------------------------------------

zscore_file  <- paths[["layer0_zscore_stack"]]
pca_file     <- paths[["layer0_pca_features"]]
metrics_file <- paths[["layer0_meanshift_metrics"]]
seg_file     <- paths[["layer0_segments"]]

## -------------------------------------------------------------------
## 2b) tmp output folder (NON-productive)
## -------------------------------------------------------------------

tmp_root <- here::here("data", "tmp", "layer0_segments", "tmp_meanshift")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Link OTB (fixed local install)
## -------------------------------------------------------------------

otb <- link2GI::linkOTB(searchLocation = "~/apps/otb911")

## -------------------------------------------------------------------
## 4) Z-score standardisation (BandMathX) — Self-describing API
## -------------------------------------------------------------------
# Technical: BandMathX expressions contain ';' and '()' and MUST be shell-quoted.
# This is handled inside run_otb_bandmathx_zscore() (see metrics-fun.R).

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
# PCA is used here as a compact, de-correlated feature space for segmentation.
# (It is not a “process model”, only a feature transform.)

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
# Parameter meaning (operational):
# - spatialr: spatial bandwidth / neighborhood radius (in pixels) for mean shift
# - ranger:   range bandwidth in feature space (tolerance of similarity)
# - minsize:  minimum region size in the internal region merging step (OTB)
#
# Scale circularity (comment only):
# - These are NOT direct “physical scales”.
# - The segmentation scale influences what is “seen” as structure and would
#   later be used to justify the segmentation; therefore we screen by stability
#   (ARI under local perturbations), not by interpretive process claims.

# --- build scales: CHM-based if available, else fallback manual grid ----------
scales <- build_scales_A(
  master_10m_raster = pca_file,
  chm_1m = NULL  # später: Pfad zum CHM/DOM/DSM
)

# If build_scales_A returned fallback with ranger already filled, keep it.
# If ranger is NA (CHM path case), derive ranger from PCA feature space
if (any(!is.finite(scales$ranger))) {
  rvec <- derive_ranger_from_pca(pca_file, n = 50000L, k = 10L, seed = 1L)
  scales <- assign_ranger_to_scales(scales, ranger_candidates = rvec)
}

metrics_list <- vector("list", nrow(scales))

for (i in seq_len(nrow(scales))) {
  
  s <- scales[i, ]
  message("Processing scale: ", s$scale_id)
  
  # THIS is the adaptive perturbation used downstream
  pert <- make_adaptive_perturb(s$spatialr, s$ranger, s$minsize, K = 8L)
  
  res <- compute_ari_prev(
    otblink     = otb,
    feat_raster = pca_file,
    out_dir_tmp = tmp_root,
    spatialr    = s$spatialr,
    ranger      = s$ranger,
    minsize     = s$minsize,
    perturb     = pert,
    sample_fact = 4L,
    seed        = 1L,
    tilesize    = 1024L,
    ram         = 32768,
    quiet       = FALSE,
    verbose     = TRUE
  )
  
  metrics_list[[i]] <- tibble::tibble(
    scale_id = s$scale_id,
    spatialr = s$spatialr,
    ranger   = s$ranger,
    minsize  = s$minsize,
    ari_prev = res$ari_prev,
    ari_sd   = res$ari_sd
  )
}

metrics_df <- dplyr::bind_rows(metrics_list)
stopifnot(all(c("ari_prev", "ari_sd") %in% names(metrics_df)))

metrics_df <- dplyr::mutate(metrics_df, score = ari_prev - 0.5 * ari_sd)
metrics_df <- dplyr::arrange(metrics_df, dplyr::desc(score))

print(metrics_df, n = Inf)

best <- metrics_df[1, ]
best

readr::write_csv(metrics_df, metrics_file)


## -------------------------------------------------------------------
## 7) Select best scale (max ARI_prev)
## -------------------------------------------------------------------
best_id <- metrics_df %>%
  arrange(desc(score)) %>%
  slice(1) %>%
  pull(scale_id)

message("Best scale selected: ", best_id)

best_params <- scales %>% filter(scale_id == best_id)

best_spatialr <- best_params$spatialr[[1]]
best_ranger   <- best_params$ranger[[1]]
best_minsize  <- best_params$minsize[[1]]


#---------------------------------------------------------------
## 8) Final segmentation: labels -> size distribution -> cutoff -> mask -> polygonize
## -------------------------------------------------------------------

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

# (8.2) Label-size distribution (pixels per label) and jump-based cutoff
# Science intent: we target “representative homogeneous areas” and therefore
# remove the small-label regime below the detected size jump (knee). This is
# a pragmatic size-based postcondition, not a claim about physical process scale.

st <- label_size_stats(label_file)

message("Label-size quantiles (pixels per label):")
print(st$quantile)
message("Shares of small labels (by count):")
print(st$shares)

min_px <- choose_min_px_jump(
  st,
  min_px_floor = 1L,   # never below 1 px
  smooth_k     = 101L  # running median window (odd); stabilizes knee detection
)
message("Chosen min_px from jump-analysis (mask labels with n_px < min_px): ", min_px)

# (8.3) Mask segments smaller than cutoff (set to NA)
label_file_clean <- label_file

if (is.finite(min_px) && min_px > 1L) {
  
  lab_clean <- mask_small_segments(label_file, min_px = min_px)
  
  label_file_clean <- file.path(
    tmp_root,
    sprintf("layer0_labels_%s_masked_minpx%d.tif", best_id, min_px)
  )
  terra::writeRaster(lab_clean, label_file_clean, overwrite = TRUE)
  stopifnot(file.exists(label_file_clean))
  
  message("Masked label raster written to (tmp): ", label_file_clean)
  
} else {
  
  message("No masking applied (min_px <= 1).")
}

# (8.4) Polygonize cleaned labels and write productive GPKG
# Technical: polygonize AFTER masking, so the productive segmentation obeys the
# postcondition (no small segments from the low-size regime).
lab_r <- terra::rast(label_file_clean)

pol <- terra::as.polygons(lab_r, dissolve = TRUE, values = TRUE, na.rm = TRUE)

# The label attribute name depends on the raster layer name; make it stable.
lab_name <- names(pol)[1]
names(pol)[names(pol) == lab_name] <- "segment_id"

seg_sf <- sf::st_as_sf(pol)
seg_sf$segment_id <- as.integer(seg_sf$segment_id)

# Write to productive path from outputs.tsv (overwrite if exists).
if (file.exists(seg_file)) {
  try(sf::st_delete_dsn(seg_file), silent = TRUE)
}
sf::st_write(seg_sf, seg_file, delete_dsn = TRUE, quiet = TRUE)

message("Final segmentation written to (productive):")
message(seg_file)

