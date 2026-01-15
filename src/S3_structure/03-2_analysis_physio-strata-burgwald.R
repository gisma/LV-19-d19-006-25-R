#!/usr/bin/env Rscript

############################################################
# 03-0_layer0_segmentation_otb.R
#
# Layer 0:
#   - Build a collinearity-robust feature image from CDSE predictor stacks
#   - Run OTB LargeScaleMeanShift in multiple scales
#   - Output: segmentation polygons (GPKG) per scale
#
# Requires:
#   - 00_setup-burgwald.R has been sourced (here::here(), tmp dirs, etc.)
#   - data/pred_stack_<year>.tif exists (CDSE workflow)
#
# Notes:
#   - OTB app: LargeScaleMeanShift (otbcli_LargeScaleMeanShift)
#     Parameters: -in, -spatialr, -ranger, -minsize, -mode vector, -mode.vector.out
#     (OTB Cookbook) https://www.orfeo-toolbox.org/CookBook/Applications/app_LargeScaleMeanShift.html
############################################################

library(here)
library(fs)
library(terra)

message("Project root: ", here::here())

# ---------------------------------------------------------
# 0) Paths
# ---------------------------------------------------------
year <- 2018
pred_stack_file <- here::here("data", paste0("pred_stack_", year, ".tif"))

out_root <- here::here("data", "processed", "layer0_segments")
fs::dir_create(out_root)

feat_dir <- file.path(out_root, "features")
seg_dir  <- file.path(out_root, "segments")
fs::dir_create(feat_dir)
fs::dir_create(seg_dir)

stopifnot(file.exists(pred_stack_file))

# ---------------------------------------------------------
# 1) Load predictor stack
# ---------------------------------------------------------
x <- terra::rast(pred_stack_file)

# Optional: keep only a sane subset first (if your stacks contain lots of near-duplicates)
# Adjust to your actual band names. Keep a "mixed" set: RED/NIR + 2-4 indices max.
# If names don't match, this block is skipped and we use all layers.
keep_try <- c("B02","B03","B04","B08","B11","B12","kNDVI","EVI","SAVI","NDVI")
if (all(keep_try %in% names(x))) {
  x <- x[[keep_try]]
}

message("Layers used: ", paste(names(x), collapse = ", "))

# ---------------------------------------------------------
# 2) Build a collinearity-robust feature image
#    Choose ONE of two strategies:
#      A) correlation pruning (keep real bands/indices, drop redundant)
#      B) PCA compression (write first N components)
#
# Default here: PCA (most robust against collinearity)
# ---------------------------------------------------------

# ---- helper: sample pixels safely
sample_pixels <- function(r, n = 50000, seed = 123) {
  set.seed(seed)
  # spatSample returns data.frame (values only) when as.df=TRUE
  d <- terra::spatSample(r, size = n, method = "random", na.rm = TRUE, as.df = TRUE)
  # drop rows with any NA (should already be na.rm=TRUE, but keep safe)
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d
}

# ---- strategy B: PCA compression
make_pca_features <- function(r, ncomp = 6, n_samp = 50000) {
  d <- sample_pixels(r, n = n_samp)
  
  # scale features -> avoids band range imbalance (OTB warns about this)
  d_scaled <- scale(d)
  
  # PCA in feature space
  p <- stats::prcomp(d_scaled, center = FALSE, scale. = FALSE)
  
  ncomp <- min(ncomp, ncol(p$rotation))
  rot <- p$rotation[, 1:ncomp, drop = FALSE]
  
  # apply PCA to raster: for each cell, (x - mu)/sd %*% rot
  # get mu/sd from the sampled scaling attributes:
  mu <- attr(d_scaled, "scaled:center")
  sd <- attr(d_scaled, "scaled:scale")
  
  # terra::app on multi-layer rasters is streamed (good for big data)
  f_pca <- function(v) {
    # v: matrix [ncell x nbands]
    # scale:
    v <- sweep(v, 2, mu, "-")
    v <- sweep(v, 2, sd, "/")
    # project:
    v %*% rot
  }
  
  pcs <- terra::app(r, fun = f_pca)
  names(pcs) <- paste0("PC", seq_len(ncomp))
  pcs
}

# Build PCA features (robust default)
feat <- make_pca_features(x, ncomp = 6, n_samp = 50000)

feat_file <- file.path(feat_dir, paste0("features_pca_y", year, "_pc6.tif"))
terra::writeRaster(feat, feat_file, overwrite = TRUE)
message("Feature image written: ", feat_file)

# ---------------------------------------------------------
# 3) Run OTB LargeScaleMeanShift multiskale
# ---------------------------------------------------------

# Detect OTB CLI binary
otb_bin <- Sys.which("otbcli_LargeScaleMeanShift")
if (otb_bin == "") {
  stop("otbcli_LargeScaleMeanShift not found in PATH. Load OTB environment / module first.")
}

# If you want OTB to use more RAM:
# (OTB parameter is also -ram, but env var sometimes helps wrappers)
Sys.setenv(OTB_MAX_RAM_HINT = "4096")

run_lsms <- function(infile, spatialr, ranger, minsize, out_vector, ram = 4096,
                     tilesizex = 500, tilesizey = 500) {
  
  args <- c(
    "-in", infile,
    "-spatialr", as.integer(spatialr),
    "-ranger",   as.numeric(ranger),
    "-minsize",  as.integer(minsize),
    "-tilesizex", as.integer(tilesizex),
    "-tilesizey", as.integer(tilesizey),
    "-mode", "vector",
    "-mode.vector.out", out_vector,
    "-ram", as.integer(ram)
  )
  
  message("OTB: ", paste(c("otbcli_LargeScaleMeanShift", args), collapse = " "))
  res <- system2(otb_bin, args = args, stdout = TRUE, stderr = TRUE)
  
  # If needed, print tail of logs for debugging:
  if (length(res)) message(paste(tail(res, 20), collapse = "\n"))
  
  invisible(out_vector)
}

# IMPORTANT:
# - After PCA, features are roughly standardized; ranger can be small.
# - If you instead feed raw reflectances, ranger must be much larger.
#
# Spatial radius is in pixels. At 10 m resolution:
#   spatialr=3  -> ~30 m
#   spatialr=6  -> ~60 m
#   spatialr=12 -> ~120 m
#
# minsize is in pixels -> controls smallest object kept.
#   minsize=50 at 10 m -> 50 pixels ~ 5000 m² (not exact area, but order)
#
# Define a multiscale parameter set:
param_grid <- data.frame(
  scale_id  = c("s30m", "s60m", "s120m"),
  spatialr  = c(3, 6, 12),
  ranger    = c(0.10, 0.12, 0.15),
  minsize   = c(80, 120, 200),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(param_grid))) {
  p <- param_grid[i, ]
  
  out_vec <- file.path(
    seg_dir,
    paste0("segments_meanshift_", p$scale_id, "_y", year, ".gpkg")
  )
  
  run_lsms(
    infile     = feat_file,
    spatialr   = p$spatialr,
    ranger     = p$ranger,
    minsize    = p$minsize,
    out_vector = out_vec,
    ram        = 4096
  )
  
  message("Segmentation written: ", out_vec)
}

message("Layer 0 done.")
