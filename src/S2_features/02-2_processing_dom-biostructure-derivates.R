#!/usr/bin/env Rscript

############################################################
# Script:  02-2_processing_dom-biostructure-derivates.R
# Project: Burgwald
#
# Purpose:
#   Produce canonical S2 biostructure rasters required by:
#     04-4_signatures_biostructure_metrics.R
#
# Concept (what "biostructure" means here)
# ----------------------------------------
# This script derives simple, robust canopy-structure summaries from a DEM (ground)
# and a DOM/DSM (surface). It computes a *dCHM* (surface height above ground) and
# aggregates it from 1 m to 10 m blocks.
#
# Output layers (10 m):
#   - chm_mean_10m        : mean canopy height (m) in each 10x10 m block
#   - chm_p95_10m         : 95th percentile of canopy height (m) per block
#   - chm_sd_10m          : within-block canopy height variability (sd, m)
#   - canopy_fraction_10m : fraction of 1 m pixels above canopy_thr_m per block
#   - biostructure_stack_10m : 4-layer stack (the layers above) for convenient downstream use
#
# Downstream implications (why these rasters matter)
# -------------------------------------------------
# Downstream scripts consume these rasters as *predictor features / signatures*
# on the same 10 m grid as other S2_features.
# Typical uses:
#   - CHM mean/p95: "how tall is vegetation in this cell?"
#   - CHM sd:       "how heterogeneous is the canopy structure?"
#   - canopy_frac:  "how much of the cell is 'canopy' vs open/low vegetation?"
#
# Important design rules (enforced by this script)
# ------------------------------------------------
# - No extra productive artifacts / keys: only paths[] outputs are written.
# - Only strict block aggregation is used (no smoothing, no segmentation).
# - dCHM negatives are clamped to 0 explicitly (interpretation: no negative canopy height).
# - Aggregation requires an integer ratio (1 m -> 10 m); otherwise stop() to prevent
#   silent misalignment or non-block aggregation.
#
# NOTE about the common "writeRaster no cell values" error
# -------------------------------------------------------
# If you see:
#   Error: [writeRaster] there are no cell values
# it typically means at least one raster has no values loaded in memory for writing
# (e.g. it is a virtual object, fully NA, or not properly computed due to upstream NA).
# In this script, each layer is explicitly computed from dCHM_1m and then written.
# If this error occurs, check whether dem_1m/dom_1m have valid cell values and overlap.
############################################################

suppressPackageStartupMessages({
  library(here)
  library(terra)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))

# -------------------------------------------------------------------
# Productive inputs (canonical S1)
# -------------------------------------------------------------------
# dem_file: ground elevation model at 1 m (DGM)
# dom_file: surface elevation model at 1 m (DOM/DSM)
dem_file <- paths[["aoi_dgm"]]
dom_file <- paths[["aoi_dom"]]

# hard requirement: both files must exist
stopifnot(file.exists(dem_file), file.exists(dom_file))

# load rasters (single-layer expected)
dem_1m <- terra::rast(dem_file)
dom_1m <- terra::rast(dom_file)

# -------------------------------------------------------------------
# Productive outputs (canonical S2_features / structure)
# -------------------------------------------------------------------
# Each output is a *single GeoTIFF* (10 m grid), plus one 4-layer stack (GeoTIFF).
# No additional outputs are created.
out_chm_mean_10m <- paths[["chm_mean_10m"]]
out_chm_p95_10m  <- paths[["chm_p95_10m"]]
out_chm_sd_10m   <- paths[["chm_sd_10m"]]
out_canfrac_10m  <- paths[["canopy_fraction_10m"]]
out_3dstruc_10m  <- paths[["biostructure_stack_10m"]]

# ensure output directories exist
dir.create(dirname(out_chm_mean_10m), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_chm_p95_10m),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_chm_sd_10m),   recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_canfrac_10m),  recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Parameters (explicit)
# -------------------------------------------------------------------
# target_res_m:
#   aggregation target grid size (meters). This script assumes a 1 m base grid
#   and performs strict factor aggregation to 10 m blocks.
#
# canopy_thr_m:
#   threshold applied on the 1 m dCHM to define "canopy presence".
#   canopy_fraction_10m is the fraction of 1 m pixels with dCHM > canopy_thr_m.
target_res_m <- 10
canopy_thr_m <- 2.0   # canopy pixel threshold on dCHM (meters)

# -------------------------------------------------------------------
# 1) Harmonise grids (DOM onto DEM grid if needed)
# -------------------------------------------------------------------
# hard rule: compute dCHM on a consistent 1 m grid.
# If grids differ, we resample DOM to DEM grid.
#
# Rationale:
#   dCHM = DOM - DEM is only meaningful if both rasters align cell-by-cell.
#   If they differ in CRS/resolution/extent, subtraction would be wrong.
same_grid <- identical(terra::crs(dem_1m), terra::crs(dom_1m)) &&
  all.equal(terra::res(dem_1m), terra::res(dom_1m)) &&
  terra::ext(dem_1m) == terra::ext(dom_1m)

if (!same_grid) {
  # Resample DOM to DEM grid (bilinear: continuous elevation surface)
  dom_1m <- terra::resample(dom_1m, dem_1m, method = "bilinear")
}

# -------------------------------------------------------------------
# 2) dCHM (DOM - DEM), clamp negatives to 0 (explicit choice)
# -------------------------------------------------------------------
# dchm_1m represents height above ground (meters).
# Negative values can occur due to:
#   - interpolation/resampling artifacts,
#   - slight vertical biases between DEM and DOM,
#   - boundary effects,
#   - data gaps.
#
# Explicit policy: clamp negative values to 0 (no negative canopy height).
dchm_1m <- dom_1m - dem_1m
dchm_1m <- terra::clamp(dchm_1m, lower = 0, values = TRUE)

# -------------------------------------------------------------------
# 3) Aggregate 1 m -> 10 m as strict block stats (no segmentation)
# -------------------------------------------------------------------
# Strict block aggregation requires an integer factor:
#   fact = target_res_m / res(dchm_1m)
# For 1 m base and target_res_m=10, fact=10 (10x10 = 100 cells per block).
res_m <- terra::res(dchm_1m)[1]
ratio <- target_res_m / res_m
fact  <- as.integer(round(ratio))

# safeguard: do NOT allow non-integer ratios, otherwise aggregation is not strict blocks
if (abs(ratio - fact) > 1e-6 || fact < 2) {
  stop("Resolution ratio not integer enough for strict block aggregation: res(dCHM) = ",
       res_m, " m; target_res_m = ", target_res_m)
}

# mean canopy height per 10 m cell (na.rm=TRUE handles nodata within blocks)
chm_mean_10m <- terra::aggregate(dchm_1m, fact = fact, fun = mean, na.rm = TRUE)

# p95 (explicit quantile):
# - uses finite values only
# - returns NA if the block has no finite values
p95_fun <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = 0.95, na.rm = TRUE, names = FALSE, type = 7))
}
chm_p95_10m <- terra::aggregate(dchm_1m, fact = fact, fun = p95_fun)

# sd within each 10 m block:
# - uses finite values only
# - returns NA if <2 finite values (sd undefined)
sd_fun <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x)
}
chm_sd_10m <- terra::aggregate(dchm_1m, fact = fact, fun = sd_fun)

# canopy fraction:
# - compute 1 m boolean mask: canopy present if dCHM > canopy_thr_m
# - aggregate by mean => fraction in [0..1]
can_mask <- dchm_1m > canopy_thr_m
canopy_fraction_10m <- terra::aggregate(can_mask, fact = fact, fun = mean, na.rm = TRUE)

# -------------------------------------------------------------------
# 4) Write canonical outputs
# -------------------------------------------------------------------
# Each layer is written as a single-band raster (GeoTIFF).
# The combined 4-layer biostructure stack is written as a multi-band GeoTIFF.
terra::writeRaster(chm_mean_10m,        out_chm_mean_10m, overwrite = TRUE)
terra::writeRaster(chm_p95_10m,         out_chm_p95_10m,  overwrite = TRUE)
terra::writeRaster(chm_sd_10m,          out_chm_sd_10m,   overwrite = TRUE)
terra::writeRaster(canopy_fraction_10m, out_canfrac_10m,  overwrite = TRUE)

# Build a 4-layer stack for downstream convenience:
# Order is explicit and names are set for stable downstream layer referencing.
r = c(chm_mean_10m,chm_p95_10m,chm_sd_10m, canopy_fraction_10m)
names(r) = c("chm_mean_10m","chm_p95_10m","chm_sd_10m", "canopy_frac_10m")
terra::writeRaster(r, out_3dstruc_10m,  overwrite = TRUE)

message("Wrote: ", out_chm_mean_10m)
message("Wrote: ", out_chm_p95_10m)
message("Wrote: ", out_chm_sd_10m)
message("Wrote: ", out_canfrac_10m)
message("Wrote: ", out_3dstruc_10m)
