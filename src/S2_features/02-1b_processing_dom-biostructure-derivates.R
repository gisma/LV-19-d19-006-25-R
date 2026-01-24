#!/usr/bin/env Rscript

############################################################
# Script:  02-2_processing_biostructure_chm-derivates.R
# Project: Burgwald
#
# Purpose:
#   Produce canonical S2 biostructure rasters required by:
#     04-4_signatures_biostructure_metrics.R
#
# Canonical outputs (via paths[] / outputs.tsv):
#   - chm_mean_10m         (tif)
#   - chm_p95_10m          (tif)
#   - chm_sd_10m           (tif)
#   - canopy_fraction_10m  (tif)
#
# Inputs (canonical S1):
#   - aoi_dgm  (DEM, 1 m)
#   - aoi_dom  (DOM/DSM, 1 m)
#
# Design rules:
#   - No extra productive artifacts / keys.
#   - All outputs are explicit file paths from paths[].
#   - Intermediate dCHM is computed in-memory (or tmp only).
############################################################

suppressPackageStartupMessages({
  library(here)
  library(terra)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))

# -------------------------------------------------------------------
# Productive inputs (canonical S1)
# -------------------------------------------------------------------
dem_file <- paths[["aoi_dgm"]]
dom_file <- paths[["aoi_dom"]]

stopifnot(file.exists(dem_file), file.exists(dom_file))

dem_1m <- terra::rast(dem_file)
dom_1m <- terra::rast(dom_file)

# -------------------------------------------------------------------
# Productive outputs (canonical S2_features / structure)
# -------------------------------------------------------------------
out_chm_mean_10m <- paths[["chm_mean_10m"]]
out_chm_p95_10m  <- paths[["chm_p95_10m"]]
out_chm_sd_10m   <- paths[["chm_sd_10m"]]
out_canfrac_10m  <- paths[["canopy_fraction_10m"]]
out_3dstruc_10m  <- paths[["biostructure_stack_10m"]]

dir.create(dirname(out_chm_mean_10m), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_chm_p95_10m),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_chm_sd_10m),   recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_canfrac_10m),  recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Parameters (explicit)
# -------------------------------------------------------------------
target_res_m <- 10
canopy_thr_m <- 2.0   # canopy pixel threshold on dCHM (meters)

# -------------------------------------------------------------------
# 1) Harmonise grids (DOM onto DEM grid if needed)
# -------------------------------------------------------------------
# hard rule: compute dCHM on a consistent 1 m grid.
# If grids differ, we resample DOM to DEM grid.
same_grid <- identical(terra::crs(dem_1m), terra::crs(dom_1m)) &&
  all.equal(terra::res(dem_1m), terra::res(dom_1m)) &&
  terra::ext(dem_1m) == terra::ext(dom_1m)

if (!same_grid) {
  dom_1m <- terra::resample(dom_1m, dem_1m, method = "bilinear")
}

# -------------------------------------------------------------------
# 2) dCHM (DOM - DEM), clamp negatives to 0 (explicit choice)
# -------------------------------------------------------------------
dchm_1m <- dom_1m - dem_1m
dchm_1m <- terra::clamp(dchm_1m, lower = 0, values = TRUE)

# -------------------------------------------------------------------
# 3) Aggregate 1 m -> 10 m as strict block stats (no segmentation)
# -------------------------------------------------------------------
res_m <- terra::res(dchm_1m)[1]
ratio <- target_res_m / res_m
fact  <- as.integer(round(ratio))

if (abs(ratio - fact) > 1e-6 || fact < 2) {
  stop("Resolution ratio not integer enough for strict block aggregation: res(dCHM) = ",
       res_m, " m; target_res_m = ", target_res_m)
}

# mean
chm_mean_10m <- terra::aggregate(dchm_1m, fact = fact, fun = mean, na.rm = TRUE)

# p95 (explicit quantile)
p95_fun <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = 0.95, na.rm = TRUE, names = FALSE, type = 7))
}
chm_p95_10m <- terra::aggregate(dchm_1m, fact = fact, fun = p95_fun)

# sd
sd_fun <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x)
}
chm_sd_10m <- terra::aggregate(dchm_1m, fact = fact, fun = sd_fun)

# canopy fraction: fraction of 1 m pixels above threshold within each 10 m block
can_mask <- dchm_1m > canopy_thr_m
canopy_fraction_10m <- terra::aggregate(can_mask, fact = fact, fun = mean, na.rm = TRUE)

# -------------------------------------------------------------------
# 4) Write canonical outputs
# -------------------------------------------------------------------
terra::writeRaster(chm_mean_10m,        out_chm_mean_10m, overwrite = TRUE)
terra::writeRaster(chm_p95_10m,         out_chm_p95_10m,  overwrite = TRUE)
terra::writeRaster(chm_sd_10m,          out_chm_sd_10m,   overwrite = TRUE)
terra::writeRaster(canopy_fraction_10m, out_canfrac_10m,  overwrite = TRUE)

r = c(chm_mean_10m,chm_p95_10m,chm_sd_10m, canopy_fraction_10m)
names(r) = c("chm_mean_10m","chm_p95_10m","chm_sd_10m", "canopy_frac_10m")
terra::writeRaster(r, out_3dstruc_10m,  overwrite = TRUE)

message("Wrote: ", out_chm_mean_10m)
message("Wrote: ", out_chm_p95_10m)
message("Wrote: ", out_chm_sd_10m)
message("Wrote: ", out_canfrac_10m)
message("Wrote: ", out_3dstruc_10m)
