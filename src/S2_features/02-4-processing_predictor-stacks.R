#!/usr/bin/env Rscript

############################################################
# Script:   02-4-processing_predictor-stacks.R
# Project:  Burgwald
#
# Purpose
# -------
# Build a *RF-ready* canonical S2 predictor stack by:
#   (1) loading the existing Sentinel predictor stack (template grid),
#   (2) appending registered S2 feature rasters (bio / hydro / topo),
#   (3) applying robust variable reduction on continuous bands using caret:
#         - remove zero-variance predictors
#         - remove exact linear combinations
#         - remove highly correlated predictors (|r| > 0.95)
#   (4) writing the result to a NEW canonical registry product with suffix "_rf".
#
# Additional S2 product
# ---------------------
# Persist RF training points (CLC labels + RF predictor values) to:
#   paths[["s2_lc_rf_trainpts_2021"]] (GPKG)
#
# Contract
# --------
# - Baseline stack remains untouched.
# - RF stack is written to a separate registry key (..._rf).
# - IO exclusively via `paths[...]`.
# - Deterministic grid alignment.
# - Categorical layers are NOT reduced.
# - Training points are derived from CLC + RF stack only (no S3 semantics).
############################################################

suppressPackageStartupMessages({
  library(here)
  library(terra)
  library(caret)
  library(sf)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))

## -------------------------------------------------------------------
## 1) Canonical IO (registry truth)
## -------------------------------------------------------------------

# Baseline predictor stack (template grid)
in_file  <- paths[["s2_pred_stack_2021"]]

# RF-ready predictor stack (new product, must exist in outputs.tsv)
out_file <- paths[["s2_pred_stack_2021_rf"]]

stopifnot(file.exists(in_file))

## -------------------------------------------------------------------
## 2) Canonical feature inputs to append
## -------------------------------------------------------------------

extra_feature_keys <- c(
  # "chm_mean_10m",
  # "chm_p95_10m",
  # "chm_sd_10m",
  # "canopy_fraction_10m",
  # "watershed_id_10m",
  # "strahler_order_10m",
  # "flowacc_10m",
  # "dist_to_stream_10m",
  # "relief_stack_10m"
)

missing_keys <- setdiff(extra_feature_keys, names(paths))
if (length(missing_keys) > 0) {
  stop("Missing required paths[] keys: ", paste(missing_keys, collapse = ", "))
}

missing_files <- extra_feature_keys[!file.exists(paths[extra_feature_keys])]
if (length(missing_files) > 0) {
  stop("Missing required files for keys: ", paste(missing_files, collapse = ", "))
}

# Categorical / ordinal predictors (never reduced)
categorical_keys <- c("watershed_id_10m", "strahler_order_10m")

## -------------------------------------------------------------------
## 3) Reduction parameters
## -------------------------------------------------------------------

corr_cutoff <- 0.95
sample_n    <- 5000L
seed        <- 1L

quiet <- FALSE
msg <- function(...) if (!quiet) message(...)

## -------------------------------------------------------------------
## 4) Helpers
## -------------------------------------------------------------------

align_to_template <- function(r, tpl, method = c("bilinear", "near")) {
  method <- match.arg(method)
  
  if (!terra::same.crs(r, tpl)) {
    r <- terra::project(r, tpl)
  }
  
  if (!terra::compareGeom(
    r, tpl, stopOnError = FALSE,
    crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE
  )) {
    r <- terra::resample(r, tpl, method = method)
  }
  
  r
}

## -------------------------------------------------------------------
## 5) Build RF predictor stack
## -------------------------------------------------------------------

run_if_missing(out_file, {
  
  msg("Loading template predictor stack: ", in_file)
  tpl <- terra::rast(in_file)
  
  msg("Appending registered S2 feature rasters ...")
  
  extra_list <- lapply(extra_feature_keys, function(k) {
    r <- terra::rast(paths[[k]])
    
    method <- if (k %in% categorical_keys) "near" else "bilinear"
    r <- align_to_template(r, tpl, method = method)
    
    # deterministic band naming
    if (terra::nlyr(r) == 1) {
      names(r) <- k
    } else {
      names(r) <- paste0(k, "__", names(r))
    }
    
    r
  })
  
  extra   <- if (length(extra_list) > 0) do.call(c, extra_list) else NULL
  stk_all <- if (!is.null(extra)) c(tpl, extra) else tpl
  
  msg("Stack before reduction: ", terra::nlyr(stk_all), " layers")
  
  ## -----------------------------------------------------------------
  ## Split continuous vs categorical predictors
  ## -----------------------------------------------------------------
  
  all_names <- names(stk_all)
  
  is_cat <- rep(FALSE, length(all_names))
  for (k in categorical_keys) {
    is_cat <- is_cat | (all_names == k) |
      startsWith(all_names, paste0(k, "__"))
  }
  
  cont_names <- all_names[!is_cat]
  cat_names  <- all_names[ is_cat]
  
  cont <- stk_all[[cont_names]]
  cat  <- if (length(cat_names) > 0) stk_all[[cat_names]] else NULL
  
  ## -----------------------------------------------------------------
  ## Variable reduction on continuous predictors (caret)
  ## -----------------------------------------------------------------
  
  set.seed(as.integer(seed))
  
  msg("Sampling pixels for variable reduction (n = ", sample_n, ") ...")
  
  pts <- terra::spatSample(
    cont,
    size      = as.integer(sample_n),
    method    = "random",
    na.rm     = TRUE,
    as.points = TRUE
  )
  
  if (is.null(pts) || nrow(pts) == 0) {
    stop("Sampling produced zero points.")
  }
  
  smp <- terra::extract(cont, pts)
  smp <- smp[, -1, drop = FALSE]
  smp <- as.data.frame(smp)
  
  if (nrow(smp) == 0) {
    stop("Sampling / extraction produced zero rows (check NA coverage / alignment).")
  }
  
  ## ---- 1) Drop zero-variance predictors ---------------------------
  
  nzv <- caret::nearZeroVar(smp, saveMetrics = TRUE)
  drop_nzv <- rownames(nzv)[nzv$zeroVar]
  
  if (length(drop_nzv) > 0) {
    msg("Dropping zero-variance bands: ", paste(drop_nzv, collapse = ", "))
    smp  <- smp[, !names(smp) %in% drop_nzv, drop = FALSE]
    cont <- cont[[setdiff(names(cont), drop_nzv)]]
  }
  
  ## ---- 2) Drop exact linear combinations --------------------------
  
  if (ncol(smp) >= 2) {
    msg("Detecting linear combinations ...")
    lc <- caret::findLinearCombos(smp)
    
    if (!is.null(lc$remove) && length(lc$remove) > 0) {
      drop_lc <- colnames(smp)[lc$remove]
      msg("Dropping linear-combo bands: ", paste(drop_lc, collapse = ", "))
      smp  <- smp[, -lc$remove, drop = FALSE]
      cont <- cont[[setdiff(names(cont), drop_lc)]]
    }
  }
  
  ## ---- 3) Drop highly correlated predictors -----------------------
  
  if (ncol(smp) >= 2) {
    msg("Computing correlation matrix ...")
    cor_mat <- stats::cor(smp, use = "pairwise.complete.obs", method = "pearson")
    
    msg("Applying correlation cutoff |r| > ", corr_cutoff)
    drop_cor_idx <- caret::findCorrelation(
      cor_mat,
      cutoff = corr_cutoff,
      names  = FALSE,
      exact  = TRUE
    )
    
    if (length(drop_cor_idx) > 0) {
      drop_cor <- colnames(smp)[drop_cor_idx]
      msg("Dropping correlated bands: ", paste(drop_cor, collapse = ", "))
      cont <- cont[[setdiff(names(cont), drop_cor)]]
    } else {
      msg("No correlated bands removed.")
    }
  }
  
  msg("Continuous predictors after reduction: ", terra::nlyr(cont))
  
  ## -----------------------------------------------------------------
  ## Reassemble final stack
  ## -----------------------------------------------------------------
  
  stk_final <- if (!is.null(cat)) c(cont, cat) else cont
  msg("Final RF stack: ", terra::nlyr(stk_final), " layers")
  
  terra::writeRaster(stk_final, out_file, overwrite = TRUE)
  msg("Wrote RF predictor stack: ", out_file)
  
  TRUE
})

## -------------------------------------------------------------------
## 6) S2 product: RF training points (CLC labels + RF predictor values)
## -------------------------------------------------------------------

out_trainpts_gpkg <- paths[["s2_lc_rf_trainpts_2021"]]

run_if_missing(out_trainpts_gpkg, {
  
  source(here::here("src", "r-libs", "helper.R"))
  
  # predictors = RF-ready stack product
  pred <- terra::rast(out_file)
  stopifnot(inherits(pred, "SpatRaster"))
  
  # CLC label raster (S1 product)
  clc_file <- paths[["aoi_clc"]]
  stopifnot(file.exists(clc_file))
  
  msg("Loading CLC label raster: ", clc_file)
  clc <- terra::rast(clc_file)
  if (terra::nlyr(clc) != 1) stop("aoi_clc must be a single-layer categorical raster.")
  
  # strict alignment to predictor grid (NEAR)
  clc <- align_to_template(clc, pred, method = "near")
  
  if (!terra::compareGeom(pred, clc, stopOnError = FALSE,
                          crs=TRUE, ext=TRUE, rowcol=TRUE, res=TRUE)) {
    stop("CLC is not aligned to RF predictor stack (geom mismatch after alignment).")
  }
  
  msg("Sampling RF training points from CLC (stratified) ...")
  
  pts_sf <- sample_training_points_from_clc(
    clc                = clc,
    n_per_class        = 500L,
    aoi                = NULL,
    predictors         = pred,
    min_dist_m         = NULL,
    seed               = as.integer(seed),
    drop_na_predictors = TRUE
  )
  
  if (nrow(pts_sf) == 0) stop("Training point sampling returned 0 points.")
  
  dir.create(dirname(out_trainpts_gpkg), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(pts_sf, out_trainpts_gpkg, delete_dsn = TRUE, quiet = TRUE)
  
  msg(sprintf("Wrote training points: %s (n=%d)", out_trainpts_gpkg, nrow(pts_sf)))
  
  TRUE
})

msg("Done.")
