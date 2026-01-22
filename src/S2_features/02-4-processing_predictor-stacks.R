#!/usr/bin/env Rscript

############################################################
# Script:   02-5-processing_predictor-stacks.R
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
# - Training points derived from *static valid mask* (CLC + predictors),
#   balanced by CLC class without per-class raster masks.
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
##     FAST balanced sampling via static valid mask + cell indices
## -------------------------------------------------------------------

out_trainpts_gpkg <- paths[["s2_lc_rf_trainpts_2021"]]

run_if_missing(out_trainpts_gpkg, {
  
  # predictors = RF-ready stack product
  pred <- terra::rast(out_file)
  stopifnot(inherits(pred, "SpatRaster"))
  
  # CLC label raster (S1 product)
  clc_file <- paths[["aoi_clc"]]
  stopifnot(file.exists(clc_file))
  
  msg("Loading CLC label raster: ", clc_file)
  clc <- terra::rast(clc_file)
  if (terra::nlyr(clc) != 1) stop("aoi_clc must be a single-layer categorical raster.")
  
  # align CLC strictly to predictor grid (NEAR)
  clc <- align_to_template(clc, pred, method = "near")
  
  if (!terra::compareGeom(pred, clc, stopOnError = FALSE,
                          crs=TRUE, ext=TRUE, rowcol=TRUE, res=TRUE)) {
    stop("CLC is not aligned to RF predictor stack (geom mismatch after alignment).")
  }
  
  ## ---- static valid mask: CLC not NA + all predictors finite (once) ----
  msg("Building static valid mask (CLC + predictors) ...")
  
  ok_pred <- terra::app(pred, \(...) as.integer(all(is.finite(c(...)))))
  ok      <- (ok_pred == 1) & !is.na(clc)
  
  cells_ok <- terra::which(ok, cells = TRUE)
  if (length(cells_ok) == 0) stop("Valid mask has zero cells (ok-mask empty).")
  
  ## ---- read class labels for ok-cells (vectorised) -----------------------
  clc_vals <- terra::values(clc, cells_ok, mat = FALSE)
  ok_finite <- is.finite(clc_vals)
  cells_ok  <- cells_ok[ok_finite]
  clc_vals  <- clc_vals[ok_finite]
  
  if (length(cells_ok) == 0) stop("After filtering non-finite CLC labels, no valid cells remain.")
  
  classes <- sort(unique(as.integer(clc_vals)))
  if (length(classes) == 0) stop("No classes found in valid mask.")
  
  msg("Valid classes in AOI mask: ", length(classes))
  
  ## ---- balanced sampling by class via indices (no raster masks) ----------
  set.seed(as.integer(seed))
  
  n_per_class <- 500L  # balanced target per class
  idx_by_class <- split(seq_along(clc_vals), clc_vals)
  
  sampled_cells <- integer(0)
  sampled_class <- integer(0)
  
  for (cls in names(idx_by_class)) {
    idx <- idx_by_class[[cls]]
    if (length(idx) == 0) next
    
    take <- min(as.integer(n_per_class), length(idx))
    pick <- sample(idx, size = take, replace = FALSE)
    
    sampled_cells <- c(sampled_cells, cells_ok[pick])
    sampled_class <- c(sampled_class, rep.int(as.integer(cls), take))
  }
  
  if (length(sampled_cells) == 0) {
    stop("Balanced sampling returned zero points (no class had any available cells).")
  }
  
  msg(sprintf("Sampled points: %d (target ~ %d)", length(sampled_cells), length(classes) * n_per_class))
  
  ## ---- build points from sampled cells -----------------------------------
  xy <- terra::xyFromCell(pred, sampled_cells)
  pts <- terra::vect(as.data.frame(xy), geom = c("x", "y"), crs = terra::crs(pred))
  pts$class <- sampled_class
  
  ## ---- extract predictors once (fast) ------------------------------------
  x <- terra::extract(pred, pts)
  x <- x[, -1, drop = FALSE]  # drop ID
  
  train_df <- data.frame(class = pts$class, x, check.names = FALSE)
  
  # drop rows with any NA/NaN/Inf predictors
  ok_row <- stats::complete.cases(train_df) & is.finite(train_df$class)
  train_df <- train_df[ok_row, , drop = FALSE]
  pts      <- pts[ok_row]
  
  if (nrow(train_df) == 0) stop("Training table empty after predictor NA filtering (unexpected).")
  
  # attach predictor columns to points
  for (nm in names(train_df)[-1]) pts[[nm]] <- train_df[[nm]]
  
  ## ---- write GPKG ---------------------------------------------------------
  pts_sf <- sf::st_as_sf(pts)
  
  dir.create(dirname(out_trainpts_gpkg), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(pts_sf, out_trainpts_gpkg, delete_dsn = TRUE, quiet = TRUE)
  
  msg(sprintf("Wrote training points: %s (n=%d)", out_trainpts_gpkg, nrow(pts_sf)))
  
  TRUE
})

msg("Done.")
