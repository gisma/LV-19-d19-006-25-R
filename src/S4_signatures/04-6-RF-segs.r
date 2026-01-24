#!/usr/bin/env Rscript

# =============================================================================
# Script:   03-1-RF-segs.R
# Project:  Burgwald Decision Stack
#
# GOAL
# ----
# Segment-based landcover classification using the preserved 040-style workflow:
#   - Segment labels transferred from point labels (hard filter: n_points >= 3, purity >= 0.8)
#   - Segment signatures from predictor stack via polygon extraction (segment signatures; multiple stats)
#   - sperrorest::partition_kmeans() on segment centroids -> seg_df$partition
#   - CAST::CreateSpacetimeFolds(spacevar="partition", k=5)
#   - Step-1: preliminary hyperparameter tuning (random search; tuneLength=5; ntree=200)
#   - Step-2: CAST::ffs with mtry fixed to Step-1 bestTune; ntree=200
#   - Step-3: final hyperparameter sweep (ntree in {200,400,600} × mtry in 2..8) + resamples comparison
#   - Step-4: final refit with best hyperparams
#   - Prediction: segment-wise probs + class label + class_id (and optional AOA/DI/LPD in segment feature space)
#
# NOTES
# -----
# - This is segment-based. Training units are segments, not points/pixels.
# - Outputs are written with a "_segmentbased" suffix derived from the canonical 03-3 paths.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(terra)
  library(dplyr)
  library(tibble)
  library(readr)
  library(caret)
  library(CAST)
  library(sperrorest)
  library(randomForest)
  library(doParallel)
  library(parallel)
  library(exactextractr)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))

quiet <- FALSE
msg <- function(...) if (!quiet) message(...)

# -----------------------------------------------------------------------------
# 0) Canonical IO (paths[] only) + derived "_segmentbased" outputs
# -----------------------------------------------------------------------------

pred_stack_file  <- paths[["s2_pred_stack_2021_rf"]]
seg_file         <- paths[["layer0_segments"]]
train_seg_rds <-   paths[["train_seg_rds"]]

# canonical outputs from 03-3
out_class_file0 <- paths[["s4_lc_rf_2021"]]
out_prob_file0  <- paths[["s4_lc_rf_prob_2021"]]
out_model_rds0  <- paths[["s4_lc_rf_model_2021"]]
out_cv_csv0     <- paths[["s4_lc_rf_cv_metrics_2021"]]

add_suffix <- function(path, suffix = "_segmentbased") {
  sub("(\\.[A-Za-z0-9]+)$", paste0(suffix, "\\1"), path)
}


out_class_file <- add_suffix(out_class_file0)
out_prob_file  <- add_suffix(out_prob_file0)
out_model_rds  <- add_suffix(out_model_rds0)
out_cv_csv     <- add_suffix(out_cv_csv0)

dir.create(dirname(out_class_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_prob_file),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_model_rds),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_cv_csv),     recursive = TRUE, showWarnings = FALSE)

# vector output for segment predictions
out_seg_gpkg <- add_suffix(sub("\\.tif$", ".gpkg", out_class_file0))

# -----------------------------------------------------------------------------
# 1) Parameters (keep 040 intent)
# -----------------------------------------------------------------------------

seed <- 1
set.seed(seed)

k_folds        <- 5L
ffs_ntree      <- 200L
random_ntree   <- 200L
random_tuneLen <- 5L

ntree_grid <- c(200L, 400L, 600L)
mtry_max   <- 8L

# prediction parallelism
n_cores_predict <- 10L
cor_cutoff <- 0.95

# lon/lat handling for partitioning: "reproject" (default) or "allow"
partition_longlat_mode <- "reproject"

# Keep the 040 subsample behaviour (50% by default)
subsample_p <- 0.1

# -----------------------------------------------------------------------------
# 2) Load predictor stack
# -----------------------------------------------------------------------------

msg("Loading predictor stack: ", pred_stack_file)
pred_stack <- terra::rast(pred_stack_file)
stopifnot(terra::nlyr(pred_stack) >= 2)

# -----------------------------------------------------------------------------
# 3) Load segments + training points
# -----------------------------------------------------------------------------

msg("Loading segments: ", seg_file)
segments <- sf::st_read(seg_file, quiet = TRUE)
stopifnot("segment_id" %in% names(segments))
if (any(!sf::st_is_valid(segments))) segments <- sf::st_make_valid(segments)

# -----------------------------------------------------------------------------
# 3) Training table from upstream (CLC-supervised segment labeling)
# -----------------------------------------------------------------------------

# Upstream script 02-6 writes a segment-level training table (RDS) that already
# contains:
#   segment_id, class, class_id, purity (+ predictor signatures)
# This 040 script keeps the original CV/tuning/ffs/training/prediction workflow,
# but it does NOT redo any label-transfer / class_id mapping / CORINE handling.

seg_df_up <- readRDS(train_seg_rds)

req_up_cols <- c("segment_id", "class")
missing_up <- setdiff(req_up_cols, names(seg_df_up))
if (length(missing_up) > 0) stop("Upstream training table missing columns: ", paste(missing_up, collapse=", "))

# -----------------------------------------------------------------------------
# 5) Segment signatures from predictor stack (segment signatures; multiple stats)
# -----------------------------------------------------------------------------

msg("Extracting segment signatures from predictor stack (mean/stdev/min/max/cv/quantiles) ...")

# Ensure CRS alignment between segments and raster
segments_m <- sf::st_transform(segments, terra::crs(pred_stack))

ext_stats <- exactextractr::exact_extract(
  pred_stack,
  segments_m,
  fun = c("mean", "stdev", "min", "max", "coefficient_of_variation", "quantile"),
  quantiles = c(0.25, 0.50, 0.75)
)

stopifnot(nrow(ext_stats) == nrow(segments_m))

# normalise names: stat.band -> stat_band
names(ext_stats) <- gsub("\\.", "_", names(ext_stats))

seg_feat <- tibble::as_tibble(ext_stats)
seg_feat$segment_id <- segments_m$segment_id

# feature contract for RF
pred_names <- setdiff(names(seg_feat), "segment_id")
stopifnot(length(pred_names) >= 2)

# -----------------------------------------------------------------------------
# 5b) Build segment training table (UPSTREAM labels + predictor signatures)
# ----------------------------------train_seg_rds -------------------------------------------
# IMPORTANT:
# - seg_df_up already contains the target label `class` (and usually `class_id`,
#   `purity`, ...). We do NOT recompute any label transfer here.
# - We *attach* predictor signatures from the canonical predictor stack to ensure
#   prediction uses the exact same feature source as in the 040 workflow.

seg_df <- seg_df_up %>%
  dplyr::select(dplyr::any_of(c("segment_id", "class", "class_id", "purity", "n_classes", "n_points"))) %>%
  dplyr::left_join(seg_feat, by = "segment_id") %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(pred_names), ~ is.finite(.x)))


if (nrow(seg_df) == 0) stop("Segment training table empty after finite filtering of predictors.")

seg_df$class <- droplevels(as.factor(seg_df$class))

msg(sprintf("Segment training table: n=%d  predictors=%d  class_levels=%d",
            nrow(seg_df), length(pred_names), nlevels(seg_df$class)))

# -----------------------------------------------------------------------------
# 5c) Caret standard predictor pruning (NZV / linear combos / correlation)
# -----------------------------------------------------------------------------

# NOTE: This is applied on the segment-level training matrix (not on rasters).
# The same reduced predictor set is then used throughout tuning/FFS/training AND
# for prediction on all segments.

X0 <- seg_df[, pred_names, drop = FALSE]

# 1) Near-zero variance predictors
nzv_idx <- caret::nearZeroVar(X0)
nzv_drop <- if (length(nzv_idx) > 0) pred_names[nzv_idx] else character(0)

# 2) Linear combinations
X1 <- X0[, setdiff(pred_names, nzv_drop), drop = FALSE]
lc <- caret::findLinearCombos(X1)
lc_drop <- if (!is.null(lc$remove) && length(lc$remove) > 0) colnames(X1)[lc$remove] else character(0)

# 3) High correlation pruning
X2 <- X1[, setdiff(colnames(X1), lc_drop), drop = FALSE]
corr_drop <- character(0)
if (ncol(X2) >= 2) {
  cm <- suppressWarnings(stats::cor(X2, use = "pairwise.complete.obs"))
  corr_idx <- caret::findCorrelation(cm, cutoff = cor_cutoff, names = FALSE, verbose = FALSE)
  if (length(corr_idx) > 0) corr_drop <- colnames(X2)[corr_idx]
}

drop_vars <- unique(c(nzv_drop, lc_drop, corr_drop))

if (length(drop_vars) > 0) {
  msg(sprintf(
    "Pruning predictors: drop=%d (NZV=%d, LC=%d, CORR=%d @ %.2f)",
    length(drop_vars), length(nzv_drop), length(lc_drop), length(corr_drop), cor_cutoff
  ))
  pred_names <- setdiff(pred_names, drop_vars)
  # keep seg_feat/seg_df consistent for downstream selects
  seg_feat <- seg_feat %>% dplyr::select(-dplyr::any_of(drop_vars))
  seg_df   <- seg_df   %>% dplyr::select(-dplyr::any_of(drop_vars))
} else {
  msg("Pruning predictors: drop=0 (no NZV/LC/CORR removals)")
}

stopifnot(length(pred_names) >= 2)

# -----------------------------------------------------------------------------
# 6) Spatial partitions (sperrorest) on segment centroids
# -----------------------------------------------------------------------------

# Centroids derived from the segment geometry (no external objects).
# Use the same CRS as segments_m (aligned to predictor stack) for consistency.
seg_cent_sf <- sf::st_centroid(segments_m)

# join centroid coords to seg_df
cent_xy <- sf::st_coordinates(seg_cent_sf)

cent_df <- tibble::tibble(
  segment_id = seg_cent_sf$segment_id,
  x = cent_xy[,1],
  y = cent_xy[,2]
)

seg_df <- seg_df %>% left_join(cent_df, by="segment_id")

if (any(!is.finite(seg_df$x)) || any(!is.finite(seg_df$y))) stop("Centroid coordinates contain non-finite values.")

# lon/lat guard for partitioning
if (sf::st_is_longlat(seg_cent_sf)) {
  if (identical(partition_longlat_mode, "allow")) {
    msg("WARNING: Segments are lon/lat; partitioning will run in degrees (allow mode).")
  } else {
    msg("Segments are lon/lat. Reprojecting to metric CRS for sperrorest partitioning ...")
    bb <- sf::st_bbox(seg_cent_sf)
    lon_mid <- (bb["xmin"] + bb["xmax"]) / 2
    lat_mid <- (bb["ymin"] + bb["ymax"]) / 2
    utm_zone <- floor((lon_mid + 180) / 6) + 1
    epsg <- if (lat_mid >= 0) 32600 + utm_zone else 32700 + utm_zone
    msg("Using UTM EPSG:", epsg)

    # reproject centroid coords and segment features for partitioning only
    seg_cent_sf2 <- sf::st_transform(seg_cent_sf, epsg)
    xy2 <- sf::st_coordinates(seg_cent_sf2)
    seg_df <- seg_df %>%
      select(-x, -y) %>%
      left_join(tibble(segment_id = seg_cent_sf2$segment_id, x = xy2[,1], y = xy2[,2]), by="segment_id")
  }
}

xy_df <- seg_df[, c("x","y")]
names(xy_df) <- c("X", "Y")

x_partition <- sperrorest::partition_kmeans(
  xy_df,
  coords = c("X", "Y"),
  nfold = k_folds,
  seed1 = 1,
  return_factor = TRUE
)

seg_df$partition <- as.vector(x_partition[[1]])
if (anyNA(seg_df$partition)) stop("sperrorest partition returned NA assignments.")

# -----------------------------------------------------------------------------
# 7) Subsample (keep 040 behaviour)
# -----------------------------------------------------------------------------

if (is.finite(subsample_p) && subsample_p > 0 && subsample_p < 1) {
  msg("Subsampling segments with p=", subsample_p, " (stratified by partition) ...")
  set.seed(100)
  idx <- caret::createDataPartition(seg_df$partition, p = subsample_p, list = FALSE)
  seg_df_small <- seg_df[idx, , drop = FALSE]
} else {
  seg_df_small <- seg_df
}

seg_df_small$class <- droplevels(seg_df_small$class)

# optional safety: gleiche Levels für class_id-map etc. später

msg("Segments used for model selection/training: ", nrow(seg_df_small))

# -----------------------------------------------------------------------------
# 8) CV folds from CAST using spacevar=partition (040 style)
# -----------------------------------------------------------------------------

folds <- CAST::CreateSpacetimeFolds(
  seg_df_small,
  spacevar = "partition",
  k = k_folds
)

ctrl <- caret::trainControl(
  method          = "cv",
  index           = folds$index,
  indexOut        = folds$indexOut,
  savePredictions = "final",
  classProbs      = TRUE,
  summaryFunction = caret::multiClassSummary,
  allowParallel   = TRUE
)

# -----------------------------------------------------------------------------
# 9) Step-1: preliminary random tuning (040 style; adapted to classification)
# -----------------------------------------------------------------------------

msg("Step-1: preliminary random tuning (tuneLength=", random_tuneLen, ", ntree=", random_ntree, ") ...")

cl1 <- parallel::makeCluster(16)
doParallel::registerDoParallel(cl1)
set.seed(seed)

rf_random <- caret::train(
  x         = seg_df_small[, pred_names, drop = FALSE],
  y         = seg_df_small$class,
  method    = "rf",
  metric    = "Kappa",
  trControl = ctrl,
  tuneLength = random_tuneLen,
  ntree     = random_ntree
)

parallel::stopCluster(cl1)

msg("Best mtry from Step-1: ", rf_random$bestTune$mtry)

# -----------------------------------------------------------------------------
# 10) Step-2: FFS (040 style)
# -----------------------------------------------------------------------------

msg("Step-2: CAST::ffs with mtry fixed to Step-1 bestTune; ntree=", ffs_ntree, " ...")

cl2 <- parallel::makeCluster(16)
doParallel::registerDoParallel(cl2)
set.seed(seed)

rf_ffs <- CAST::ffs(
  predictors = seg_df_small[, pred_names, drop = FALSE],
  response   = seg_df_small$class,
  method     = "rf",
  metric     = "Kappa",
  trControl  = ctrl,
  tuneGrid   = data.frame(mtry = rf_random$bestTune$mtry),
  ntree      = ffs_ntree,
  verbose    = TRUE
)

parallel::stopCluster(cl2)

sel_preds <- as.character(rf_ffs$selectedvars)
sel_preds <- sel_preds[sel_preds %in% pred_names]
if (length(sel_preds) == 0) stop("FFS selected no predictors.")

msg("Selected predictors: ", paste(sel_preds, collapse = ", "))

# -----------------------------------------------------------------------------
# 11) Step-3: final sweep over ntree × mtry (040 style) using selected predictors
# -----------------------------------------------------------------------------

msg("Step-3: hyperparameter sweep over ntree × mtry (selected predictors only) ...")

mtry_hi <- min(mtry_max, length(sel_preds))
mtry_vals <- seq.int(2L, mtry_hi)

models <- list()
grid_log <- list()

cl3 <- parallel::makeCluster(16)
doParallel::registerDoParallel(cl3)

for (nt in ntree_grid) {
  for (mt in mtry_vals) {
    set.seed(seed)
    rf_tmp <- caret::train(
      x         = seg_df_small[, sel_preds, drop = FALSE],
      y         = seg_df_small$class,
      method    = "rf",
      metric    = "Kappa",
      trControl = ctrl,
      tuneGrid  = data.frame(mtry = mt),
      ntree     = as.integer(nt)
    )
    key <- sprintf("ntree_%d__mtry_%d", as.integer(nt), as.integer(mt))
    models[[key]] <- rf_tmp
    grid_log[[key]] <- rf_tmp$results
    msg("  trained: ", key)
  }
}

parallel::stopCluster(cl3)

# Compare via resamples (as in 040)
resamp <- caret::resamples(models)
res_sum <- summary(resamp, metric = "Kappa")

# choose best by max mean Kappa across resamples
kappa_means <- res_sum$statistics$Kappa[, "Mean"]
best_key <- names(which.max(kappa_means))
msg("Best model from sweep: ", best_key, " (Mean Kappa=", max(kappa_means), ")")

rf_best <- models[[best_key]]
best_mtry <- rf_best$bestTune$mtry

# also parse ntree from key
best_ntree <- as.integer(
  sub("ntree_(\\d+)__mtry_\\d+", "\\1", best_key)
)


# -----------------------------------------------------------------------------
# 12) Step-4: final refit (040 style)
# -----------------------------------------------------------------------------

msg("Step-4: final refit with best hyperparams (mtry=", best_mtry, ", ntree=", best_ntree, ") ...")

cl4 <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl4)
set.seed(seed)

rf_prod <- caret::train(
  x         = seg_df_small[, sel_preds, drop = FALSE],
  y         = seg_df_small$class,
  method    = "rf",
  metric    = "Kappa",
  trControl = ctrl,
  tuneGrid  = data.frame(mtry = best_mtry),
  ntree     = best_ntree
)

parallel::stopCluster(cl4)

saveRDS(list(model = rf_prod, ffs = rf_ffs, random = rf_random, sweep_summary = res_sum),
        out_model_rds)

cv_res <- rf_prod$results %>% as_tibble() %>% arrange(desc(Kappa))
readr::write_csv(cv_res, out_cv_csv)

# -----------------------------------------------------------------------------
# 13) Prediction on ALL segments (segment-wise probs + class/class_id)
# -----------------------------------------------------------------------------

msg("Predicting on all segments (segment-wise) ...")

# Build features for ALL segments (already extracted into seg_feat)
seg_all <- segments %>% select(segment_id)
seg_all_df <- seg_feat %>% select(segment_id, all_of(sel_preds)) %>%
  filter(if_all(all_of(sel_preds), ~ is.finite(.x)))

# predict probabilities
clp <- parallel::makeCluster(n_cores_predict)
doParallel::registerDoParallel(clp)
set.seed(seed)

pp <- predict(rf_prod, newdata = seg_all_df[, sel_preds, drop=FALSE], type = "prob")

parallel::stopCluster(clp)

# derive class label + class_id
cls_levels <- colnames(pp)
mx_idx <- max.col(pp, ties.method = "first")
pred_class <- cls_levels[mx_idx]
pred_pmax  <- pp[cbind(seq_len(nrow(pp)), mx_idx)]

class_map <- seg_df_small %>% distinct(class, class_id) %>% mutate(class = as.character(class))
pred_class_id <- class_map$class_id[match(pred_class, class_map$class)]

pred_tab <- tibble(
  segment_id = seg_all_df$segment_id,
  pred_class = pred_class,
  pred_class_id = pred_class_id,
  pred_pmax = as.numeric(pred_pmax)
)

# attach probs (optional) as separate columns
pp_df <- as_tibble(pp)
names(pp_df) <- paste0("p_", names(pp_df))
pred_tab <- bind_cols(pred_tab, pp_df)

# join to segments and write
seg_out <- segments %>% left_join(pred_tab, by="segment_id")
sf::st_write(seg_out, out_seg_gpkg, delete_dsn = TRUE, quiet = TRUE)

# -----------------------------------------------------------------------------
# 14) Optional: rasterize segment class_id to canonical out_class_file
# -----------------------------------------------------------------------------

msg("Rasterizing segment predictions to class_id raster ...")
tpl <- pred_stack[[1]]
seg_v2 <- terra::vect(seg_out)

class_r <- terra::rasterize(seg_v2, tpl, field = "pred_class_id", touches = TRUE)
terra::writeRaster(class_r, out_class_file, overwrite = TRUE)

# Write probability max raster (single band) as a convenience, and keep multi-class probs only in GPKG
pmax_r <- terra::rasterize(seg_v2, tpl, field = "pred_pmax", touches = TRUE)
terra::writeRaster(pmax_r, out_prob_file, overwrite = TRUE)

# -----------------------------------------------------------------------------
# 15) Optional: AOA/DI/LPD in segment feature space (not raster space)
# -----------------------------------------------------------------------------

msg("Computing CAST::aoa in segment feature space (DI/LPD) ...")
AOA_seg <- CAST::aoa(seg_all_df[, sel_preds, drop = FALSE], rf_prod, LPD = TRUE, verbose = TRUE)

# Attach DI/LPD/AOA to segment output if present
try({
  if (!is.null(AOA_seg$DI)) seg_out$DI  <- AOA_seg$DI
  if (!is.null(AOA_seg$LPD)) seg_out$LPD <- AOA_seg$LPD
  if (!is.null(AOA_seg$AOA)) seg_out$AOA <- AOA_seg$AOA
  sf::st_write(seg_out, out_seg_gpkg, delete_dsn = TRUE, quiet = TRUE)
}, silent = TRUE)

msg("DONE.")

