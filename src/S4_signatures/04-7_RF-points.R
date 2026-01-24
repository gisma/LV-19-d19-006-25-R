#!/usr/bin/env Rscript

# =============================================================================
# Script:   03-3_classification_rf_cast_pixelwise.R
# Project:  Burgwald Decision Stack
#
# PURPOSE
# -------
# Pixel-wise landcover classification using Random Forest (caret + randomForest),
# including:
#   - spatial CV via sperrorest blocks (kmeans)
#   - CAST forward feature selection (ffs) with pairwise start (mtry=2)
#   - final RF fit on selected predictors
#   - pixel-wise probability + hard class_id raster
#   - CAST AOA (with DI / LPD) + errorProfiles for DI/LPD
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
})

source(here::here("src", "_core", "01-setup-burgwald.R"))

quiet <- FALSE
msg <- function(...) if (!quiet) message(...)

# -----------------------------------------------------------------------------
# 0) Strict IO (paths[] only)
# -----------------------------------------------------------------------------

pred_stack_file  <- paths[["s2_pred_stack_2021_rf"]]
train_pts_file   <- paths[["s2_lc_rf_trainpts_2021"]]
seg_file         <- paths[["layer0_segments"]]

stopifnot(file.exists(pred_stack_file))
stopifnot(file.exists(train_pts_file))
stopifnot(file.exists(seg_file))

out_class_file  <- paths[["s2_lc_rf_2021"]]
out_prob_file   <- paths[["s2_lc_rf_prob_2021"]]
out_model_rds   <- paths[["s2_lc_rf_model_2021"]]
out_cv_csv      <- paths[["s2_lc_rf_cv_metrics_2021"]]

dir.create(dirname(out_class_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_prob_file),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_model_rds),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_cv_csv),     recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1) Parameters
# -----------------------------------------------------------------------------

set.seed(1)

k_folds    <- 5L
ffs_ntree  <- 500L
ffs_mtry   <- 2L          # pairwise start => mtry=2 in vignette logic
ffs_verbose <- TRUE       # set TRUE if CAST prints progress internally

# -----------------------------------------------------------------------------
# 2) Load predictor stack
# -----------------------------------------------------------------------------

msg("Loading predictor stack: ", pred_stack_file)
pred_stack <- terra::rast(pred_stack_file)
pred_names <- names(pred_stack)

p <- terra::nlyr(pred_stack)
stopifnot(is.finite(p) && p >= 2)
msg("Predictor bands (p): ", p)

# -----------------------------------------------------------------------------
# 3) Load segments (for join / segment_id)
# -----------------------------------------------------------------------------

msg("Loading segments: ", seg_file)
segments <- sf::st_read(seg_file, quiet = TRUE)
stopifnot("segment_id" %in% names(segments))
if (any(!sf::st_is_valid(segments))) segments <- sf::st_make_valid(segments)

# -----------------------------------------------------------------------------
# 4) Load training points, validate required columns, join segment_id
# -----------------------------------------------------------------------------

msg("Loading training points: ", train_pts_file)
train_pts <- sf::st_read(train_pts_file, quiet = TRUE)
if (nrow(train_pts) == 0) stop("Training point file is empty: ", train_pts_file)

req_cols <- c("class", "class_id", "code", "class_name")
missing_req <- setdiff(req_cols, names(train_pts))
if (length(missing_req) > 0) {
  stop("Training point file missing required columns: ", paste(missing_req, collapse = ", "))
}

missing_pred <- setdiff(pred_names, names(train_pts))
if (length(missing_pred) > 0) {
  stop("Training points missing predictor columns: ", paste(missing_pred, collapse = ", "))
}

# CRS align for join
if (sf::st_crs(train_pts)$wkt != sf::st_crs(segments)$wkt) {
  train_pts <- sf::st_transform(train_pts, sf::st_crs(segments))
}
if (any(!sf::st_is_valid(train_pts))) train_pts <- sf::st_make_valid(train_pts)

msg("Joining training points -> segment_id (within) ...")
train_pts_join <- sf::st_join(
  train_pts,
  segments[, c("segment_id")],
  join = sf::st_within,
  left = FALSE
)
if (nrow(train_pts_join) == 0) {
  stop("No training points fall within segments. Check CRS / alignment / coverage.")
}

if (sf::st_is_longlat(train_pts_join)) {
  stop("Training points are in lon/lat degrees. sperrorest partitioning requires a projected CRS (meters). Reproject your segments/train points to UTM (or another metric CRS) before running this script.")
}

xy <- sf::st_coordinates(train_pts_join)

train_df <- train_pts_join %>%
  mutate(x = xy[,1], y = xy[,2]) %>%
  sf::st_drop_geometry() %>%
  mutate(
    class      = as.factor(class),
    class_id   = as.integer(class_id),
    code       = as.integer(code),
    class_name = as.character(class_name)
  ) %>%
  select(x, y, segment_id, class, class_id, code, class_name, all_of(pred_names)) %>%
  filter(
    is.finite(class_id),
    is.finite(x),
    is.finite(y),
    !is.na(class_name),
    !is.na(class),
    if_all(all_of(pred_names), ~ is.finite(.x))
  )

if (nrow(train_df) == 0) stop("Training table empty after finite filtering.")

train_df$class <- droplevels(as.factor(train_df$class))

msg(sprintf("Training table: n=%d  predictors=%d  class_levels=%d",
            nrow(train_df), length(pred_names), nlevels(train_df$class)))

# -----------------------------------------------------------------------------
# 5) Spatial CV folds (sperrorest blocks; leave-one-block-out)
# -----------------------------------------------------------------------------

msg("Building spatial CV folds (sperrorest::partition_kmeans; nfold=", k_folds, ") ...")

xy_df <- train_df[, c("x", "y")]
names(xy_df) <- c("X", "Y")

# returns list (one element per repetition). We use repetition=1.
x_partition <- sperrorest::partition_kmeans(
  xy_df,
  coords = c("X", "Y"),
  nfold = k_folds,
  seed1 = 1,
  return_factor = TRUE
)

fold_fac <- as.vector(x_partition[[1]])
if (anyNA(fold_fac)) stop("sperrorest partition returned NA fold assignments.")

fold_ids <- sort(unique(fold_fac))
all_idx  <- seq_len(nrow(train_df))

# leave-one-block-out over the sperrorest blocks
fold_indexOut <- lapply(fold_ids, function(k) which(fold_fac == k))
fold_index    <- lapply(fold_indexOut, function(test_idx) setdiff(all_idx, test_idx))

names(fold_index)    <- paste0("SP_", fold_ids)
names(fold_indexOut) <- paste0("SP_", fold_ids)

ctrl <- caret::trainControl(
  method          = "cv",
  index           = fold_index,
  indexOut        = fold_indexOut,
  savePredictions = "final",
  classProbs      = TRUE,
  summaryFunction = caret::multiClassSummary,
  allowParallel   = TRUE
)

# -----------------------------------------------------------------------------
# 6) CAST Forward Feature Selection (ffs)
# -----------------------------------------------------------------------------

msg("--------------------------------------------------")
msg("CAST::ffs (Forward Feature Selection)")
msg("Predictors available: ", length(pred_names))
msg("Initial 2-var combinations: ", choose(length(pred_names), 2))
msg("mtry fixed at 2 for pairwise start; ntree=", ffs_ntree)
msg("--------------------------------------------------")

set.seed(1)
ffs_fit <- CAST::ffs(
  predictors = train_df[, pred_names, drop = FALSE],
  response   = train_df$class,
  method     = "rf",
  tuneGrid   = data.frame(mtry = ffs_mtry),
  ntree      = ffs_ntree,
  trControl  = ctrl,
  verbose    = ffs_verbose
)

msg("FFS finished.")
print(ffs_fit)

sel_preds <- as.character(ffs_fit$selectedvars)
sel_preds <- sel_preds[sel_preds %in% pred_names]

pred_stack_sel <- pred_stack[[sel_preds]]

# -----------------------------------------------------------------------------
# 7) Final RF training on selected predictors (caret + randomForest)
# -----------------------------------------------------------------------------

msg("Training final RF (randomForest) on selected predictors ...")

if (length(sel_preds) < 1) stop("No predictors selected by FFS.")
if (length(sel_preds) < 2) msg("Note: only 1 predictor selected; mtry will be forced to 1.")

mtry_grid <- data.frame(
  mtry = unique(pmax(1L, pmin(length(sel_preds), c(1L, 2L, 4L, floor(sqrt(length(sel_preds)))))))
)

rf_fit <- caret::train(
  x         = train_df[, sel_preds, drop = FALSE],
  y         = train_df$class,
  method    = "rf",
  trControl = ctrl,
  tuneGrid  = mtry_grid,
  metric    = "Kappa",
  ntree     = 500L
)

saveRDS(rf_fit, out_model_rds)

cv_res <- rf_fit$results %>%
  as_tibble() %>%
  arrange(desc(Kappa))

readr::write_csv(cv_res, out_cv_csv)

print(rf_fit)
print(cv_res, n = Inf)

# -----------------------------------------------------------------------------
# 8) Predict pixel-wise class probabilities and hard class_id map
# -----------------------------------------------------------------------------

predict_raster_probs <- function(rast_stack, model, filename, pred_names, overwrite = TRUE) {
  levs <- levels(model$trainingData$.outcome)[[1]]
  ncls <- length(levs)
  
  tpl <- rast_stack[[1]]
  out <- terra::rast(tpl, nlyrs = ncls)
  names(out) <- paste0("p_", levs)
  
  blks <- terra::blocks(rast_stack)
  if (is.null(blks$row) || is.null(blks$nrows)) {
    stop("terra::blocks() did not return $row/$nrows as expected.")
  }
  
  out <- terra::writeStart(out, filename = filename, overwrite = overwrite)
  
  for (i in seq_along(blks$row)) {
    v <- terra::values(rast_stack, row = blks$row[i], nrows = blks$nrows[i], mat = TRUE)
    
    ok <- apply(v, 1, function(x) all(is.finite(x)))
    pred_block <- matrix(NA_real_, nrow = nrow(v), ncol = ncls)
    
    if (any(ok)) {
      df <- as.data.frame(v[ok, , drop = FALSE])
      colnames(df) <- pred_names
      
      pp <- predict(model, newdata = df, type = "prob")
      
      # ensure column order matches levs
      pp <- pp[, levs, drop = FALSE]
      pred_block[ok, ] <- as.matrix(pp)
    }
    
    terra::writeValues(out, pred_block, blks$row[i])
    if (i %% 10 == 0) msg("  block ", i, "/", length(blks$row))
  }
  
  out <- terra::writeStop(out)
  out
}

msg("Predicting pixel-wise probabilities -> ", out_prob_file)
prob_r <- predict_raster_probs(pred_stack_sel, rf_fit, out_prob_file, sel_preds, overwrite = TRUE)

# hard class map: argmax over probs -> class label -> class_id mapping
msg("Building hard class_id raster -> ", out_class_file)

# factor levels from training
levs <- levels(rf_fit$trainingData$.outcome)[[1]]

# map class label -> class_id (from training table)
class_map <- train_df %>%
  distinct(class, class_id) %>%
  mutate(class = as.character(class))

# ensure 1-1 mapping
if (any(duplicated(class_map$class))) {
  stop("class -> class_id mapping is not unique in training data.")
}

# argmax index raster
mx <- terra::which.max(prob_r)
mx_vals <- terra::values(mx, mat = FALSE)

# convert to class label then class_id
cls_lab <- levs[mx_vals]
cls_id <- class_map$class_id[match(cls_lab, class_map$class)]

class_r <- terra::rast(mx)
terra::values(class_r) <- cls_id
names(class_r) <- "class_id"

terra::writeRaster(class_r, out_class_file, overwrite = TRUE)

# -----------------------------------------------------------------------------
# 9) CAST AOA (DI/LPD) + error profiles
# -----------------------------------------------------------------------------

msg("Computing CAST::aoa (DI/LPD) ...")
AOA <- CAST::aoa(pred_stack_sel, rf_fit, LPD = TRUE, verbose = TRUE)

# errorProfiles for DI / LPD
msg("Computing CAST::errorProfiles for DI / LPD ...")
DI_error_model  <- CAST::errorProfiles(rf_fit, AOA, variable = "DI")
LPD_error_model <- CAST::errorProfiles(rf_fit, AOA, variable = "LPD")

# save AOA layers next to outputs (keep paths[] discipline? -> use out_prob_file dir)
aoa_dir <- file.path(dirname(out_prob_file), "aoa")
dir.create(aoa_dir, recursive = TRUE, showWarnings = FALSE)

terra::writeRaster(AOA$AOA, file.path(aoa_dir, "aoa_mask.tif"), overwrite = TRUE)
terra::writeRaster(AOA$DI,  file.path(aoa_dir, "aoa_DI.tif"),   overwrite = TRUE)
terra::writeRaster(AOA$LPD, file.path(aoa_dir, "aoa_LPD.tif"),  overwrite = TRUE)

saveRDS(list(DI = DI_error_model, LPD = LPD_error_model), file.path(aoa_dir, "aoa_errorProfiles.rds"))

msg("DONE.")
