#!/usr/bin/env Rscript

# =============================================================================
# Script:   03-3_classification_rf_cast_pixelwise.R
# Project:  Burgwald Decision Stack
#
# PURPOSE
# -------
# Pixel-wise landcover classification (S2 output) using Random Forest (RF),
# trained from a persisted S2 training product:
#
#   paths[["s2_lc_rf_trainpts_2021"]]  (GPKG; class + predictor columns)
#
# ARCHITECTURE NOTE
# -----------------
# - Predictors are consumed as one canonical S2 product: s2_pred_stack_2021_rf
# - Training points are generated in S2 by 02-5-processing_predictor-stacks.R
# - This script performs model training, spatial CV (CAST), and prediction.
#
# CAST ALIGNMENT
# --------------
# Spatial CV folds are grouped by segment_id to reduce spatial leakage.
#
# INPUTS (canonical via paths[])
# -----------------------------
# Predictors (RF-ready stack produced in S2):
#   paths[["s2_pred_stack_2021_rf"]]
#
# Training points (S2; persisted):
#   paths[["s2_lc_rf_trainpts_2021"]]
#
# Segmentation polygons (S3; only for CAST grouping by segment_id):
#   paths[["layer0_segments"]]
#
# OUTPUTS (must exist in outputs.tsv / paths[])
# ---------------------------------------------
# Pixel class map (S2):
#   paths[["s2_lc_rf_2021"]]
#
# Pixel probability stack (S2):
#   paths[["s2_lc_rf_prob_2021"]]
#
# Model artefacts:
#   paths[["s2_lc_rf_model_2021"]]
#   paths[["s2_lc_rf_cv_metrics_2021"]]
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
})

source(here::here("src", "_core", "01-setup-burgwald.R"))

# -----------------------------------------------------------------------------
# 0) Strict IO (paths[] only)
# -----------------------------------------------------------------------------

pred_stack_file  <- paths[["s2_pred_stack_2021_rf"]]
train_pts_file   <- paths[["s2_lc_rf_trainpts_2021"]]
seg_file         <- paths[["layer0_segments"]]

stopifnot(file.exists(pred_stack_file))
stopifnot(file.exists(train_pts_file))
stopifnot(file.exists(seg_file))

# Outputs
out_class_file  <- paths[["s2_lc_rf_2021"]]
out_prob_file   <- paths[["s2_lc_rf_prob_2021"]]
out_model_rds   <- paths[["s2_lc_rf_model_2021"]]
out_cv_csv      <- paths[["s2_lc_rf_cv_metrics_2021"]]

dir.create(dirname(out_class_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_prob_file),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_model_rds),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_cv_csv),     recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1) Parameters (explicit)
# -----------------------------------------------------------------------------

set.seed(1)

# Spatial CV design (CAST)
k_folds <- 5L

# ranger tuning constants (mtry is derived from p later, see section 3)
rf_splitrule_vals     <- "gini"
rf_min_node_size_vals <- c(1, 5, 10)

# -----------------------------------------------------------------------------
# 3) Load data + derive canonical RF tuning grid from predictor stack
# -----------------------------------------------------------------------------

pred_stack <- terra::rast(pred_stack_file)
segments   <- sf::st_read(seg_file, quiet = TRUE)
stopifnot("segment_id" %in% names(segments))

if (any(!sf::st_is_valid(segments))) segments <- sf::st_make_valid(segments)

pred_names <- names(pred_stack)

p <- terra::nlyr(pred_stack)
stopifnot(is.finite(p) && p >= 2)

mtry_raw <- unique(as.integer(round(c(
  sqrt(p),
  p / 3,
  0.10 * p,
  0.20 * p,
  0.30 * p,
  0.50 * p,
  0.70 * p
))))
mtry_vals <- sort(unique(mtry_raw[mtry_raw >= 1 & mtry_raw <= p]))
if (length(mtry_vals) < 3) {
  mtry_vals <- sort(unique(pmax(1L, pmin(p, as.integer(c(1, floor(p / 2), p))))))
}

rf_grid <- expand.grid(
  mtry          = mtry_vals,
  splitrule     = rf_splitrule_vals,
  min.node.size = rf_min_node_size_vals
)

message("Predictor bands (p): ", p)
message("Derived mtry grid: ", paste(mtry_vals, collapse = ", "))

# -----------------------------------------------------------------------------
# 4) Load persisted S2 training points, attach segment_id (CAST grouping)
# -----------------------------------------------------------------------------

train_pts <- sf::st_read(train_pts_file, quiet = TRUE)
if (nrow(train_pts) == 0) stop("Training point file is empty: ", train_pts_file)
if (!("class" %in% names(train_pts))) stop("Training point file has no 'class' column: ", train_pts_file)

# predictor columns must exist in training points
missing <- setdiff(pred_names, names(train_pts))
if (length(missing) > 0) {
  stop("Training points missing predictor columns: ", paste(missing, collapse = ", "))
}

# CRS align for join
if (sf::st_crs(train_pts)$wkt != sf::st_crs(segments)$wkt) {
  train_pts <- sf::st_transform(train_pts, sf::st_crs(segments))
}
if (any(!sf::st_is_valid(train_pts))) train_pts <- sf::st_make_valid(train_pts)

# attach segment_id for CAST fold grouping
train_pts_join <- sf::st_join(
  train_pts,
  segments[, c("segment_id")],
  join = sf::st_within,
  left = FALSE
)
if (nrow(train_pts_join) == 0) {
  stop("No training points fall within segments. Check CRS / alignment / coverage.")
}

train_df <- train_pts_join %>%
  sf::st_drop_geometry() %>%
  mutate(class = as.character(class)) %>%
  select(segment_id, class, all_of(pred_names)) %>%
  filter(if_all(all_of(pred_names), ~ is.finite(.x)))

if (nrow(train_df) == 0) stop("Training table empty after finite filtering.")
train_df$class <- factor(train_df$class)

message(sprintf("Training table loaded: n=%d  p=%d predictors", nrow(train_df), length(pred_names)))

# -----------------------------------------------------------------------------
# 5) Spatial CV folds (CAST)
# -----------------------------------------------------------------------------

folds <- CAST::CreateSpacetimeFolds(
  x        = train_df,
  spacevar = "segment_id",
  timevar  = NA,
  k        = k_folds,
  class    = "class",
  seed     = 1
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
# 6) Train RF (ranger) with caret + CAST folds
# -----------------------------------------------------------------------------

rf_fit <- caret::train(
  x          = train_df[, pred_names],
  y          = train_df$class,
  method     = "ranger",
  trControl  = ctrl,
  tuneGrid   = rf_grid,
  importance = "impurity",
  metric     = "ROC"
)

saveRDS(rf_fit, out_model_rds)

cv_res <- rf_fit$results %>%
  as_tibble() %>%
  arrange(desc(ROC))

readr::write_csv(cv_res, out_cv_csv)

print(rf_fit)
print(cv_res, n = Inf)

# -----------------------------------------------------------------------------
# 7) Predict pixel-wise class probabilities and hard class map
# -----------------------------------------------------------------------------

predict_raster_probs <- function(rast_stack, model, filename, overwrite = TRUE) {
  levs <- levels(model$trainingData$.outcome)[[1]]
  ncls <- length(levs)
  
  out <- rast_stack[[1]]
  out <- terra::rast(out)
  out <- terra::setValues(out, NA_real_)
  out <- terra::wrap(out)
  out <- out[[rep(1, ncls)]]
  names(out) <- paste0("p_", levs)
  
  bs <- terra::blockSize(rast_stack)
  out <- terra::writeStart(out, filename = filename, overwrite = overwrite)
  
  for (i in seq_len(bs$n)) {
    v <- terra::values(rast_stack, row = bs$row[i], nrows = bs$nrows[i], mat = TRUE)
    ok <- apply(v, 1, function(x) all(is.finite(x)))
    
    pred_block <- matrix(NA_real_, nrow = nrow(v), ncol = ncls)
    if (any(ok)) {
      df <- as.data.frame(v[ok, , drop = FALSE])
      colnames(df) <- pred_names
      pp <- predict(model, newdata = df, type = "prob")
      pred_block[ok, ] <- as.matrix(pp)
    }
    
    out <- terra::writeValues(out, pred_block, bs$row[i])
  }
  
  terra::writeStop(out)
  terra::rast(filename)
}

# 7.1 Probabilities
prob_r <- predict_raster_probs(pred_stack, rf_fit, out_prob_file, overwrite = TRUE)

# 7.2 Hard class (argmax of probabilities)
prob_mat_to_class <- function(p) {
  if (!is.matrix(p)) p <- as.matrix(p)
  apply(p, 1, function(x) if (all(is.na(x))) NA_integer_ else which.max(x))
}

bs <- terra::blockSize(prob_r)
class_r <- prob_r[[1]]
class_r <- terra::setValues(class_r, NA_integer_)
class_r <- terra::writeStart(class_r, filename = out_class_file, overwrite = TRUE)

for (i in seq_len(bs$n)) {
  pmat <- terra::values(prob_r, row = bs$row[i], nrows = bs$nrows[i], mat = TRUE)
  cls  <- prob_mat_to_class(pmat)
  class_r <- terra::writeValues(class_r, cls, bs$row[i])
}
class_r <- terra::writeStop(class_r)

names(class_r) <- "class_id"

message("Predictor stack used:    ", pred_stack_file)
message("Training points used:    ", train_pts_file)
message("Wrote probability stack: ", out_prob_file)
message("Wrote class raster:      ", out_class_file)
