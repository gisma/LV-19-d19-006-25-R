#!/usr/bin/env Rscript

# =============================================================================
# Script:   03-3_classification_rf_cast_pixelwise.R
# Project:  Burgwald Decision Stack
#
# PURPOSE
# -------
# Pixel-wise landcover classification (S2 output) using a Random Forest (RF),
# trained from point-based ground truth (GT) that is *propagated to segments*
# (S3) to obtain high-quality training masks.
#
# Key idea (your design):
# - You want a pixel model (S2 classes), but you want segments as an annotation
#   carrier to speed up training data creation and improve boundary quality.
# - Segments are NOT the prediction object; they are only used to define
#   where GT labels are trusted (label propagation / weak supervision).
#
# CAST alignment:
# - We use spatial cross-validation folds (CAST) to avoid overly optimistic
#   performance due to spatial autocorrelation and leakage.
# - Crucially: folds are created on the *segment level* (or spatial blocks),
#   so pixels from the same segment cannot end up in train and test at once.
#
# INPUTS (canonical via paths[])
# -----------------------------
# - Predictor raster stack (S2 features):
#     paths[["s2_pred_stack_2021"]]   (multi-band GeoTIFF; 10 m)
#
# - Segmentation polygons (S3):
#     paths[["layer0_segments"]]      (GPKG; segment_id + geometry)
#
# - Ground truth points (S1 labels):
#     paths[["gt_points_lc"]]         (GPKG; class + geometry)
#   Required attribute columns in GT:
#     - class  (factor or character; categorical land-cover class label
#               used as training target, e.g. "deciduous", "mixed",
#               "coniferous", "grassland", "cropland")
# OUTPUTS (define in outputs.tsv and map to paths[])
# -------------------------------------------------
# - Pixel class map (S2 label raster):
#     paths[["s2_lc_rf_2021"]]
#
# - Optional: class probability stack (S2):
#     paths[["s2_lc_rf_prob_2021"]]
#
# - Model report (csv/rds):
#     paths[["s2_lc_rf_model_2021"]]
#     paths[["s2_lc_rf_cv_metrics_2021"]]
#
# NOTES
# -----
# - This script is written as an operational, reproducible control script
#   in the CAST spirit (caret + spatial CV folds).
# - It does NOT introduce new abstractions; everything is explicit.
#
# REFERENCES
# ----------
# CAST package introduction and spatial CV concept:
#   https://cran.r-project.org/web/packages/CAST/vignettes/cast01-CAST-intro.html
# CreateSpacetimeFolds documentation (used here for spatial folds):
#   https://www.rdocumentation.org/packages/CAST/topics/CreateSpacetimeFolds
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
pred_stack_file <- paths[["s2_pred_stack_2021"]]
seg_file        <- paths[["layer0_segments"]]
gt_file         <- paths[["gt_points_lc"]]

stopifnot(file.exists(pred_stack_file))
stopifnot(file.exists(seg_file))
stopifnot(file.exists(gt_file))

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

# Vegetation threshold (if you later want rules using CHM; kept explicit)
veg_thr_m <- 2.0

# Segment-label propagation rules:
# - min_points_per_segment: discard segments with too few GT points
# - purity_thr: require dominant class fraction >= purity_thr
min_points_per_segment <- 2L
purity_thr             <- 0.70

# Pixel sampling design:
# - max_pixels_per_segment: limit to avoid huge segments dominating
# - class_balance: if TRUE, downsample per class after extraction
max_pixels_per_segment <- 400L
class_balance          <- TRUE
target_n_per_class     <- 20000L

# Spatial CV design (CAST):
# - k folds
# - spacevar: we use a spatial grouping variable (segment_id or spatial blocks)
k_folds <- 5L

# Random Forest tuning grid (ranger via caret):
# - mtry tuned; splitrule fixed; min.node.size tuned
# You can extend this grid if you want more exhaustive tuning.
rf_grid <- expand.grid(
  mtry          = c(2, 4, 6, 8),          # adapt to number of predictors
  splitrule     = "gini",
  min.node.size = c(1, 5, 10)
)

# -----------------------------------------------------------------------------
# 2) Load data
# -----------------------------------------------------------------------------
pred_stack <- terra::rast(pred_stack_file)         # 10 m predictors
segments   <- sf::st_read(seg_file, quiet = TRUE) # S3 polygons
gt_points  <- sf::st_read(gt_file, quiet = TRUE)  # S1 points with class

stopifnot("segment_id" %in% names(segments))
stopifnot("class"      %in% names(gt_points))

# Harmonize CRS
if (sf::st_crs(gt_points)$wkt != sf::st_crs(segments)$wkt) {
  gt_points <- sf::st_transform(gt_points, sf::st_crs(segments))
}

# Optional: ensure valid geometries
if (any(!sf::st_is_valid(segments)))  segments  <- sf::st_make_valid(segments)
if (any(!sf::st_is_valid(gt_points))) gt_points <- sf::st_make_valid(gt_points)

# -----------------------------------------------------------------------------
# 3) Point-in-polygon: assign GT points to segments (segment_id)
# -----------------------------------------------------------------------------
# This is the key step: points define labels; segments define trusted areas.
gt_join <- sf::st_join(
  gt_points,
  segments[, c("segment_id")],
  join = sf::st_within,
  left = FALSE
)

if (nrow(gt_join) == 0) {
  stop("No GT points fall within segments. Check CRS/alignment or GT coverage.")
}

# -----------------------------------------------------------------------------
# 4) Propagate GT to segment labels (explicit rules)
# -----------------------------------------------------------------------------
# For each segment:
#   - count points per class
#   - determine dominant class + purity
#   - accept label only if min_points and purity threshold satisfied
seg_labels <- gt_join %>%
  st_drop_geometry() %>%
  mutate(class = as.character(class)) %>%
  count(segment_id, class, name = "n") %>%
  group_by(segment_id) %>%
  mutate(n_total = sum(n),
         frac    = n / n_total) %>%
  arrange(segment_id, desc(frac), desc(n)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    accept = (n_total >= min_points_per_segment) & (frac >= purity_thr)
  )

# Keep only accepted training segments
train_segments <- seg_labels %>%
  filter(accept) %>%
  select(segment_id, class, n_total, frac)

if (nrow(train_segments) == 0) {
  stop("No segments passed GT propagation rules. Lower purity_thr or min_points_per_segment.")
}

# Attach labels to segment geometry (only training subset)
train_pol <- segments %>%
  inner_join(train_segments, by = "segment_id")

# -----------------------------------------------------------------------------
# 5) Pixel sampling from labelled segments (training data on pixel level)
# -----------------------------------------------------------------------------
# We sample pixel values from pred_stack inside each labelled segment polygon.
# Important: training rows are pixels; target is the propagated segment class.

# Convert sf polygons to terra vector for extraction
train_vect <- terra::vect(train_pol)

# Extract raster values for all pixels under training polygons
# terra::extract returns a data.frame with polygon ID mapping.
# We request cell=TRUE to later allow deduplication / balancing if needed.
ext <- terra::extract(
  pred_stack,
  train_vect,
  cells = TRUE
)

# ext includes an ID column mapping to polygons in train_vect
# We map that back to segment_id + class.
id_map <- train_pol %>%
  mutate(.poly_id = seq_len(n())) %>%
  st_drop_geometry() %>%
  select(.poly_id, segment_id, class)

ext <- ext %>%
  as_tibble() %>%
  rename(.poly_id = ID)

# Join labels
train_df <- ext %>%
  inner_join(id_map, by = ".poly_id") %>%
  select(segment_id, class, cell, everything(), -.poly_id)

# Remove rows with NA predictors
pred_names <- names(pred_stack)
train_df <- train_df %>%
  filter(if_all(all_of(pred_names), ~ is.finite(.x)))

if (nrow(train_df) == 0) {
  stop("Training extraction produced 0 valid pixel rows (all NA?).")
}

# Limit pixels per segment to avoid dominance by huge polygons
train_df <- train_df %>%
  group_by(segment_id) %>%
  slice_sample(n = min(n(), max_pixels_per_segment)) %>%
  ungroup()

# Optional: class balancing (downsample to target_n_per_class)
if (class_balance) {
  train_df <- train_df %>%
    group_by(class) %>%
    slice_sample(n = min(n(), target_n_per_class)) %>%
    ungroup()
}

train_df$class <- factor(train_df$class)

# -----------------------------------------------------------------------------
# 6) Spatial CV folds (CAST) to avoid leakage
# -----------------------------------------------------------------------------
# We must not split pixels randomly:
# - pixels from the same segment are extremely correlated
# - segments are spatially autocorrelated
#
# Here: use segment_id as the grouping variable for folds.
# CAST's CreateSpacetimeFolds can create folds by group variables.
folds <- CAST::CreateSpacetimeFolds(
  x        = train_df,
  spacevar = "segment_id",
  timevar  = NA,
  k        = k_folds,
  class    = "class",
  seed     = 1
)

# folds$index and folds$indexOut are directly usable in caret::trainControl
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
# 7) Train RF (ranger) with caret + CAST folds
# -----------------------------------------------------------------------------
# Use ranger backend because it is fast and supports probability prediction.
rf_fit <- caret::train(
  x          = train_df[, pred_names],
  y          = train_df$class,
  method     = "ranger",
  trControl  = ctrl,
  tuneGrid   = rf_grid,
  importance = "impurity",
  metric     = "ROC"  # multiClassSummary provides ROC as mean of one-vs-all AUCs
)

# Persist model object
saveRDS(rf_fit, out_model_rds)

# Save CV predictions + metrics summary
cv_res <- rf_fit$results %>%
  as_tibble() %>%
  arrange(desc(ROC))

readr::write_csv(cv_res, out_cv_csv)

print(rf_fit)
print(cv_res, n = Inf)

# -----------------------------------------------------------------------------
# 8) Predict pixel-wise class map (S2 output)
# -----------------------------------------------------------------------------
# Predict class probabilities on the full predictor raster.
# We do not use terra::predict(method="caret") because caret objects are not
# always handled robustly; instead we use predict() and chunk manually.

# Helper: chunked prediction to avoid memory blowups
predict_raster_probs <- function(rast_stack, model, filename, overwrite = TRUE) {
  # Prepare output
  levs <- levels(model$trainingData$.outcome)[[1]]
  ncls <- length(levs)
  
  # Create empty SpatRaster with ncls layers
  out <- rast_stack[[1]]
  out <- terra::rast(out)
  out <- terra::setValues(out, NA_real_)
  out <- terra::wrap(out)
  out <- out[[rep(1, ncls)]]
  names(out) <- paste0("p_", levs)
  
  # Write in chunks
  bs <- terra::blockSize(rast_stack)
  out <- terra::writeStart(out, filename = filename, overwrite = overwrite)
  
  for (i in seq_len(bs$n)) {
    v <- terra::values(rast_stack, row = bs$row[i], nrows = bs$nrows[i], mat = TRUE)
    
    # v is matrix: (n_cells_in_block x n_predictors)
    # Handle all-NA rows quickly
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

# 8.1 Probabilities
prob_r <- predict_raster_probs(pred_stack, rf_fit, out_prob_file, overwrite = TRUE)

# 8.2 Hard class (argmax over probabilities)
# Convert probability stack to class raster (factor levels)
prob_mat_to_class <- function(p) {
  # p: numeric matrix (n x k)
  if (!is.matrix(p)) p <- as.matrix(p)
  m <- apply(p, 1, function(x) if (all(is.na(x))) NA_integer_ else which.max(x))
  m
}

bs <- terra::blockSize(prob_r)
class_r <- prob_r[[1]]
class_r <- terra::setValues(class_r, NA_integer_)
class_r <- terra::writeStart(class_r, filename = out_class_file, overwrite = TRUE)

levs <- sub("^p_", "", names(prob_r))

for (i in seq_len(bs$n)) {
  p <- terra::values(prob_r, row = bs$row[i], nrows = bs$nrows[i], mat = TRUE)
  cls <- prob_mat_to_class(p)
  class_r <- terra::writeValues(class_r, cls, bs$row[i])
}
class_r <- terra::writeStop(class_r)

# Store class as RAT / levels (optional)
# (terra factor handling can be fragile across GDAL; keep minimal)
names(class_r) <- "class_id"

message("Wrote probability stack: ", out_prob_file)
message("Wrote class raster:      ", out_class_file)

# -----------------------------------------------------------------------------
# 9) Optional: Area of Applicability (AOA) (CAST)
# -----------------------------------------------------------------------------
# AOA helps detect where the model extrapolates beyond training feature space.
# This is useful for mapping trustworthiness of the pixel classification.
#
# NOTE: AOA can be expensive for big rasters. Use only if needed.
#
# Example (commented by default):
#
# train_features <- train_df[, pred_names]
# aoa <- CAST::aoa(
#   newdata   = pred_stack,      # can be SpatRaster
#   model     = rf_fit,
#   training  = train_features,
#   response  = train_df$class
# )
# terra::writeRaster(aoa$AOA, paths[["s2_lc_rf_aoa_2021"]], overwrite = TRUE)
# terra::writeRaster(aoa$DI,  paths[["s2_lc_rf_di_2021"]],  overwrite = TRUE)
#
# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------
