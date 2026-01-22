#!/usr/bin/env Rscript

# =============================================================================
# Script:   03-3_classification_rf_cast_pixelwise.R
# Project:  Burgwald Decision Stack
#
# PURPOSE
# -------
# Pixel-wise landcover classification using Random Forest (caret + ranger),
# including:
#   - spatial CV grouped by segment_id (CAST)
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
ffs_verbose <- TRUE      # set TRUE if CAST prints progress internally

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
# 3) Load segments (CAST grouping)
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

train_df <- train_pts_join %>%
  sf::st_drop_geometry() %>%
  mutate(
    class      = as.factor(class),          # keep factor from 02 (now fixed)
    class_id   = as.integer(class_id),
    code       = as.integer(code),
    class_name = as.character(class_name)
  ) %>%
  select(segment_id, class, class_id, code, class_name, all_of(pred_names)) %>%
  filter(
    is.finite(class_id),
    !is.na(class_name),
    !is.na(class),
    if_all(all_of(pred_names), ~ is.finite(.x))
  )

if (nrow(train_df) == 0) stop("Training table empty after finite filtering.")

# ensure factor outcome (no per-row unique suffixes; 02 now guarantees class is sane)
train_df$class <- droplevels(as.factor(train_df$class))

msg(sprintf("Training table: n=%d  predictors=%d  class_levels=%d",
            nrow(train_df), length(pred_names), nlevels(train_df$class)))

# -----------------------------------------------------------------------------
# 4.1) Build level -> class_id lookup (for writing class_id raster)
# -----------------------------------------------------------------------------

class_lookup <- train_df %>%
  distinct(class, class_id, code, class_name)

ambig <- class_lookup %>%
  count(class) %>%
  filter(n != 1)

if (nrow(ambig) > 0) {
  stop("Ambiguous class mapping (class -> class_id not unique). Fix training points: ",
       paste(ambig$class, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 5) Spatial CV folds (CAST)
# -----------------------------------------------------------------------------

msg("Building spatial CV folds (segment_id, k=", k_folds, ") ...")

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
# 6) CAST Forward Feature Selection (ffs)
# -----------------------------------------------------------------------------

msg("--------------------------------------------------")
msg("CAST::ffs (Forward Feature Selection)")
msg("Predictors available: ", length(pred_names))
msg("Initial 2-var combinations: ", choose(length(pred_names), 2))
msg("mtry fixed at 2 for pairwise start; ntree=", ffs_ntree)
msg("--------------------------------------------------")

# IMPORTANT: CAST::ffs uses 'predictors' + 'response' (not x/y, not names)
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

sel_preds <- ffs_fit$optVariables
if (is.null(sel_preds) || length(sel_preds) < 2) {
  stop("FFS did not return >=2 selected predictors. optVariables is empty/too small.")
}

msg("Selected predictors (n=", length(sel_preds), "): ", paste(sel_preds, collapse = ", "))

# subset stack to selected predictors for prediction + AOA
pred_stack_sel <- pred_stack[[sel_preds]]

# -----------------------------------------------------------------------------
# 7) Final RF training on selected predictors (ranger)
# -----------------------------------------------------------------------------

msg("Training final RF (ranger) on selected predictors ...")

rf_fit <- caret::train(
  x          = train_df[, sel_preds, drop = FALSE],
  y          = train_df$class,
  method     = "ranger",
  trControl  = ctrl,
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
# 8) Predict pixel-wise class probabilities and hard class_id map
# -----------------------------------------------------------------------------

predict_raster_probs <- function(rast_stack, model, filename, pred_names, overwrite = TRUE) {
  levs <- levels(model$trainingData$.outcome)[[1]]
  ncls <- length(levs)
  
  out <- rast_stack[[1]]
  out <- terra::rast(out)
  out <- terra::setValues(out, NA_real_)
  out <- terra::wrap(out)
  out <- out[[rep(1, ncls)]]
  names(out) <- paste0("p_", levs)
  
  bs  <- terra::blockSize(rast_stack)
  out <- terra::writeStart(out, filename = filename, overwrite = overwrite)
  
  for (i in seq_len(bs$n)) {
    v  <- terra::values(rast_stack, row = bs$row[i], nrows = bs$nrows[i], mat = TRUE)
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

msg("Predicting probability stack ...")
prob_r <- predict_raster_probs(
  pred_stack_sel,
  rf_fit,
  out_prob_file,
  sel_preds,
  overwrite = TRUE
)

# hard class_id
msg("Writing hard class_id raster (argmax over probs) ...")

levs <- levels(rf_fit$trainingData$.outcome)[[1]]
class_id_by_level <- class_lookup$class_id[match(levs, class_lookup$class)]
if (any(!is.finite(class_id_by_level))) stop("Lookup failed: some model levels missing in class_lookup.")

prob_mat_to_class_id <- function(pmat, class_id_by_level) {
  if (!is.matrix(pmat)) pmat <- as.matrix(pmat)
  apply(pmat, 1, function(x) {
    if (all(is.na(x))) NA_integer_
    else class_id_by_level[which.max(x)]
  })
}

bs <- terra::blockSize(prob_r)
class_r <- prob_r[[1]]
class_r <- terra::setValues(class_r, NA_integer_)
class_r <- terra::writeStart(class_r, filename = out_class_file, overwrite = TRUE)

for (i in seq_len(bs$n)) {
  pmat <- terra::values(prob_r, row = bs$row[i], nrows = bs$nrows[i], mat = TRUE)
  cls  <- prob_mat_to_class_id(pmat, class_id_by_level)
  class_r <- terra::writeValues(class_r, cls, bs$row[i])
}
class_r <- terra::writeStop(class_r)
names(class_r) <- "class_id"

# -----------------------------------------------------------------------------
# 9) AOA + Error Profiles (DI / LPD)
# -----------------------------------------------------------------------------

msg("Computing AOA (LPD=TRUE) ...")
AOA <- CAST::aoa(
  pred_stack_sel,
  rf_fit,
  LPD = TRUE,
  verbose = FALSE
)
msg("AOA finished.")

msg("Computing error profiles (DI / LPD) ...")
DI_error_model <- CAST::errorProfiles(rf_fit, AOA, variable = "DI")
LPD_error_model <- CAST::errorProfiles(rf_fit, AOA, variable = "LPD")
msg("Error profiling finished.")

# -----------------------------------------------------------------------------
# 10) Summary
# -----------------------------------------------------------------------------

msg("Predictor stack used:    ", pred_stack_file)
msg("Training points used:    ", train_pts_file)
msg("Selected predictors:     ", paste(sel_preds, collapse = ", "))
msg("Wrote probability stack: ", out_prob_file)
msg("Wrote class raster:      ", out_class_file)
msg("Wrote model:             ", out_model_rds)
msg("Wrote CV metrics:        ", out_cv_csv)
msg("Done.")

# Optional: return bundle when sourced interactively (harmless under Rscript)
invisible(list(
  ffs = ffs_fit,
  model = rf_fit,
  AOA = AOA,
  DI_error_model = DI_error_model,
  LPD_error_model = LPD_error_model,
  prob = prob_r,
  class = class_r
))
