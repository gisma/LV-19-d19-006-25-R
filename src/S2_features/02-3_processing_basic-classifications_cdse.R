#!/usr/bin/env Rscript

############################################################
# Script:   02-3_processing_basic-classifications_cdse.R
# Author:   [Your Name]
# Project:  Burgwald
#
# Purpose:
# --------
# Perform classification for ONE CDSE-based predictor stack (single time slice)
# using:
#   * k-means (unsupervised, exploratory)
#   * Maximum Likelihood Classification (MLC) via RStoolbox::superClass(model="mlc")
#   * Random Forest via caret (method="rf")
# and save RF confusion matrix.
#
# Contract rule (project core):
# ----------------------------
# ALL filenames/paths come from metadata/outputs.tsv via `paths` created in
# src/_core/01-setup-burgwald.R.
#
# Reads (S1):
#   - cdse_pred_stack (tif)
#   - training_polygons (gpkg)
# Writes (S2):
#   - classification_kmeans (tif)
#   - classification_mlc (tif)
#   - classification_rf (tif)
#   - classification_rf_model (rds)
#   - classification_rf_confusion (rds)
############################################################

# -----------------------------
# 0) Setup + packages
# -----------------------------
source(here::here("src", "_core", "01-setup-burgwald.R"))

# -----------------------------
# 1) Contracts (outputs.tsv -> paths)
# -----------------------------
pred_stack_file <- paths[["cdse_pred_stack"]]
train_poly_file <- paths[["training_polygons"]]

out_kmeans   <- paths[["classification_kmeans"]]
out_mlc      <- paths[["classification_mlc"]]
out_rf       <- paths[["classification_rf"]]
out_rf_model <- paths[["classification_rf_model"]]
out_rf_conf  <- paths[["classification_rf_confusion"]]

class_field <- "class"

# -----------------------------
# 2) Load predictor stack + training polygons
# -----------------------------
# RStoolbox expects Raster*; exactextractr works fine with terra SpatRaster.
pred_stack_t <- terra::rast(pred_stack_file)
pred_stack_r <- raster::stack(pred_stack_file)

train_areas <- sf::st_read(train_poly_file, quiet = TRUE)

# -----------------------------
# 3) Unsupervised k-means (RStoolbox)
# -----------------------------
set.seed(123)

kmeans_result <- RStoolbox::unsuperClass(
  img       = pred_stack_r,
  nClasses  = 5,
  norm      = TRUE,
  algorithm = "MacQueen"
)

raster::writeRaster(
  kmeans_result$map,
  filename  = out_kmeans,
  overwrite = TRUE
)

# -----------------------------
# 4) Training data extraction (exactextractr)
# -----------------------------
tDF_list <- exactextractr::exact_extract(
  pred_stack_t,
  train_areas,
  force_df      = TRUE,
  include_cell  = TRUE,
  include_xy    = TRUE,
  full_colnames = TRUE,
  include_cols  = class_field
)

tDF <- dplyr::bind_rows(tDF_list)

tDF <- tDF[stats::complete.cases(tDF), ]
tDF[[class_field]] <- as.factor(tDF[[class_field]])

# predictor columns (exclude bookkeeping)
drop_cols <- c("cell", "x", "y", class_field, "coverage_fraction")
predictor_cols <- setdiff(colnames(tDF), drop_cols)

# -----------------------------
# 5) Train/test split
# -----------------------------
set.seed(123)
idx <- caret::createDataPartition(
  y    = tDF[[class_field]],
  p    = 0.7,
  list = FALSE
)

trainDat <- tDF[idx,  c(class_field, predictor_cols, "x", "y")]
testDat  <- tDF[-idx, c(class_field, predictor_cols, "x", "y")]

# -----------------------------
# 6) MLC (RStoolbox::superClass)
# -----------------------------
sp_trainDat <- trainDat
sp_testDat  <- testDat

sp::coordinates(sp_trainDat) <- ~ x + y
sp::coordinates(sp_testDat)  <- ~ x + y

# keep it as in your previous codebase: proj4 from raster stack
sp::proj4string(sp_trainDat) <- sp::CRS(raster::crs(pred_stack_r))
sp::proj4string(sp_testDat)  <- sp::CRS(raster::crs(pred_stack_r))

mlc_model <- RStoolbox::superClass(
  img         = pred_stack_r,
  trainData   = sp_trainDat,
  valData     = sp_testDat,
  responseCol = class_field,
  model       = "mlc",
  tuneLength  = 1,
  verbose     = TRUE
)

raster::writeRaster(
  mlc_model$map,
  filename  = out_mlc,
  overwrite = TRUE
)

# -----------------------------
# 7) Random Forest (caret)
# -----------------------------
ctrl <- caret::trainControl(
  method          = "cv",
  number          = 10,
  savePredictions = TRUE
)

set.seed(123)
rf_model <- caret::train(
  x          = trainDat[, predictor_cols, drop = FALSE],
  y          = trainDat[[class_field]],
  method     = "rf",
  metric     = "Kappa",
  trControl  = ctrl,
  importance = TRUE
)

# Predict RF on raster (use raster::predict with caret model)
rf_map <- raster::predict(
  object    = pred_stack_r,
  model     = rf_model,
  type      = "raw",
  progress  = "text",
  na.rm     = TRUE
)

raster::writeRaster(
  rf_map,
  filename  = out_rf,
  overwrite = TRUE
)

# -----------------------------
# 8) RF confusion matrix + save model
# -----------------------------
rf_pred_test <- stats::predict(rf_model, newdata = testDat[, predictor_cols, drop = FALSE])

cm_rf <- caret::confusionMatrix(
  data      = rf_pred_test,
  reference = testDat[[class_field]]
)

saveRDS(rf_model, out_rf_model)
saveRDS(cm_rf,    out_rf_conf)
