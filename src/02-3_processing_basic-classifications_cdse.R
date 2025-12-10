#!/usr/bin/env Rscript

############################################################
# Script:   02-3_processing_basic-classifications_cdse.R
# Author:   [Your Name]
# Project:  [Your Project Name]
#
# Purpose:
# --------
# - Perform classification for ONE CDSE-based predictor stack
#   (one time slice, e.g. 2018-06-29) using:
#     * k-means (unsupervised, exploratory)
#     * Maximum Likelihood Classification (MLC, via RStoolbox::superClass)
#     * Random Forest (via caret + randomForest)
#   and compute a confusion matrix for RF.
#
# Assumptions:
# - 01_setup-burgwald.R has already been run in this project.
# - here::here() points to the project root.
# - A CDSE-based predictor stack exists, written by your CDSE script:
#       here("data", "pred_stack_2018.tif")
#   (adjust pred_stack_file below if you use a different name).
# - A training polygon layer exists with a "class" column, e.g.:
#       here("data", "processed", "train_areas_burgwald.gpkg")
#
# NO gdalcubes, NO time-series cube here – only single-date classification.
############################################################

# -----------------------------
# 0) Packages
# -----------------------------
library(here)
library(terra)
library(sf)
library(dplyr)
library(exactextractr)
library(caret)
library(RStoolbox)
library(randomForest)
library(sp)  # for Spatial* objects used by superClass()

message("Project root (here): ", here::here())

# -----------------------------
# 1) Input paths (adjust only if needed)
# -----------------------------

# CDSE-derived predictor stack for a single time slice (e.g. best-summer day 2018)
pred_stack_file <- here::here("data", "pred_stack_2018.tif")

# Training polygons with a "class" field (e.g. clearcut / other)
# If your training file has a different name, change ONLY this path.
train_poly_file <- here::here("data", "processed", "train_areas_burgwald.gpkg")

# Name of the class attribute in the training polygons
class_field <- "class"

# -----------------------------
# 2) Load raster stack + training polygons
# -----------------------------
if (!file.exists(pred_stack_file)) {
  stop("Predictor stack not found: ", pred_stack_file,
       "\nMake sure your CDSE script has written pred_stack_2018.tif into data/.")
}

if (!file.exists(train_poly_file)) {
  stop("Training polygons not found: ", train_poly_file,
       "\nAdapt train_poly_file or export your training data to this path.")
}

pred_stack   <- terra::rast(pred_stack_file)
train_areas  <- sf::st_read(train_poly_file, quiet = TRUE)

message("Loaded predictor stack with ", nlyr(pred_stack), " layers.")
message("Loaded training polygons: ", nrow(train_areas), " features.")

# Quick sanity check that class field exists
if (!(class_field %in% names(train_areas))) {
  stop("Field '", class_field, "' not found in training polygons. Available fields: ",
       paste(names(train_areas), collapse = ", "))
}

# -----------------------------
# 3) Unsupervised k-means (exploratory)
# -----------------------------

# By default, use ALL layers in pred_stack for clustering.
# If you want to restrict, you can specify explicit layer names, e.g.:
# cluster_layers <- pred_stack[[c("B04", "B08", "EVI", "kNDVI", "SAVI")]]
cluster_layers <- pred_stack

set.seed(123)  # reproducibility

kmeans_result <- RStoolbox::unsuperClass(
  img       = cluster_layers,
  nClasses  = 5,
  norm      = TRUE,
  algorithm = "MacQueen"
)

kmeans_map <- kmeans_result$map

# Optional: write to disk for inspection
terra::writeRaster(
  kmeans_map,
  here::here("data", "classification_kmeans_2018.tif"),
  overwrite = TRUE
)

message("Unsupervised k-means classification written to data/classification_kmeans_2018.tif")

# -----------------------------
# 4) Extract training data from CDSE stack
# -----------------------------

# exact_extract returns a list of data.frames (one per polygon). We include:
# - all layers from pred_stack
# - x/y coordinates
# - cell indices
# - the polygon class attribute
tDF_list <- exactextractr::exact_extract(
  pred_stack,
  train_areas,
  force_df       = TRUE,
  include_cell   = TRUE,
  include_xy     = TRUE,
  full_colnames  = TRUE,
  include_cols   = class_field
)

tDF <- dplyr::bind_rows(tDF_list)

# Drop rows with any NA in predictors or missing class
tDF <- tDF[complete.cases(tDF), ]

# Ensure class is a factor
tDF[[class_field]] <- as.factor(tDF[[class_field]])

message("Extracted ", nrow(tDF), " training pixels with class information.")
print(summary(tDF[[class_field]]))

# -----------------------------
# 5) Train/test split (caret)
# -----------------------------

set.seed(123)

# Identify predictor columns (exclude cell/x/y/class/coverage_fraction)
drop_cols <- c("cell", "x", "y", class_field, "coverage_fraction")
predictor_cols <- setdiff(colnames(tDF), drop_cols)

idx <- caret::createDataPartition(
  y = tDF[[class_field]],
  p = 0.7,        # 70% training, 30% test
  list = FALSE
)

trainDat <- tDF[idx, c(class_field, predictor_cols, "x", "y")]
testDat  <- tDF[-idx, c(class_field, predictor_cols, "x", "y")]

# -----------------------------
# 6) Maximum Likelihood Classification (MLC) via superClass()
# -----------------------------

# superClass still expects Spatial* objects for training data
sp_trainDat <- trainDat
sp_testDat  <- testDat

# Convert to SpatialPointsDataFrame
sp::coordinates(sp_trainDat) <- ~ x + y
sp::coordinates(sp_testDat)  <- ~ x + y

# Set CRS from the CDSE stack
sp::proj4string(sp_trainDat) <- terra::crs(pred_stack, proj = TRUE)
sp::proj4string(sp_testDat)  <- terra::crs(pred_stack, proj = TRUE)

# Fit MLC only for this ONE time slice
mlc_model <- RStoolbox::superClass(
  img           = pred_stack,
  trainData     = sp_trainDat,
  valData       = sp_testDat,
  responseCol   = class_field,
  model         = "mlc",
  tuneLength    = 1,
  trainPartition = 0.7,
  verbose       = TRUE
)

mlc_map <- mlc_model$map

terra::writeRaster(
  mlc_map,
  here::here("data", "classification_mlc_2018.tif"),
  overwrite = TRUE
)

message("MLC classification written to data/classification_mlc_2018.tif")

# -----------------------------
# 7) Random Forest classification (caret + randomForest)
# -----------------------------

ctrl <- caret::trainControl(
  method          = "cv",
  number          = 10,
  savePredictions = TRUE
)

set.seed(123)

rf_model <- caret::train(
  x         = trainDat[, predictor_cols],
  y         = trainDat[[class_field]],
  method    = "rf",
  metric    = "Kappa",
  trControl = ctrl,
  importance = TRUE
)

print(rf_model)

# Predict RF on the CDSE stack for this one time slice
rf_map <- terra::predict(pred_stack, rf_model, na.rm = TRUE)

terra::writeRaster(
  rf_map,
  here::here("data", "classification_rf_2018.tif"),
  overwrite = TRUE
)

message("RF classification written to data/classification_rf_2018.tif")

# -----------------------------
# 8) Confusion matrix for RF (test set)
# -------------------------------------------------------------------
#  Confusion matrix ONLY for Random Forest classification
#
#  IMPORTANT:
#  - We compute a confusion matrix for the Random Forest model because:
#       • RF is trained using caret::train(), which explicitly separates
#         training and test data and expects external accuracy assessment.
#
#  - We DO NOT compute a confusion matrix for the Maximum Likelihood (MLC)
#    classification because:
#       • RStoolbox::superClass() already performs its own internal
#         train/validation split and produces accuracy statistics internally.
#       • Running an additional confusion matrix outside superClass()
#         would mix evaluation strategies and lead to inconsistent metrics.
#
#  → Therefore: RF gets an external confusion matrix, MLC does not.
# -------------------------------------------------------------------

# Predict classes for the independent RF test data
rf_pred_test <- predict(rf_model, newdata = testDat)

# Compute confusion matrix for Random Forest
cm_rf <- caret::confusionMatrix(
  data      = rf_pred_test,   # predicted classes
  reference = testDat$class   # reference labels
)

cm_rf
essage("Random Forest confusion matrix (test data):")
print(cm_rf)

# Optionally save the model + confusion matrix
saveRDS(rf_model, here::here("data", "rf_model_2018_cdse.rds"))
saveRDS(cm_rf,    here::here("data", "rf_confusion_2018_cdse.rds"))

message("Done: k-means, MLC, and RF classification for ONE CDSE time slice (pred_stack_2018).")
