#!/usr/bin/env Rscript

############################################################
# Script:   02-5-processing_predictor-stacks.R  (SPATSTAT + CARET + LABEL)
# Project:  Burgwald
#
# Purpose
# -------
#   (1) load predictor stack (terra)
#   (2) load CLC raster and align to predictor grid (near)
#   (3) build AOI polygon from predictor extent
#   (4) sample points in AOI with minimum distance (spatstat::rSSI)
#   (5) map points -> raster cells
#   (6) extract predictors + CLC class
#   (7) map class_id -> code/name via clc_legend, create caret-safe label factor
#   (8) variable reduction (caret)
#   (9) mirror reduction to predictor stack (canonical write)
#  (10) export training points as sf (always overwrite)
############################################################

suppressPackageStartupMessages({
  library(here)
  library(terra)
  library(sf)
  library(caret)
  library(spatstat.geom)
  library(spatstat.random)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))

quiet <- FALSE
msg <- function(...) if (!quiet) message(...)

## -------------------------------------------------------------------
## 1) Canonical IO
## -------------------------------------------------------------------

in_file   <- paths[["s2_pred_stack_2021"]]
out_train <- paths[["s2_lc_rf_trainpts_2021"]]

stopifnot(!is.null(in_file), file.exists(in_file))
stopifnot(!is.null(out_train))

clc_key <- "aoi_clc"
stopifnot(clc_key %in% names(paths))
stopifnot(file.exists(paths[[clc_key]]))

# canonical reduced stack output (must exist; do not invent keys)
out_stack_key <- "s2_pred_stack_2021_rf"
if (!out_stack_key %in% names(paths)) {
  stop("Missing required paths[] key for canonical reduced stack: ", out_stack_key)
}
out_stack_file <- paths[[out_stack_key]]

## -------------------------------------------------------------------
## 2) Sampling + Reduction Parameters  (CENTRAL)
## -------------------------------------------------------------------

# sampling
target_n  <- 10000L
mindist_m <- 10
seed      <- 1L
crs_m     <- 25832   # metric CRS for distance

# variable reduction
corr_cutoff <- 0.95

set.seed(as.integer(seed))

## -------------------------------------------------------------------
## 3) Helpers
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

write_training <- function(sf_obj, out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(sf_obj, out_path, delete_dsn = TRUE, quiet = TRUE)
  invisible(TRUE)
}

## -------------------------------------------------------------------
## 4) Main
## -------------------------------------------------------------------

msg("Loading predictor stack ...")
pred <- terra::rast(in_file)

msg("Loading CLC and aligning to predictor grid (near) ...")
clc  <- align_to_template(terra::rast(paths[[clc_key]]), pred, method = "near")

# --------------------------------------------------------------------
# AOI polygon from raster extent
# --------------------------------------------------------------------
msg("Building AOI polygon from raster extent ...")
aoi_v <- terra::as.polygons(terra::ext(pred), crs = terra::crs(pred))
aoi   <- sf::st_as_sf(aoi_v)

# --------------------------------------------------------------------
# Project everything to metric CRS (for mindist)
# --------------------------------------------------------------------
msg("Projecting to EPSG:", crs_m, " ...")
pred_m <- terra::project(pred, paste0("EPSG:", crs_m))
clc_m  <- terra::project(clc,  paste0("EPSG:", crs_m))
aoi_m  <- sf::st_transform(aoi, crs_m)

# --------------------------------------------------------------------
# Sampling (spatstat rSSI)
# --------------------------------------------------------------------
msg("Sampling points with rSSI: target = ", target_n,
    ", mindist = ", mindist_m, " m")

win <- spatstat.geom::as.owin(sf::st_geometry(aoi_m))

pts_ppp <- spatstat.random::rSSI(
  r   = as.numeric(mindist_m),
  n   = as.integer(target_n),
  win = win
)

xy <- cbind(pts_ppp$x, pts_ppp$y)
colnames(xy) <- c("x", "y")

if (nrow(xy) == 0)
  stop("rSSI returned zero points (AOI too small or mindist too large).")

msg("Generated ", nrow(xy), " sampling points.")

# --------------------------------------------------------------------
# Map points -> raster cells
# --------------------------------------------------------------------
cells <- terra::cellFromXY(pred_m, xy)

ok <- is.finite(cells)
cells <- as.integer(cells[ok])
xy    <- xy[ok, , drop = FALSE]

if (length(cells) == 0)
  stop("No points mapped to raster cells.")

# --------------------------------------------------------------------
# Fast extraction
# --------------------------------------------------------------------
msg("Extracting predictor values ...")
X <- terra::values(pred_m, mat = TRUE)[cells, , drop = FALSE]

msg("Extracting CLC labels ...")
y <- terra::values(clc_m, mat = FALSE)[cells]

# --------------------------------------------------------------------
# CLC legend mapping + caret-safe label factor  (FIXED)
# --------------------------------------------------------------------
if (!exists("clc_legend")) stop("Object 'clc_legend' not found in this script/session.")

req_cols <- c("class_id", "code", "name")
if (!all(req_cols %in% names(clc_legend))) {
  stop("clc_legend must contain columns: ", paste(req_cols, collapse = ", "))
}

class_id   <- as.integer(y)
m          <- match(class_id, clc_legend$class_id)

code       <- clc_legend$code[m]
class_name <- clc_legend$name[m]

# FIX: build labels at CLASS LEVEL (not per-row unique suffixes)
# - one caret-safe label per class_id
# - suffixes only if there are actual collisions between different classes
map <- unique(data.frame(class_id = class_id, class_name = class_name, stringsAsFactors = FALSE))
map <- map[is.finite(map$class_id) & !is.na(map$class_name), , drop = FALSE]
map <- map[order(map$class_id), , drop = FALSE]

map$class_lab <- make.names(map$class_name, unique = TRUE)

class_lab <- map$class_lab[match(class_id, map$class_id)]

train_df <- data.frame(
  x = xy[,1],
  y = xy[,2],
  class_id   = class_id,
  code       = code,
  class_name = class_name,
  class      = factor(class_lab, levels = map$class_lab),
  X,
  check.names = FALSE
)

# NA filter: after building all columns, drop rows with missing class mapping or predictors
train_df <- train_df[
  is.finite(train_df$class_id) &
    !is.na(train_df$class_name) &
    !is.na(train_df$class) &
    stats::complete.cases(train_df),
  , drop = FALSE
]

if (nrow(train_df) == 0)
  stop("All samples dropped after NA filtering (CLC mapping / predictor NA).")

msg("Training rows before reduction: ", nrow(train_df))
msg("Predictors before reduction: ", ncol(train_df) - 6)
msg("Class levels (n): ", nlevels(train_df$class))

# --------------------------------------------------------------------
# Variable reduction (caret) on predictor columns only
# --------------------------------------------------------------------

pred_names <- setdiff(names(train_df), c("x", "y", "class", "class_id", "code", "class_name"))
smp <- train_df[, pred_names, drop = FALSE]

# --- near-zero variance ---
nzv <- caret::nearZeroVar(smp, saveMetrics = TRUE)
drop_nzv <- rownames(nzv)[nzv$zeroVar]

if (length(drop_nzv) > 0) {
  msg("Dropping near-zero variance: ", paste(drop_nzv, collapse = ", "))
  smp <- smp[, !names(smp) %in% drop_nzv, drop = FALSE]
}

# --- linear combinations ---
if (ncol(smp) >= 2) {
  lc <- caret::findLinearCombos(smp)
  if (!is.null(lc$remove) && length(lc$remove) > 0) {
    drop_lc <- colnames(smp)[lc$remove]
    msg("Dropping linear combos: ", paste(drop_lc, collapse = ", "))
    smp <- smp[, -lc$remove, drop = FALSE]
  }
}

# --- high correlation ---
if (ncol(smp) >= 2) {
  cor_mat <- stats::cor(smp, use = "pairwise.complete.obs")
  drop_cor_idx <- caret::findCorrelation(
    cor_mat,
    cutoff = corr_cutoff,
    names  = FALSE,
    exact  = TRUE
  )
  
  if (length(drop_cor_idx) > 0) {
    drop_cor <- colnames(smp)[drop_cor_idx]
    msg("Dropping correlated predictors: ", paste(drop_cor, collapse = ", "))
    smp <- smp[, -drop_cor_idx, drop = FALSE]
  }
}

keep_pred <- names(smp)
msg("Predictors after reduction: ", length(keep_pred))

# rebuild training table (keep metadata + reduced predictors)
keep_cols <- c("x", "y", "class_id", "code", "class_name", "class", keep_pred)
train_df  <- train_df[, keep_cols, drop = FALSE]

# mirror reduction to predictor stack (canonical grid/CRS)
pred_red <- pred[[keep_pred]]

msg("Writing reduced predictor stack (canonical): ", out_stack_file)
dir.create(dirname(out_stack_file), recursive = TRUE, showWarnings = FALSE)
terra::writeRaster(pred_red, out_stack_file, overwrite = TRUE)

# --------------------------------------------------------------------
# Convert to sf + write
# --------------------------------------------------------------------
train_sf <- sf::st_as_sf(
  train_df,
  coords = c("x", "y"),
  crs    = sf::st_crs(aoi_m),
  remove = FALSE
)

msg("Writing training dataset: ", out_train)
write_training(train_sf, out_train)

msg("Done.")
