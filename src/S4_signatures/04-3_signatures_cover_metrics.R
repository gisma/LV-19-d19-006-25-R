#!/usr/bin/env Rscript
############################################################
# Script:  04-4_attributes_cover_metrics.R
# Project: Burgwald
#
# Purpose:
#   Build S4_signatures cover metrics per *stable* segment:
#     - per-class area fractions from a landcover/classification raster
#     - forest_fraction (sum over a configurable class set)
#
#   - NO diversity metrics
#   - NO clustering
#   - NO candidates
#   - NO decision logic
#
# Inputs (via paths[...] registry):
#   - layer0_segments            (S3_structure, gpkg)
#   - classification_rf          (S2_features, tif)  preferred
#       fallback: aoi_clc        (S1_observation, tif)
#
# Output (via paths[...] registry):
#   - layer0_attr_cover_metrics  (S4_signatures, gpkg)
#
# Output columns (segment-level):
#   segment_id
#   frac_class_<VALUE>    (one column per class value found)
#   forest_fraction       (sum of selected forest classes)
#
# Notes:
#   - Forest classes must be configured explicitly below.
############################################################

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(tibble)
  library(exactextractr)
  library(purrr)
  library(tidyr)
})

## ---------------------------------------------------------
## 0) Setup / Registry
## ---------------------------------------------------------
if (!exists("paths")) {
  if (file.exists("00-setup-burgwald.R")) {
    source("00-setup-burgwald.R")
  } else {
    stop("Object 'paths' not found and 00-setup-burgwald.R not present in working directory.")
  }
}
stopifnot(is.list(paths))

req_keys <- c("layer0_segments", "layer0_attr_cover_metrics")
missing_keys <- setdiff(req_keys, names(paths))
if (length(missing_keys) > 0) {
  stop("Missing required registry keys in 'paths': ", paste(missing_keys, collapse = ", "))
}

seg_file <- paths[["layer0_segments"]]
out_file <- paths[["layer0_attr_cover_metrics"]]

# Prefer RF classification; fall back to CLC
lc_key <- NULL
if ("classification_rf" %in% names(paths) && file.exists(paths[["classification_rf"]])) {
  lc_key <- "classification_rf"
} else if ("aoi_clc" %in% names(paths) && file.exists(paths[["aoi_clc"]])) {
  lc_key <- "aoi_clc"
} else {
  stop("No landcover raster found. Expected paths[['classification_rf']] or paths[['aoi_clc']].")
}
lc_file <- paths[[lc_key]]

message("Input segments:  ", seg_file)
message("Input landcover: ", lc_file, " (key: ", lc_key, ")")
message("Output:          ", out_file)

## ---------------------------------------------------------
## 1) Read inputs
## ---------------------------------------------------------
segments_sf <- sf::read_sf(seg_file)
if (!"segment_id" %in% names(segments_sf)) {
  stop("Segments file does not contain 'segment_id': ", seg_file)
}
lc <- terra::rast(lc_file)

## ---------------------------------------------------------
## 2) CRS harmonisation
## ---------------------------------------------------------
crs_lc <- terra::crs(lc, proj = TRUE)
if (is.na(crs_lc) || crs_lc == "") stop("Landcover raster has no CRS: ", lc_file)

if (is.na(sf::st_crs(segments_sf))) stop("Segments have no CRS: ", seg_file)

if (!identical(sf::st_crs(segments_sf)$wkt, crs_lc)) {
  segments_sf <- sf::st_transform(segments_sf, crs_lc)
}

if (any(!sf::st_is_valid(segments_sf))) {
  message("Fixing invalid segment geometries with st_make_valid() ...")
  segments_sf <- sf::st_make_valid(segments_sf)
}

# Ensure integer-coded classes
terra::values(lc) <- round(terra::values(lc))

## ---------------------------------------------------------
## 3) Configure forest class set (EDIT THIS)
## ---------------------------------------------------------
# IMPORTANT: Set this to match YOUR legend (RF class IDs or CLC codes).
forest_classes <- c(23, 24, 25)

## ---------------------------------------------------------
## 4) Helper: class fractions per polygon
## ---------------------------------------------------------
compute_class_fracs <- function(values) {
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(tibble(class = integer(0), frac = numeric(0)))
  }
  tab <- table(values)
  tibble(
    class = as.integer(names(tab)),
    frac  = as.numeric(tab) / sum(tab)
  )
}

message("Computing class fractions per segment ...")

frac_list <- exactextractr::exact_extract(
  lc,
  segments_sf,
  fun = function(values, coverage_fraction) {
    compute_class_fracs(values)
  },
  progress = TRUE
)

frac_long_df <- purrr::imap_dfr(frac_list, function(tbl, idx) {
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
  tbl %>% mutate(segment_id = segments_sf$segment_id[idx])
})

if (nrow(frac_long_df) == 0) {
  stop("No class fractions computed. Check overlap and NA coverage in landcover raster.")
}

## ---------------------------------------------------------
## 5) Pivot to wide: frac_class_<value>
## ---------------------------------------------------------
frac_wide <- frac_long_df %>%
  tidyr::pivot_wider(
    id_cols      = segment_id,
    names_from   = class,
    values_from  = frac,
    names_prefix = "frac_class_",
    values_fill  = 0
  )

frac_cols <- grep("^frac_class_", names(frac_wide), value = TRUE)

# forest_fraction: sum of selected forest classes (if present)
forest_cols <- paste0("frac_class_", forest_classes)
forest_cols_present <- intersect(forest_cols, names(frac_wide))

if (length(forest_cols_present) == 0) {
  message("WARNING: none of the configured forest_classes exist in raster values. forest_fraction will be 0.")
  frac_wide$forest_fraction <- 0
} else {
  frac_wide$forest_fraction <- rowSums(frac_wide[, forest_cols_present, drop = FALSE])
}

cover_df <- frac_wide %>%
  dplyr::select(segment_id, dplyr::all_of(frac_cols), forest_fraction)

## ---------------------------------------------------------
## 6) Attach geometry + write product
## ---------------------------------------------------------
out_sf <- segments_sf %>%
  dplyr::select(segment_id, geometry) %>%
  dplyr::left_join(cover_df, by = "segment_id")

out_dir <- dirname(out_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sf::st_write(out_sf, out_file, delete_dsn = TRUE, quiet = TRUE)

message("Wrote: ", out_file)
message("Forest classes configured: ", paste(forest_classes, collapse = ", "))
message("Fraction columns: ", length(frac_cols))
