#!/usr/bin/env Rscript
############################################################
# Script:  04-2_signatures_physio_metrics.R
# Project: Burgwald
#
# Purpose:
#   Build S4_signatures physiographic metrics per *stable* segment.
#   - NO clustering
#   - NO candidates
#   - NO decision logic
#
# Inputs (via paths[...] registry):
#   - layer0_segments            (S3_structure, gpkg)
#   - aoi_dgm                    (S1_observation, tif)  OR
#     optionally relief_stack_10m (S2_features, tif) if you want to use it later
#   - wind_means_summary         (S2_features, rds) OPTIONAL (for windwardness)
#
# Output (via paths[...] registry):
#   - layer0_attr_physio_metrics (S4_signatures, gpkg)
#
# Output columns (segment-level):
#   segment_id
#   elev_mean
#   slope_mean_deg
#   aspect_mean_deg
#   southness_mean   (scaled to [0,1], 0=north, 1=south)
#   windward_summer  (optional; mean cos(aspect - mean_wd_summer), [-1..1])
#   windward_winter  (optional; mean cos(aspect - mean_wd_winter), [-1..1])
############################################################

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(tibble)
  library(exactextractr)
})

## ---------------------------------------------------------
## 0) Setup / Registry
## ---------------------------------------------------------

source(here::here("src","_core","01-setup-burgwald.R"))

seg_file <- paths[["layer0_segments"]]
dem_file <- paths[["aoi_dgm"]]
out_file <- paths[["layer0_attr_physio_metrics"]]

message("Input segments: ", seg_file)
message("Input DEM:      ", dem_file)
message("Output:         ", out_file)

## ---------------------------------------------------------
## 1) Read inputs
## ---------------------------------------------------------
segments_sf <- sf::read_sf(seg_file)
if (!"segment_id" %in% names(segments_sf)) {
  stop("Segments file does not contain 'segment_id': ", seg_file)
}

dem <- terra::rast(dem_file)

## ---------------------------------------------------------
## 2) CRS harmonisation
## ---------------------------------------------------------
crs_dem <- terra::crs(dem, proj = TRUE)
if (is.na(crs_dem) || crs_dem == "") stop("DEM has no CRS: ", dem_file)

if (is.na(sf::st_crs(segments_sf))) {
  stop("Segments have no CRS: ", seg_file)
}

if (!identical(sf::st_crs(segments_sf)$wkt, crs_dem)) {
  segments_sf <- sf::st_transform(segments_sf, crs_dem)
}

# Optional: fix invalid geometries if present (keep minimal)
# (Only touch if needed)
if (any(!sf::st_is_valid(segments_sf))) {
  message("Fixing invalid segment geometries with st_make_valid() ...")
  segments_sf <- sf::st_make_valid(segments_sf)
}

## ---------------------------------------------------------
## 3) Terrain derivatives
## ---------------------------------------------------------
# Slope / aspect from DEM
slope_deg  <- terra::terrain(dem, v = "slope",  unit = "degrees")
aspect_deg <- terra::terrain(dem, v = "aspect", unit = "degrees")

# Southness: cos(aspect - 180°) scaled to [0,1]
aspect_rad  <- aspect_deg * pi / 180
southness   <- cos(aspect_rad - pi)     # [-1..1]
southness01 <- (southness + 1) / 2      # [0..1]
names(southness01) <- "southness"

## ---------------------------------------------------------
## 4) Optional: windwardness using wind_means_summary
## ---------------------------------------------------------
has_wind <- ("wind_means_summary" %in% names(paths)) && file.exists(paths[["wind_means_summary"]])

windward_summer_r <- NULL
windward_winter_r <- NULL

if (has_wind) {
  wind_means_summary <- readRDS(paths[["wind_means_summary"]])
  
  if (!all(c("period", "mean_wd") %in% names(wind_means_summary))) {
    message("wind_means_summary exists but lacks columns {period, mean_wd}. Skipping windwardness.")
    has_wind <- FALSE
  } else {
    wd_summer <- wind_means_summary$mean_wd[wind_means_summary$period == "Summer"]
    wd_winter <- wind_means_summary$mean_wd[wind_means_summary$period == "Winter"]
    
    if (length(wd_summer) != 1L || length(wd_winter) != 1L) {
      message("wind_means_summary does not have exactly one Summer and one Winter mean_wd. Skipping windwardness.")
      has_wind <- FALSE
    } else {
      # windwardness: cos((aspect - mean_wd) * pi/180) in [-1..1]
      windward_summer_r <- cos((aspect_deg - wd_summer) * pi / 180)
      windward_winter_r <- cos((aspect_deg - wd_winter) * pi / 180)
      names(windward_summer_r) <- "windward_summer"
      names(windward_winter_r) <- "windward_winter"
      message("Windwardness enabled (Summer/Winter).")
    }
  }
} else {
  message("wind_means_summary not available. Windwardness will be NA.")
}

## ---------------------------------------------------------
## 5) Zonal extraction per segment
## ---------------------------------------------------------
# exact_extract returns numeric vectors aligned with polygon order
elev_mean         <- exactextractr::exact_extract(dem,        segments_sf, "mean")
slope_mean_deg    <- exactextractr::exact_extract(slope_deg,  segments_sf, "mean")
aspect_mean_deg   <- exactextractr::exact_extract(aspect_deg, segments_sf, "mean")
southness_mean    <- exactextractr::exact_extract(southness01,segments_sf, "mean")

windward_summer <- rep(NA_real_, length(elev_mean))
windward_winter <- rep(NA_real_, length(elev_mean))

if (has_wind) {
  windward_summer <- exactextractr::exact_extract(windward_summer_r, segments_sf, "mean")
  windward_winter <- exactextractr::exact_extract(windward_winter_r, segments_sf, "mean")
}

physio_df <- tibble(
  segment_id       = segments_sf$segment_id,
  elev_mean        = as.numeric(elev_mean),
  slope_mean_deg   = as.numeric(slope_mean_deg),
  aspect_mean_deg  = as.numeric(aspect_mean_deg),
  southness_mean   = as.numeric(southness_mean),
  windward_summer  = as.numeric(windward_summer),
  windward_winter  = as.numeric(windward_winter)
)

## ---------------------------------------------------------
## 6) Attach to geometry + write product
## ---------------------------------------------------------
segments_sf <- sf::read_sf(seg_file)

if (!identical(sf::st_crs(segments_sf)$wkt, crs_dem)) {
  segments_sf <- sf::st_transform(segments_sf, crs_dem)
}
# sicherstellen, dass es sf ist
stopifnot(inherits(segments_sf, "sf"))

# Geometriename robust holen
geom_col <- attr(segments_sf, "sf_column")
stopifnot(is.character(geom_col), geom_col %in% names(segments_sf))

out_sf <- segments_sf %>%
  dplyr::select(segment_id, dplyr::all_of(geom_col)) %>%
  dplyr::left_join(physio_df, by = "segment_id") %>%
  sf::st_as_sf()

# Ensure output directory exists
out_dir <- dirname(out_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Overwrite existing layer/file
sf::st_write(out_sf, out_file, delete_dsn = TRUE, quiet = TRUE)

message("Wrote: ", out_file)
message("Columns: ", paste(names(out_sf), collapse = ", "))
