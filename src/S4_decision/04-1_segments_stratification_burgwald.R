#!/usr/bin/env Rscript
############################################################
# 04-1_segments_stratification_burgwald.R  (S4)
#
# Purpose (Layer-conform):
#   Consume artifacts from other layers and produce segment-based
#   strata + station-candidate sets.
#
# Inputs (artifacts):
#   (S2) Segments (GPKG): sf polygons with unique segment_id
#   (S1/AOI) DEM (GeoTIFF)
#   (S1/AOI) Landcover (GeoTIFF, categorical; CLC-like)
#   (Wind layer) wind_means_summary (RDS):
#       columns: period ∈ {"Summer","Winter"}, mean_wd (0..360 deg)
#
# Outputs (artifacts):
#   - segments_full_sf.gpkg:
#       segments with IT + physio + windwardness + strata ids
#   - station_it_candidates.gpkg:
#       representative segment centroids per IT stratum
#   - station_physio_candidates.gpkg:
#       representative segment centroids per Physio stratum
#
# Notes:
#   - This keeps IT + Physio in one S4 step for pragmatism,
#     but consumes only upstream artifacts (no implicit objects).
#   - Windwardness is computed as cos(aspect - mean_wind_dir).
############################################################

suppressPackageStartupMessages({
  library(here)
  library(terra)
  library(sf)
  library(landscapemetrics)
  library(exactextractr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})


source(here::here("src","_core", "00-setup-burgwald.R"))
## ---------------------------------------------------------
## 0) Config (edit here once; keep paths stable)
## ---------------------------------------------------------

# Where to read the segment polygons produced by S2 (OTB MeanShift)
segments_file <- here::here(
  "data", "processed", "layer0_segments", "segments",
  "segments_meanshift_y2018.gpkg"
)

# AOI inputs (your existing convention)
aoi_root <- here::here("data", "raw", "AOI_Burgwald")

dem_file <- file.path(aoi_root, "dem", "dem_dgm1_burgwald.tif")
clc_file <- file.path(aoi_root, "clc", "clc5_2018_burgwald.tif")

# Wind layer artifact (produced by your DWD wind script)
wind_means_file <- here::here(
  "data", "processed", "dwd_wind",
  "wind_means_summary.rds"
)

# Output directory (S4 products)
out_dir <- here::here("data", "processed", "layer4_segments_strata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_segments_full <- file.path(out_dir, "segments_full_sf.gpkg")
out_station_it    <- file.path(out_dir, "station_it_candidates.gpkg")
out_station_phys  <- file.path(out_dir, "station_physio_candidates.gpkg")

# Class codes for "forest" in your LC raster (adapt to your legend)
forest_codes <- c(23, 24, 25)

# Strata controls
k_it   <- 6
k_phys <- 6

# Number of candidates per stratum
n_per_it_stratum    <- 3
n_per_physio_stratum <- 3

set.seed(123)

message("Project root: ", here::here())
message("Segments:      ", segments_file)
message("DEM:           ", dem_file)
message("Landcover:     ", clc_file)
message("Wind means:    ", wind_means_file)
message("Output dir:    ", out_dir)

## ---------------------------------------------------------
## 1) Read artifacts
## ---------------------------------------------------------
if (!file.exists(segments_file)) stop("Missing segments_file: ", segments_file)
if (!file.exists(dem_file))     stop("Missing dem_file: ", dem_file)
if (!file.exists(clc_file))     stop("Missing clc_file: ", clc_file)
if (!file.exists(wind_means_file)) stop("Missing wind_means_file: ", wind_means_file)

segments_sf <- sf::read_sf(segments_file)
dem_burgwald <- terra::rast(dem_file)
lc_burgwald  <- terra::rast(clc_file)

wind_means_summary <- readRDS(wind_means_file)

## Validate wind means
req_cols <- c("period", "mean_wd")
if (!all(req_cols %in% names(wind_means_summary))) {
  stop("wind_means_summary must contain: ", paste(req_cols, collapse = ", "))
}

wind_dir_summer <- wind_means_summary$mean_wd[wind_means_summary$period == "Summer"]
wind_dir_winter <- wind_means_summary$mean_wd[wind_means_summary$period == "Winter"]

if (length(wind_dir_summer) != 1L || length(wind_dir_winter) != 1L) {
  stop("wind_means_summary must have exactly one row for Summer and Winter.")
}

## Ensure segment_id exists
if (!"segment_id" %in% names(segments_sf)) {
  warning("segments_sf has no 'segment_id' → creating from row_number().")
  segments_sf$segment_id <- dplyr::row_number()
}

## ---------------------------------------------------------
## 2) Harmonise CRS / grids
## ---------------------------------------------------------
crs_dem <- terra::crs(dem_burgwald)

# Reproject segments to DEM CRS if needed
if (!identical(sf::st_crs(segments_sf)$wkt, crs_dem)) {
  segments_sf <- sf::st_transform(segments_sf, crs_dem)
}

# Reproject landcover to DEM grid (categorical)
lc_burgwald_utm <- terra::project(lc_burgwald, dem_burgwald, method = "near")
terra::values(lc_burgwald_utm) <- round(terra::values(lc_burgwald_utm))

## ---------------------------------------------------------
## 3) IT metrics per segment (entropy + relmutinf)
## ---------------------------------------------------------
# Each segment is treated as a mini-landscape.
lsm_seg_it <- landscapemetrics::sample_lsm(
  landscape = lc_burgwald_utm,
  y         = segments_sf,
  what      = c("lsm_l_ent", "lsm_l_relmutinf"),
  level     = "landscape"
)

it_seg_wide <- lsm_seg_it %>%
  dplyr::select(segment_id = plot_id, metric, value) %>%
  tidyr::pivot_wider(names_from = metric, values_from = value) %>%
  dplyr::rename(
    ent = lsm_l_ent,
    relmutinf = lsm_l_relmutinf
  )

segments_it_sf <- segments_sf %>%
  dplyr::left_join(it_seg_wide, by = "segment_id")

## ---------------------------------------------------------
## 4) Physiographic metrics per segment + windwardness
## ---------------------------------------------------------
# Terrain derivatives
slope_deg  <- terra::terrain(dem_burgwald, v = "slope",  unit = "degrees")
aspect_deg <- terra::terrain(dem_burgwald, v = "aspect", unit = "degrees")

# Southness scaled to [0,1]
aspect_rad  <- aspect_deg * pi / 180
southness   <- cos(aspect_rad - pi)     # cos(aspect - 180°)
southness01 <- (southness + 1) / 2

# Forest fraction from LC
forest_binary <- terra::ifel(lc_burgwald_utm %in% forest_codes, 1, 0)

# Windwardness rasters (Summer/Winter)
windward_summer_r <- cos((aspect_deg - wind_dir_summer) * pi / 180)
windward_winter_r <- cos((aspect_deg - wind_dir_winter) * pi / 180)
names(windward_summer_r) <- "windward_summer"
names(windward_winter_r) <- "windward_winter"

# Segment-wise extraction (means)
elev_vals        <- exactextractr::exact_extract(dem_burgwald,         segments_it_sf, "mean")
slope_vals       <- exactextractr::exact_extract(slope_deg,           segments_it_sf, "mean")
south_vals       <- exactextractr::exact_extract(southness01,         segments_it_sf, "mean")
forest_vals      <- exactextractr::exact_extract(forest_binary,       segments_it_sf, "mean")
wind_summer_vals <- exactextractr::exact_extract(windward_summer_r,   segments_it_sf, "mean")
wind_winter_vals <- exactextractr::exact_extract(windward_winter_r,   segments_it_sf, "mean")

physio_seg_df <- tibble::tibble(
  segment_id      = segments_it_sf$segment_id,
  elev_mean       = as.numeric(elev_vals),
  slope_mean_deg  = as.numeric(slope_vals),
  southness_mean  = as.numeric(south_vals),
  forest_fraction = as.numeric(forest_vals),
  windward_summer = as.numeric(wind_summer_vals),
  windward_winter = as.numeric(wind_winter_vals)
)

segments_full_sf <- segments_it_sf %>%
  dplyr::left_join(physio_seg_df, by = "segment_id")

## ---------------------------------------------------------
## 5) IT strata (pattern-based)
## ---------------------------------------------------------
it_clust_df <- segments_full_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(segment_id, ent, relmutinf) %>%
  dplyr::filter(!is.na(ent), !is.na(relmutinf))

it_scaled <- it_clust_df %>%
  dplyr::mutate(
    ent_z       = as.numeric(scale(ent)),
    relmutinf_z = as.numeric(scale(relmutinf))
  )

km_it <- stats::kmeans(
  it_scaled[, c("ent_z", "relmutinf_z")],
  centers = k_it,
  nstart  = 50
)

it_scaled$it_stratum_id <- km_it$cluster

segments_full_sf <- segments_full_sf %>%
  dplyr::left_join(
    it_scaled %>% dplyr::select(segment_id, it_stratum_id),
    by = "segment_id"
  )

## ---------------------------------------------------------
## 6) Physio strata (process-based)
## ---------------------------------------------------------
# Core physio features for clustering:
# elev_mean, slope_mean_deg, southness_mean, forest_fraction
# Windwardness stays as descriptive attribute (not in clustering).
physio_clust_df <- segments_full_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(
    segment_id,
    elev_mean, slope_mean_deg, southness_mean, forest_fraction
  ) %>%
  dplyr::filter(
    !is.na(elev_mean),
    !is.na(slope_mean_deg),
    !is.na(southness_mean),
    !is.na(forest_fraction)
  )

physio_scaled <- physio_clust_df %>%
  dplyr::mutate(
    elev_z       = as.numeric(scale(elev_mean)),
    slope_z      = as.numeric(scale(slope_mean_deg)),
    southness_z  = as.numeric(scale(southness_mean)),
    forestfrac_z = as.numeric(scale(forest_fraction))
  )

km_phys <- stats::kmeans(
  physio_scaled[, c("elev_z", "slope_z", "southness_z", "forestfrac_z")],
  centers = k_phys,
  nstart  = 50
)

physio_scaled$physio_stratum_id <- km_phys$cluster

segments_full_sf <- segments_full_sf %>%
  dplyr::left_join(
    physio_scaled %>% dplyr::select(segment_id, physio_stratum_id),
    by = "segment_id"
  )

## ---------------------------------------------------------
## 7) Representative segments as station candidates
## ---------------------------------------------------------

### 7.1 IT candidates
it_centroids <- it_scaled %>%
  dplyr::group_by(it_stratum_id) %>%
  dplyr::summarise(
    ent_z_mean       = mean(ent_z, na.rm = TRUE),
    relmutinf_z_mean = mean(relmutinf_z, na.rm = TRUE),
    .groups = "drop"
  )

it_scaled2 <- it_scaled %>%
  dplyr::left_join(it_centroids, by = "it_stratum_id") %>%
  dplyr::mutate(
    dist_to_it_center = sqrt((ent_z - ent_z_mean)^2 + (relmutinf_z - relmutinf_z_mean)^2)
  )

rep_it_segments <- it_scaled2 %>%
  dplyr::group_by(it_stratum_id) %>%
  dplyr::slice_min(order_by = dist_to_it_center, n = n_per_it_stratum, with_ties = FALSE)_
