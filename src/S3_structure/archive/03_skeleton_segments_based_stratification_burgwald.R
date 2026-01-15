#!/usr/bin/env Rscript
############################################################
# 03_xx_segments_based_stratification_burgwald.R
#
# Assumption:
#   You already have a segmentation of the Burgwald region:
#     segments_sf (sf polygons, one row per segment),
#     with a unique segment_id.
#
# Goal:
#   1) Compute IT metrics (entropy, relative mutual information)
#      per segment.
#   2) Compute physiographic metrics per segment:
#        - elevation
#        - slope
#        - southness
#        - forest_fraction
#        - windwardness (Summer / Winter)
#   3) Derive:
#        - IT strata  (pattern-based)
#        - Physio strata (process-based)
#   4) Select representative segments as station candidates
#      for both strata systems, using segment centroids.
#
# Layers:
#   - Base units: segmentation (segments_sf).
#   - Layer 1: IT-based stratification.
#   - Layer 2: Physiographic stratification.
#   - Layer 3 (not in this script): info-gain & hydrological
#     coupling on candidate sets.
############################################################

## ---------------------------------------------------------
## 0) Packages
## ---------------------------------------------------------
library(here)
library(terra)
library(sf)
library(landscapemetrics)
library(exactextractr)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

message("Project root: ", here::here())

## ---------------------------------------------------------
## 0.1) Inputs and assumptions
## ---------------------------------------------------------
# You must provide:
#   - segments_sf: sf polygons with column 'segment_id'
#   - dem_burgwald: DEM raster
#   - lc_burgwald: land-cover raster (CLC-like)
#   - wind_means_summary: data.frame with columns:
#       period ∈ { "Summer", "Winter" }
#       mean_wd (mean wind direction in degrees 0–360)

# Example: read segments from GPKG (adapt path/name)
segments_file <- here::here("data", "processed", "burgwald_segments_meanshift.gpkg")
segments_sf   <- sf::read_sf(segments_file)

# DEM and land cover (adapt to your project paths)
aoi_root  <- here::here("data", "raw", "AOI_Burgwald")
dem_file  <- file.path(aoi_root, "dem", "dem_dgm1_burgwald.tif")
clc_file  <- file.path(aoi_root, "clc", "clc5_2018_burgwald.tif")

dem_burgwald <- terra::rast(dem_file)
lc_burgwald  <- terra::rast(clc_file)

# Wind means (must be created beforehand by your DWD wind script)
#   wind_means_summary:
#     period  = "Summer" / "Winter"
#     mean_wd = mean wind direction in degrees
if (!exists("wind_means_summary")) {
  stop("wind_means_summary not found. Run wind analysis first.")
}
if (!all(c("period", "mean_wd") %in% names(wind_means_summary))) {
  stop("wind_means_summary must contain 'period' and 'mean_wd'.")
}
wind_dir_summer <- wind_means_summary$mean_wd[wind_means_summary$period == "Summer"]
wind_dir_winter <- wind_means_summary$mean_wd[wind_means_summary$period == "Winter"]
if (length(wind_dir_summer) != 1L || length(wind_dir_winter) != 1L) {
  stop("wind_means_summary must have exactly one 'Summer' and one 'Winter' row.")
}

## ---------------------------------------------------------
## 1) Harmonise CRS
## ---------------------------------------------------------
# We adopt DEM CRS as reference for all raster/segment operations.

crs_dem <- terra::crs(dem_burgwald)

# Reproject segments to DEM CRS (if needed)
if (!identical(sf::st_crs(segments_sf)$wkt, crs_dem)) {
  segments_sf <- sf::st_transform(segments_sf, crs_dem)
}

# Reproject land cover to DEM grid (nearest neighbour, categorical)
lc_burgwald_utm <- terra::project(lc_burgwald, dem_burgwald, method = "near")
terra::values(lc_burgwald_utm) <- round(terra::values(lc_burgwald_utm))

## ---------------------------------------------------------
## 2) IT metrics per segment (entropy + rel. mutual information)
## ---------------------------------------------------------
# Idea:
#   Each segment is treated as a "mini-landscape".
#   sample_lsm() computes IT metrics on these segments.
#   Output is linked by 'plot_id' which we map to segment_id.

# Ensure segments_sf has 'segment_id'
if (!"segment_id" %in% names(segments_sf)) {
  segments_sf$segment_id <- dplyr::row_number()
}

# Landscapemetrics: sample IT metrics on land cover
lsm_seg_it <- sample_lsm(
  landscape = lc_burgwald_utm,
  y         = segments_sf,
  what      = c("lsm_l_ent", "lsm_l_relmutinf"),
  level     = "landscape"
)

# Wide table: segment_id | ent | relmutinf
it_seg_wide <- lsm_seg_it %>%
  dplyr::select(segment_id = plot_id, metric, value) %>%
  tidyr::pivot_wider(names_from = metric, values_from = value)

# Attach to segments
segments_it_sf <- segments_sf %>%
  dplyr::left_join(it_seg_wide, by = "segment_id")

## ---------------------------------------------------------
## 3) Physiographic metrics per segment
## ---------------------------------------------------------
# Compute:
#   - elev_mean
#   - slope_mean_deg
#   - southness_mean
#   - forest_fraction
#   - windward_summer
#   - windward_winter

# 3.1 Terrain derivatives: slope, aspect, southness
slope_deg  <- terra::terrain(dem_burgwald, v = "slope",  unit = "degrees")
aspect_deg <- terra::terrain(dem_burgwald, v = "aspect", unit = "degrees")

aspect_rad  <- aspect_deg * pi / 180
southness   <- cos(aspect_rad - pi)         # cos(aspect - 180°)
southness01 <- (southness + 1) / 2          # scale to [0,1]

# 3.2 Forest binary from land cover
#     Adapt class codes to your CLC legend
forest_binary <- terra::ifel(lc_burgwald_utm %in% c(23, 24, 25), 1, 0)

# 3.3 Windwardness rasters (Summer / Winter)
#     cos((aspect - mean_wind_dir) * π/180):
#       +1 → fully windward (Luv)
#       0  → side
#       -1 → leeward (Lee)

windward_summer_r <- cos((aspect_deg - wind_dir_summer) * pi / 180)
windward_winter_r <- cos((aspect_deg - wind_dir_winter) * pi / 180)
names(windward_summer_r) <- "windward_summer"
names(windward_winter_r) <- "windward_winter"

# 3.4 Exact extraction per segment
# Note:
#   exact_extract(…, polygons) returns a numeric vector if fun="mean".
#   We coerce to numeric to avoid data.frame quirks.

elev_vals <- exactextractr::exact_extract(dem_burgwald,       segments_it_sf, "mean")
slope_vals <- exactextractr::exact_extract(slope_deg,         segments_it_sf, "mean")
south_vals <- exactextractr::exact_extract(southness01,       segments_it_sf, "mean")
forest_vals <- exactextractr::exact_extract(forest_binary,    segments_it_sf, "mean")
wind_summer_vals <- exactextractr::exact_extract(windward_summer_r, segments_it_sf, "mean")
wind_winter_vals <- exactextractr::exact_extract(windward_winter_r, segments_it_sf, "mean")

# Assemble physio table
physio_seg_df <- tibble::tibble(
  segment_id        = segments_it_sf$segment_id,
  elev_mean         = as.numeric(elev_vals),
  slope_mean_deg    = as.numeric(slope_vals),
  southness_mean    = as.numeric(south_vals),
  forest_fraction   = as.numeric(forest_vals),
  windward_summer   = as.numeric(wind_summer_vals),
  windward_winter   = as.numeric(wind_winter_vals)
)

# Attach physio metrics
segments_full_sf <- segments_it_sf %>%
  dplyr::left_join(physio_seg_df, by = "segment_id")

## ---------------------------------------------------------
## 4) IT-based strata (pattern stratification)
## ---------------------------------------------------------
# Feature space: ent, relmutinf (+ optional forest_fraction).
# No DEM-derivatives here.

it_clust_df <- segments_full_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(segment_id, ent, relmutinf) %>%
  dplyr::filter(!is.na(ent), !is.na(relmutinf))

it_scaled <- it_clust_df %>%
  dplyr::mutate(
    ent_z       = as.numeric(scale(ent)),
    relmutinf_z = as.numeric(scale(relmutinf))
  )

k_it <- 6  # number of IT strata (tune as needed)

set.seed(123)
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
## 5) Physio-based strata (process stratification)
## ---------------------------------------------------------
# Core physio features:
#   elev_mean, slope_mean_deg, southness_mean, forest_fraction
# Windwardness is kept as descriptive attribute but not in clustering.

physio_clust_df <- segments_full_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(
    segment_id,
    elev_mean,
    slope_mean_deg,
    southness_mean,
    forest_fraction
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

k_phys <- 6  # number of physio strata (tune as needed)

set.seed(123)
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
## 6) Representative segments as station candidates
## ---------------------------------------------------------
# For each strata system, select segments closest to the centroid
# in the relevant feature space, then use polygon centroids as
# station candidate locations.
#
# IT strata:
#   - based on ent_z, relmutinf_z
# Physio strata:
#   - based on elev_z, slope_z, southness_z, forestfrac_z

### 6.1 IT strata candidates --------------------------------

# Stratum centroids (IT space)
it_centroids <- it_scaled %>%
  dplyr::group_by(it_stratum_id) %>%
  dplyr::summarise(
    ent_z_mean       = mean(ent_z,       na.rm = TRUE),
    relmutinf_z_mean = mean(relmutinf_z, na.rm = TRUE),
    .groups = "drop"
  )

# Distance of each segment to IT centroid
it_scaled <- it_scaled %>%
  dplyr::left_join(it_centroids, by = "it_stratum_id") %>%
  dplyr::mutate(
    dist_to_it_center = sqrt(
      (ent_z       - ent_z_mean)^2 +
        (relmutinf_z - relmutinf_z_mean)^2
    )
  )

n_per_it_stratum <- 3  # number of candidates per IT stratum

rep_it_segments <- it_scaled %>%
  dplyr::group_by(it_stratum_id) %>%
  dplyr::slice_min(
    order_by  = dist_to_it_center,
    n         = n_per_it_stratum,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

rep_it_sf <- segments_full_sf %>%
  dplyr::inner_join(
    rep_it_segments %>% dplyr::select(segment_id, it_stratum_id),
    by = c("segment_id", "it_stratum_id")
  )

station_it_candidates_sf <- rep_it_sf %>%
  sf::st_centroid() %>%
  dplyr::mutate(
    station_it_id = dplyr::row_number()
  ) %>%
  dplyr::select(
    station_it_id,
    segment_id,
    it_stratum_id,
    ent,
    relmutinf
    # geometry kept implicitly
  )

### 6.2 Physio strata candidates ---------------------------

# Merge physio_scaled with physio_stratum_id and z-means
physio_centroids <- physio_scaled %>%
  dplyr::group_by(physio_stratum_id) %>%
  dplyr::summarise(
    elev_z_mean       = mean(elev_z,       na.rm = TRUE),
    slope_z_mean      = mean(slope_z,      na.rm = TRUE),
    southness_z_mean  = mean(southness_z,  na.rm = TRUE),
    forestfrac_z_mean = mean(forestfrac_z, na.rm = TRUE),
    .groups = "drop"
  )

physio_scaled <- physio_scaled %>%
  dplyr::left_join(physio_centroids, by = "physio_stratum_id") %>%
  dplyr::mutate(
    dist_to_physio_center = sqrt(
      (elev_z       - elev_z_mean)^2 +
        (slope_z      - slope_z_mean)^2 +
        (southness_z  - southness_z_mean)^2 +
        (forestfrac_z - forestfrac_z_mean)^2
    )
  )

n_per_physio_stratum <- 3  # number of candidates per physio stratum

rep_physio_segments <- physio_scaled %>%
  dplyr::group_by(physio_stratum_id) %>%
  dplyr::slice_min(
    order_by  = dist_to_physio_center,
    n         = n_per_physio_stratum,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

rep_physio_sf <- segments_full_sf %>%
  dplyr::inner_join(
    rep_physio_segments %>%
      dplyr::select(segment_id, physio_stratum_id),
    by = c("segment_id", "physio_stratum_id")
  )

station_physio_candidates_sf <- rep_physio_sf %>%
  sf::st_centroid() %>%
  dplyr::mutate(
    station_physio_id = dplyr::row_number()
  ) %>%
  dplyr::select(
    station_physio_id,
    segment_id,
    physio_stratum_id,
    elev_mean,
    slope_mean_deg,
    southness_mean,
    forest_fraction,
    windward_summer,
    windward_winter
    # geometry kept implicitly
  )

## ---------------------------------------------------------
## 7) Object overview
## ---------------------------------------------------------

message("Segments with strata (first rows, no geometry):")
print(
  segments_full_sf %>%
    sf::st_drop_geometry() %>%
    dplyr::select(
      segment_id,
      ent,
      relmutinf,
      elev_mean,
      slope_mean_deg,
      southness_mean,
      forest_fraction,
      windward_summer,
      windward_winter,
      it_stratum_id,
      physio_stratum_id
    ) %>%
    head()
)

message("\nIT-based station candidates:")
print(
  station_it_candidates_sf %>%
    sf::st_drop_geometry()
)

message("\nPhysio-based station candidates:")
print(
  station_physio_candidates_sf %>%
    sf::st_drop_geometry()
)

############################################################
# Summary:
#   - segments_full_sf:
#       base segmentation with IT + physio + windwardness +
#       it_stratum_id + physio_stratum_id.
#   - station_it_candidates_sf:
#       representative segments (centroids) per IT stratum.
#   - station_physio_candidates_sf:
#       representative segments (centroids) per physio stratum.
#
# Next step:
#   - Combine both candidate sets.
#   - Apply information-gain optimisation (e.g. kriging variance,
#     radar–gauge residuals) and hydrological coupling (P–Q coherence)
#     on this joint candidate pool.
############################################################
