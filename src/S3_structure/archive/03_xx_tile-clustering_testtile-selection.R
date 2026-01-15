#!/usr/bin/env Rscript
############################################################
# 03_xx_tile-clustering_testtile-selection.R
#
# Goal
# ----
# 1. Build a coarse grid ("tiles") over the AOI.
# 2. Compute a small, robust feature vector per tile:
#      - Information-theoretic metrics from land cover:
#          * ent (entropy H(x))
#          * relmutinf (relative mutual information U)
#      - Simple physiographic metrics:
#          * elev_mean
#          * forest_fraction
# 3. Cluster tiles in this feature space → "landscape types".
# 4. Select 1–N representative tiles per type as *test tiles*
#    for segmentation / supercell / meanshift calibration.
#
# This is deliberately lean:
#   - Kachel-Clustering (not pixel-based)
#   - No segmentation here
#   - No optimisation yet
############################################################

# ---------------------------------------------------------
# 0) Packages
# ---------------------------------------------------------
library(here)
library(sf)
library(terra)
library(landscapemetrics)
library(exactextractr)
library(dplyr)
library(tidyr)
library(tibble)

message("Project root: ", here::here())

# ---------------------------------------------------------
# 1) Paths & input data
# ---------------------------------------------------------
# Adjust paths if needed – this assumes your existing Burgwald setup.

aoi_file  <- here::here("data", "processed", "aoi_burgwald.gpkg")
clc_file  <- here::here("data", "raw", "AOI_Burgwald", "clc", "clc5_2018_burgwald.tif")
dem_file  <- here::here("data", "raw", "AOI_Burgwald", "dem", "dem_dgm1_burgwald.tif")

# Read AOI, land cover, DEM
aoi_burgwald <- sf::read_sf(aoi_file)
lc_burgwald  <- terra::rast(clc_file)
dem_burgwald <- terra::rast(dem_file)

# ---------------------------------------------------------
# 2) Harmonise CRS and clip to AOI
# ---------------------------------------------------------
# Use DEM CRS as metric reference (e.g. EPSG:25832).

crs_dem  <- terra::crs(dem_burgwald)
aoi_utm  <- sf::st_transform(aoi_burgwald, crs_dem)

# Clip DEM to AOI
dem_burgwald <- dem_burgwald |>
  terra::crop(terra::vect(aoi_utm)) |>
  terra::mask(terra::vect(aoi_utm))

# Reproject land cover to DEM grid (nearest neighbour, categorical)
lc_burgwald_utm <- terra::project(lc_burgwald, dem_burgwald, method = "near")
plot(lc_burgwald_utm)

# ---------------------------------------------------------
# 3) Build tile grid (Kacheln) over AOI
# ---------------------------------------------------------
# IMPORTANT: Tile size controls scale of "landscape type".
# Example: 500 m → 0.5 km tiles

tile_size <- 500  # [m]

tiles_utm <- sf::st_make_grid(
  aoi_utm,
  cellsize = tile_size,
  what     = "polygons"
) |>
  sf::st_as_sf() |>
  dplyr::mutate(tile_id = dplyr::row_number())

# For landscapemetrics::sample_lsm(), tiles must be in raster CRS
tiles_lc_crs <- sf::st_transform(tiles_utm, sf::st_crs(lc_burgwald_utm))

# ---------------------------------------------------------
# 4) Information-theoretic metrics per tile (H, U)
# ---------------------------------------------------------
# Each tile is treated as a small "landscape".
# sample_lsm() computes landscape-level metrics:
#   - lsm_l_ent       -> ent
#   - lsm_l_relmutinf -> relmutinf

lsm_tiles <- landscapemetrics::sample_lsm(
  landscape = lc_burgwald_utm,
  y         = tiles_lc_crs,
  what      = c("lsm_l_ent", "lsm_l_relmutinf"),
  level     = "landscape"
)

# Wide table: tile_id | ent | relmutinf
it_tiles_wide <- lsm_tiles |>
  dplyr::select(tile_id = plot_id, metric, value) |>
  tidyr::pivot_wider(names_from = metric, values_from = value)

# ---------------------------------------------------------
# 5) Simple physiographic metrics per tile
# ---------------------------------------------------------
# Here: elev_mean and forest_fraction (CLC codes 23,24,25 as forest).
# Use exactextractr for speed / robustness on large rasters.

# 5.1 Elevation mean
elev_vals <- exactextractr::exact_extract(
  dem_burgwald,
  tiles_utm,
  "mean",
  progress = TRUE
)
elev_vals <- as.numeric(elev_vals)

# 5.2 Forest fraction from land cover
forest_binary <- terra::ifel(lc_burgwald_utm %in% c(23, 24, 25), 1, 0)

forest_vals <- exactextractr::exact_extract(
  forest_binary,
  tiles_utm,
  "mean",        # mean of binary → proportion of forest
  progress = TRUE
)
forest_vals <- as.numeric(forest_vals)

physio_tiles <- tibble::tibble(
  tile_id         = tiles_utm$tile_id,
  elev_mean       = elev_vals,
  forest_fraction = forest_vals
)

# ---------------------------------------------------------
# 6) Combine features into one tile feature table
# ---------------------------------------------------------

tile_features <- tiles_utm |>
  dplyr::left_join(it_tiles_wide,   by = "tile_id") |>
  dplyr::left_join(physio_tiles,    by = "tile_id")

# Optional: sanity check
message("Tile feature summary (no geometry):")
print(
  tile_features |>
    sf::st_drop_geometry() |>
    dplyr::select(tile_id, ent, relmutinf, elev_mean, forest_fraction) |>
    summary()
)

# ---------------------------------------------------------
# 7) Prepare feature matrix for clustering
# ---------------------------------------------------------
# You decide which features to use:
#   - ent, relmutinf     → pattern / heterogeneity
#   - elev_mean          → orographic control
#   - forest_fraction    → structural / interception control
#
# All are scaled (z-scores) to avoid dominance by units.

clust_df <- tile_features |>
  sf::st_drop_geometry() |>
  dplyr::select(
    tile_id,
    ent,
    relmutinf,
    elev_mean,
    forest_fraction
  ) |>
  dplyr::filter(
    !is.na(ent),
    !is.na(relmutinf),
    !is.na(elev_mean),
    !is.na(forest_fraction)
  )

clust_scaled <- clust_df |>
  dplyr::mutate(
    ent_z       = as.numeric(scale(ent)),
    relmutinf_z = as.numeric(scale(relmutinf)),
    elev_z      = as.numeric(scale(elev_mean)),
    forest_z    = as.numeric(scale(forest_fraction))
  )

# ---------------------------------------------------------
# 8) Kachel-Clustering → "Landscape types"
# ---------------------------------------------------------
# k_tiles = number of tile-clusters (= landscape types).
# Typical range: 4–8. Adjust based on AOI size / heterogeneity.

k_tiles <- 8

set.seed(123)
km_tiles <- stats::kmeans(
  clust_scaled[, c("ent_z", "relmutinf_z", "elev_z", "forest_z")],
  centers = k_tiles,
  nstart  = 50
)

clust_scaled$tile_cluster_id <- km_tiles$cluster

# Attach cluster IDs to tile geometries
# starting point: tiles_sf has tile_id + geometry
# clust_scaled has tile_id, ent, relmutinf, elev_mean, forest_fraction, tile_cluster_id, ...

tiles_cluster_sf <- tiles_utm |>
  dplyr::left_join(
    clust_scaled |>
      dplyr::select(
        tile_id,
        tile_cluster_id,
        ent,
        relmutinf,
        elev_mean,
        forest_fraction
      ),
    by = "tile_id"
  )


# ---------------------------------------------------------
# 9) Select representative test tiles per cluster
# ---------------------------------------------------------
# Strategy:
#   - Compute cluster centroids in z-space.
#   - For each tile: distance to its cluster centroid.
#   - Select n_rep_per_cluster tiles with minimal distance
#     → "typical" tiles for that landscape type.

n_rep_per_cluster <- 5  # 1–2 is usually enough

# 9.1 Cluster centroids in z-space
cluster_centroids <- clust_scaled |>
  dplyr::group_by(tile_cluster_id) |>
  dplyr::summarise(
    ent_z_mean       = mean(ent_z,       na.rm = TRUE),
    relmutinf_z_mean = mean(relmutinf_z, na.rm = TRUE),
    elev_z_mean      = mean(elev_z,      na.rm = TRUE),
    forest_z_mean    = mean(forest_z,    na.rm = TRUE),
    .groups = "drop"
  )

clust_scaled <- clust_scaled |>
  dplyr::left_join(cluster_centroids, by = "tile_cluster_id") |>
  dplyr::mutate(
    dist_to_center = sqrt(
      (ent_z       - ent_z_mean)^2 +
        (relmutinf_z - relmutinf_z_mean)^2 +
        (elev_z      - elev_z_mean)^2 +
        (forest_z    - forest_z_mean)^2
    )
  )

# 9.3 Representative tiles per cluster
rep_tiles <- clust_scaled |>
  dplyr::group_by(tile_cluster_id) |>
  dplyr::slice_min(
    order_by  = dist_to_center,
    n         = n_rep_per_cluster,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

# 9.4 Attach geometries
rep_tiles_sf <- tiles_cluster_sf |>
  dplyr::right_join(
    rep_tiles |>
      dplyr::select(tile_id),
    by = "tile_id"
  )

# ---------------------------------------------------------
# 10) Outputs & interpretation
# ---------------------------------------------------------

message("Tile clusters (first rows, no geometry):")
print(
  tiles_cluster_sf |>
    sf::st_drop_geometry() |>
    dplyr::select(
      tile_id,
      ent,
      relmutinf,
      elev_mean,
      forest_fraction,
      tile_cluster_id
    ) |>
    head()
)

message("\nRepresentative test tiles per cluster:")
print(
  rep_tiles_sf |>
    sf::st_drop_geometry() |>
    dplyr::select(
      tile_id,
      tile_cluster_id,
      ent,
      relmutinf,
      elev_mean,
      forest_fraction
    )
)

# Optional interactive checks (comment out if run in batch)
library(mapview)
mapview::mapview(tiles_cluster_sf, zcol = "tile_cluster_id") +
  mapview::mapview(rep_tiles_sf, col.regions = "black")

############################################################
# Objects of interest:
#
#   tiles_cluster_sf
#     - sf grid of tiles with:
#         * ent, relmutinf
#         * elev_mean, forest_fraction
#         * tile_cluster_id  (landscape-type ID)
#
#   rep_tiles_sf
#     - subset of tiles_cluster_sf:
#         * 1–n tiles per cluster
#         * geometries = test tiles for segmentation calibration
#
# Typische Nutzung:
#   - Für jede tile_cluster_id:
#       * rep_tiles_sf %>% filter(tile_cluster_id == k)
#       * Kachel-Extent holen, Sentinel-Stack / CHM croppen
#       * MeanShift / Supercells / SLIC etc. mit
#         unterschiedlichen Parametern testen.
#
#   - Parameterwahl dann nicht anhand eines einzigen Ausschnitts,
#     sondern robust über alle Landschaftstypen.
############################################################
