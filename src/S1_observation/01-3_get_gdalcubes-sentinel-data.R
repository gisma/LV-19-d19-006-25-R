#!/usr/bin/env Rscript

############################################################
# Script: 01-3_get-gdalcubes-sentinel-data.R
# Author:   [Your Name]
# Project:  [Your Project Name]
#
# Purpose:
# --------
#   Retrieve and preprocess Sentinel-2 L2A cloud-optimized scenes
#   for the Burgwald AOI using a complete STAC → gdalcubes workflow.
#
#   This script performs the following tasks:
#     1) Query the AWS Earth-Search STAC API for Sentinel-2 L2A COGs
#        over the Burgwald region (matching the AOI defined in
#        src/00-setup-burgwald.R).
#
#     2) Filter scenes by:
#          - spatial extent (AOI bounding box in EPSG:4326)
#          - temporal window (2018-06 to 2022-09)
#          - metadata (cloud cover < 5%)
#
#     3) Convert the STAC search result into a gdalcubes
#        image collection, including all relevant Sentinel-2 bands
#        and the Scene Classification Layer (SCL).
#
#     4) Define a gdalcubes cube view:
#          - Projection: EPSG:32632 (UTM zone 32N)
#          - 10 m spatial resolution
#          - Monthly temporal aggregation (median)
#          - Bilinear spatial resampling
#
#     5) Build and mask the raster cube
#        (masking SCL classes 3, 8, 9 = cloud shadow / cloud).
#
#     6) Write the time-aggregated data cube as a NetCDF file:
#          data/burgwald_2018_2022_all.nc
#
#     7) Example analysis:
#          - Compute per-pixel kNDVI using apply_pixel()
#          - Aggregate over time via reduce_time("mean(kNDVI)")
#          - Plot the mean kNDVI distribution
#
#
# Important Notes:
#   - All directories and AOI objects (aoi_burgwald_wgs,
#     root_folder, etc.) are created in 00-setup-burgwald.R.
#     This script assumes that the setup has been run first.
#
#   - The DWD meteorological data pipeline is handled
#     separately in:
#         src/01-data-retrieval.R
#     and is intentionally NOT included here.
#     Satellite and climate data remain separate pipelines.
#
#   - This script strictly avoids setwd() and relies entirely
#     on here::here() for reproducible project structures.
#
#   - The STAC → gdalcubes workflow keeps data *lazy* until the
#     final write_ncdf() step. No local download of individual
#     Sentinel-2 scenes occurs outside the cube-building process.
#
############################################################


# sourcing of the setup and specific used funtions
source("src/00-setup-burgwald.R")

out_file <- file.path(root_folder, "data/burgwald_2018_2022_all.nc")

## ----stac_burgwald_gdalcubes----------------------------------------------

# 1) Build a STAC client for the Earth-search API hosted by Element84 on AWS.
#    This catalog exposes Sentinel-2 L2A COGs via a STAC-compliant interface.
s <- rstac::stac("https://earth-search.aws.element84.com/v0")

# 2) Derive the Burgwald bounding box from the AOI polygon aoi_burgwald_wgs.
#    aoi_burgwald_wgs is assumed to be defined in your setup script
#    (e.g. 00-setup-burgwald.R) and uses EPSG:4326 (lon/lat).
bbox_bw <- sf::st_bbox(aoi_burgwald_wgs)

# 3) Search for Sentinel-2 L2A COG scenes over the Burgwald AOI.
#    - collections: "sentinel-s2-l2a-cogs" (cloud-optimized Sentinel-2)
#    - bbox: uses the WGS84 bounding box of the AOI
#    - datetime: restricts to 2018-06 through 2022-09 (focus on summers)
#    - limit: maximum number of STAC items returned
items <- s |>
  rstac::stac_search(
    collections = "sentinel-s2-l2a-cogs",
    bbox        = c(
      bbox_bw["xmin"],
      bbox_bw["ymin"],
      bbox_bw["xmax"],
      bbox_bw["ymax"]
    ),
    datetime    = "2018-06-01/2022-09-01",  # same summer-focused window
    limit       = 600
  ) |>
  rstac::post_request()

# Print basic STAC search response for inspection
items


# 4) Build a gdalcubes STAC image collection from the returned STAC items.
#    This is a *logical* collection: it does not download any data yet, but
#    describes which assets (bands) are available and where they are stored.
#    We explicitly select relevant Sentinel-2 bands plus SCL (scene classification).
s2_collection <- gdalcubes::stac_image_collection(
  items$features,
  asset_names = c(
    "B01","B02","B03","B04","B05","B06",
    "B07","B08","B8A","B09","B11","SCL"
  ),
  # property_filter: only keep scenes with < 5% cloud cover.
  # x is a list of STAC item properties; we access the eo:cloud_cover field.
  property_filter = function(x) {
    x[["eo:cloud_cover"]] < 5
  }
)

# Inspect the resulting image collection (metadata, number of images, etc.)
s2_collection


## ----burgwald_cubeview-----------------------------------------------------

# 5) Define a gdalcubes "cube view" in EPSG:32632 (UTM zone 32N).
#    This step decides:
#    - spatial reference system (srs)
#    - spatial extent (in UTM coordinates)
#    - spatial resolution (dx, dy)
#    - temporal extent (t0, t1)
#    - temporal resolution (dt)
#    - aggregation (how multiple scenes in a time step are reduced)
#    - resampling (spatial interpolation method)

#    a) Transform Burgwald bbox from WGS84 to UTM 32N.
#       We create an sfc from the bbox, reproject it, and compute its new bbox.
sf::st_as_sfc(bbox_bw) |>
  sf::st_transform("EPSG:32632") |>
  sf::st_bbox() -> bbox_utm

#    b) Create the cube_view:
#       - srs: UTM 32N
#       - extent: slightly expanded UTM bounding box (+/- 10 m buffer)
#       - dx, dy: 10 m × 10 m pixels (Sentinel-2 native resolution)
#       - dt: "P1M" → aggregate to monthly time steps
# ISO-8601 Duration (super compact):
# P = period (date), T = time separator
# Date units:  Y = years, M = months, W = weeks, D = days
# Time units:  H = hours, M = minutes, S = seconds
# Examples:
#   P1M      = 1 month
#   P1Y      = 1 year
#   P3W      = 3 weeks
#   P10D     = 10 days
#   PT6H     = 6 hours
#   PT30M    = 30 minutes
#   PT30S    = 30 seconds
#       - aggregation: "median" over all scenes per month
#       - resampling: "bilinear" for spatial interpolation
v <- gdalcubes::cube_view(
  srs   = "EPSG:32632",
  extent = list(
    t0    = "2018-06",
    t1    = "2022-09",
    left  = bbox_utm["xmin"] - 10,
    right = bbox_utm["xmax"] + 10,
    bottom= bbox_utm["ymin"] - 10,
    top   = bbox_utm["ymax"] + 10
  ),
  dx          = 10,          # 10 m pixel size in x-direction
  dy          = 10,          # 10 m pixel size in y-direction
  dt          = "P1M",       # P1M = 1 month 
  aggregation = "median",    # temporal reducer: median of all scenes per month
  resampling  = "bilinear"   # spatial resampling of original COGs
)

# Inspect the cube view definition
v


## ----burgwald_get_data_write_ncdf------------------------------------------

# 6) Build the actual raster cube from the STAC collection and cube view,
#    apply a cloud mask based on SCL, and write the result to a NetCDF file.
#
#    - image_mask("SCL", values = c(3, 8, 9)) creates a mask where SCL equals
#      one of the specified classes and sets those pixels to NA.
#      Common interpretation for SCL:
#        3 = Cloud shadow
#        8 = Medium probability cloud
#        9 = High probability cloud
#
#    - gdalcubes_options(parallel = 16) uses 16 CPU threads (if available).
#    - ncdf_compression_level = 5 sets moderate NetCDF compression.

s2_mask <- gdalcubes::image_mask("SCL", values = c(3, 8, 9))

gdalcubes::gdalcubes_options(
  parallel               = 16,
  ncdf_compression_level = 5
)

# raster_cube(...) creates a virtual data cube over the COGs.
# write_ncdf(...) is the *action* that actually triggers reading from the
# cloud and writing the aggregated cube to disk as a NetCDF file.
cube <- gdalcubes::raster_cube(s2_collection, v, mask = s2_mask)
write_nc_if_missing(cube, out_file)




# 7) Open the NetCDF cube and compute mean kNDVI over time.
#
#    kNDVI (as used here) is defined as:
#      kNDVI = tanh( ((B08 - B04) / (B08 + B04))^2 )
#
#    where:
#      - B08 = NIR band
#      - B04 = red band
#
#    Steps:
#      a) ncdf_cube(...) connects the NetCDF file back into gdalcubes.
#      b) apply_pixel(...) computes kNDVI per pixel and time step.
#      c) reduce_time("mean(kNDVI)") collapses the time dimension by computing
#         the temporal mean per pixel.
#      d) plot(...) visualizes the resulting mean-kNDVI map.
gdalcubes::ncdf_cube(
  file.path(root_folder, "data/burgwald_2018_2022_all.nc")
) |>
  gdalcubes::apply_pixel(
    "tanh(((B08-B04)/(B08+B04))^2)",  # per-pixel kNDVI expression
    "kNDVI"                           # output band name
  ) |>
  gdalcubes::reduce_time("mean(kNDVI)") |>
  plot(
    key.pos  = 1,           # legend position (bottom)
    col      = ndvi.col(11),# NDVI-like color palette (user-defined function)
    nbreaks  = 12           # number of color breaks
  )
