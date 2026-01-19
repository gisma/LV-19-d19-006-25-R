#!/usr/bin/env Rscript

############################################################
# Script: 01-3_get_gdalcubes-sentinel-data.R
#
# Purpose:
#   Retrieve and preprocess Sentinel-2 L2A cloud-optimized scenes
#   for the Burgwald AOI using a STAC -> gdalcubes workflow.
#
# Canonical output:
#   NetCDF cube is written to the canonical S2 path defined by outputs.tsv:
#     key = s2_gdalcubes_cube_2018_2022_all
#
# Requirements:
#   - src/_core/01-setup-burgwald.R provides:
#       * aoi_burgwald_wgs, paths, root_folder, ndvi.col, write_nc_if_missing()
############################################################

source("src/_core/01-setup-burgwald.R")

# Canonical S2 output path (must exist as key in metadata/outputs.tsv)
out_file <- paths[["s2_gdalcubes_cube_2018_2022_all"]]

## ----stac_burgwald_gdalcubes----------------------------------------------

# 1) Build a STAC client for the Earth-search API hosted by Element84 on AWS.
s <- rstac::stac("https://earth-search.aws.element84.com/v0")

# 2) Derive the Burgwald bounding box from the AOI polygon (EPSG:4326).
bbox_bw <- sf::st_bbox(aoi_burgwald_wgs)

# 3) Search for Sentinel-2 L2A COG scenes over the Burgwald AOI.
items <- s |>
  rstac::stac_search(
    collections = "sentinel-s2-l2a-cogs",
    bbox        = c(
      bbox_bw["xmin"],
      bbox_bw["ymin"],
      bbox_bw["xmax"],
      bbox_bw["ymax"]
    ),
    datetime    = "2018-06-01/2022-09-01",
    limit       = 600
  ) |>
  rstac::post_request()

items

# 4) Build a gdalcubes STAC image collection from the returned STAC items.
s2_collection <- gdalcubes::stac_image_collection(
  items$features,
  asset_names = c(
    "B01","B02","B03","B04","B05","B06",
    "B07","B08","B8A","B09","B11","SCL"
  ),
  property_filter = function(x) {
    x[["eo:cloud_cover"]] < 5
  }
)

s2_collection

## ----burgwald_cubeview-----------------------------------------------------

# 5) Define cube view in EPSG:32632 (UTM zone 32N).
sf::st_as_sfc(bbox_bw) |>
  sf::st_transform("EPSG:32632") |>
  sf::st_bbox() -> bbox_utm

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
  dx          = 10,
  dy          = 10,
  dt          = "P1M",
  aggregation = "median",
  resampling  = "bilinear"
)

v

## ----burgwald_get_data_write_ncdf------------------------------------------

# 6) Build the raster cube, apply cloud mask, and write to NetCDF.
s2_mask <- gdalcubes::image_mask("SCL", values = c(3, 8, 9))

gdalcubes::gdalcubes_options(
  parallel               = 16,
  ncdf_compression_level = 5
)

cube <- gdalcubes::raster_cube(s2_collection, v, mask = s2_mask)
write_nc_if_missing(cube, out_file)

kndvi_mean <- gdalcubes::ncdf_cube(out_file) |>
  gdalcubes::apply_pixel(
    "tanh(((B08-B04)/(B08+B04))^2)",
    "kNDVI"
  ) |>
  gdalcubes::reduce_time("mean(kNDVI)")

kndvi_out <- paths[["s2_kndvi_mean_2018_2022"]]


# write GeoTIFF
gdalcubes::write_tif(kndvi_mean, kndvi_out)
