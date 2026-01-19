#!/usr/bin/env Rscript

############################################################
# Script: 01-1_get_base-geodata.R
#
# Products (via paths):
#   - aoi_dgm
#   - aoi_clc
#   - osm_by_key
############################################################

source("src/_core/01-setup-burgwald.R")


############################################################
# 1) DGM1 DEM → mosaic and clip to AOI
############################################################

dem_out_file <- paths[["aoi_dgm"]]

run_if_missing(dem_out_file, {
  
  today <- format(Sys.Date(), "%Y%m%d")
  
  base_url_wf <- paste0(
    "https://gds.hessen.de/downloadcenter/", today,
    "/3D-Daten/Digitales%20Gel%C3%A4ndemodell%20(DGM1)/Landkreis%20Waldeck-Frankenberg/"
  )
  
  base_url_mr <- paste0(
    "https://gds.hessen.de/downloadcenter/", today,
    "/3D-Daten/Digitales%20Gelände­modell%20(DGM1)/Landkreis%20Marburg-Biedenkopf/"
  )
  
  dgm1_urls <- c(
    burgwald     = paste0(base_url_wf, "Burgwald%20-%20DGM1.zip"),
    gemuenden    = paste0(base_url_wf, "Gem%C3%BCnden%20(Wohra)%20-%20DGM1.zip"),
    rosenthal    = paste0(base_url_wf, "Rosenthal%20-%20DGM1.zip"),
    muenchhausen = paste0(base_url_mr, "M%C3%BCnchhausen%20-%20DGM1.zip"),
    rauschenberg = paste0(base_url_mr, "Rauschenberg%20-%20DGM1.zip"),
    coelbe       = paste0(base_url_mr, "C%C3%B6lbe%20-%20DGM1.zip"),
    lahntal      = paste0(base_url_mr, "Lahntal%20-%20DGM1.zip"),
    wohra        = paste0(base_url_mr, "Wohratal%20-%20DGM1.zip"),
    wetter       = paste0(base_url_mr, "Wetter%20(Hessen)%20-%20DGM1.zip"),
    haina        = paste0(base_url_wf, "Haina%20(Kloster)%20-%20DGM1.zip")
  )
  
  all_tif_files <- character(0)
  
  for (nm in names(dgm1_urls)) {
    
    url_i      <- dgm1_urls[[nm]]
    zip_file_i <- here::here("data", "raw", paste0("dgm1_", nm, ".zip"))
    unzip_dir  <- here::here("data", "raw", paste0("dgm1_", nm))
    
    download_if_missing(url_i, zip_file_i, mode = "wb")
    
    if (!fs::dir_exists(unzip_dir)) {
      fs::dir_create(unzip_dir)
      unzip(zip_file_i, exdir = unzip_dir)
    }
    
    tif_i <- dir(
      unzip_dir,
      pattern = "\\.tif$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    all_tif_files <- c(all_tif_files, tif_i)
  }
  
  dem_list <- lapply(all_tif_files, function(f) {
    r <- terra::rast(f)
    terra::crs(r) <- "EPSG:25832"
    r
  })
  
  dem_merged <- do.call(terra::mosaic, c(dem_list, fun = "min"))
  
  aoi_dem_sf <- sf::st_transform(aoi_burgwald_wgs, terra::crs(dem_merged))
  aoi_dem_v  <- terra::vect(aoi_dem_sf)
  
  dem_burgwald <- dem_merged |>
    terra::crop(aoi_dem_v) |>
    terra::mask(aoi_dem_v)
  
  terra::writeRaster(dem_burgwald, dem_out_file, overwrite = TRUE)
})


############################################################
# 2) CLC5 2018 → clip to AOI
############################################################

clc_out_file <- paths[["aoi_clc"]]

run_if_missing(clc_out_file, {
  
  clc_zip  <- here::here("data", "raw", "31916.zip")
  clc_root <- here::here("data", "raw", "clc5_2018_copernicus")
  unzip(clc_zip, exdir = clc_root)
  
  results_dir <- dir(clc_root, pattern = "Results", full.names = TRUE)[1]
  zipname     <- dir(results_dir, pattern = "\\.zip$", full.names = TRUE)[1]
  
  unzip(zipname, exdir = here::here("data", "raw"))
  
  data_dir <- here::here(
    "data", "raw",
    "u2018_clc2018_v2020_20u1_raster100m", "DATA"
  )
  
  clc_tif <- dir(data_dir, pattern = "\\.tif$", full.names = TRUE)[1]
  clc_rast <- terra::rast(clc_tif)
  
  aoi_clc   <- sf::st_transform(aoi_burgwald_wgs, terra::crs(clc_rast))
  aoi_clc_v <- terra::vect(aoi_clc)
  
  clc_crop <- clc_rast |>
    terra::crop(aoi_clc_v) |>
    terra::mask(aoi_clc_v)
  
  terra::writeRaster(clc_crop, clc_out_file, overwrite = TRUE)
})


############################################################
# 3) OpenStreetMap layers by key
############################################################

osm_dir <- paths[["osm_by_key"]]

run_if_missing(osm_dir, {
  Sys.setenv(OSM_TIMEOUT = 300)
  osm_by_key <- get_osm_burgwald_by_key(
    aoi_wgs84 = aoi_burgwald_wgs,
    out_dir   = osm_dir,
    keys      = c("highway", "landuse", "natural", "waterway",
                  "building", "railway"),
    write     = TRUE
  )
  
  mapview::mapview(osm_by_key$highway$lines)
  mapview::mapview(osm_by_key$landuse$polygons)
})
