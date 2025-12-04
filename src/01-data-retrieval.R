############################################################
# 01-data-retrieval.R
#
# Purpose:
#   This script prepares all *base geodata* for the Burgwald AOI:
#
#   1) Topography:
#        - Download and mosaic DGM1 (1 m) for the surrounding counties
#        - Clip the DEM to the Burgwald AOI
#
#   2) Land cover:
#        - Extract and crop CLC5 2018 (Copernicus 100 m raster)
#        - Store a clipped land-cover raster for the AOI
#
#   3) Vector context data (OpenStreetMap):
#        - Download OSM feature groups (highway, landuse, natural, waterway,
#          building, railway) intersecting the AOI
#        - Save one file per OSM key under data/raw/AOI_Burgwald/osm_by_key/
#
#   4) Meteorological forcing (DWD):
#        - Download and preprocess DWD station data for the Burgwald region:
#            * Hourly wind (FF, DD)
#            * Hourly precipitation (R1)
#            * Sub-hourly precipitation (10-min and 5-min resolution)
#        - Time period is controlled by start_date / end_date (default: last 2 years)
#
# Output:
#   All base layers are stored under:
#     data/raw/AOI_Burgwald/
#
#   plus processed DWD station data under:
#     data/raw/dwd-stations/
#     data/processed/dwd/
#
# Notes:
#   - Satellite data (Sentinel / CDSE / gdalcubes) are fetched and processed
#     in *separate* scripts; this file only prepares DEM, land cover,
#     OSM context layers and DWD meteorological data.
#
# Requirements:
#   - 00-setup-burgwald.R must be executed first
#       → loads packages
#       → defines aoi_burgwald_wgs, aoi_root, run_if_missing(), download_if_missing()
#       → sets global DWD paths (path_dwd_raw, path_dwd_processed)
############################################################


# sourcing of the setup and specific used funtions
source("src/00-setup-burgwald.R")


# statdate/ enddate determines the meteo data download period

start_date <- Sys.Date() - 365 * 2
end_date   <- Sys.Date()

start_dt <- as.POSIXct(paste0(start_date, " 00:00:00"), tz = "UTC")
end_dt   <- as.POSIXct(paste0(end_date,   " 23:59:59"), tz = "UTC")



############################################################
# 1) DGM1 processing (1 m) → mosaic + clip
#    Output: data/raw/AOI_Burgwald/dem/dem_dgm1_burgwald.tif
############################################################

dem_out_file <- file.path(aoi_root, "dem", "dem_dgm1_burgwald.tif")

run_if_missing(dem_out_file, {
  
  # HVBG provides date-specific links, so URLs are built dynamically.
  today <- format(Sys.Date(), "%Y%m%d")
  
  base_url_wf <- paste0(
    "https://gds.hessen.de/downloadcenter/", today,
    "/3D-Daten/Digitales%20Gel%C3%A4ndemodell%20(DGM1)/Landkreis%20Waldeck-Frankenberg/"
  )
  
  base_url_mr <- paste0(
    "https://gds.hessen.de/downloadcenter/", today,
    "/3D-Daten/Digitales%20Gelände­modell%20(DGM1)/Landkreis%20Marburg-Biedenkopf/"
  )
  
  # DGM1 ZIPs covering the Burgwald region (Hesse).
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
  
  # ----------------------------------------------
  # Download all DGM1 ZIP files + extract all TIFs
  # ----------------------------------------------
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
  
  # Create a mosaic of all DGM tiles
  dem_list <- lapply(all_tif_files, function(f) {
    r <- terra::rast(f)
    terra::crs(r) <- "EPSG:25832"  # ETRS89 / UTM 32N
    r
  })
  
  dem_merged <- do.call(terra::mosaic, c(dem_list, fun = "min"))
  
  # Clip mosaic to AOI
  aoi_dem_sf <- sf::st_transform(aoi_burgwald_wgs, terra::crs(dem_merged))
  aoi_dem_v  <- terra::vect(aoi_dem_sf)
  
  dem_burgwald <- dem_merged |>
    terra::crop(aoi_dem_v) |>
    terra::mask(aoi_dem_v)
  
  terra::writeRaster(dem_burgwald, dem_out_file, overwrite = TRUE)
  message("✓ DGM1 saved to: ", dem_out_file)
})


############################################################
# 2) CLC5 2018 – extract Copernicus raster → clip to AOI
#    Output: data/raw/AOI_Burgwald/clc/clc5_2018_burgwald.tif
############################################################

clc_out_file <- file.path(aoi_root, "clc", "clc5_2018_burgwald.tif")

run_if_missing(clc_out_file, {
  
  # Main CLC ZIP containing subdirectories and an inner ZIP
  clc_zip  <- here::here("data", "raw", "31916.zip")
  clc_root <- here::here("data", "raw", "clc5_2018_copernicus")
  unzip(clc_zip, exdir = clc_root)
  
  # Find "Results/" inside Copernicus archive
  results_dir <- dir(clc_root, pattern = "Results", full.names = TRUE)[1]
  zipname     <- dir(results_dir, pattern = "\\.zip$", full.names = TRUE)[1]
  
  # Extract raster100m package into data/raw/
  unzip(zipname, exdir = here::here("data", "raw"))
  
  # Path inside the Copernicus structure
  data_dir <- here::here(
    "data", "raw",
    "u2018_clc2018_v2020_20u1_raster100m", "DATA"
  )
  
  clc_tif <- dir(data_dir, pattern = "\\.tif$", full.names = TRUE)[1]
  
  clc_rast <- terra::rast(clc_tif)
  
  # Clip to AOI
  aoi_clc   <- sf::st_transform(aoi_burgwald_wgs, terra::crs(clc_rast))
  aoi_clc_v <- terra::vect(aoi_clc)
  
  clc_crop <- clc_rast |>
    terra::crop(aoi_clc_v) |>
    terra::mask(aoi_clc_v)
  
  terra::writeRaster(clc_crop, clc_out_file, overwrite = TRUE)
  message("✓ CLC5 raster saved to: ", clc_out_file)
})


############################################################
# 3) OSM download by key (roads, landuse, natural, …)
#    Output: one file per OSM key under:
#            data/raw/AOI_Burgwald/osm_by_key/
############################################################

osm_dir <- file.path(aoi_root, "osm_by_key")

osm_by_key <- get_osm_burgwald_by_key(
  aoi_wgs84 = aoi_burgwald_wgs,
  out_dir   = osm_dir,
  keys      = c("highway", "landuse", "natural", "waterway",
                "building", "railway"),
  write     = TRUE
)

# Example visualization
mapview::mapview(osm_by_key$highway$lines)
mapview::mapview(osm_by_key$landuse$polygons)

FD_bw <- burgwald_get_hourly_dwd(
  var        = "wind",           # << NICHT "precipitation"
  params     = c("FF", "DD"),    # Windgeschwindigkeit + Richtung
  start_date = start_date,
  end_date   = end_date
)

  rr_60min_bw <- burgwald_get_hourly_dwd(
    var        = "precipitation",
    params     = "R1",
    start_date = start_date,
    end_date   = end_date
  )

rr_10min <- burgwald_get_subhourly_precip(
  resolution = "10min",
  start_date = start_date,
  end_date   = end_date
)

rr_5min <- burgwald_get_subhourly_precip(
  resolution = "5min",
  start_date = start_date,
  end_date   = end_date
)

