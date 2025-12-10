#!/usr/bin/env Rscript

############################################################
# Script:  01-2-get_CDSE-sentinel-data.R
# Author:   [Your Name]
# Project:  [Your Project Name]
#
# Purpose:
# --------
## - Query Sentinel-2 L2A scenes from the Copernicus Data Space
##   Ecosystem (CDSE) over the Burgwald AOI for 2018–2022.
## - For each year, identify the **clearest summer day**
##   (June–August, minimum mean cloud cover).
## - For each selected date, request:
##     * 12 raw Sentinel-2 bands via `RawBands.js`
##     * kNDVI via `kndvi.js`
##     * SAVI via `savi.js`
##     * EVI via `evi.js`
## - Assemble all layers into a **predictor stack** and save it as
##   `pred_stack_<year>.tif` in the project `data/` directory.
##
## Inputs / assumptions:
## - A here()-based project layout with:
##     * `root_folder` defined (e.g. in 00-setup-burgwald.R)
##     * `aoi_burgwald_wgs` (sf polygon in EPSG:4326)
##     * `burgwald_bbox` (xmin/xmax/ymin/ymax)
## - Environment variables:
##     * `CDSE_ID`
##     * `CDSE_SECRET`
##   for CDSE OAuth authentication.
## - Local JavaScript scripts in `scripts/`:
##     * `kndvi.js`
##     * `savi.js`
##     * `evi.js`
##
## Output:
## - One multi-layer GeoTIFF per year:
##     `data/pred_stack_<year>.tif`
##   containing 12 S2 bands + EVI + kNDVI + SAVI.
##
## Notes:
## - The scientific / spatial behaviour is **not** modified here;
##   only documentation and English comments are added.
## ============================================================


## ------------------------------------------------------------
## 0) Setup: packages, project path, helper functions ----
## ------------------------------------------------------------

# pacman – convenience package to install & load packages
# (install pacman if necessary, then load it)
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)

# Install CDSE from GitHub if it is not yet available.
# CDSE = R client for Copernicus Data Space Ecosystem.
# (only runs once per machine / environment)
if (!requireNamespace("CDSE", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("zivankaraman/CDSE")
}

# Load all required packages in one go (installing if missing).
# This gives you spatial, remote sensing, ML and helper packages.
pacman::p_load(
  mapview, mapedit, tmap, tmaptools,           # interactive and thematic maps
  raster, terra, tidyterra, stars, gdalcubes, sf, # raster & vector data, cubes
  RStoolbox, exactextractr,                   # RS utilities, exact raster extraction
  randomForest, ranger, e1071, caret,         # ML / classification
  dplyr, ggplot2, tidyr,                      # tidyverse core
  link2GI, rstac, OpenStreetMap,              # GIS bindings, STAC, OSM tiles
  colorspace, ows4R, httr,                    # palettes, OGC services, HTTP
  here, CDSE, lubridate                       # project paths, CDSE client, dates
)



# Project setup: you assume a "here"-based project layout.
# 00-setup-burgwald.R is expected to define:
# - aoi_burgwald_wgs (sf polygon in EPSG:4326)
# - burgwald_bbox    (bbox vector xmin/xmax/ymin/ymax)
# Both are used as spatial constraints for the STAC queries.
source(here::here(root_folder,"src", "00-setup-burgwald.R"))
source(here::here(root_folder,"src", "01-fun-data-retrieval.R"))

## ------------------------------------------------------------
## 1) Geometries / CLC check (optional visual sanity check) ----
## ------------------------------------------------------------

# Simple ggplot showing the Burgwald AOI in red.
# This is only for visual verification and has no effect on later steps.
bg <- ggplot() +
  geom_sf(
    data      = sf::st_sf(geometry = aoi_burgwald_wgs),
    fill      = NA,
    color     = "red",
    linewidth = 0.8
  ) +
  coord_sf(
    xlim = c(burgwald_bbox["xmin"], burgwald_bbox["xmax"]),
    ylim = c(burgwald_bbox["ymin"], burgwald_bbox["ymax"])
  ) +
  theme_bw()

print(bg)

# Optional: load CORINE Land Cover raster for the Burgwald,
# if it exists in your project structure.
# This is a pure sanity check / inspection step.
aoi_root     <- here::here("data", "raw", "AOI_Burgwald")
clc_out_file <- file.path(aoi_root, "clc", "clc5_2018_burgwald.tif")
if (file.exists(clc_out_file)) {
  corine_burgwald <- terra::rast(clc_out_file)
  # mapview_clc() presumably is defined in helper-func.R
  mapview_clc(corine_burgwald)
}

## ------------------------------------------------------------
## 2) STAC client for Copernicus Data Space Ecosystem ----
## ------------------------------------------------------------

# Build a STAC client pointing to the CDSE STAC endpoint.
# This is the entry point for Sentinel-2 item discovery.
cdse_stac <- rstac::stac("https://stac.dataspace.copernicus.eu/v1/")

# Build a slightly buffered AOI (in degrees).
# buffer_deg = 0.0 means effectively no buffer; you can increase if needed.
aoi_burgwald_wgs_buf <- sf::st_buffer(aoi_burgwald_wgs, buffer_deg)
bbox_bw <- sf::st_bbox(aoi_burgwald_wgs_buf)

# Define the time interval for which you want to search Sentinel-2 imagery.
# Here: 2018-01-01 to 2022-12-31 (UTC).
start <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
end   <- as.POSIXct("2022-12-31 23:59:59", tz = "UTC")

# STAC datetime string in RFC3339 interval format: "start/end".
datetime_range <- paste0(
  format(start, "%Y-%m-%dT%H:%M:%SZ"), "/",
  format(end,   "%Y-%m-%dT%H:%M:%SZ")
)

## ------------------------------------------------------------
## 3) Search Sentinel-2 L2A over Burgwald and build daily stats ----
## ------------------------------------------------------------

# STAC search: Sentinel-2 L2A items intersecting your AOI bbox
# during the specified time range, limit to max 500 items.
s2_search_range <- cdse_stac |>
  rstac::stac_search(
    collections = "sentinel-2-l2a",
    bbox        = as.numeric(bbox_bw),
    datetime    = datetime_range,
    limit       = 500
  ) |>
  rstac::post_request()

# Fetch all items (may involve paging).
s2_items_range <- rstac::items_fetch(s2_search_range)
# Convert STAC items to a tibble, with one row per scene.
tbl_range      <- rstac::items_as_tibble(s2_items_range)

# Build a per-day summary:
# - date: Date extracted from datetime
# - year: year part of date
# - cloud: eo:cloud_cover property
# Then compute mean cloud coverage per date and year.
day_diag <- tbl_range |>
  dplyr::mutate(
    date  = as.Date(datetime),
    year  = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m")),
    cloud = as.numeric(`eo:cloud_cover`)
  ) |>
  dplyr::group_by(year, date) |>
  dplyr::summarise(
    n_scenes   = dplyr::n(),                 # number of scenes that day
    mean_cloud = mean(cloud, na.rm = TRUE),  # average cloud cover
    .groups    = "drop"
  )

# Filter only the summer months (June–August).
day_summer <- day_diag |>
  dplyr::filter(month(date) %in% 6:8)

# For each year in the summer subset, select the date with the
# smallest mean_cloud (i.e., the "cleanest" summer day per year).
best_summer_days <- day_summer |>
  dplyr::group_by(year) |>
  dplyr::slice_min(mean_cloud, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

print(best_summer_days)
# At this point you have one row per year, containing:
# - year
# - date (best summer day)
# - n_scenes
# - mean_cloud

## ------------------------------------------------------------
## 4) CDSE OAuth client (ID/Secret from environment variables) ----
## ------------------------------------------------------------

# You expect the environment variables CDSE_ID and CDSE_SECRET
# to be set (e.g. in ~/.Renviron or exported in the shell).
id     <- Sys.getenv("CDSE_ID")
secret <- Sys.getenv("CDSE_SECRET")

if (id == "" || secret == "") {
  # Hard safety guard: stop early if credentials are missing.
  stop("CDSE_ID und/oder CDSE_SECRET sind nicht gesetzt (Umgebungsvariablen).")
}

# Create an OAuth client object using the CDSE package.
# This encapsulates client credentials and token handling.
OAuthClient <- CDSE::GetOAuthClient(id = id, secret = secret)

## ------------------------------------------------------------
## 5) JavaScript processing scripts for RawBands, kNDVI, SAVI, EVI ----
## ------------------------------------------------------------

# RawBands.js is shipped inside the CDSE package (in its inst/scripts/).
# It outputs 12 Sentinel-2 bands as separate layers.
script_file_raw <- system.file("scripts", "RawBands.js", package = "CDSE")
if (script_file_raw == "") {
  stop("RawBands.js wurde im CDSE-Paket nicht gefunden.")
}

# kNDVI, SAVI, EVI are assumed to be local JS files in your project
# (these scripts define how the EOxHub processing service computes the index).
script_file_kndvi <- here::here("scripts", "kndvi.js")
script_file_savi  <- here::here("scripts", "savi.js")
script_file_evi   <- here::here("scripts", "evi.js")

# Safety guard: fail fast if any JS file is missing.
for (f in c(script_file_kndvi, script_file_savi, script_file_evi)) {
  if (!file.exists(f)) {
    stop("JS-Script nicht gefunden: ", f)
  }
}

## ------------------------------------------------------------
## 6) Helper function: build predictor stack for ONE date ----
## ------------------------------------------------------------

# This function:
#   - takes a single date, an AOI, the OAuth client, and 4 JS scripts
#   - calls CDSE::GetImage 4 times:
#       1) RawBands.js  (12 bands)
#       2) kndvi.js     (1 band: kNDVI)
#       3) savi.js      (1 band: SAVI)
#       4) evi.js       (1 band: EVI)
#   - checks that each result is a terra::SpatRaster
#   - renames all bands to meaningful names
#   - returns a stacked SpatRaster with all predictors.
get_pred_stack_for_date <- function(
    the_date,
    aoi,
    client,
    script_file_raw,
    script_file_kndvi,
    script_file_savi,
    script_file_evi
) {
  # Normalize date to character "YYYY-MM-DD".
  day_str <- format(as.Date(the_date), "%Y-%m-%d")
  message("==> Hole Daten für ", day_str)
  
  # ---- 6.1 Raw bands (12 bands via RawBands.js)
  raw <- CDSE::GetImage(
    aoi              = aoi,                   # sf polygon in WGS84
    time_range       = day_str,              # one-day window
    script           = script_file_raw,      # RawBands.js
    collection       = "sentinel-2-l2a",     # S2 L2A collection
    format           = "image/tiff",         # GeoTIFF
    mosaicking_order = "leastCC",            # least cloud cover
    resolution       = 10,                   # 10 m resolution
    client           = client                # OAuth client
  )
  
  # Guard: ensure we got a SpatRaster.
  if (!inherits(raw, "SpatRaster")) {
    stop("RawBands-GetImage() hat kein SpatRaster zurückgegeben (Datum: ", day_str, ").")
  }
  
  # If we have 12 bands, assign S2 band names.
  if (terra::nlyr(raw) == 12L) {
    names(raw) <- c(
      "B01", "B02", "B03", "B04", "B05", "B06",
      "B07", "B08", "B8A", "B09", "B11", "B12"
    )
  }
  
  # ---- 6.2 kNDVI image
  kndvi <- CDSE::GetImage(
    aoi              = aoi,
    time_range       = day_str,
    script           = script_file_kndvi,
    collection       = "sentinel-2-l2a",
    format           = "image/tiff",
    mosaicking_order = "leastCC",
    resolution       = 10,
    mask             = TRUE,    # apply cloud/snow mask as defined in JS
    client           = client
  )
  if (!inherits(kndvi, "SpatRaster")) {
    stop("kNDVI-GetImage() hat kein SpatRaster zurückgegeben (Datum: ", day_str, ").")
  }
  names(kndvi)[1] <- "kNDVI"
  
  # ---- 6.3 SAVI image
  savi <- CDSE::GetImage(
    aoi              = aoi,
    time_range       = day_str,
    script           = script_file_savi,
    collection       = "sentinel-2-l2a",
    format           = "image/tiff",
    mosaicking_order = "leastCC",
    resolution       = 10,
    mask             = TRUE,
    client           = client
  )
  if (!inherits(savi, "SpatRaster")) {
    stop("SAVI-GetImage() hat kein SpatRaster zurückgegeben (Datum: ", day_str, ").")
  }
  names(savi)[1] <- "SAVI"
  
  # ---- 6.4 EVI image
  evi <- CDSE::GetImage(
    aoi              = aoi,
    time_range       = day_str,
    script           = script_file_evi,
    collection       = "sentinel-2-l2a",
    format           = "image/tiff",
    mosaicking_order = "leastCC",
    resolution       = 10,
    mask             = TRUE,
    client           = client
  )
  if (!inherits(evi, "SpatRaster")) {
    stop("EVI-GetImage() hat kein SpatRaster zurückgegeben (Datum: ", day_str, ").")
  }
  names(evi)[1] <- "EVI"
  
  # ---- 6.5 Build the final predictor stack
  # Concatenate:
  #   - all 12 raw spectral bands
  #   - one-band EVI layer
  #   - one-band kNDVI layer
  #   - one-band SAVI layer
  pred_stack <- c(
    raw,
    evi[[1]],
    kndvi[[1]],
    savi[[1]]
  )
  
  pred_stack
}

## ------------------------------------------------------------
## 7) Loop: for each year, download one summer predictor stack ----
## ------------------------------------------------------------

# Choice of AOI to be passed to GetImage: here, the buffered AOI.
# You could also use aoi_burgwald_wgs directly if desired.
aoi_for_download <- aoi_burgwald_wgs_buf   # could also be aoi_burgwald_wgs

# Output directory where the yearly stacks will be stored as GeoTIFFs.
out_dir <- here::here("data")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Iterate over each row in best_summer_days.
# For each year:
#   - get the "best" summer date
#   - call get_pred_stack_for_date()
#   - write the result as pred_stack_<year>.tif
for (i in seq_len(nrow(best_summer_days))) {
  yr <- best_summer_days$year[i]
  dt <- best_summer_days$date[i]
  
  message("============================================")
  message("Jahr ", yr, " – Datum ", dt, " (Sommer, min. cloud)")
  
  stack_year <- get_pred_stack_for_date(
    the_date          = dt,
    aoi               = aoi_for_download,
    client            = OAuthClient,
    script_file_raw   = script_file_raw,
    script_file_kndvi = script_file_kndvi,
    script_file_savi  = script_file_savi,
    script_file_evi   = script_file_evi
  )
  
  # Build filename like "pred_stack_2018.tif", "pred_stack_2019.tif", ...
  out_file <- file.path(out_dir, sprintf("pred_stack_%d.tif", yr))
  
  # Write the SpatRaster stack to disk (overwrite if exists).
  terra::writeRaster(
    stack_year,
    out_file,
    overwrite = TRUE
  )
  
  message(">> geschrieben: ", out_file)
}

message("Fertig: Sommer-Prädiktor-Stacks für alle Jahre in ", out_dir)
