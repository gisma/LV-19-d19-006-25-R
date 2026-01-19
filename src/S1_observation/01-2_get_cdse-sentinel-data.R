#!/usr/bin/env Rscript

#######################################################################
# Script:   01-2_get_cdse-sentinel-data.R
#
# Purpose:
#   Build yearly Sentinel-2 summer predictor stacks for Burgwald via CDSE.
#   For each year, pick the best summer day (min mean cloud cover) and
#   request a multi-layer GeoTIFF (12 S2 bands + EVI + kNDVI + SAVI).
#
# Contract logic:
#   - Inputs are taken from setup objects (AOI etc.) and canonical S-paths.
#   - Outputs are written ONLY to canonical S2 paths defined in metadata/outputs.tsv.
#
# Requirements:
#   - src/_core/01-setup-burgwald.R provides:
#       * aoi_burgwald_wgs, burgwald_bbox, paths
#   - metadata/outputs.tsv contains keys:
#       * s2_pred_stack_<year> (S2_features)
#   - Environment variables:
#       * CDSE_ID, CDSE_SECRET
#######################################################################

source("src/_core/01-setup-burgwald.R")

## ------------------------------------------------------------
## 1) Geometries / optional CLC check (S1)
## ------------------------------------------------------------

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

# Optional: load canonical S1 CLC if present
clc_s1 <- paths[["aoi_clc"]]
if (!is.null(clc_s1) && file.exists(clc_s1)) {
  corine_burgwald <- terra::rast(clc_s1)
  mapview_clc(corine_burgwald)
}

## ------------------------------------------------------------
## 2) STAC client for Copernicus Data Space Ecosystem
## ------------------------------------------------------------

cdse_stac <- rstac::stac("https://stac.dataspace.copernicus.eu/v1/")

buffer_deg <- 0.01
aoi_burgwald_wgs_buf <- sf::st_buffer(aoi_burgwald_wgs, buffer_deg)
bbox_bw <- sf::st_bbox(aoi_burgwald_wgs_buf)

start <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
end   <- as.POSIXct("2022-12-31 23:59:59", tz = "UTC")

datetime_range <- paste0(
  format(start, "%Y-%m-%dT%H:%M:%SZ"), "/",
  format(end,   "%Y-%m-%dT%H:%M:%SZ")
)

## ------------------------------------------------------------
## 3) Search Sentinel-2 L2A over Burgwald and build daily stats
## ------------------------------------------------------------

s2_search_range <- cdse_stac |>
  rstac::stac_search(
    collections = "sentinel-2-l2a",
    bbox        = as.numeric(bbox_bw),
    datetime    = datetime_range,
    limit       = 500
  ) |>
  rstac::post_request()

s2_items_range <- rstac::items_fetch(s2_search_range)
tbl_range      <- rstac::items_as_tibble(s2_items_range)

day_diag <- tbl_range |>
  dplyr::mutate(
    date  = as.Date(datetime),
    year  = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m")),
    cloud = as.numeric(`eo:cloud_cover`)
  ) |>
  dplyr::group_by(year, date) |>
  dplyr::summarise(
    n_scenes   = dplyr::n(),
    mean_cloud = mean(cloud, na.rm = TRUE),
    .groups    = "drop"
  )

day_summer <- day_diag |>
  dplyr::filter(lubridate::month(date) %in% 6:8)

best_summer_days <- day_summer |>
  dplyr::group_by(year) |>
  dplyr::slice_min(mean_cloud, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

print(best_summer_days)

## ------------------------------------------------------------
## 4) CDSE OAuth client (ID/Secret from env)
## ------------------------------------------------------------

id     <- Sys.getenv("CDSE_ID")
secret <- Sys.getenv("CDSE_SECRET")
if (id == "" || secret == "") {
  stop("CDSE_ID and/or CDSE_SECRET are not set (environment variables).")
}
OAuthClient <- CDSE::GetOAuthClient(id = id, secret = secret)

## ------------------------------------------------------------
## 5) JavaScript processing scripts
## ------------------------------------------------------------

script_file_raw <- system.file("scripts", "RawBands.js", package = "CDSE")
if (script_file_raw == "") stop("RawBands.js was not found in the CDSE package.")

script_file_kndvi <- here::here("src","tools", "kndvi.js")
script_file_savi  <- here::here("src","tools", "savi.js")
script_file_evi   <- here::here("src", "tools","evi.js")

for (f in c(script_file_kndvi, script_file_savi, script_file_evi)) {
  if (!file.exists(f)) stop("JS script not found: ", f)
}

## ------------------------------------------------------------
## 6) Helper: build predictor stack for ONE date
## ------------------------------------------------------------

get_pred_stack_for_date <- function(
    the_date,
    aoi,
    client,
    script_file_raw,
    script_file_kndvi,
    script_file_savi,
    script_file_evi
) {
  day_str <- format(as.Date(the_date), "%Y-%m-%d")
  message("==> Downloading data for ", day_str)
  
  raw <- CDSE::GetImage(
    aoi              = aoi,
    time_range       = day_str,
    script           = script_file_raw,
    collection       = "sentinel-2-l2a",
    format           = "image/tiff",
    mosaicking_order = "leastCC",
    resolution       = 10,
    client           = client
  )
  if (!inherits(raw, "SpatRaster")) {
    stop("RawBands GetImage() did not return a SpatRaster (date: ", day_str, ").")
  }
  if (terra::nlyr(raw) == 12L) {
    names(raw) <- c("B01","B02","B03","B04","B05","B06","B07","B08","B8A","B09","B11","B12")
  }
  
  kndvi <- CDSE::GetImage(
    aoi              = aoi,
    time_range       = day_str,
    script           = script_file_kndvi,
    collection       = "sentinel-2-l2a",
    format           = "image/tiff",
    mosaicking_order = "leastCC",
    resolution       = 10,
    mask             = TRUE,
    client           = client
  )
  if (!inherits(kndvi, "SpatRaster")) {
    stop("kNDVI GetImage() did not return a SpatRaster (date: ", day_str, ").")
  }
  names(kndvi)[1] <- "kNDVI"
  
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
    stop("SAVI GetImage() did not return a SpatRaster (date: ", day_str, ").")
  }
  names(savi)[1] <- "SAVI"
  
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
    stop("EVI GetImage() did not return a SpatRaster (date: ", day_str, ").")
  }
  names(evi)[1] <- "EVI"
  
  c(raw, evi[[1]], kndvi[[1]], savi[[1]])
}

## ------------------------------------------------------------
## 7) Loop: one summer predictor stack per year -> canonical S2 output
## ------------------------------------------------------------

aoi_for_download <- aoi_burgwald_wgs_buf

for (i in seq_len(nrow(best_summer_days))) {
  yr <- best_summer_days$year[i]
  dt <- best_summer_days$date[i]
  
  message("============================================")
  message("Year ", yr, " – date ", dt, " (summer, minimum mean cloud cover)")
  
  stack_year <- get_pred_stack_for_date(
    the_date          = dt,
    aoi               = aoi_for_download,
    client            = OAuthClient,
    script_file_raw   = script_file_raw,
    script_file_kndvi = script_file_kndvi,
    script_file_savi  = script_file_savi,
    script_file_evi   = script_file_evi
  )
  
  key <- sprintf("s2_pred_stack_%d", yr)
  out_file <- paths[[key]]
  if (is.null(out_file)) {
    stop("Missing outputs.tsv entry for key: ", key)
  }
  
  terra::writeRaster(stack_year, out_file, overwrite = TRUE)
  message(">> Written predictor stack: ", out_file)
}

message("Done.")
