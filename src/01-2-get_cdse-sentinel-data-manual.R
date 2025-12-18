#!/usr/bin/env Rscript

############################################################
# Script:  01-2-get_s2_cdse_manual.R
# Author:  [Your Name]
# Project: [Your Project Name]
#
# Purpose:
# --------
# - Demonstrate how to build Sentinel-2 predictor stacks for the
#   Burgwald AOI *without* the CDSE R wrapper.
# - All remote calls use the CDSE Process API directly (httr + JSON).
# - Spectral indices (kNDVI, SAVI, EVI) are computed in R.
#
# Workflow:
# ---------
# 1) Use STAC (rstac) to:
#      - query Sentinel-2 L2A scenes (2018–2022) over Burgwald
#      - find the "best" summer day (Jun–Aug, min mean cloud) per year
# 2) For each selected date:
#      - call Process API once to get 12 raw bands (B01–B12)
#      - compute kNDVI, SAVI, EVI in R
#      - write pred_stack_<year>.tif with 12 bands + 3 indices
#
# Assumptions:
# ------------
# - here()-based project layout
# - src/00-setup-burgwald.R defines:
#     * aoi_burgwald_wgs (sf polygon, EPSG:4326)
#     * burgwald_bbox    (xmin/xmax/ymin/ymax)
# - Environment variables:
#     * CDSE_ID
#     * CDSE_SECRET
############################################################


## ------------------------------------------------------------
## 0) Packages & setup
## ------------------------------------------------------------

needed <- c(
  "terra", "sf", "rstac", "httr", "jsonlite",
  "here", "dplyr", "lubridate", "ggplot2"
)

inst <- needed[!needed %in% rownames(installed.packages())]
if (length(inst) > 0) install.packages(inst)


# Project-specific setup: AOI etc.
# Expects aoi_burgwald_wgs and burgwald_bbox to be defined there.
source(here::here("src", "00-setup-burgwald.R"))

# Optional: quick geometry sanity check
bg <- ggplot() +
  geom_sf(
    data      = st_sf(geometry = aoi_burgwald_wgs),
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


## ------------------------------------------------------------
## 1) Pure R implementations of kNDVI, SAVI, EVI
## ------------------------------------------------------------

# All functions assume reflectances in [0,1]

kndvi_from_red_nir <- function(red, nir) {
  x <- (nir - red) / (nir + red)
  tanh(x^2)
}

savi_from_red_nir <- function(red, nir, L = 0.5) {
  (nir - red) * (1 + L) / (nir + red + L)
}

evi_from_red_nir_blue <- function(red, nir, blue,
                                  G  = 2.5,
                                  C1 = 6.0,
                                  C2 = 7.5,
                                  L  = 1.0) {
  num <- nir - red
  den <- nir + C1 * red - C2 * blue + L
  G * num / den
}

# terra-helper: add indices to a SpatRaster stack
add_veg_indices_to_stack <- function(
    stack,
    nir_band  = "B08",
    red_band  = "B04",
    blue_band = "B02",
    savi_L    = 0.5,
    evi_G     = 2.5,
    evi_C1    = 6.0,
    evi_C2    = 7.5,
    evi_L     = 1.0
) {
  stopifnot(inherits(stack, "SpatRaster"))
  
  needed <- c(nir_band, red_band, blue_band)
  missing <- setdiff(needed, names(stack))
  if (length(missing) > 0) {
    stop("Missing bands in stack: ", paste(missing, collapse = ", "))
  }
  
  nir  <- stack[[nir_band]]
  red  <- stack[[red_band]]
  blue <- stack[[blue_band]]
  
  # kNDVI
  kndvi <- (nir - red) / (nir + red)
  kndvi <- tanh(kndvi^2)
  names(kndvi) <- "kNDVI"
  
  # SAVI
  savi <- (nir - red) * (1 + savi_L) / (nir + red + savi_L)
  names(savi) <- "SAVI"
  
  # EVI
  num <- nir - red
  den <- nir + evi_C1 * red - evi_C2 * blue + evi_L
  evi <- evi_G * num / den
  names(evi) <- "EVI"
  
  c(stack, kndvi, savi, evi)
}


## ------------------------------------------------------------
## 2) CDSE OAuth + Process API (manual)
## ------------------------------------------------------------

cdse_get_access_token <- function(client_id, client_secret) {
  resp <- httr::POST(
    url  = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token",
    body = list(
      grant_type    = "client_credentials",
      client_id     = client_id,
      client_secret = client_secret
    ),
    encode = "form"
  )
  
  # Debug-Ausgabe: Status + Body
  message("Token status: ", httr::status_code(resp))
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  #message("Token response body:\n", txt)
  
  httr::stop_for_status(resp)
  
  tok <- jsonlite::fromJSON(txt)$access_token
  
  if (is.null(tok) || !nzchar(tok)) {
    stop("Failed to obtain access token from CDSE (no access_token field).")
  }
  tok
}


# Evalscript for 12 raw Sentinel-2 bands as reflectance
evalscript_raw <- "
//VERSION=3
function setup() {
  return {
    input: [{
      bands: [
        'B01','B02','B03','B04','B05','B06',
        'B07','B08','B8A','B09','B11','B12'
      ],
      units: 'REFLECTANCE'
    }],
    output: {
      bands: 12,
      sampleType: 'FLOAT32'
    }
  };
}

function evaluatePixel(s) {
  return [
    s.B01, s.B02, s.B03, s.B04, s.B05, s.B06,
    s.B07, s.B08, s.B8A, s.B09, s.B11, s.B12
  ];
}
"

cdse_get_image_raw <- function(the_date,
                               aoi,
                               access_token,
                               evalscript = evalscript_raw,
                               width      = 1024,
                               height     = 1024) {
  stopifnot(inherits(aoi, c("sf", "sfc")))
  if (inherits(aoi, "sf")) aoi <- st_geometry(aoi)
  
  if (st_crs(aoi)$epsg != 4326) {
    aoi <- st_transform(aoi, 4326)
  }
  
  day_str <- format(as.Date(the_date), "%Y-%m-%d")
  message("==> Requesting raw bands for ", day_str)
  
  bb <- st_bbox(aoi)
  
  body <- list(
    input = list(
      bounds = list(
        properties = list(
          crs = "http://www.opengis.net/def/crs/EPSG/0/4326"
        ),
        bbox = as.numeric(bb)
      ),
      data = list(list(
        type = "sentinel-2-l2a",
        dataFilter = list(
          timeRange = list(
            from = paste0(day_str, "T00:00:00Z"),
            to   = paste0(day_str, "T23:59:59Z")
          ),
          mosaickingOrder = "leastCC"
        )
      ))
    ),
    output = list(
      width  = as.integer(width),
      height = as.integer(height),
      responses = list(list(
        identifier = "default",
        format = list(type = "image/tiff")
      ))
    ),
    evalscript = evalscript
  )
  
  resp <- httr::POST(
    url = "https://sh.dataspace.copernicus.eu/api/v1/process",
    httr::add_headers(
      Authorization = paste("Bearer", access_token),
      "Content-Type" = "application/json",
      Accept         = "image/tiff"     # <- HIER muss es stehen!
    ),
    body = jsonlite::toJSON(body, auto_unbox = TRUE)
  )
  
  status <- httr::status_code(resp)
  # txt    <- httr::content(resp, as = "text", encoding = "UTF-8")
  # 
  # message("CDSE status: ", status)
  # message("CDSE response (first 2000 chars):\n",
  #         substr(txt, 1, 2000))
  # 
  # # KEIN stop_for_status() hier!
  # if (status >= 400) {
  #   stop("CDSE request failed with status ", status)
  # }
  # 
  httr::stop_for_status(resp)
  
  tf <- tempfile(fileext = ".tif")
  writeBin(httr::content(resp, as = "raw"), tf)
  terra::rast(tf)
}


## ------------------------------------------------------------
## 3) Helper: get predictor stack for ONE date (manual)
## ------------------------------------------------------------

get_pred_stack_for_date_manual <- function(the_date,
                                           aoi,
                                           access_token) {
  raw_stack <- cdse_get_image_raw(
    the_date     = the_date,
    aoi          = aoi,
    access_token = access_token,
    evalscript   = evalscript_raw
  )
  
  # Benenne Bänder, falls nicht gesetzt
  if (nlyr(raw_stack) == 12L) {
    names(raw_stack) <- c(
      "B01","B02","B03","B04","B05","B06",
      "B07","B08","B8A","B09","B11","B12"
    )
  }
  
  add_veg_indices_to_stack(
    stack     = raw_stack,
    nir_band  = "B08",
    red_band  = "B04",
    blue_band = "B02"
  )
}


## ------------------------------------------------------------
## 4) STAC: find best summer day per year (2018–2022)
## ------------------------------------------------------------

cdse_stac <- rstac::stac("https://stac.dataspace.copernicus.eu/v1/")

buffer_deg = 0.01 #
aoi_burgwald_wgs_buf <- st_buffer(aoi_burgwald_wgs, buffer_deg)
bbox_bw <- st_bbox(aoi_burgwald_wgs_buf)

start <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
end   <- as.POSIXct("2022-12-31 23:59:59", tz = "UTC")

datetime_range <- paste0(
  format(start, "%Y-%m-%dT%H:%M:%SZ"), "/",
  format(end,   "%Y-%m-%dT%H:%M:%SZ")
)

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
## 5) Main loop: build pred_stack_<year>.tif (manual CDSE)
## ------------------------------------------------------------

id     <- Sys.getenv("CDSE_ID")
secret <- Sys.getenv("CDSE_SECRET")
if (id == "" || secret == "") {
  stop("CDSE_ID and/or CDSE_SECRET not set as environment variables.")
}

access_token <- cdse_get_access_token(id, secret)

aoi_for_download <- aoi_burgwald_wgs_buf
out_dir <- here::here("data")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (i in seq_len(nrow(best_summer_days))) {
  yr <- best_summer_days$year[i]
  dt <- best_summer_days$date[i]
  
  message("============================================")
  message("Year ", yr, " – date ", dt, " (summer, min cloud)")
  
  stack_year <- get_pred_stack_for_date_manual(
    the_date     = dt,
    aoi          = aoi_for_download,
    access_token = access_token
  )
  
  out_file <- file.path(out_dir, sprintf("pred_stack_%d.tif", yr))
  
  terra::writeRaster(
    stack_year,
    out_file,
    overwrite = TRUE
  )
  
  message(">> written: ", out_file)
}

message("Done: manual CDSE Process API predictor stacks in ", out_dir)
