#!/usr/bin/env Rscript

# 00_clim_main_RSK_Hessen.R
# Purpose: Interpolate DWD KL daily precipitation (RSK) over AOI = Hessen.
# Notes:
# - newpipe is fixed TRUE (kl_daily_* pipeline).
# - No municipality extraction here. Later extraction happens on watershed polygons.

# --- Setup (packages + folders + envrmt paths) -----------------------------
source("src/_core/01-setup-burgwald_RSK.R")

# --- Functions required for newpipe ---------------------------------------
source("src/functions/fun_clim_core_RSK.R")

# ---- Parameters -----------------------------------------------------------
epsg <- 3035
res  <- 500

startDate <- "2003-01-01"
endDate   <- "2024-12-31"

minStations <- 150
param <- "RSK"
newpipe <- TRUE

# ---- DEM (expects dem.rds under data/raw/dem.rds or adjust here) ----------
# If you already have dem.rds (stars) in data/raw/, keep as-is.
dem_path <- file.path(envrmt$path_data_lev0, "dem.rds")
if (!file.exists(dem_path)) {
  stop("DEM missing: ", dem_path, " (provide dem.rds as stars object in EPSG:3035 or add your DEM preparation here).")
}
dem <- readRDS(dem_path)
dem <- stars::st_warp(dem, crs = epsg)

# ---- AOI: Hessen (federal state boundary) --------------------------------
aoi_hessen <- geodata::gadm(country = "DEU", level = 1, path = envrmt$path_data_lev0)
aoi_hessen <- sf::st_as_sf(aoi_hessen)
aoi_hessen <- aoi_hessen[aoi_hessen$NAME_1 == "Hessen", ]
aoi_hessen <- sf::st_transform(aoi_hessen, crs = epsg)

# crop + mask DEM to Hessen
dem <- stars::st_crop(dem, sf::st_bbox(aoi_hessen))
dem <- stars::st_mask(dem, aoi_hessen)

# ---- DWD KL daily: build dt_core + sf for RSK -----------------------------
root_raw  <- file.path(envrmt$path_CDC_KL, "kl_daily_raw")
out_dir   <- file.path(envrmt$path_CDC_KL, "kl_daily_extracted")
meta_file <- file.path(root_raw, "KL_Tageswerte_Beschreibung_Stationen.txt")

dt_core <- kl_daily_core_prepare(
  start_date = startDate,
  end_date   = endDate,
  max_missing_frac = 0.25,
  delta_t = 7L,
  keep_cols = c("STATIONS_ID","MESS_DATUM","QN_4","RSK","RSKF","eor"),
  root_raw   = root_raw,
  out_dir    = out_dir,
  meta_file  = meta_file,
  do_download = FALSE,
  force_unzip = FALSE,
  force_merge = TRUE,
  auto_resolve_duplicates = TRUE,
  write_merged_csv = TRUE
)

dt_core_rds <- file.path(out_dir, sprintf("dt_core_%s_%s.rds",
  format(as.Date(startDate), "%Y%m%d"),
  format(as.Date(endDate),   "%Y%m%d")
))
saveRDS(dt_core, dt_core_rds)

cVar.sf <- kl_daily_extract_sf(dt = dt_core, param = "RSK")

# ---- Dates to process -----------------------------------------------------
dat_list_all <- sort(unique(as.character(cVar.sf$MESS_DATUM)))

target_dir <- file.path(envrmt$path_data_lev1, param)
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

tifs_existing <- list.files(target_dir, pattern = "\\.tif$", full.names = FALSE)
dates_done <- character(0)
if (length(tifs_existing) > 0) {
  dates_done <- sub(paste0("_", param, "\\.tif$"), "", basename(tifs_existing))
}

dat_list_window <- dat_list_all[
  as.Date(dat_list_all) >= as.Date(startDate) &
  as.Date(dat_list_all) <= as.Date(endDate)
]

dat_list_todo <- setdiff(dat_list_window, dates_done)

message(sprintf("[KED] %s: already: %d | todo: %d", param, length(dates_done), length(dat_list_todo)))

if (length(dat_list_todo) == 0) quit(save = "no", status = 0)

# ---- Interpolation loop ---------------------------------------------------
invisible(
  pbmcapply::pbmclapply(seq_along(dat_list_todo), function(n) {

    currentDate <- dat_list_todo[n]
    cd <- substr(currentDate, 1, 10)
    target_file <- file.path(target_dir, paste0(cd, "_", param, ".tif"))

    if (file.exists(target_file)) return(NULL)

    day_sf <- cVar.sf[cVar.sf$MESS_DATUM == as.Date(currentDate), ]

    dat <- sanitize_climate_param(day_sf, param, date = currentDate)

    geom_col <- intersect(c("geometry", "geom"), names(dat))
    if (length(geom_col) == 1 && geom_col != "geometry") {
      names(dat)[names(dat) == geom_col] <- "geometry"
    }
    sf::st_geometry(dat) <- "geometry"
    dat <- dat[, c("Stationshoehe", "tmp", "geometry")]

    dat$tmp <- as.numeric(dat$tmp)
    names(dat)[names(dat) == "tmp"] <- param

    data <- tidyr::drop_na(dat)

    if (sf::st_crs(data) != sf::st_crs(dem)) {
      data <- sf::st_transform(data, crs = epsg)
    }

    if (sum(!is.na(data[[param]])) <= minStations) {
      message(sprintf("Skip %s (n=%d < minStations=%d)", cd, sum(!is.na(data[[param]])), minStations))
      return(NULL)
    }

    data <- dplyr::distinct(data, geometry, .keep_all = TRUE)
    data <- sf::st_transform(data, sf::st_crs(dem))

    vm.auto <- automap::autofitVariogram(
      as.formula(paste(param, "~1")),
      input_data = data
    )

    pred <- gstat::krige(
      formula     = as.formula(paste(param, "~Stationshoehe")),
      locations   = data,
      newdata     = dem,
      model       = vm.auto$var_model,
      debug.level = -1
    )

    stars::write_stars(
      pred,
      target_file,
      NA_value  = -9999,
      overwrite = TRUE
    )

    rm(pred, data, dat, day_sf)
    gc()

    NULL
  }, mc.cores = max(1L, parallel::detectCores() - 1L))
)

