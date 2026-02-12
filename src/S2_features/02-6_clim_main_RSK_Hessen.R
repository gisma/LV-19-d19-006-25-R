#!/usr/bin/env Rscript

# 00_clim_main_RSK_Hessen.R
# Purpose: Interpolate DWD KL daily precipitation (RSK) over AOI = Hessen.
# Notes:
# - newpipe is fixed TRUE (kl_daily_* pipeline).
# - No municipality extraction here. Later extraction happens on watershed polygons.

root_folder= here::here()

# --- Setup (packages + folders + envrmt paths) -----------------------------
source(file.path(root_folder, "src", "_core", "01-setup-burgwald.R"))

# --- Functions required for newpipe ---------------------------------------

source(file.path(root_folder, "src", "r-libs/", "fun_clim_data.R"))
# ---- Parameters -----------------------------------------------------------
epsg <- 3035
res  <- 500

start_date <- "2022-01-01"
end_date   <- "2024-12-31"

minStations <- 150
param <- "RSK"
newpipe <- TRUE

# ---- Default parameters ----
crs   = raster::crs("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
epsg  = 3035
res   = 500             # grid resolution for the prediction DEM
type      = "historical"   # or "recent" (last ~500 days)
var       = "kl"           # DWD climate variable group
reso      = "daily"
downloadDEM = FALSE        # set TRUE only once
getClimate  = TRUE         # download new climate data
calc_commu  = TRUE         # extract results per municipality
calc_bl     = TRUE        # optionally calculate one state only
bl          = "Hessen"     # name of state

PM          = FALSE        # specific handling for air pressure
param       = "RSK"        # one of: c("RSK","SDK","NM","UPM","TXK","TNK","TMK","TGK","VPM","PM ")
newpipe =TRUE

# ---- Prepare auxiliary data ----
if (downloadDEM) {
  dem_DEU <- prepare_DEM(envrmt = envrmt, bl = "DEU", crs = crs, res = res)
}
dem <- readRDS(file.path(envrmt$path_data_lev0, "dem.rds"))
dem <- stars::st_warp(dem, crs = epsg)

# ---- Climate variable loop (main prediction) ----
# Use one single parameter for now (override cVar)
cVar <- param

# ---------------------------------
# 0) Basis: Zeitfenster + Pfade
# ---------------------------------
# start_date <- as.Date("2003-01-01")
# end_date   <- as.Date("2024-12-31")

root_raw <- file.path(envrmt$path_CDC_KL, "kl_daily_raw")          # oder dein Pfad
out_dir  <- file.path(envrmt$path_CDC_KL, "kl_daily_extracted")    # oder dein Pfad
meta_file <- file.path(root_raw, "KL_Tageswerte_Beschreibung_Stationen.txt")
dt_core <- kl_daily_core_prepare(
  start_date = start_date,
  end_date   = end_date,
  max_missing_frac = 0.25,
  delta_t = 7L,
  keep_cols = c(
    "STATIONS_ID","MESS_DATUM","QN_4",
    "RSK","RSKF","SDK","SHK_TAG","NM","VPM","PM","TMK","UPM","TXK","TNK","TGK","eor"
  ),
  root_raw   = root_raw,
  out_dir    = out_dir,
  meta_file  = meta_file,
  do_download = TRUE,     # TRUE wenn ZIPs fehlen und geladen werden sollen
  force_unzip = TRUE,
  force_merge = TRUE,
  auto_resolve_duplicates = TRUE,
  write_merged_csv = TRUE
)

# optional: für schnelle Wiederverwendung einfrieren
dt_core_rds <- file.path(out_dir, sprintf("dt_core_%s_%s.rds",
                                          format(as.Date(start_date), "%Y%m%d"),
                                          format(as.Date(end_date),   "%Y%m%d"))
)
saveRDS(dt_core, dt_core_rds)




if(newpipe) {
  
  sf_PM  <- kl_daily_extract_sf(dt = dt_core, param = "PM")
  # sf_TMK <- kl_daily_extract_sf(dt = dt_core, param = "TMK")
  # sf_TXK <- kl_daily_extract_sf(dt = dt_core, param = "TXK")
  
  # oder aus RDS:
  dt_core <- readRDS(dt_core_rds)
  cVar.sf   <- kl_daily_extract_sf(dt = dt_core, param = "PM")
  
  
  
  # cVar.sf <- ex_clim_hard(startDate, endDate,
  #                reso = "daily", var = "kl",
  #              type_hist = "historical", type_recent = "recent",
  #              param = cVar,
  #              max_missing_frac = 0.25,
  #              delta_t = 7L,
  #              base_dir = envrmt$path_CDC_KL,
  #              do_download = FALSE,
  #              force_unzip = FALSE,
  #              auto_resolve_duplicates = TRUE,
  #           write_gpkg = F)
  #   
  
  
} else {
  # Download filtered and preprocessed daily climate data
  cVar.sf <- ex_clim_new(
    startDate = startDate,
    endDate   = endDate,
    reso      = reso,
    var       = var,
    type      = type,
    param     = param
  )
}


# ---- Extract all unique available dates ----
dat_list_all <- sort(unique(as.character(cVar.sf$MESS_DATUM)))

# Zielordner einmal definieren
target_dir  <- file.path(envrmt$path_data_lev1, cVar)
if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)

# Existierende TIFFs einlesen und Datum herausziehen (Format: YYYY-MM-DD_PARAM.tif)
tifs_existing <- list.files(target_dir, pattern = "\\.tif$", full.names = FALSE)
dates_done <- character(0)
if (length(tifs_existing) > 0) {
  dates_done <- sub(paste0("_", cVar, "\\.tif$"), "", basename(tifs_existing))
}

# Nur Zielzeitraum behalten (Schnittmenge mit start/end)
dat_list_window <- dat_list_all[
  as.Date(dat_list_all) >= as.Date(startDate) &
    as.Date(dat_list_all) <= as.Date(endDate)
]

# Fehlende Tage bestimmen
dat_list_todo <- setdiff(dat_list_window, dates_done)

message(sprintf("[KED] %s: bereits vorhanden: %d | zu berechnen: %d",
                cVar, length(dates_done), length(dat_list_todo)))



# Kurzschluss: falls nichts zu tun, spätere Schritte überspringen
if (length(dat_list_todo) > 0) {
  
  # ---- Loop über fehlende Tage (parallel mit Fortschrittsbalken) ----
  invisible(
    pbmcapply::pbmclapply(seq_along(dat_list_todo), function(n) {
      currentDate <- dat_list_todo[n]
      cd <- substr(currentDate, 1, 10)
      target_file <- file.path(target_dir, paste0(cd, "_", cVar, ".tif"))
      
      # Sicherheitsnetz – falls parallel zwei Jobs denselben Tag erwischen
      if (file.exists(target_file)) {
        message("Skipping existing file (race-safe): ", target_file)
        return(NULL)
      }
      
      # Tageswerte filtern
      #cVar.sf.day <- cVar.sf[as.character(cVar.sf$MESS_DATUM) == as.Date(currentDate), ]
      cVar.sf.day <- cVar.sf[cVar.sf$MESS_DATUM == as.Date(currentDate), ]
      
      # Einheiten-/Range-Check
      dat <- sanitize_climate_param(cVar.sf.day, cVar, date = currentDate)
      
      # Geometriespalte harmonisieren
      geom_col <- intersect(c("geometry", "geom"), names(dat))
      if (length(geom_col) == 1 && geom_col != "geometry") {
        names(dat)[names(dat) == geom_col] <- "geometry"
      }
      sf::st_geometry(dat) <- "geometry"
      dat <- dat[, c("Stationshoehe", "tmp", "geometry")]
      
      dat$tmp <- as.numeric(dat$tmp)
      names(dat)[names(dat) == "tmp"] <- cVar
      
      data <- tidyr::drop_na(dat)
      
      # CRS an DEM angleichen
      if (sf::st_crs(data) != sf::st_crs(dem)) {cVar.sf
        data <- sf::st_transform(data, crs = epsg)
      }
      
      # Mindestanzahl prüfen
      if (sum(!is.na(data[[cVar]])) > minStations) {
        
        # Dubletten raus, CRS sicherstellen
        data <- dplyr::distinct(data, geometry, .keep_all = TRUE)
        data <- sf::st_transform(data, sf::st_crs(dem))
        
        # Variogramm automatisch
        vm.auto <- automap::autofitVariogram(
          as.formula(paste(cVar, "~1")),
          input_data = data
        )
        
        # KED (Stationshöhe als Drift)
        pred <- gstat::krige(
          formula     = as.formula(paste(cVar, "~Stationshoehe")),
          locations   = data,
          newdata     = dem,
          model       = vm.auto$var_model,
          debug.level = -1
        )
        
        # Schreiben
        stars::write_stars(
          pred,
          target_file,
          NA_value  = -9999,
          overwrite = TRUE
        )
        
        rm(pred); gc()
        
      } else {
        # Fallback: leeres Raster
        stars::write_stars(
          dem * 0 - 9999,
          target_file,
          overwrite = TRUE
        )
      }
      
      message("Written: ", target_file)
      return(target_file)
      
    }, mc.cores = 14, mc.allow.recursive = TRUE)
  )
  
} else {
  message(sprintf("[KED] %s: Keine offenen Tage im Zeitraum %s bis %s.",
                  cVar, startDate, endDate))
}

# ---- COG-Konvertierung via GNU parallel (falls TIFFs existieren) ----
tifs <- list.files(target_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(tifs) > 0)

cores <- max(1L, parallel::detectCores() - 1L)

run_one_cog <- function(f) {
  out_tmp <- sub("\\.tif(f)?$", "_cog.tif", f, ignore.case = TRUE)
  if (identical(out_tmp, f)) out_tmp <- paste0(f, "_cog.tif")  # fallback
  
  args <- c(
    f, out_tmp,
    "-of", "COG",
    "-co", "COMPRESS=ZSTD",
    "-co", "LEVEL=7",
    "-co", "PREDICTOR=2",
    "-co", "BLOCKSIZE=256"
  )
  
  # call gdal_translate without shell quoting issues
  res <- system2("gdal_translate", args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(res, "status")
  if (is.null(status)) status <- 0L
  
  if (status != 0L) {
    return(list(file = f, ok = FALSE, status = status, log = res))
  }
  
  ok_mv <- file.rename(out_tmp, f)
  if (!ok_mv) {
    return(list(file = f, ok = FALSE, status = 1L, log = c(res, "file.rename() failed")))
  }
  
  list(file = f, ok = TRUE, status = 0L, log = res)
}

# Linux: mclapply nutzt Fork; auf Windows -> lapply
is_windows <- identical(tolower(Sys.info()[["sysname"]]), "windows")

results <- if (is_windows) {
  lapply(tifs, run_one_cog)
} else {
  parallel::mclapply(tifs, run_one_cog, mc.cores = cores)
}

# Quick summary
ok <- vapply(results, `[[`, logical(1), "ok")
cat("OK:", sum(ok), " / FAIL:", sum(!ok), "\n")

# Print first failure (real error message)
if (any(!ok)) {
  i <- which(!ok)[1]
  cat("\nFirst failure:", results[[i]]$file, "\nExit:", results[[i]]$status, "\n")
  cat(paste(results[[i]]$log, collapse = "\n"), "\n")
}


# ---- Final step: extract to municipality level ----
# (Assumes community-level processing script is modular)
#for (cVar in c("RSK", "SDK", "NM", "UPM", "TXK", "TNK", "TMK", "TGK", "VPM", "PM ")) {
extract_climate_statistics(
  cVar = param,
  envrmt = envrmt,
  calc_commu = TRUE,
  gemeinden_sf_3035 = st_read(file.path(envrmt$path_data_lev0, "gemeinden_DE_3035.gpkg")),
  dig = 1
)
#}
