#' Create AOI geometry with a buffer in kilometers
#'
#' @param bbox Named numeric vector xmin, xmax, ymin, ymax (WGS84).
#' @param buffer_km Numeric; buffer distance in kilometers.
#'
#' @return An `sfc` polygon in EPSG:4326.
aoi_with_buffer <- function(bbox, buffer_km = 50) {
  # Build AOI polygon in WGS84 from bbox
  aoi <- sf::st_as_sfc(sf::st_bbox(bbox, crs = 4326))
  
  # Transform to metric CRS, apply buffer in meters, transform back to WGS84
  aoi_buf <- aoi |>
    sf::st_transform(3857) |>
    sf::st_buffer(buffer_km * 1000) |>
    sf::st_transform(4326)
  
  aoi_buf
}


#' Get OpenStreetMap (OSM) data for an AOI by key
#'
#' Downloads and clips OSM features for a given area of interest (AOI) in
#' WGS84 (EPSG:4326) for one or more OSM keys (e.g. `"highway"`, `"building"`).
#' Optionally writes the results to GeoPackage files (one per key).
#'
#' @param aoi_wgs84 An `sf` or `sfc` object representing the area of interest
#'   as a (multi-)polygon in CRS EPSG:4326 (WGS84).
#' @param keys Character vector of OSM keys to query (e.g. `"highway"`,
#'   `"landuse"`, `"natural"`, `"waterway"`, `"building"`, `"railway"`).
#' @param out_dir Optional output directory where GeoPackage files will be
#'   written. If `NULL`, no files are written unless `write = FALSE`.
#' @param write Logical; if `TRUE`, write each key's data into a GeoPackage
#'   file in `out_dir`. If `FALSE`, only return the data in memory.
#'
#' @return A named list with one element per key in `keys`. Each element is
#'   itself a list with three components:
#'   \describe{
#'     \item{points}{`sf` object of OSM points clipped to the AOI, or `NULL`.}
#'     \item{lines}{`sf` object of OSM lines clipped to the AOI, or `NULL`.}
#'     \item{polygons}{`sf` object of OSM polygons clipped to the AOI,
#'       or `NULL`.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Example AOI in EPSG:4326
#' aoi <- st_as_sfc(st_bbox(c(xmin = 8.9, xmax = 9.1,
#'                            ymin = 50.9, ymax = 51.1),
#'                          crs = 4326))
#'
#' osm_data <- get_osm_burgwald_by_key(
#'   aoi_wgs84 = aoi,
#'   keys      = c("highway", "building"),
#'   out_dir   = "data/osm",
#'   write     = TRUE
#' )
#' }
get_osm_burgwald_by_key <- function(aoi_wgs84,
                                    keys   = c("highway",
                                               "landuse",
                                               "natural",
                                               "waterway",
                                               "building",
                                               "railway"),
                                    out_dir = NULL,
                                    write   = TRUE) {
  
  # expects: aoi_wgs84 as sf/sfc POLYGON in EPSG:4326
  stopifnot(inherits(aoi_wgs84, c("sf", "sfc")))
  
  # If AOI is an sf object, extract geometry; otherwise use as-is (sfc)
  if (inherits(aoi_wgs84, "sf")) {
    aoi <- sf::st_geometry(aoi_wgs84)
  } else {
    aoi <- aoi_wgs84
  }
  
  # Ensure AOI has a CRS and is in EPSG:4326 (WGS84)
  if (is.na(sf::st_crs(aoi))) stop("AOI has no CRS.")
  if (sf::st_crs(aoi)$epsg != 4326) {
    stop("AOI must be in EPSG:4326 (WGS84).")
  }
  
  # Create output directory if requested and it does not yet exist
  if (!is.null(out_dir) && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Compute bounding box of AOI to use in OSM query
  bbox <- sf::st_bbox(aoi)
  
  # helper: clip sf object to AOI
  clip_to_aoi <- function(x, aoi) {
    # Return NULL if object is empty or NULL
    if (is.null(x) || nrow(x) == 0) return(NULL)
    
    # Ensure input is in WGS84 (4326) before intersection
    x <- sf::st_transform(x, 4326)
    
    # Spatially intersect with AOI (suppressing topology warnings)
    x <- suppressWarnings(sf::st_intersection(x, aoi))
    
    # If intersection yields no features, return NULL
    if (nrow(x) == 0) return(NULL)
    x
  }
  
  # Initialize result list to hold data per key
  res <- list()
  
  # Loop over requested OSM keys
  for (k in keys) {
    message("Fetching OSM for key: ", k)
    
    # Build OSM query for the AOI bbox and given key
    q  <- osmdata::opq(bbox = bbox)
    q  <- osmdata::add_osm_feature(q, key = k)
    
    # Download OSM data as sf objects
    od <- osmdata::osmdata_sf(q)
    
    # Clip points, lines and polygons to AOI
    pts   <- clip_to_aoi(od$osm_points,   aoi)
    lines <- clip_to_aoi(od$osm_lines,    aoi)
    polys <- clip_to_aoi(od$osm_polygons, aoi)
    
    # optional schreiben: ein GPKG pro Key, bis zu drei Layer
    if (write && !is.null(out_dir)) {
      gpkg_file <- file.path(out_dir, paste0("osm_", k, "_burgwald.gpkg"))
      
      # Track whether we've already written the first layer
      # erstes vorhandenes Layer überschreibt Datei
      first_written <- FALSE
      
      # Write points layer (overwrite file on first write)
      if (!is.null(pts)) {
        sf::st_write(pts, gpkg_file,
                     layer      = "points",
                     delete_dsn = TRUE,
                     quiet      = TRUE)
        first_written <- TRUE
      }
      
      # Append lines layer if available
      if (!is.null(lines)) {
        sf::st_write(lines, gpkg_file,
                     layer  = "lines",
                     append = first_written,
                     quiet  = TRUE)
        first_written <- TRUE
      }
      
      # Append polygons layer if available
      if (!is.null(polys)) {
        sf::st_write(polys, gpkg_file,
                     layer  = "polygons",
                     append = first_written,
                     quiet  = TRUE)
      }
      
      message("  -> written: ", gpkg_file)
    }
    
    # Store clipped sf objects in result list under current key
    res[[k]] <- list(
      points   = pts,
      lines    = lines,
      polygons = polys
    )
  }
  
  return(res)
}

# -------------------------------------------------------
#' Download a file only if it does not already exist
#'
#' Convenience wrapper around [download.file()] that first checks if
#' `destfile` exists. If the file is already present, the download is skipped.
#'
#' @param url Character scalar; URL to download from.
#' @param destfile Character scalar; local file path to save the downloaded
#'   file.
#' @param mode Character scalar passed to [download.file()] (e.g. `"wb"` for
#'   binary mode on Windows).
#'
#' @return Invisibly returns `NULL`. Called for its side-effect of downloading
#'   the file (if missing).
#'
#' @examples
#' \dontrun{
#' download_if_missing(
#'   url      = "https://example.com/data.tif",
#'   destfile = "data/raw/data.tif"
#' )
#' }
# -------------------------------------------------------
download_if_missing <- function(url, destfile, mode = "wb") {
  # If destfile does not exist, download it; otherwise skip
  if (!file.exists(destfile)) {
    message("Downloading:\n  ", url, "\n  -> ", destfile)
    download.file(url, destfile = destfile, mode = mode, quiet = FALSE)
  } else {
    message("File already exists, skipping download:\n  ", destfile)
  }
}

# -------------------------------------------------------
#' Run an expression only if target files are missing
#'
#' Evaluates an expression (typically a code block in `{ ... }`) only if
#' at least one of the specified target files does not yet exist. If all
#' target files exist, the expression is skipped.
#'
#' @param targets Character vector of file paths (or a single path) that
#'   should exist after the expression has been run.
#' @param expr An expression containing the work to be done, usually
#'   written as `{ ... }`. The expression is forced only if at least one
#'   target does not exist.
#'
#' @return Invisibly returns `TRUE` if `expr` was executed, `FALSE` if the
#'   step was skipped because all targets already existed.
#'
#' @examples
#' \dontrun{
#' run_if_missing(
#'   targets = c("data/processed.rds", "logs/preprocessing.log"),
#'   expr = {
#'     # expensive computation here
#'     # saveRDS(result, "data/processed.rds")
#'     # writeLines("done", "logs/preprocessing.log")
#'   }
#' )
#' }
# -------------------------------------------------------
run_if_missing <- function(targets, expr) {
  # Ensure targets is a character vector
  targets <- as.character(targets)
  
  # If all target files exist, skip execution
  if (all(file.exists(targets))) {
    message(
      "All target files already exist, skipping step:\n  ",
      paste(targets, collapse = "\n  ")
    )
    return(invisible(FALSE))
  }
  
  # Evaluate expression only if at least one target is missing
  force(expr)
  invisible(TRUE)
}


#' Get hourly DWD station data (precipitation / wind) for Burgwald AOI + buffer
#'
#' This function downloads and reads hourly DWD station data from the CDC
#' archive (`observations_germany/climate/hourly/...`) for either:
#' - `var = "precipitation"`: hourly precipitation (R1 / RR)
#' - `var = "wind"`: hourly wind speed and direction (F/D, FF/DD)
#'
#' It combines both `historical` and `recent` data, filters stations to those
#' inside the Burgwald AOI plus a user-defined buffer, and restricts the time
#' range to `start_date`–`end_date`. Results are returned as an `sf` object and
#' can optionally be written to CSV and RDS.
#'
#' There is an internal mapping from "nice" parameter names (RR, FF, DD) to the
#' actual column names in the DWD product files (R1, F, D). See Details.
#'
#' @details
#' For `var = "precipitation"` (folder `hourly/precipitation`), the product
#' files `produkt_rr_stunde_*.txt` contain, among others:
#' - `R1`: hourly precipitation (mm)
#'
#' Accepted `params`:
#' - `"RR"`: legacy name, mapped internally to `"R1"`
#' - `"R1"`: actual column name in the product file
#'
#' For `var = "wind"` (folder `hourly/wind`), the product files
#' `produkt_ff_stunde_*.txt` contain, among others:
#' - `F`: wind speed (m/s)
#' - `D`: wind direction (degrees)
#'
#' Accepted `params`:
#' - `"FF"`: mapped to `"F"`
#' - `"DD"`: mapped to `"D"`
#' - `"F"` / `"D"`: can be used directly
#'
#' If none of the user-specified `params` can be mapped to a real column name,
#' the function stops with an error.
#'
#' @param var Character, DWD subfolder under `hourly`, e.g. `"precipitation"`
#'   or `"wind"`.
#' @param params Character vector of target columns, e.g. `"RR"` or
#'   `c("FF", "DD")`. These are mapped to the actual column names used in the
#'   DWD product files.
#' @param start_date,end_date Date or character in `"YYYY-MM-DD"` format.
#'   Defines the time window for filtering (start inclusive, end exclusive of
#'   `end_date + 1`).
#' @param aoi_wgs AOI polygon as `sf` or `sfc` in EPSG:4326 (WGS84). By
#'   default this is `aoi_burgwald_wgs` from your setup.
#' @param buffer_km Numeric; buffer distance in kilometers around the AOI to
#'   include additional stations in the surroundings.
#' @param write_csv,write_rds Logical; if `TRUE`, write CSV (no geometry) and
#'   RDS (with geometry) to the project directories.
#'
#' @return An `sf` object with:
#' - station metadata (e.g. station height, name, Bundesland)
#' - column `STATIONS_ID`
#' - column `datetime` (POSIXct, UTC)
#' - one or more climate variables as specified in `params`
#'
#' Additionally, if `write_csv` / `write_rds = TRUE`, files are written under
#' `data/raw/dwd-stations` and `data/processed/dwd-stations`.
#'
#' @examples
#' \dontrun{
#' # 2 years of hourly precipitation (R1 / RR) in Burgwald AOI + 50 km buffer
#' rr_bw <- burgwald_get_hourly_dwd(
#'   var        = "precipitation",
#'   params     = "RR",   # mapped internally to "R1"
#'   start_date = Sys.Date() - 365 * 2,
#'   end_date   = Sys.Date()
#' )
#'
#' # 2 years of hourly wind (speed + direction) for the same AOI + buffer
#' wind_bw <- burgwald_get_hourly_dwd(
#'   var        = "wind",
#'   params     = c("FF", "DD"),  # mapped to c("F", "D")
#'   start_date = Sys.Date() - 365 * 2,
#'   end_date   = Sys.Date()
#' )
#' }
burgwald_get_hourly_dwd <- function(var        = c("precipitation", "wind"),
                                    params     = "R1",
                                    start_date = Sys.Date() - 365 * 2,
                                    end_date   = Sys.Date(),
                                    aoi_wgs    = aoi_burgwald_wgs,
                                    buffer_km  = 50,
                                    write_csv  = TRUE,
                                    write_rds  = TRUE) {
  
  var <- match.arg(var)
  
  if (var == "precipitation") {
    param_map <- c("RR" = "R1", "R1" = "R1")
  } else {
    param_map <- c("FF" = "F", "DD" = "D", "F" = "F", "D" = "D")
  }
  
  params <- unname(param_map[params])
  params <- params[!is.na(params)]
  if (!length(params)) stop("No valid parameters specified for var = '", var, "'.")
  
  raw_dir  <- here::here("data", "raw",       "dwd-stations")
  proc_dir <- here::here("data", "processed", "dwd-stations")
  fs::dir_create(raw_dir,  recurse = TRUE)
  fs::dir_create(proc_dir, recurse = TRUE)
  
  stopifnot(inherits(aoi_wgs, c("sf", "sfc")))
  aoi_geom <- if (inherits(aoi_wgs, "sf")) sf::st_geometry(aoi_wgs) else aoi_wgs
  if (is.na(sf::st_crs(aoi_geom))) stop("AOI has no CRS.")
  if (sf::st_crs(aoi_geom)$epsg != 4326) stop("AOI must be in EPSG:4326 (WGS84).")
  
  aoi_buf <- aoi_geom |>
    sf::st_transform(3857) |>
    sf::st_buffer(buffer_km * 1000) |>
    sf::st_transform(4326)
  
  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)
  
  base_http <- paste0(
    "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/",
    var, "/"
  )
  
  types <- c("historical", "recent")
  
  meta_df  <- NULL
  all_urls <- character(0)
  
  for (ty in types) {
    
    type_url <- paste0(base_http, ty, "/")
    
    file_listing <- rvest::read_html(type_url) |>
      rvest::html_elements("a") |>
      rvest::html_attr("href")
    
    if (is.null(meta_df)) {
      desc_file <- grep("_Stundenwerte_Beschreibung_Stationen", file_listing, value = TRUE)
      if (!length(desc_file)) stop("No station description found for var = '", var, "'.")
      
      desc_url <- paste0(type_url, desc_file[1])
      
      meta_df <- utils::read.fwf(
        desc_url,
        widths = c(6, 9, 9, 15, 12, 10, 41, 41, 4),
        skip = 2,
        col.names = c(
          "Stations_id", "von_datum", "bis_datum", "Stationshoehe",
          "geoBreite", "geoLaenge", "Stationsname", "Bundesland", "Abgabe"
        ),
        strip.white = TRUE,
        fileEncoding = "Latin1"
      )
    }
    
    zip_files <- grep("^stundenwerte_.*\\.zip$", file_listing, value = TRUE)
    if (length(zip_files)) {
      all_urls <- c(all_urls, paste0(type_url, zip_files))
    }
  }
  
  if (is.null(meta_df)) stop("Failed to read station metadata.")
  if (!length(all_urls)) stop("No ZIP URLs found (historical + recent).")
  
  meta_df$geoBreite <- as.numeric(meta_df$geoBreite)
  meta_df$geoLaenge <- as.numeric(meta_df$geoLaenge)
  
  stations_sf <- sf::st_as_sf(meta_df, coords = c("geoLaenge", "geoBreite"), crs = 4326)
  stations_aoi <- suppressWarnings(sf::st_join(stations_sf, sf::st_as_sf(aoi_buf), join = sf::st_within))
  stations_aoi <- stations_aoi[!is.na(stations_aoi$Abgabe), ]
  if (!nrow(stations_aoi)) stop("No DWD stations within AOI+buffer.")
  
  station_ids <- sprintf("%05d", as.integer(stations_aoi$Stations_id))
  
  zip_paths <- file.path(raw_dir, basename(all_urls))
  for (i in seq_along(all_urls)) {
    if (!file.exists(zip_paths[i])) {
      utils::download.file(all_urls[i], zip_paths[i], mode = "wb", quiet = TRUE)
    }
  }
  zip_paths <- zip_paths[file.exists(zip_paths)]
  if (!length(zip_paths)) stop("No ZIP files available after download.")
  
  data_list <- list()
  unzip_dir <- file.path(raw_dir, paste0("unzipped_hourly_", var))
  fs::dir_create(unzip_dir, recurse = TRUE)
  
  for (zp in zip_paths) {
    txt_files <- utils::unzip(zp, exdir = unzip_dir)
    txt_files <- txt_files[grepl("^produkt_.*\\.txt$", basename(txt_files))]
    if (!length(txt_files)) next
    
    for (tf in txt_files) {
      df <- tryCatch(
        utils::read.table(
          tf, sep = ";", header = TRUE,
          stringsAsFactors = FALSE,
          na.strings = c("-999", "-9999", "-999.0", "-9999.0")
        ),
        error = function(e) NULL
      )
      if (is.null(df)) next
      if (!"STATIONS_ID" %in% names(df) || !"MESS_DATUM" %in% names(df)) next
      
      df$STATIONS_ID <- sprintf("%05d", as.integer(df$STATIONS_ID))
      df <- df[df$STATIONS_ID %in% station_ids, ]
      if (!nrow(df)) next
      
      df$datetime <- as.POSIXct(strptime(as.character(df$MESS_DATUM), "%Y%m%d%H", tz = "UTC"))
      df <- df[df$datetime >= as.POSIXct(start_date) & df$datetime < as.POSIXct(end_date + 1), ]
      if (!nrow(df)) next
      
      keep_params <- intersect(params, names(df))
      if (!length(keep_params)) next
      
      df <- df[, c("STATIONS_ID", "MESS_DATUM", "datetime", keep_params), drop = FALSE]
      data_list[[length(data_list) + 1]] <- df
    }
  }
  
  if (!length(data_list)) {
    stop("No data found for params = ", paste(params, collapse = ", "),
         " in the requested time range.")
  }
  
  dat <- data.table::rbindlist(data_list, fill = TRUE)
  
  stations_aoi$STATIONS_ID <- sprintf("%05d", as.integer(stations_aoi$Stations_id))
  stations_used <- stations_aoi[stations_aoi$STATIONS_ID %in% unique(dat$STATIONS_ID), ]
  
  merged_sf <- merge(stations_used, dat, by = "STATIONS_ID")
  
  # --- filenames -----------------------------------------------------------
  stub <- paste0(
    "burgwald_hourly_", var, "_",
    paste(params, collapse = "-"), "_",
    format(start_date, "%Y%m%d"), "_",
    format(end_date, "%Y%m%d")
  )
  
  if (write_csv) {
    csv_file <- file.path(raw_dir, paste0(stub, ".csv"))
    data.table::fwrite(sf::st_drop_geometry(merged_sf), csv_file, dec = ".")
  }
  
  if (write_rds) {
    # canonical if available, else fallback to processed/
    out_key <- if (var == "wind") "dwd_hourly_wind" else "dwd_hourly_precip"
    if (exists("paths", inherits = TRUE) && !is.null(paths[[out_key]])) {
      rds_file <- paths[[out_key]]
      fs::dir_create(dirname(rds_file), recurse = TRUE)
      saveRDS(merged_sf, rds_file)
    } else {
      rds_file <- file.path(proc_dir, paste0(stub, ".rds"))
      saveRDS(merged_sf, rds_file)
    }
  }
  
  merged_sf
}


#' Burgwald: DWD subhourly precipitation (10 min / 5 min) for AOI + buffer
#'
#' Downloads DWD subhourly station data (CDC, observations_germany/climate)
#' for precipitation at 10-minute or 5-minute resolution, combines
#' historical and recent data, filters stations to the Burgwald AOI + buffer,
#' and optionally writes CSV and RDS output.
#'
#' The function currently covers:
#' - resolution = "10min": directory `10_minutes/precipitation`,
#'   main precipitation column is `RWS_10` (10-minute sum).
#' - resolution = "5min": directory `5_minutes/precipitation`,
#'   main precipitation column is `RS_05` (5-minute sum).
#'
#' All timestamps are converted to POSIXct (`datetime`, UTC) and filtered
#' to the requested time window. Only stations whose locations fall inside
#' the AOI + buffer are considered.
#'
#' Data are stored under:
#' - `here("data", "raw", "dwd-stations")`  (CSV and unzipped txt)
#' - `here("data", "processed", "dwd-stations")` (RDS with sf)
#'
#' @param resolution Character, one of `"10min"` or `"5min"`.
#' @param start_date,end_date Date or `"YYYY-MM-DD"`; time window for filtering.
#' @param aoi_wgs AOI polygon (sf or sfc) in EPSG:4326. Default: `aoi_burgwald_wgs`.
#' @param buffer_km Numeric buffer distance (km) around AOI. Default: 50.
#' @param write_csv Logical; if `TRUE`, write a flat CSV (no geometry).
#' @param write_rds Logical; if `TRUE`, write an `sf` object as RDS.
#'
#' @return An `sf` object containing:
#'   - station metadata (id, name, height, etc.)
#'   - geometry of stations (EPSG:4326)
#'   - `MESS_DATUM` and `datetime` (POSIXct, UTC)
#'   - the precipitation column (`RWS_10` or `RS_05`)
burgwald_get_subhourly_precip <- function(resolution = c("10min", "5min"),
                                          start_date = Sys.Date() - 365 * 2,
                                          end_date   = Sys.Date(),
                                          aoi_wgs    = aoi_burgwald_wgs,
                                          buffer_km  = 50,
                                          write_csv  = TRUE,
                                          write_rds  = TRUE) {
  
  resolution <- match.arg(resolution)
  
  # Map resolution to DWD directory and main precip column
  if (resolution == "10min") {
    res_dir     <- "10_minutes"
    precip_col  <- "RWS_10"
    desc_file   <- "zehn_min_rr_Beschreibung_Stationen.txt"
    zip_pattern <- "^10minutenwerte_nieder_.*\\.zip$"
    out_key     <- "dwd_precip_10min"
  } else {
    res_dir     <- "5_minutes"
    precip_col  <- "RS_05"
    desc_file   <- "5min_rr_Beschreibung_Stationen.txt"
    zip_pattern <- "^5minutenwerte_nieder_.*\\.zip$"
    out_key     <- "dwd_precip_5min"
  }
  
  # Base HTTP directory for this resolution and precipitation
  base_http <- paste0(
    "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/",
    res_dir, "/precipitation/"
  )
  
  # RAW provider cache
  raw_dir <- here::here("data", "raw", "providers", "dwd")
  fs::dir_create(raw_dir, recurse = TRUE)
  
  # --- AOI and buffer handling ---------------------------------------------
  stopifnot(inherits(aoi_wgs, c("sf", "sfc")))
  aoi_geom <- if (inherits(aoi_wgs, "sf")) sf::st_geometry(aoi_wgs) else aoi_wgs
  if (is.na(sf::st_crs(aoi_geom))) stop("AOI has no CRS.")
  if (sf::st_crs(aoi_geom)$epsg != 4326) stop("AOI must be in EPSG:4326 (WGS84).")
  
  aoi_buf <- aoi_geom |>
    sf::st_transform(3857) |>
    sf::st_buffer(buffer_km * 1000) |>
    sf::st_transform(4326)
  
  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)
  
  # --- 1) Read station metadata --------------------------------------------
  desc_url <- paste0(base_http, "historical/", desc_file)
  
  meta_df <- utils::read.fwf(
    desc_url,
    widths = c(6, 9, 9, 15, 12, 10, 41, 41, 4),
    skip   = 2,
    col.names = c(
      "Stations_id", "von_datum", "bis_datum", "Stationshoehe",
      "geoBreite", "geoLaenge", "Stationsname", "Bundesland", "Abgabe"
    ),
    strip.white  = TRUE,
    fileEncoding = "Latin1"
  )
  
  meta_df$geoBreite <- as.numeric(meta_df$geoBreite)
  meta_df$geoLaenge <- as.numeric(meta_df$geoLaenge)
  
  stations_sf <- sf::st_as_sf(meta_df, coords = c("geoLaenge", "geoBreite"), crs = 4326)
  
  inside <- sf::st_intersects(stations_sf, aoi_buf, sparse = FALSE)[, 1]
  stations_aoi <- stations_sf[inside, ]
  if (!nrow(stations_aoi)) stop("No stations inside AOI+buffer for this resolution.")
  
  stations_aoi$Stations_id <- sprintf("%05d", as.integer(as.character(stations_aoi$Stations_id)))
  station_ids <- unique(stations_aoi$Stations_id)
  
  # --- 2) Collect ZIP URLs (historical + recent) ---------------------------
  all_urls <- character(0)
  
  if (resolution == "10min") {
    hist_url <- paste0(base_http, "historical/")
    fl_hist  <- rvest::read_html(hist_url) |>
      rvest::html_elements("a") |>
      rvest::html_attr("href")
    zip_hist <- grep(zip_pattern, fl_hist, value = TRUE)
    if (length(zip_hist)) all_urls <- c(all_urls, paste0(hist_url, zip_hist))
  } else {
    hist_root <- paste0(base_http, "historical/")
    fl_root   <- rvest::read_html(hist_root) |>
      rvest::html_elements("a") |>
      rvest::html_attr("href")
    
    year_dirs <- fl_root[grepl("^[0-9]{4}/$", fl_root)]
    if (length(year_dirs)) {
      years_all <- as.integer(sub("/$", "", year_dirs))
      year_min  <- as.integer(format(start_date, "%Y")) - 1L
      year_max  <- as.integer(format(end_date,   "%Y")) + 1L
      years_use <- years_all[years_all >= year_min & years_all <= year_max]
      
      for (y in years_use) {
        y_url <- paste0(hist_root, y, "/")
        fl_y  <- rvest::read_html(y_url) |>
          rvest::html_elements("a") |>
          rvest::html_attr("href")
        zip_y <- grep(zip_pattern, fl_y, value = TRUE)
        if (length(zip_y)) all_urls <- c(all_urls, paste0(y_url, zip_y))
      }
    } else {
      zip_hist <- grep(zip_pattern, fl_root, value = TRUE)
      if (length(zip_hist)) all_urls <- c(all_urls, paste0(hist_root, zip_hist))
    }
  }
  
  recent_url <- paste0(base_http, "recent/")
  fl_recent  <- rvest::read_html(recent_url) |>
    rvest::html_elements("a") |>
    rvest::html_attr("href")
  zip_recent <- grep(zip_pattern, fl_recent, value = TRUE)
  if (length(zip_recent)) all_urls <- c(all_urls, paste0(recent_url, zip_recent))
  
  if (!length(all_urls)) stop("No ZIP URLs found for this resolution (historical + recent).")
  
  zip_station_ids <- sub("^.*_(\\d{5})_.*\\.zip$", "\\1", basename(all_urls))
  keep      <- zip_station_ids %in% station_ids
  urls_keep <- all_urls[keep]
  
  if (!length(urls_keep)) stop("No ZIP files for stations inside AOI+buffer.")
  
  # --- 3) Download ZIPs if missing -----------------------------------------
  zip_paths <- file.path(raw_dir, basename(urls_keep))
  for (i in seq_along(urls_keep)) {
    if (!file.exists(zip_paths[i])) {
      utils::download.file(urls_keep[i], destfile = zip_paths[i], mode = "wb", quiet = FALSE)
    }
  }
  
  zip_paths <- zip_paths[file.exists(zip_paths)]
  if (!length(zip_paths)) stop("No ZIP files available locally.")
  
  # --- 4) Unzip and read produkt_*.txt -------------------------------------
  data_list <- list()
  unzip_dir <- file.path(raw_dir, paste0("unzipped_", resolution, "_precip"))
  fs::dir_create(unzip_dir, recurse = TRUE)
  
  for (zp in zip_paths) {
    
    txt_files <- utils::unzip(zp, exdir = unzip_dir)
    txt_files <- txt_files[grepl("^produkt_.*\\.txt$", basename(txt_files))]
    if (!length(txt_files)) next
    
    for (tf in txt_files) {
      
      df <- tryCatch(
        utils::read.table(
          tf,
          sep = ";",
          header = TRUE,
          stringsAsFactors = FALSE,
          na.strings = c("-999", "-9999", "-999.0", "-9999.0")
        ),
        error = function(e) NULL
      )
      if (is.null(df)) next
      if (!"STATIONS_ID" %in% names(df) || !"MESS_DATUM" %in% names(df)) next
      if (!precip_col %in% names(df)) next
      
      df$STATIONS_ID <- sprintf("%05d", as.integer(as.character(df$STATIONS_ID)))
      df <- df[df$STATIONS_ID %in% station_ids, ]
      if (!nrow(df)) next
      
      ts_chr <- as.character(df$MESS_DATUM)
      ts_chr <- ts_chr[!is.na(ts_chr)]
      if (!length(ts_chr)) next
      
      fmt <- if (nchar(ts_chr[1]) == 12L) "%Y%m%d%H%M" else "%Y%m%d%H"
      
      df$datetime <- as.POSIXct(strptime(as.character(df$MESS_DATUM), fmt, tz = "UTC"))
      
      df <- df[df$datetime >= as.POSIXct(start_date) &
                 df$datetime <  as.POSIXct(end_date + 1), ]
      if (!nrow(df)) next
      
      df <- df[, c("STATIONS_ID", "MESS_DATUM", "datetime", precip_col), drop = FALSE]
      data_list[[length(data_list) + 1]] <- df
    }
  }
  
  if (!length(data_list)) {
    stop("No data for column ", precip_col, " in the requested period for resolution=", resolution, ".")
  }
  
  dat <- data.table::rbindlist(data_list, fill = TRUE)
  
  # --- 5) Merge with station sf --------------------------------------------
  stations_aoi$STATIONS_ID <- sprintf("%05d", as.integer(as.character(stations_aoi$Stations_id)))
  stations_used <- stations_aoi[stations_aoi$STATIONS_ID %in% unique(dat$STATIONS_ID), ]
  
  merged_sf <- merge(stations_used, dat, by = "STATIONS_ID")
  
  # --- 6) Optional writes ---------------------------------------------------
  stub <- paste0(
    "burgwald_", resolution, "_precip_",
    format(start_date, "%Y%m%d"), "_",
    format(end_date,   "%Y%m%d")
  )
  
  if (write_csv) {
    csv_file <- file.path(raw_dir, paste0(stub, ".csv"))
    data.table::fwrite(sf::st_drop_geometry(merged_sf), csv_file, dec = ".")
  }
  
  if (write_rds) {
    rds_file <- paths[[out_key]]
    fs::dir_create(dirname(rds_file), recurse = TRUE)
    saveRDS(merged_sf, rds_file)
  }
  
  merged_sf
}



#' Burgwald: RADOLAN RW (recent hourly radar) for AOI + buffer
#'
#' Downloads and reads recent hourly RADOLAN RW radar composites from
#' the DWD gridbase (`hourly/radolan/recent/bin` via rdwd), converts them
#' to raster, reprojects to lat/lon, clips to the Burgwald AOI + buffer,
#' and optionally writes GeoTIFFs and an RDS stack.
#'
#' This uses `rdwd::indexFTP()` + `rdwd::dataDWD()` + `rdwd::readDWD()` +
#' `rdwd::projectRasterDWD()`, following the examples in the rdwd
#' "Raster data" and `readDWD.radar()` documentation.
#'
#' Note:
#' - This only covers the period available under
#'   `hourly/radolan/recent/bin` on DWD's gridbase.
#' - It is *not* a full 2-year reproc/historical fetcher (that would
#'   need additional logic for the tarballs under reproc/historical).
#'
#' @param start_datetime,end_datetime POSIXct or character (e.g. "2023-01-01 00:00").
#'   Time range in UTC for selecting RADOLAN RW files. Start is inclusive,
#'   end is exclusive of `end_datetime + 1 hour`.
#' @param aoi_wgs AOI polygon (sf/sfc) in EPSG:4326. Default: `aoi_burgwald_wgs`.
#' @param buffer_km Buffer distance (km) around AOI to clip a slightly larger area.
#' @param write_rasters Logical; if TRUE, write cropped GeoTIFFs to
#'   `data/raw/radolan-rw`.
#' @param write_rds Logical; if TRUE, write a SpatRaster stack (terra) to
#'   `data/processed/radolan-rw/burgwald_radolan_rw_stack.rds`.
#'
#' @return A `terra::SpatRaster` with one layer per RADOLAN time step,
#'   cropped to AOI + buffer, and a time vector in `terra::time(x)`.
burgwald_get_radolan_rw_recent <- function(start_datetime,
                                           end_datetime,
                                           aoi_wgs       = aoi_burgwald_wgs,
                                           buffer_km     = 50,
                                           write_rasters = TRUE,
                                           write_rds     = TRUE) {
  
  # ---- Basic sanity checks -------------------------------------------------
  if (missing(start_datetime) || missing(end_datetime)) {
    stop("Please provide start_datetime and end_datetime (POSIXct or character).")
  }
  
  # Normalise to POSIXct (UTC)
  start_datetime <- as.POSIXct(start_datetime, tz = "UTC")
  end_datetime   <- as.POSIXct(end_datetime,   tz = "UTC")
  if (is.na(start_datetime) || is.na(end_datetime)) {
    stop("start_datetime / end_datetime could not be converted to POSIXct.")
  }
  if (end_datetime <= start_datetime) {
    stop("end_datetime must be > start_datetime.")
  }
  
  message("RADOLAN RW time range: ", start_datetime, " – ", end_datetime)
  
  # ---- AOI + buffer in WGS84 ----------------------------------------------
  stopifnot(inherits(aoi_wgs, c("sf", "sfc")))
  if (inherits(aoi_wgs, "sf")) {
    aoi_geom <- sf::st_geometry(aoi_wgs)
  } else {
    aoi_geom <- aoi_wgs
  }
  if (is.na(sf::st_crs(aoi_geom))) stop("AOI has no CRS.")
  if (sf::st_crs(aoi_geom)$epsg != 4326) {
    stop("AOI must be in EPSG:4326 (WGS84).")
  }
  
  # Buffer AOI (metric CRS) and back to 4326
  aoi_buf <- aoi_geom |>
    sf::st_transform(3857) |>
    sf::st_buffer(buffer_km * 1000) |>
    sf::st_transform(4326)
  
  bbox <- sf::st_bbox(aoi_buf)
  message("AOI+buffer bbox (WGS84): ", paste(round(bbox, 3), collapse = ", "))
  
  # ---- RAW provider cache path --------------------------------------------
  raw_dir <- here::here("data", "raw", "providers", "dwd", "radolan")
  
  # ---- 1) List available recent RW files via rdwd -------------------------
  if (!requireNamespace("rdwd", quietly = TRUE)) {
    stop("Package 'rdwd' is required for RADOLAN handling.")
  }
  
  gridbase <- rdwd::gridbase
  
  message("Querying index for hourly/radolan/recent/bin at gridbase ...")
  rw_links <- rdwd::indexFTP(
    "hourly/radolan/recent/bin",
    base = gridbase,
    dir  = raw_dir
  )
  
  if (!length(rw_links)) {
    stop("No RADOLAN RW links found under hourly/radolan/recent/bin.")
  }
  
  fn   <- basename(rw_links)
  ts_s <- sub("^.*-(\\d{10})-.*$", "\\1", fn)   # extract 10-digit time code
  ok   <- grepl("^[0-9]{10}$", ts_s)
  
  rw_links <- rw_links[ok]
  ts_s     <- ts_s[ok]
  
  # Convert to POSIXct: YYMMDDHHMM (UTC)
  ts_posix <- as.POSIXct(strptime(ts_s, "%y%m%d%H%M", tz = "UTC"))
  
  # Filter to requested time range
  keep <- ts_posix >= start_datetime & ts_posix <= end_datetime
  rw_links <- rw_links[keep]
  ts_posix <- ts_posix[keep]
  
  if (!length(rw_links)) {
    stop("No RADOLAN RW files in requested time range (recent).")
  }
  
  # Sort by time
  ord      <- order(ts_posix)
  rw_links <- rw_links[ord]
  ts_posix <- ts_posix[ord]
  
  message("Number of RADOLAN RW files selected: ", length(rw_links))
  
  # ---- 2) Download files (if not present) ---------------------------------
  gz_paths <- file.path(raw_dir, basename(rw_links))
  
  for (i in seq_along(rw_links)) {
    if (!file.exists(gz_paths[i])) {
      message("Download RADOLAN RW: ", basename(gz_paths[i]))
      rdwd::dataDWD(
        rw_links[i],
        base   = gridbase,
        joinbf = TRUE,
        dir    = raw_dir,
        read   = FALSE
      )
    } else {
      message("Already present: ", basename(gz_paths[i]))
    }
  }
  
  gz_paths <- gz_paths[file.exists(gz_paths)]
  if (!length(gz_paths)) {
    stop("No local RADOLAN RW files found after download step.")
  }
  
  # ---- 3) Read, project, crop/mask to AOI buffer --------------------------
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for raster handling.")
  }
  
  rasters_list <- vector("list", length(gz_paths))
  
  for (i in seq_along(gz_paths)) {
    message("Reading RADOLAN RW file: ", basename(gz_paths[i]))
    
    rad <- rdwd::readDWD(gz_paths[i])
    r   <- rad$dat
    
    r_ll <- rdwd::projectRasterDWD(
      r,
      proj     = "radolan",
      extent   = "radolan",
      adjust05 = FALSE
    )
    
    aoi_vect <- terra::vect(aoi_buf)
    r_crop   <- terra::crop(r_ll, aoi_vect)
    r_mask   <- terra::mask(r_crop, aoi_vect)
    
    terra::time(r_mask) <- ts_posix[i]
    
    if (write_rasters) {
      tstamp  <- format(ts_posix[i], "%Y%m%d%H%M")
      tifname <- file.path(raw_dir, paste0("radolan_rw_burgwald_", tstamp, ".tif"))
      terra::writeRaster(
        r_mask,
        tifname,
        overwrite = TRUE,
        filetype  = "GTiff",
        gdal      = c("COMPRESS=DEFLATE", "TILED=YES")
      )
    }
    
    rasters_list[[i]] <- r_mask
  }
  
  # ---- 4) Stack all rasters into one SpatRaster ---------------------------
  r_stack <- terra::rast(rasters_list)
  terra::time(r_stack) <- ts_posix
  
  # ---- 5) Optional: save stack as canonical RDS (S1) -----------------------
  if (write_rds) {
    out_rds <- paths[["radolan_rw_recent"]]
    saveRDS(r_stack, out_rds)
    message("RADOLAN RW stack RDS: ", out_rds)
  }
  
  r_stack
}







# CORINE CLC (Raster mit Werten 1–44 oder 2–41) → mapview mit Textlegende
mapview_clc <- function(r,
                        legend = clc_legend,
                        layer_name = "CORINE Land Cover 2018") {
  
  if (!inherits(r, "SpatRaster")) {
    stop("Input 'r' must be a terra SpatRaster.")
  }
  
  # vorhandene Klassen im Raster
  vals <- sort(unique(terra::values(r)))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) stop("Raster has no non-NA values.")
  
  # nur die Klassen, die im Ausschnitt vorkommen
  leg_sub <- subset(legend, class_id %in% vals)
  if (nrow(leg_sub) == 0) {
    stop("No matching classes in legend for values: ",
         paste(vals, collapse = ", "))
  }
  
  # Raster -> Polygone je Klasse
  r_pol <- terra::as.polygons(r, dissolve = TRUE, values = TRUE, na.rm = TRUE)
  r_sf  <- sf::st_as_sf(r_pol)
  
  # Werte-Spalte auf "class_id" normieren
  value_col <- names(r)[1]               # z.B. "U2018_CLC2018_V2020_20u1"
  names(r_sf)[names(r_sf) == value_col] <- "class_id"
  
  # Legendeninformationen anhängen
  r_sf <- merge(
    r_sf,
    leg_sub[, c("class_id", "name", "color")],
    by = "class_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  # Farbpalette nach Klassenname
  pal <- leg_sub$color
  names(pal) <- leg_sub$name
  
  mapview::mapviewOptions(fgb = FALSE)
  
  mapview::mapview(
    r_sf,
    zcol        = "name",   # kategorial nach CLC-Namen
    col.regions = pal,
    legend      = TRUE,
    layer.name  = layer_name
  )
}

write_nc_if_missing <- function(cube, path) {
  if (file.exists(path)) {
    message("NetCDF existiert bereits, überspringe: ", path)
    return(invisible(FALSE))
  }
  gdalcubes::write_ncdf(cube, path, overwrite = TRUE)
  invisible(TRUE)
}
