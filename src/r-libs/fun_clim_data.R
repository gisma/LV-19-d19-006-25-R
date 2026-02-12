# get_cl_data <- function(l)
# {
# 
#   urls <- selectDWD(id = lki[l], res="daily", var="kl", per=type)#, outvec=TRUE)
#   clims <- dataDWD(urls, varnames=FALSE, dir=envrmt$path_GhcnDaily)
#   climdata <- clims[c("STATIONS_ID","MESS_DATUM","RSK","SDK","NM","VPM","PM","TMK","UPM","TXK","TNK","TGK")]
#   return(climdata)
# 
# }
# 
# get_climdata <- function(l,var="pressure",res ="hourly",col="P")
# {
#   
#   urls <- selectDWD(id = lki[l], res=res, var=var, per=type)#, outvec=TRUE)
#   clims <- dataDWD(urls, varnames=FALSE, dir=envrmt$path_GhcnDaily)
#   climdata <- clims[c("STATIONS_ID","MESS_DATUM",col)]
#   return(climdata)
#   
# }

#'
#' Partly the code is taken from
#' https://bookdown.org/brry/rdwd/use-case-get-all-hourly-rainfall-data-20142016.html
#' 
read_cl_param <- function(file, fread=TRUE,sd=NULL,ed = NULL, par = NULL, ...)
{
  climdata <- readDWD(file, fread=fread)
  climdata <- climdata[climdata$MESS_DATUM > as.POSIXct(as.Date(sd)) & 
                         climdata$MESS_DATUM < as.POSIXct(as.Date(ed))    , ]
  climdata <- climdata[ , c("STATIONS_ID","MESS_DATUM", par)]
  climdata$MESS_DATUM <- as.Date(climdata$MESS_DATUM) # might save some memory space...
  return(climdata)
}
#'
#' Partly the code is taken from
#' https://bookdown.org/brry/rdwd/use-case-get-all-hourly-rainfall-data-20142016.html
#' 
#' 
#' 
ex_clim = function(startDate=NULL,endDate=NULL,reso=NULL,var=NULL,type=NULL,param=NULL){
  message("::: get climate data :::")
  # Select daily climate data:
  data("metaIndex")
  myIndex <- metaIndex[
    metaIndex$von_datum < as.Date(startDate) &
      metaIndex$bis_datum > as.Date(endDate) & metaIndex$hasfile   ,  ]
  data("fileIndex")
  links <- fileIndex[
    suppressWarnings(as.numeric(fileIndex$id)) %in% myIndex$Stations_id &
      fileIndex$res==reso &
      fileIndex$var==var &
      fileIndex$per==type         , "path" ]
  localfiles <- dataDWD(
    links,
    joinbf = TRUE,
    sleep = 0.2,
    read = FALSE,
    dir = envrmt$path_CDC_KL
  )
  localfiles <- dataDWD(links, joinbf=TRUE, sleep=0.2, read=FALSE,dir=envrmt$path_CDC_KL)
  localfiles = localfiles[file.exists(localfiles)]
  
  # extract station ids 
  matrix_of_params <- parallel::mclapply( seq_along(localfiles), function(n){
    #for (n in seq_along(1:3)){
    read_cl_param(localfiles[n],sd = startDate,ed=endDate,par = param,)
  }, mc.cores = 16, mc.allow.recursive = TRUE)
  var_all = data.table::rbindlist(matrix_of_params) 
  
  # Transform into spatial object:
  msf <- sf::st_as_sf(myIndex, coords=c("geoLaenge", "geoBreite"), crs=4326)
  stations = msf[msf$Stations_id %in%  unique(var_all$STATIONS_ID),]
  stations = stations[stations$res == reso & stations$var==var & stations$per == type & stations$hasfile==TRUE   , ]
  names(stations)[1] = "STATIONS_ID"
  
  merge_cl = merge(stations,var_all)
  # transform to crs
  cVar.sf <- st_transform(merge_cl, crs)
  saveRDS(cVar.sf,paste0(envrmt$path_data_lev0,"/daily_climate_",type,"_",param,".rds"))  

  gc()
  
  return(cVar.sf )
  # actually this means to extract hourly data in this case PM 
  # it is hard wired so far, have a look at get_climdata
  #   if (PM) {
  #     pressureLK <- pbapply::pblapply(1:length(lki), get_climdata)
  #     pressureLK_all = do.call(rbind,pressureLK)
  #     merge_PM = merge(stations,pressureLK_all)
  #     cVar_PM.sf <- st_transform(merge_PM, crs)
  #     saveRDS(cVar_PM.sf,paste0(envrmt$path_data_lev0,"/hourly_PM_",type,".rds"))
  #      # calculate daily mean from hourly data  
  #     cVar_PM.sf$Date <- as.Date(cVar_PM.sf$MESS_DATUM, format = "%Y%m%d%H")
  #     # Group by STATIONS_ID and Date, then calculate mean of P
  #     cVar_PM.sf <- cVar_PM.sf %>%
  #       group_by(STATIONS_ID,Stationshoehe,Date) %>%
  #       summarise(PM = mean(P))
  #     names(cVar_PM.sf) = c("STATIONS_ID", "Stationshoehe", "MESS_DATUM", "PM","geometry")
  #     saveRDS(daily_means,paste0(envrmt$path_data_lev0,"/daily_PM_",type,".rds"))
  #     
  # }
}


#' Modern DWD Climate Data Extractor (robuste Version mit Cache, Filter und Export)
#'
#' @param startDate Startdatum (z. B. "2010-01-01")
#' @param endDate Enddatum (z. B. "2020-12-31")
#' @param reso Auflösung (z. B. "daily")
#' @param var DWD Variable (z. B. "kl")
#' @param type Zeitraumtyp (z. B. "historical")
#' @param param Zielparameter (z. B. "TMK")
#' @return sf-Objekt mit Messwerten
#' @export
ex_clim_new <- function(startDate, endDate, reso = "daily", var = "kl", type = "historical", param = NULL) {

  
  message("::: Suche verfügbare DWD-Stationen :::")
  outname <- paste0(envrmt$path_CDC_KL, "/", gsub("-", "", startDate), "_", gsub("-", "", endDate), "_", param, ".gpkg")
  if (!file.exists(outname)) {
    
  # Zielverzeichnis sicherstellen
  if (!dir.exists(envrmt$path_CDC_KL)) dir.create(envrmt$path_CDC_KL, recursive = TRUE)
  
  # Pfad zur gecachten Datei
  filename <- paste0(gsub("-", "", startDate), "_", gsub("-", "", endDate), "_", param, ".gpkg")
  result_path <- file.path(envrmt$path_CDC_KL, filename)
  if (file.exists(result_path)) {
    message("✔️ Vorhandene Datei gefunden: ", filename)
    return(st_read(result_path, quiet = TRUE))
  }
  
  # --- ab hier unveränderter Code ---
  
  # Basis-URLs
  base_url_http <- paste0("https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/", reso, "/", var, "/", type, "/")
  base_url <- paste0("ftp://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/", reso, "/", var, "/", type, "/")
  
  # Datei-Listing
  file_listing <- read_html(base_url_http) |>
    html_elements("a") |>
    html_attr("href")
  
  zip_files <- grep("^tageswerte_.*\\.zip$", file_listing, value = TRUE)
  if (length(zip_files) == 0) stop("⚠️ Keine ZIP-Dateien gefunden.")
  
  # Metadaten laden
  station_file <- grep("Beschreibung_Stationen", file_listing, value = TRUE)
  if (length(station_file) == 0) stop("⚠️ Keine Stationsbeschreibung gefunden.")
  
  station_url <- paste0(base_url, station_file)
  meta_df <- read.fwf(
    station_url,
    widths = c(6, 9, 9, 15, 12, 10, 41, 41, 4),
    skip = 2,
    col.names = c("Stations_id", "von_datum", "bis_datum", "Stationshoehe",
                  "geoBreite", "geoLaenge", "Stationsname", "Bundesland", "Abgabe"),
    strip.white = TRUE,
    fileEncoding = "Latin1"
  )
  
  # ZIP-Dateien herunterladen und entpacken
  for (zf in zip_files) {
    zip_url <- paste0(base_url, zf)
    zip_path <- file.path(envrmt$path_CDC_KL, basename(zf))
    
    if (!file.exists(zip_path)) {
      download.file(zip_url, destfile = zip_path, mode = "wb", quiet = TRUE)
    } else {
      message(paste("✔️ Vorhanden:", basename(zf)))
    }
    
    exdir <- file.path(envrmt$path_CDC_KL, tools::file_path_sans_ext(basename(zf)))
    if (!dir.exists(exdir)) dir.create(exdir, recursive = TRUE)
    unzip(zip_path, exdir = exdir)
  }
  
  # 🔁 produkt_klima_tag-Dateien aus allen Unterordnern
  unzipped_dirs <- list.dirs(envrmt$path_CDC_KL, recursive = TRUE, full.names = TRUE)
  produkt_files <- list.files(
    path = unzipped_dirs,
    pattern = "^produkt_klima_tag_.*\\.txt$",
    full.names = TRUE,
    recursive = FALSE
  )
  
  message("📂 Lese und filtere alle produkt_klima_tag_*.txt-Dateien – das kann einige Minuten dauern ...")
  
  
  data_list <- list()
  for (txt_file_full in produkt_files) {
    df <- tryCatch(
      read.table(txt_file_full, sep = ";", header = TRUE, na.strings = c("-999", "-999.0")),
      error = function(e) return(NULL)
    )
    
    if (!is.null(df) && param %in% names(df)) {
      df$MESS_DATUM <- as.Date(as.character(df$MESS_DATUM), format = "%Y%m%d")
      df <- df[df$MESS_DATUM >= as.Date(startDate) & df$MESS_DATUM <= as.Date(endDate), ]
      df <- df[, c("STATIONS_ID", "MESS_DATUM", param)]
      data_list[[length(data_list) + 1]] <- df
    }
  }
  
  var_df <- rbindlist(data_list)
  if (nrow(var_df) == 0) stop("⚠️ Kein Parameterwert in gewünschtem Zeitraum gefunden.")
  
  # Spatial Join
  # Spatial Join
  meta_df$geoBreite <- as.numeric(meta_df$geoBreite)
  meta_df$geoLaenge <- as.numeric(meta_df$geoLaenge)
  
  # Vereinheitliche Namen und Typen
  names(meta_df)[names(meta_df) == "Stations_id"] <- "STATIONS_ID"
  meta_df$STATIONS_ID <- as.numeric(meta_df$STATIONS_ID)
  var_df$STATIONS_ID <- as.numeric(var_df$STATIONS_ID)
  
  # Spatial DataFrame
  stations_sf <- st_as_sf(meta_df, coords = c("geoLaenge", "geoBreite"), crs = 4326)
  
  # Merge Klimawerte mit Metadaten
  merge_sf <- merge(stations_sf,var_df, by = "STATIONS_ID")
  
  message("::: Entferne unvollständige Einträge :::")
  n_before <- nrow(merge_sf)

  merge_sf <- merge_sf[complete.cases(st_drop_geometry(merge_sf)), ]
  n_after <- nrow(merge_sf)
  message(paste("⚠️ Gefiltert:", n_before - n_after, "Einträge mit NA"))
  message(paste("✅ Verbleibende Einträge:", n_after))
  message(paste("📡 Anzahl eindeutiger Stationen:", length(unique(merge_sf$STATIONS_ID))))
  
  # Speichern, wenn noch nicht vorhanden
  outname <- paste0(envrmt$path_CDC_KL, "/", gsub("-", "", startDate), "_", gsub("-", "", endDate), "_", param, ".gpkg")
  if (!file.exists(outname)) {
    message(paste("💾 Datei wird gespeichert unter:", outname))
    st_write(merge_sf, outname, delete_dsn = TRUE, quiet = TRUE)
  } else {
    message(paste("📁 Datei existiert bereits:", outname))
  }
  
  return(st_as_sf(merge_sf))
  }
  merge_sf = st_read(outname)
  message(paste("📁 Existierende Datei \n", outname, " \n wird verwendet"))
  return(merge_sf)
}

#' Klassifiziere Kriging-Unsicherheit in Qualitätsstufen und speichere Raster
#'
#' @param stars_obj stars-Objekt mit mindestens dem Layer "var1.var"
#' @param export_path Zeichenkette: vollständiger Dateipfad für das Output-Raster (optional)
#' @param compress Option für GeoTIFF-Kompression (default: "LZW")
#' @return stars-Objekt mit Qualitätsklassen als Faktor
#' @export
classify_kriging_uncertainty <- function(stars_obj, export_path = NULL, compress = "LZW") {
  require(stars)
  
  # Prüfen, ob Layer "var1.var" vorhanden ist
  if (!"var1.var" %in% names(stars_obj)) {
    stop("Layer 'var1.var' fehlt im stars-Objekt")
  }
  
  # Berechne die Standardabweichung σ = sqrt(var)
  sigma_raster <- sqrt(stars_obj[["var1.var"]])
  
  # Klassifiziere Unsicherheiten in vier Kategorien:
  # hoch:        σ < 1
  # mittel:      1 ≤ σ < 2
  # niedrig:     2 ≤ σ < 4
  # sehr niedrig: σ ≥ 4
  quality_class <- cut(
    sigma_raster,
    breaks = c(0, 1, 2, 4, Inf),
    labels = c("hoch", "mittel", "niedrig", "sehr niedrig"),
    right = FALSE
  )
  
  # Wandle in stars-Objekt um
  quality_raster <- st_as_stars(quality_class)
  names(quality_raster) <- "kriging_quality"
  
  # Optional speichern
  if (!is.null(export_path)) {
    stars::write_stars(quality_raster, export_path, overwrite = TRUE, options = paste0("COMPRESS=", compress))
  }
  
  return(quality_raster)
}
#' Bereinigt und begrenzt Klimadaten für einen gegebenen Parameter
#'
#' @param sf_day sf-Objekt mit Tageswerten eines Parameters
#' @param cVar Zeichenkette des DWD-KL-Parameters (z. B. "TXK", "RSK", ...)
#' @param date Datum im Format "YYYY-MM-DD" (nur erforderlich für SDK)
#' @return sf-Objekt mit neuer Spalte `tmp`, die die bereinigten Werte enthält
sanitize_climate_param <- function(sf_day, cVar, date = NULL) {
  x <- as.numeric(sf_day[[cVar]])
  
  # Fehldatenwert -999 zu NA
  x[x == -999] <- NA
  
  if (cVar == "SDK" && !is.null(date)) {
    # Sonnenauf- und untergang → maximale Tageslänge (in Stunden)
    dt <- suncalc::getSunlightTimes(date = as.Date(date), lat = 51.0, lon = 9.0, tz = "UTC")
    td <- dt$sunset - dt$sunrise
    maxDaylight <- ceiling(as.numeric(td))
    x[x > maxDaylight] <- maxDaylight
  }
  
  if (cVar == "PM") {
    x[x > 1060.6] <- 1060.6
    #x[x < 954.9]  <- 954.9
    x[x < 500.0]  <- 500.0
  }
  
  if (cVar == "UPM") {
    x[x > 100] <- 100
    x[x < 0]   <- 0
  }
  
  if (cVar %in% c("TXK", "TNK", "TMK")) {
    x[x > 42]    <- 42
    x[x < -46.0] <- -46.0
  }
  
  if (cVar == "NM") {
    x[x > 8.0] <- 8.0
    x[x < 0]   <- 0
  }
  
  if (cVar == "RSK") {
    x[x > 312] <- 312
    x[x < 0]   <- 0
  }
  
  sf_day$tmp <- x
  return(sf_day)
}

#' @title Prepare German Elevation and Boundary Data
#'
#' @description
#' Downloads and processes DEM and administrative boundary data for Germany.
#' The result is a ready-to-use elevation raster (`dem`) for kriging or spatial modeling, saved as `dem.rds`.
#'
#' @param envrmt A named list of project paths, e.g. from `envimaR::createEnv()`.
#' @param bl Character: name of the federal state to extract (e.g. "Hessen").
#' @param crs Coordinate reference system (e.g. EPSG code or `terra::crs` object).
#' @param res Numeric: target grid resolution in meters (e.g. 500).
#' @param downloadDEM Logical: if `TRUE`, downloads and processes SRTM DEM; if `FALSE`, loads existing file.
#'
#' @return A `stars` object representing the processed DEM.
#' @export
#'
#' @author Chris Reudenbach
#' @examples
#' \dontrun{
#' envrmt <- createEnv()
#' dem <- prepare_DEM(envrmt, bl = "Hessen", crs = 3035, res = 500, downloadDEM = TRUE)
#' }

prepare_DEM <- function(envrmt, bl, crs, res = 500, downloadDEM = TRUE) {
  message("::: Load boundary data via geodata :::")

  
  # Load administrative boundaries
  germany <- geodata::gadm(country = "DEU", level = 1, path = tempdir())
  germany.sf <- st_as_sf(germany)
  germany.sf <- st_transform(germany.sf, crs = crs)
  
  # Filter selected state
  bl_sp <- germany[germany$NAME_1 == bl, ]
  bl_sf <- st_as_sf(bl_sp)
  bl_sf <- st_transform(bl_sf, crs = crs)
  
  # Dissolve all DE states into a unified polygon
  states_special <- c("Baden-Württemberg", "Nordrhein-Westfalen", "Hessen", "Bayern",
                      "Niedersachsen", "Sachsen-Anhalt", "Rheinland-Pfalz", "Sachsen",
                      "Mecklenburg-Vorpommern", "Schleswig-Holstein", "Brandenburg",
                      "Thüringen", "Saarland", "Berlin", "Hamburg", "Bremen")
  DE.states <- germany.sf[germany.sf$NAME_1 %in% states_special, ]
  DE <- DE.states %>% group_by(NAME_1) %>% summarize()
  
  if (downloadDEM) {
    message("::: Download and prepare SRTM elevation data :::")
    de_4326 <- st_transform(DE, 4326)
    st_write(de_4326, file.path(envrmt$path_data_lev0, "de_4326.shp"), delete_dsn = TRUE)
    
    download.url <- "https://opendem.info/downloads/srtm_germany_dtm.zip"
    zipfile <- file.path(envrmt$path_data_lev0, "srtm_germany_dtm.zip")
    download.file(download.url, zipfile, mode = "wb")
    unzip(zipfile, exdir = envrmt$path_data_lev0)
    
    srtm.germany <- terra::mask(
      terra::rast(file.path(envrmt$path_data_lev0, "srtm_germany_dtm.tif")),
      de_4326
    )
    
    message("::: Create template raster :::")
    grid.DE <- expand.grid(
      x = seq(round(st_bbox(germany.sf)["xmin"]), round(st_bbox(germany.sf)["xmax"]), by = res),
      y = seq(round(st_bbox(germany.sf)["ymin"]), round(st_bbox(germany.sf)["ymax"]), by = res)
    )
    coordinates(grid.DE) <- ~x + y
    crs(grid.DE) <- crs
    template_raster <- rasterFromXYZ(grid.DE, crs = crs)
    
    srtm.germany <- terra::project(srtm.germany, crs(template_raster))
    srtm500 <- terra::resample(srtm.germany, rast(template_raster))
    srtm500 <- round(srtm500, 0)
    names(srtm500) <- "Stationshoehe"
    
    dem <- st_as_stars(srtm500)
    saveRDS(dem, file.path(envrmt$path_data_lev0, "dem.rds"))
    rm(srtm.germany, srtm500, template_raster, grid.DE)
  }
  
  dem <- readRDS(file.path(envrmt$path_data_lev0, "dem.rds"))
  return(dem)
}

# ---- Define reusable classification logic ----
classify_raster_by_cVar <- function(cVar, rast) {
  ranges <- list(
    #PM  = c(0, 954.9, 954.9, 1060.6, 99999, 1060.6),
    UPM = c(-999, 0, 0, 100, 99999, 100),
    TXK = c(-999, -46, -46, 42, 99999, 42),
    TMK = c(-999, -46, -46, 42, 99999, 42),
    TNK = c(-999, -46, -46, 42, 99999, 42),
    NM  = c(-1000, 0, 0, 8, 99999, 8),
    RSK = c(-1000, 0, 0, 312, 99999, 312)
  )
  if (cVar %in% names(ranges)) {
    m <- matrix(ranges[[cVar]], ncol = 3, byrow = TRUE)
    return(terra::classify(rast, m, include.lowest = TRUE))
  } else {
    return(rast)
  }
}


# ---- CSV export helper ----
write_stat_csv <- function(sf_df, file,dig) {
  df <- st_drop_geometry(sf_df)
  numeric_cols <- sapply(df, is.numeric)
  df[numeric_cols] <- lapply(df[numeric_cols], round, digits = dig)
  var_fin <- sjmisc::replace_columns(sf_df, df)
  data.table::fwrite(st_drop_geometry(var_fin), file = file, dec = ".")
}

correct_daytime = function(fn){
  for (f in fn){
    current = terra::rast(f)
    dt=suncalc::getSunlightTimes(date = as.Date(substr(basename(f),1,10)), lat = 51.0, lon = 9.0, tz = "UTC")
    td=dt$sunset-dt$sunrise
    maxDaylight= ceiling(as.numeric(unlist(stringr::str_split(td,"Time difference of "))))
    m <- c(-999, 0,0, maxDaylight,99999,maxDaylight)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    current <- terra::classify(current, rclmat, include.lowest=TRUE)
    names(current) = xfun::sans_ext(basename(f))
    if (calc_bl) 
    {  
      # Calculate data frame of min and max precipitation for all months
      var <- cbind(bl_sf, exactextractr::exact_extract(terra::rast(f), bl_sf, c("min", "max","count","majority","median","quantile","minority","variance","stdev","coefficient_of_variation"),quantiles = c(0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.25,0.5,0.75,0.9)))
      var$date = substr(tools::file_path_sans_ext(basename(f)),1,10)
      #c("min",	"max",	"count",	"majority",	"median",	"q10",	"q20",	"q30",	"q40",	"q60",	"q70",	"q80",	"q25",	"q50",	"q75",	"q90",	"minority","variance","stdev","coefficient_of_variation")
      vr=st_drop_geometry(var[,c("min",	"max",	"count",	"majority",	"median",	"q10",	"q20",	"q30",	"q40",	"q60",	"q70",	"q80",	"q25",	"q50",	"q75",	"q90",	"minority")]) %>%
        mutate_if(is.numeric, round, digits=dig)
      var_fin=sjmisc::replace_columns(var,vr)
      #saveRDS(var,file.path(envrmt$path_data_lev2,cVar,paste0(tools::file_path_sans_ext(basename(clim_files[i])),".rds")))
      data.table::fwrite(st_drop_geometry(var_fin),file=file.path(envrmt$path_data_lev2,bl,cVar,paste0(tools::file_path_sans_ext(basename(clim_files[i])),".csv")),dec = ".")
      current = terra::mask(raster::raster(f), bl_sf)
      current = terra::crop(current,bl_sf)
      raster::writeRaster(current,file.path(envrmt$path_data_lev1,bl,cVar,paste0(bl,basename(f))),overwrite=TRUE)
    }  else {
      writeRaster(current, f,gdal=c("COMPRESS=NONE"),overwrite=TRUE)
    }
  }
}

#' @title Extraction of Climate Raster Statistics by Community or Federal State
#'
#' @description
#' Performs zonal statistics on daily interpolated climate raster data per administrative unit
#' (either municipalities or federal states). Includes:
#' \itemize{
#'   \item Classification of physically implausible values (value capping)
#'   \item Efficient extraction of descriptive statistics using \code{exactextractr}
#'   \item Optional masking and clipping of raster data to federal states
#'   \item Export of individual \code{.csv} statistics per day and merged outputs
#' }
#'
#' Raster values are processed for a given climate variable (e.g., temperature, precipitation).
#' The workflow supports both municipality-level and state-level statistics, controlled via
#' \code{calc_commu} and \code{calc_bl}.
#'
#' @param cVar [character] name of the climate variable (e.g., "TMK", "RSK")
#' @param envrmt [list] project environment list with standardized paths
#' @param calc_commu [logical] whether to compute municipality-level stats
#' @param calc_bl [logical] whether to compute federal state (BL) stats
#' @param bl [character] short name of Bundesland (e.g., "Hessen")
#' @param gemeinden_sf_3035 [sf] polygon layer of all municipalities in EPSG:3035
#' @param bl_sf [sf] polygon layer of one federal state (if \code{calc_bl = TRUE})
#' @param common_quantiles [numeric] quantiles used for descriptive statistics
#' @param dig [integer] numeric rounding precision for each variable
#'
#' @return
#' Creates daily `.csv` files per administrative unit containing statistical summaries.
#' Also exports clipped raster files for BL-level analysis and concatenated `.out` tables
#' across all days.
#'
#' @author
#' Chris Reudenbach, \email{creuden@@gmail.com}
extract_climate_statistics <- function(cVar,
                                       envrmt,
                                       calc_commu = TRUE,
                                       calc_bl = FALSE,
                                       bl = NULL,
                                       gemeinden_sf_3035,
                                       bl_sf = NULL,
                                       common_quantiles = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.25, 0.5, 0.75, 0.9),
                                       dig = 1) {
  
  # Read climate files
  clim_files <- sort(list.files(file.path(envrmt$path_data_lev1, cVar), pattern = paste0("*", cVar, "\\.tif$"), full.names = TRUE))
  
  # Create output directories
  dir.create(file.path(envrmt$path_data_lev2, cVar), recursive = TRUE, showWarnings = FALSE)
  if (calc_bl) {
    dir.create(file.path(envrmt$path_data_lev1, bl, cVar), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(envrmt$path_data_lev2, bl, cVar), recursive = TRUE, showWarnings = FALSE)
  }
  
  # Loop over climate raster files
  parallel::mclapply(seq_along(clim_files), function(i) {
   # for (i in seq_along(clim_files)) {   
    outfile_commu <- file.path(envrmt$path_data_lev2, cVar, paste0(tools::file_path_sans_ext(basename(clim_files[i])), ".csv"))
    #if (file.exists(outfile_commu)) return(NULL)
    if (calc_commu && file.exists(outfile_commu)) return(NULL)
    if (cVar == "SDK") {
      correct_daytime(fn = clim_files[i])
     # return(NULL)
    }
    
    current <- terra::rast(clim_files[i])
    current <- classify_raster_by_cVar(cVar, current)
    names(current) <- xfun::sans_ext(basename(clim_files[i]))
    
    # --- HARD NoData handling: turn numeric fill values into real NA ---
    nd <- suppressWarnings(terra::NAflag(current))
    if (length(nd) == 0) nd <- NA
    nd <- nd[1]
    
    # if NAflag exists: map that exact value to NA
    if (is.finite(nd)) {
      current[current == nd] <- NA
    }
    
    # your explicit writeRaster NAflag
    current[current == -9999] <- NA
    dig <- ifelse(cVar %in% c("PM", "TXK", "TMK", "TNK", "NM", "RSK"), 1, 0)
    
    if (calc_bl) {
      stat_vars <- cbind(
        bl_sf,
        exactextractr::exact_extract(current, bl_sf,
                                     c("min", "max", "count", "majority", "median", "quantile", "minority", "variance", "stdev", "coefficient_of_variation"),
                                     quantiles = common_quantiles)
      )
      stat_vars$date <- substr(tools::file_path_sans_ext(basename(clim_files[i])), 1, 10)
      
      out_bl_file <- file.path(envrmt$path_data_lev2, bl, cVar, paste0(tools::file_path_sans_ext(basename(clim_files[i])), ".csv"))
      write_stat_csv(stat_vars, out_bl_file,dig)
      
      masked <- terra::mask(current, bl_sf)
      cropped <- terra::crop(masked, bl_sf)
      raster::writeRaster(cropped, file.path(envrmt$path_data_lev1, bl, cVar, paste0(bl, basename(clim_files[i]))), overwrite = TRUE)
    }
    
    if (calc_commu) {
      stat_vars <- cbind(
        gemeinden_sf_3035,
        exactextractr::exact_extract(current, gemeinden_sf_3035,
                                     c("min", "max", "count", "majority", "median", "quantile", "minority", "variance", "stdev", "coefficient_of_variation"),
                                     quantiles = common_quantiles)
      )
      stat_vars$date <- substr(tools::file_path_sans_ext(basename(clim_files[i])), 1, 10)
      write_stat_csv(stat_vars, outfile_commu,dig)
    }
    
  }, mc.cores = 16, mc.allow.recursive = TRUE)
  #} 
  # Merge CSV outputs
  if (calc_bl) {
    system(paste0("head -n 1 ", envrmt$path_data_lev2, "/", bl, cVar, "/2003-01-01_", cVar, ".csv > ",
                  envrmt$path_data_lev2, "/", bl, cVar, "/", bl, cVar, "_2003-2021.out && tail -n+2 -q ",
                  envrmt$path_data_lev2, "/", bl, cVar, "/*", cVar, ".csv >> ",
                  envrmt$path_data_lev2, "/", bl, cVar, "/", bl, cVar, "_2003-2021.out"), intern = FALSE)
  }
  if (calc_commu) {
    # system(paste0("head -n 1 ", envrmt$path_data_lev2, "/", cVar, "/2003-01-01_", cVar, ".csv > ",
    #               envrmt$path_data_lev2, "/", cVar, "/", cVar, "_2003-2021.out && tail -n+2 -q ",
    #               envrmt$path_data_lev2, "/", cVar, "/*", cVar, ".csv >> ",
    #               envrmt$path_data_lev2, "/", cVar, "/", cVar, "_2003-2021.out"), intern = FALSE)
    # Dynamisch Start- und Enddatum aus Dateinamen ableiten
    dates <- substr(tools::file_path_sans_ext(basename(clim_files)), 1, 10)
    start_date <- min(dates)
    end_date <- max(dates)
    
    # Ausgabedateinamen generieren
    merged_file_bl <- file.path(envrmt$path_data_lev2, bl, cVar, paste0(bl, cVar, "_", start_date, "-", end_date, ".out"))
    merged_file_commu <- file.path(envrmt$path_data_lev2, cVar, paste0(cVar, "_", start_date, "-", end_date, ".out"))
    
    # Zusammenführen für BL
    if (calc_bl) {
      first_file_bl <- file.path(envrmt$path_data_lev2, bl, cVar, paste0(start_date, "_", cVar, ".csv"))
      system(paste0("head -n 1 ", first_file_bl, " > ", merged_file_bl, 
                    " && tail -n+2 -q ", envrmt$path_data_lev2, "/", bl, cVar, "/*", cVar, ".csv >> ", merged_file_bl))
    }
    
    # Zusammenführen für Gemeinden
    if (calc_commu) {
      first_file_commu <- file.path(envrmt$path_data_lev2, cVar, paste0(start_date, "_", cVar, ".csv"))
      system(paste0("head -n 1 ", first_file_commu, " > ", merged_file_commu, 
                    " && tail -n+2 -q ", envrmt$path_data_lev2, "/", cVar, "/*", cVar, ".csv >> ", merged_file_commu))
    }
    
  }
}

###### missing days

# --- Fehlende Tage einlesen & Auswahl treffen -------------------------------
# benötigt: stringr, readr, dplyr, tibble, lubridate
req <- c("stringr","readr","dplyr","tibble","lubridate")
invisible(lapply(req, require, character.only = TRUE))

read_missing_days <- function(path_txt) {
  readr::read_lines(path_txt) |>
    stringr::str_trim() |>
    (\(x) x[x != ""])() |>                     # leere Zeilen raus
    (\(x) tibble::tibble(line = x))() |>
    tidyr::extract(
      line, into = c("var","date_raw"),
      regex = "^\\s*([A-Z]+)\\s*:\\s*(\\d{4}-\\d{1,2}-\\d{1,2})\\s*$",
      remove = TRUE
    ) |>
    dplyr::mutate(
      var  = stringr::str_trim(var),
      date = lubridate::ymd(date_raw, quiet = TRUE)
    ) |>
    dplyr::filter(!is.na(var), !is.na(date)) |>
    dplyr::distinct(var, date, .keep_all = TRUE) |>
    dplyr::arrange(var, date) |>
    dplyr::select(var, date)
}

# Interaktive Auswahl (mehrfach möglich) – alternativ unten "non-interactive" nutzen
select_missing <- function(miss_df,
                           allowed_vars = c("RSK","SDK","NM","UPM","TXK","TNK","TMK","TGK","VPM","PM"),
                           startDate = NULL, endDate = NULL) {
  df <- miss_df |> dplyr::filter(var %in% allowed_vars)
  if (!is.null(startDate)) df <- df |> dplyr::filter(date >= as.Date(startDate))
  if (!is.null(endDate))   df <- df |> dplyr::filter(date <= as.Date(endDate))
  
  # 1) Variablen wählen
  vars_all <- sort(unique(df$var))
  if (length(vars_all) == 0) stop("Keine passenden Einträge gefunden.")
  vars_pick <- utils::select.list(vars_all, title = "Variablen wählen (Mehrfachauswahl möglich)", multiple = TRUE)
  if (length(vars_pick) == 0) stop("Keine Variable gewählt.")
  
  # 2) Tage wählen (pro Variable)
  picks <- lapply(vars_pick, function(v) {
    dates_v <- df$date[df$var == v] |> unique() |> sort()
    labels  <- format(dates_v, "%Y-%m-%d")
    sel     <- utils::select.list(labels, title = paste("Tage wählen für", v, "(Mehrfachauswahl möglich)"), multiple = TRUE)
    tibble::tibble(var = v, date = as.Date(sel))
  })
  dplyr::bind_rows(picks) |> dplyr::arrange(var, date)
}

# Non-interactive Helfer: wähle gezielt
filter_missing <- function(miss_df, vars_to_run = NULL, dates_to_run = NULL,
                           startDate = NULL, endDate = NULL) {
  df <- miss_df
  if (!is.null(vars_to_run))  df <- df |> dplyr::filter(var %in% vars_to_run)
  if (!is.null(dates_to_run)) df <- df |> dplyr::filter(date %in% as.Date(dates_to_run))
  if (!is.null(startDate))    df <- df |> dplyr::filter(date >= as.Date(startDate))
  if (!is.null(endDate))      df <- df |> dplyr::filter(date <= as.Date(endDate))
  df |> dplyr::distinct(var, date) |> dplyr::arrange(var, date)
}


# --- Simple GeoTIFF write (no COG) + Sanity-Check ---------------------------
safe_write_tif <- function(pred_stars, target_file) {
  tmpf <- paste0(target_file, ".tmp")
  
  # 1) Schreiben als Float32 (kein COG)
  r <- terra::rast(pred_stars)  # stars -> SpatRaster (alle Bänder)
  terra::writeRaster(
    r, tmpf,
    overwrite = TRUE,
    filetype  = "GTiff",
    datatype  = "FLT4S",  # Float32
    gdal      = c("COMPRESS=DEFLATE","TILED=YES","BLOCKXSIZE=256","BLOCKYSIZE=256"),
    NAflag    = -9999
  )
  
  # 2) Rücklesen & prüfen (Band 1)
  rt <- try(terra::rast(tmpf), silent = TRUE)
  if (inherits(rt, "try-error")) stop("[WRITE] Rücklesen fehlgeschlagen: ", tmpf)
  
  v <- terra::values(rt[[1]], mat = FALSE)
  nz <- sum(is.finite(v) & v != 0, na.rm = TRUE)
  nfin <- sum(is.finite(v), na.rm = TRUE)
  
  # 3) Notbremse: alles 0 oder alles NA → Fehler auslösen
  if (nfin == 0 || nz == 0) {
    unlink(tmpf)
    stop(sprintf("[WRITE] Sanity-Check fail (nfin=%d, nz=%d). Abbruch.", nfin, nz))
  }
  
  file.rename(tmpf, target_file)
  invisible(TRUE)
}

# =====================================================================
# DWD CDC — Daily KL (historical + recent) — SCRIPT KAPSELUNG 1:1
#
# Ziel:
# - KEINE neue Semantik, KEINE Optimierung, KEIN Refactoring der Logik.
# - Das direktdownload.R wird nur in zwei Einheiten geschnitten:
#
#   (A) kl_daily_core_prepare():  entspricht dem Script bis inkl. "Output dt_legacy"
#       -> liefert dt (data.table) im exakt gleichen Zustand wie im Script
#          nach dem Merge/Filter/Meta-Merge und nach dem Entfernen von source/file.
#
#   (B) kl_daily_extract_sf():    entspricht dem Script-Block "Legacy target: sf..."
#       -> nimmt dt und baut sf_legacy (geom-Spalte "geom") für Kriging etc.
#
# Zusatzanforderung:
# - Der Extract-Wrapper kann optional die komplette Kette triggern (Download+Unzip+Core).
#   -> kl_daily_extract_sf(..., force_download_chain=TRUE, ...) ruft core_prepare().
#
# Referenz: /mnt/data/direktdownload.R (user-provided)
# =====================================================================


# ---------------------------------------------------------------------
# (A) CORE: 1:1-Kapselung von direktdownload.R
#     Download/Unzip/Read/Diagnostics/Duplicate-Resolve/Station-Filter
#     + Meta-Join (wie Script) + dt_legacy Output
# ---------------------------------------------------------------------
kl_daily_core_prepare <- function(
    # Zeitfenster
  start_date,
  end_date,
  
  # HARTE Regeln wie im direktdownload-Script
  max_missing_frac = 0.25,
  delta_t = 7L,
  
  # DWD URLs (wie direktdownload)
  base_hist   = "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/",
  base_recent = "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/recent/",
  
  # Spaltenschema (wie direktdownload)
  keep_cols = c(
    "STATIONS_ID","MESS_DATUM","QN_4",
    "RSK","RSKF","SDK","SHK_TAG","NM","VPM","PM","TMK","UPM","TXK","TNK","TGK","eor"
  ),
  
  # Pfade (wie direktdownload; nur parametrisierbar, kein anderer Kram)
  root_raw = file.path("data", "data_lev0", "kl_daily", "raw"),
  out_dir  = file.path("data", "data_lev0", "kl_daily", "extracted"),
  
  # Flags (wie direktdownload)
  do_download = FALSE,
  force_unzip = FALSE,
  force_merge = TRUE,
  auto_resolve_duplicates = TRUE,
  
  # Station-Meta-Datei (lokal; exakt diese Datei liest der Script-Teil 9.5)
  meta_file = file.path("data", "data_lev0", "kl_daily", "raw", "data/raw/KL_Tageswerte_Beschreibung_Stationen.txt"),
  
  # Output merged csv (wie direktdownload)
  write_merged_csv = TRUE
) {
  
  suppressPackageStartupMessages({
    library(data.table)
    library(rvest)
    library(xml2)
  })
  
  # --- Date coercion (minimal, aber zwingend) ---
  if (is.character(start_date)) start_date <- as.Date(start_date)
  if (is.character(end_date))   end_date   <- as.Date(end_date)
  stopifnot(inherits(start_date, "Date"), inherits(end_date, "Date"))
  if (end_date < start_date) stop("end_date < start_date")
  
  # --- Folder structure (1:1 direktdownload) ---
  raw_hist <- file.path(root_raw, "historical")
  raw_rec  <- file.path(root_raw, "recent")
  ex_hist  <- file.path(root_raw, "extracted_historical")
  ex_rec   <- file.path(root_raw, "extracted_recent")
  
  dir.create(raw_hist, showWarnings = FALSE, recursive = TRUE)
  dir.create(raw_rec,  showWarnings = FALSE, recursive = TRUE)
  dir.create(ex_hist,  showWarnings = FALSE, recursive = TRUE)
  dir.create(ex_rec,   showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
  
  n_days_total <- as.integer(end_date - start_date) + 1L
  full_dates   <- seq.Date(start_date, end_date, by = "day") # (wie Script)
  
  pat_prod <- "^produkt_klima_tag_\\d{8}_\\d{8}_\\d{5}\\.txt$"
  pat_zip  <- "\\.zip$"
  
  out_merged_csv_gz <- file.path(out_dir, sprintf(
    "kl_daily_merged_%s_%s_stationfilter_miss%02d_gap%d.csv",
    format(start_date, "%Y%m%d"),
    format(end_date, "%Y%m%d"),
    as.integer(max_missing_frac * 100),
    as.integer(delta_t)
  ))
  
  # diagnostics (wie Script)
  prod_files_csv <- file.path(out_dir, "diag_product_files_found.csv")
  qn4_global_csv <- file.path(out_dir, "diag_qn4_global_counts.csv")
  qn4_station_csv<- file.path(out_dir, "diag_qn4_by_station.csv")
  dup_keys_csv   <- file.path(out_dir, "diag_duplicates_station_date.csv")
  dup_rows_csv   <- file.path(out_dir, "diag_duplicates_rows_detail.csv")
  st_diag_csv    <- file.path(out_dir, "diag_station_window_missing_gap.csv")
  
  # PRE-FLIGHT (wie Script)
  if (!isTRUE(force_merge) && file.exists(out_merged_csv_gz)) {
    message("FINAL OUTPUT EXISTS → SKIP MERGE (set force_merge=TRUE): ", out_merged_csv_gz)
    stop("force_merge=FALSE and final output exists (script would quit here).")
  }
  
  # ============================================================
  # HELPERS (1:1 aus direktdownload)
  # ============================================================
  list_remote_files <- function(url, pattern = NULL) {
    doc <- rvest::read_html(url)
    files <- doc |> rvest::html_elements("a") |> rvest::html_attr("href")
    files <- files[!is.na(files) & files != "../"]
    if (!is.null(pattern)) files <- files[grepl(pattern, files)]
    files <- unique(files)
    data.table(url = paste0(url, files), file = files)
  }
  
  download_if_missing <- function(url, destfile) {
    if (file.exists(destfile)) return(invisible(FALSE))
    message("download zip: ", destfile)
    tryCatch({
      utils::download.file(url, destfile = destfile, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) {
      message("DOWNLOAD FAILED: ", url)
      message(conditionMessage(e))
      FALSE
    })
  }
  
  unzip_zip <- function(zipfile, exdir) {
    suppressWarnings(utils::unzip(zipfile, exdir = exdir))
    invisible(TRUE)
  }
  
  parse_prod_span <- function(files) {
    bn <- basename(files)
    m <- regexec("^produkt_klima_tag_(\\d{8})_(\\d{8})_(\\d{5})\\.txt$", bn)
    mm <- regmatches(bn, m)
    ok <- lengths(mm) == 4
    if (!all(ok)) stop("Unexpected product filenames encountered.")
    data.table(
      STATIONS_ID = vapply(mm, `[[`, character(1), 4),
      file_start  = as.Date(vapply(mm, `[[`, character(1), 2), "%Y%m%d"),
      file_end    = as.Date(vapply(mm, `[[`, character(1), 3), "%Y%m%d"),
      file        = bn
    )
  }
  
  read_one_product_file <- function(f, source_tag) {
    bn <- basename(f)
    m <- regexec("^produkt_klima_tag_(\\d{8})_(\\d{8})_(\\d{5})\\.txt$", bn)
    mm <- regmatches(bn, m)[[1]]
    if (length(mm) != 4) stop("Unexpected product filename: ", bn)
    
    file_sid <- mm[4]
    
    dt <- data.table::fread(
      f,
      sep = ";",
      na.strings = "-999",
      strip.white = TRUE,
      showProgress = FALSE
    )
    
    req <- c("STATIONS_ID", "MESS_DATUM", "QN_4")
    miss <- setdiff(req, names(dt))
    if (length(miss) > 0) {
      message("SKIP (missing cols): ", bn, " -> ", paste(miss, collapse = ", "))
      return(NULL)
    }
    
    dt[, STATIONS_ID := sprintf("%05d", as.integer(STATIONS_ID))]
    dt[, MESS_DATUM := as.Date(as.character(MESS_DATUM), format = "%Y%m%d")]
    dt[, QN_4 := as.integer(QN_4)]
    
    # HARD station-id validation: filename id must match content (wie Script)
    if (any(dt$STATIONS_ID != file_sid, na.rm = TRUE)) {
      bad <- dt[STATIONS_ID != file_sid, .N]
      message(
        "STATION-ID MISMATCH -> dropping mismatching rows (keep filename SID only): ",
        "filename SID=", file_sid, " file=", bn, " mismatching rows=", bad
      )
      dt <- dt[STATIONS_ID == file_sid]
      if (nrow(dt) == 0) return(NULL)
    }
    
    # date window filter only (wie Script)
    dt <- dt[!is.na(MESS_DATUM) & MESS_DATUM >= start_date & MESS_DATUM <= end_date]
    if (nrow(dt) == 0) return(NULL)
    
    present <- intersect(keep_cols, names(dt))
    dt <- dt[, ..present]
    for (cc in setdiff(keep_cols, names(dt))) dt[, (cc) := NA]
    setcolorder(dt, keep_cols)
    
    # provenance for duplicate diagnostics (wie Script)
    dt[, source := source_tag]
    dt[, file := bn]
    
    dt[]
  }
  
  # ============================================================
  # 1) Remote list zips (wie Script)
  # ============================================================
  z_hist <- list_remote_files(base_hist, pat_zip)
  z_rec  <- list_remote_files(base_recent, pat_zip)
  
  message("zip count historical: ", nrow(z_hist))
  message("zip count recent:     ", nrow(z_rec))
  
  z_hist[, dest := file.path(raw_hist, file)]
  z_rec[,  dest := file.path(raw_rec,  file)]
  
  # ============================================================
  # 2) Download (wie Script)
  # ============================================================
  if (isTRUE(do_download)) {
    miss_hist <- z_hist[!file.exists(dest)]
    miss_rec  <- z_rec [!file.exists(dest)]
    
    message("ZIPs missing (historical): ", nrow(miss_hist))
    message("ZIPs missing (recent):     ", nrow(miss_rec))
    
    for (i in seq_len(nrow(miss_hist))) download_if_missing(miss_hist$url[i], miss_hist$dest[i])
    for (i in seq_len(nrow(miss_rec)))  download_if_missing(miss_rec$url[i],  miss_rec$dest[i])
  } else {
    message("DOWNLOAD DISABLED (do_download=FALSE) → using local ZIPs only")
  }
  
  missing_hist <- z_hist[!file.exists(dest)]
  missing_rec  <- z_rec [!file.exists(dest)]
  if (nrow(missing_hist) + nrow(missing_rec) > 0) {
    stop(
      "Missing required ZIPs locally. Set do_download=TRUE or place zips in:\n",
      "  ", raw_hist, "  (historical)\n",
      "  ", raw_rec,  "  (recent)\n",
      "Missing counts: hist=", nrow(missing_hist), " recent=", nrow(missing_rec)
    )
  }
  
  # ============================================================
  # 3) Unzip (wie Script)
  # ============================================================
  have_prod_hist <- length(list.files(ex_hist, pattern = pat_prod, recursive = TRUE)) > 0
  have_prod_rec  <- length(list.files(ex_rec,  pattern = pat_prod, recursive = TRUE)) > 0
  
  if (isTRUE(force_unzip) || !have_prod_hist) {
    message("UNZIP historical -> ", ex_hist)
    for (zf in z_hist$dest) unzip_zip(zf, exdir = ex_hist)
  } else {
    message("UNZIP historical skipped")
  }
  
  if (isTRUE(force_unzip) || !have_prod_rec) {
    message("UNZIP recent -> ", ex_rec)
    for (zf in z_rec$dest) unzip_zip(zf, exdir = ex_rec)
  } else {
    message("UNZIP recent skipped")
  }
  
  # ============================================================
  # 4) Find product files (+ span table from filenames) (wie Script)
  # ============================================================
  prod_hist <- list.files(ex_hist, pattern = pat_prod, full.names = TRUE, recursive = TRUE)
  prod_rec  <- list.files(ex_rec,  pattern = pat_prod, full.names = TRUE, recursive = TRUE)
  
  message("product txt count hist:   ", length(prod_hist))
  message("product txt count recent: ", length(prod_rec))
  
  if (length(prod_hist) + length(prod_rec) == 0) {
    stop("No produkt_klima_tag_*.txt found under: ", ex_hist, " or ", ex_rec)
  }
  
  data.table::fwrite(
    rbind(
      data.table(path = prod_hist, source = "historical"),
      data.table(path = prod_rec,  source = "recent")
    ),
    prod_files_csv, sep = ";"
  )
  
  span_all <- data.table::rbindlist(list(parse_prod_span(prod_hist), parse_prod_span(prod_rec)), use.names = TRUE)
  st_span_files <- span_all[, .(
    min_file_start = min(file_start, na.rm = TRUE),
    max_file_end   = max(file_end,   na.rm = TRUE),
    n_files        = .N
  ), by = STATIONS_ID]
  
  # ============================================================
  # 5) Read all products (wie Script)
  # ============================================================
  res_hist <- vector("list", length(prod_hist))
  for (i in seq_along(prod_hist)) {
    if (i %% 250 == 0) message("read(hist) ", i, " / ", length(prod_hist))
    res_hist[[i]] <- read_one_product_file(prod_hist[i], "historical")
  }
  
  res_rec <- vector("list", length(prod_rec))
  for (i in seq_along(prod_rec)) {
    if (i %% 250 == 0) message("read(recent) ", i, " / ", length(prod_rec))
    res_rec[[i]] <- read_one_product_file(prod_rec[i], "recent")
  }
  
  dt <- data.table::rbindlist(c(res_hist, res_rec), use.names = TRUE, fill = TRUE)
  rm(res_hist, res_rec); gc()
  
  if (nrow(dt) == 0) stop("No rows after reading+date filter.")
  data.table::setDT(dt)  # **wichtig**: garantiert data.table
  
  # ============================================================
  # 6) QN_4 diagnostics (wie Script)
  # ============================================================
  data.table::fwrite(dt[, .N, by = QN_4][order(QN_4)], qn4_global_csv, sep=";")
  data.table::fwrite(dt[, .N, by = .(STATIONS_ID, QN_4)][order(STATIONS_ID, QN_4)], qn4_station_csv, sep=";")
  
  # ============================================================
  # 7) Duplicate check (wie Script)
  # ============================================================
  dup_keys <- dt[, .N, by = .(STATIONS_ID, MESS_DATUM)][N > 1L]
  if (nrow(dup_keys) > 0) {
    data.table::fwrite(dup_keys, dup_keys_csv, sep = ";")
    
    setkey(dt, STATIONS_ID, MESS_DATUM)
    dup_detail <- dt[dup_keys, on = .(STATIONS_ID, MESS_DATUM)]
    data.table::fwrite(dup_detail, dup_rows_csv, sep = ";")
    
    if (!isTRUE(auto_resolve_duplicates)) {
      stop("Stop: duplicates exist. See diagnostics in extracted/. Set auto_resolve_duplicates=TRUE to resolve.")
    }
    
    val_cols <- setdiff(keep_cols, c("STATIONS_ID","MESS_DATUM","QN_4","eor"))
    dt[, nnz := rowSums(!is.na(.SD)), .SDcols = val_cols]
    setorder(dt, STATIONS_ID, MESS_DATUM, -nnz, -QN_4)
    dt <- dt[, .SD[1], by = .(STATIONS_ID, MESS_DATUM)]
    dt[, nnz := NULL]
  }
  
  # ============================================================
  # 8) Station filter (wie Script)
  # ============================================================
  dt[, MESS_DATUM := as.Date(MESS_DATUM)]
  dt_dates <- unique(dt[, .(STATIONS_ID, MESS_DATUM)])
  setorder(dt_dates, STATIONS_ID, MESS_DATUM)
  
  st_cov <- dt_dates[, .(n_dates = .N, min_date = min(MESS_DATUM), max_date = max(MESS_DATUM)), by = STATIONS_ID]
  st_cov[, expected := n_days_total]
  st_cov[, coverage_frac := n_dates / expected]
  st_cov[, missing_frac  := 1 - coverage_frac]
  
  st_gap <- dt_dates[, {
    d <- MESS_DATUM
    gap_before <- as.integer(min(d) - start_date)
    gap_after  <- as.integer(end_date - max(d))
    dif <- diff(d)
    internal <- if (length(dif)) pmax(0L, as.integer(dif) - 1L) else integer()
    list(max_gap = max(c(gap_before, internal, gap_after), na.rm = TRUE))
  }, by = STATIONS_ID]
  
  st_diag <- merge(st_cov, st_span_files, by = "STATIONS_ID", all.x = TRUE)
  st_diag <- merge(st_diag, st_gap,       by = "STATIONS_ID", all.x = TRUE)
  
  data.table::fwrite(st_diag[order(-missing_frac, -max_gap)], st_diag_csv, sep=";")
  
  st_keep <- st_diag[
    min_file_start <= start_date &
      max_file_end   >= end_date &
      missing_frac   <= max_missing_frac &
      max_gap        <= delta_t,
    STATIONS_ID
  ]
  
  dt <- dt[STATIONS_ID %in% st_keep]
  if (nrow(dt) == 0) stop("After station filter no stations remain. See diag_station_window_missing_gap.csv")
  
  # ============================================================
  # 9) HARD POST-CHECKS (wie Script)
  # ============================================================
  stopifnot(all(nchar(dt$STATIONS_ID) == 5))
  stopifnot(all(grepl("^[0-9]{5}$", dt$STATIONS_ID)))
  stopifnot(dt[, anyDuplicated(paste(STATIONS_ID, MESS_DATUM))] == 0L)
  
  dt_dates2 <- unique(dt[, .(STATIONS_ID, MESS_DATUM)])
  st_gap2 <- dt_dates2[, {
    d <- sort(MESS_DATUM)
    gap_before <- as.integer(min(d) - start_date)
    gap_after  <- as.integer(end_date - max(d))
    dif <- diff(d)
    internal <- if (length(dif)) pmax(0L, as.integer(dif) - 1L) else integer()
    list(max_gap = max(c(gap_before, internal, gap_after), na.rm = TRUE),
         n_dates  = .N)
  }, by = STATIONS_ID]
  
  stopifnot(all(st_gap2$max_gap <= delta_t))
  stopifnot(all((1 - st_gap2$n_dates / n_days_total) <= max_missing_frac))
  
  # ============================================================
  # 9.5) Merge station metadata (wie Script-Block in direktdownload)
  # ============================================================
  suppressPackageStartupMessages({
    library(sf)
  })
  
  meta_file <- path.expand(meta_file)
  stopifnot(file.exists(meta_file))
  
  ln <- readLines(meta_file, warn = FALSE)
  x  <- ln[-c(1,2)]
  x  <- x[nzchar(trimws(x))]
  
  st_meta <- data.table::rbindlist(lapply(x, function(s){
    m <- regexec(
      "^\\s*([0-9]{5})\\s+([0-9]{8})\\s+([0-9]{8})\\s+(-?[0-9]+)\\s+([0-9]+\\.[0-9]+)\\s+([0-9]+\\.[0-9]+)\\s*(.*)$",
      s
    )
    z <- regmatches(s, m)[[1]]
    if (length(z) != 8) return(NULL)
    
    rest <- trimws(z[8])
    p <- trimws(unlist(strsplit(rest, "\\s{2,}", perl = TRUE)))
    p <- p[nzchar(p)]
    
    list(
      STATIONS_ID   = sprintf("%05d", as.integer(z[2])),
      von_datum     = as.Date(z[3], "%Y%m%d"),
      bis_datum     = as.Date(z[4], "%Y%m%d"),
      Stationshoehe = as.numeric(z[5]),
      geoBreite     = as.numeric(z[6]),
      geoLaenge     = as.numeric(z[7]),
      Stationsname  = if (length(p) >= 1) p[1] else NA_character_,
      Bundesland    = if (length(p) >= 2) p[2] else NA_character_,
      Abgabe        = if (length(p) >= 3) p[3] else NA_character_
    )
  }), fill = TRUE)
  
  fix_enc <- function(v) iconv(v, to = "UTF-8", sub = "")
  st_meta[, `:=`(
    Stationsname = fix_enc(Stationsname),
    Bundesland   = fix_enc(Bundesland),
    Abgabe       = fix_enc(Abgabe)
  )]
  
  st_meta_u <- st_meta[, .(
    Stationsname  = Stationsname[1],
    Bundesland    = Bundesland[1],
    Abgabe        = Abgabe[1],
    Stationshoehe = Stationshoehe[1],
    geoBreite     = geoBreite[1],
    geoLaenge     = geoLaenge[1]
  ), by = STATIONS_ID]
  
  st_span_u <- st_meta[, .(
    von_datum = as.integer(format(min(von_datum, na.rm = TRUE), "%Y%m%d")),
    bis_datum = as.integer(format(max(bis_datum, na.rm = TRUE), "%Y%m%d"))
  ), by = STATIONS_ID]
  
  setkey(st_span_u, STATIONS_ID)
  setkey(dt,        STATIONS_ID)
  dt <- st_span_u[dt]
  
  if (anyNA(dt$von_datum) || anyNA(dt$bis_datum)) {
    bad <- unique(dt[is.na(von_datum) | is.na(bis_datum), STATIONS_ID])
    stop("Missing von_datum/bis_datum for STATIONS_ID: ", paste(bad, collapse = ", "))
  }
  
  setkey(st_meta_u, STATIONS_ID)
  dt[, STATIONS_ID := sprintf("%05d", as.integer(STATIONS_ID))]
  setkey(dt, STATIONS_ID)
  dt <- st_meta_u[dt]
  
  # ============================================================
  # 10) Output dt_legacy (wie Script)
  # ============================================================
  dt[, c("source","file") := NULL]
  setorder(dt, STATIONS_ID, MESS_DATUM)
  setcolorder(dt, keep_cols)
  
  if (isTRUE(write_merged_csv)) {
    data.table::fwrite(dt, out_merged_csv_gz, sep = ";", na = "")
    message("merged saved: ", out_merged_csv_gz)
  }
  
  dt
}


# ---------------------------------------------------------------------
# (B) EXTRACT: dt -> sf_legacy (geom-Spalte "geom") — 1:1 Script
# ---------------------------------------------------------------------
kl_daily_extract_sf <- function(
    dt = NULL,
    
    # chain switch: TRUE => core wird mit do_download=TRUE ausgeführt
    force_download_chain = FALSE,
    
    # welche Param-Spalte als Wert in sf
    param = "PM",
    
    # core args (nur relevant wenn dt NULL oder force_download_chain=TRUE)
    start_date = as.Date("2003-01-01"),
    end_date   = as.Date("2024-12-31"),
    max_missing_frac = 0.25,
    delta_t = 7L,
    base_hist   = "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/",
    base_recent = "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/recent/",
    keep_cols = c(
      "STATIONS_ID","MESS_DATUM","QN_4",
      "RSK","RSKF","SDK","SHK_TAG","NM","VPM","PM","TMK","UPM","TXK","TNK","TGK","eor"
    ),
    root_raw = file.path("data", "data_lev0", "kl_daily", "raw"),
    out_dir  = file.path("data", "data_lev0", "kl_daily", "extracted"),
    meta_file = file.path("data", "data_lev0", "kl_daily", "raw", "KL_Tageswerte_Beschreibung_Stationen.txt"),
    do_download = TRUE,
    force_unzip = TRUE,
    force_merge = TRUE,
    auto_resolve_duplicates = TRUE,
    write_merged_csv = TRUE
) {
  
  suppressPackageStartupMessages({
    library(data.table)
    library(sf)
  })
  
  if (isTRUE(force_download_chain) || is.null(dt)) {
    dt <- kl_daily_core_prepare(
      start_date = start_date,
      end_date   = end_date,
      max_missing_frac = max_missing_frac,
      delta_t = delta_t,
      base_hist = base_hist,
      base_recent = base_recent,
      keep_cols = keep_cols,
      root_raw = root_raw,
      out_dir  = out_dir,
      do_download = if (isTRUE(force_download_chain)) TRUE else isTRUE(do_download),
      force_unzip = force_unzip,
      force_merge = force_merge,
      auto_resolve_duplicates = auto_resolve_duplicates,
      meta_file = meta_file,
      write_merged_csv = write_merged_csv
    )
  }
  
  data.table::setDT(dt) # garantiert data.table
  
  need <- c("STATIONS_ID","MESS_DATUM", param,
            "Stationshoehe","Stationsname","Bundesland","Abgabe",
            "geoBreite","geoLaenge")
  stopifnot(all(need %chin% names(dt)))
  
  dt_legacy <- data.table::copy(dt)
  
  dt_legacy[, STATIONS_ID   := as.numeric(STATIONS_ID)]
  dt_legacy[, Stationshoehe := as.integer(Stationshoehe)]
  dt_legacy[, MESS_DATUM    := as.Date(MESS_DATUM)]
  
  if (!("von_datum" %chin% names(dt_legacy))) dt_legacy[, von_datum := NA_integer_]
  if (!("bis_datum" %chin% names(dt_legacy))) dt_legacy[, bis_datum := NA_integer_]
  dt_legacy[, von_datum := as.integer(von_datum)]
  dt_legacy[, bis_datum := as.integer(bis_datum)]
  
  keep <- c("STATIONS_ID","von_datum","bis_datum","Stationshoehe",
            "Stationsname","Bundesland","Abgabe",
            "MESS_DATUM", param, "geoLaenge","geoBreite")
  dt_legacy <- dt_legacy[, ..keep]
  
  sf_legacy <- sf::st_as_sf(
    dt_legacy,
    coords = c("geoLaenge","geoBreite"),
    crs = 4326,
    remove = TRUE
  )
  sf_legacy$STATIONS_ID <- as.numeric(sf_legacy$STATIONS_ID)
  sf_legacy$geom <- sf::st_geometry(sf_legacy)
  sf_legacy <- sf::st_set_geometry(sf_legacy, "geom")
  if ("geometry" %in% names(sf_legacy)) sf_legacy$geometry <- NULL
  
  sf_legacy
}
