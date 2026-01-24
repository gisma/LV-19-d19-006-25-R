#!/usr/bin/env Rscript

#######################################################################
# Script:   02-3_processing_burgwald-dwd-windfield
# Project:  Burgwald
#
# CONTENT + TECHNICAL OVERVIEW (ENGLISH)
# --------------------------------------
# This script performs a wind climatology exploration for DWD hourly wind
# observations and writes all derived artefacts (tables + figures) to
# canonical output locations defined in `src/_core/outputs.tsv` via `paths`.
#
# IMPORTANT PROJECT RULES IMPLEMENTED HERE
# ----------------------------------------
# (1) INPUT RULE (your decision):
#     Raw input data are NOT managed by `outputs.tsv`.
#     The script therefore reads DWD wind ZIPs directly from:
#       data/raw/providers/dwd/**/stundenwerte_FF*.zip
#     (folder convention, not a contract path).
#
# (2) OUTPUT RULE:
#     ALL output artefacts (S2 tables and S2 figures) are written only via
#     `paths[...]` created in `src/_core/01-setup-burgwald.R`.
#
# (3) AOI RULE:
#     The Burgwald AOI is provided by the setup as an in-memory object:
#       aoi_burgwald_wgs  (typically EPSG:4326)
#
# WHAT THE SCRIPT DOES
# --------------------
# A) Load DWD hourly wind measurements from ZIPs (stundenwerte_FF*.zip).
#    - Extract and read the "produkt_ff_stunde_*" data file inside each ZIP.
#    - Extract and read the "Metadaten_Geographie_*" file inside each ZIP.
#    - Harmonize the geo-metadata column names to:
#         STATIONS_ID, geoBreite, geoLaenge, Stationsname
#      because DWD geo metadata uses different labels (e.g. "Stations_id",
#      "Geogr.Breite", "Geogr.Laenge").
#    - Join measurements with station coordinates and build an sf object.
#
# B) Compute wind roses + summaries for two station sets:
#    (1) ALL stations available in the raw ZIP collection.
#    (2) ONLY those stations located within a 20 km buffer around the AOI.
#
# C) Persist all artefacts (S2):
#    - openair-ready input tables (RDS)
#    - wind rose figures (PNG)
#    - mean wind summaries + direction frequency tables (RDS)
#    - buffer-mean time series + seasonal/annual mean table (RDS)
#
# NOTES ON UNITS AND CRS
# ----------------------
# - The AOI object is assumed WGS84 lon/lat. Distances in WGS84 degrees
#   are not meaningful, so the 20 km buffer is computed after transforming
#   to a metric CRS (EPSG:25832 / UTM32N) and then used for spatial filtering.
#
# DEPENDENCIES
# ------------
# - dplyr, sf, lubridate: data handling + spatial operations
# - openair: windRose plotting
# - base R utils: unzip, read.delim for DWD semicolon files
#
#######################################################################

## -------------------------------------------------------------------
## 0) Libraries + setup
## -------------------------------------------------------------------

library(here)
library(dplyr)
library(sf)
library(openair)
library(lubridate)

# Setup must define:
# - `paths` (from outputs.tsv; used ONLY for outputs)
# - `aoi_burgwald_wgs` (AOI polygon in WGS, used for buffering)
source(here::here("src", "_core", "01-setup-burgwald.R"))

## -------------------------------------------------------------------
## 1) RAW input: read DWD hourly wind ZIPs
## -------------------------------------------------------------------

# RAW provider root (folder convention, NOT via outputs.tsv)
# This is the only hardwired input location in the script by design.
wind_zip_root <- here::here("data", "raw", "providers", "dwd")

# We expect DWD hourly wind packages named like:
#   stundenwerte_FF_XXXXXX_YYYYMMDD_YYYYMMDD_hist.zip
# or similar. We search recursively under the provider root.
zip_files <- list.files(
  wind_zip_root,
  pattern = "^stundenwerte_FF.*\\.zip$",
  full.names = TRUE,
  recursive = TRUE
)

# Unzip workspace in project tmp; this is a working directory, not an artefact.
unzip_root <- here::here("tmp", "dwd_hourly_wind_unzip")
dir.create(unzip_root, recursive = TRUE, showWarnings = FALSE)

# DWD CDC text files in these ZIPs are typically semicolon-separated (";").
# We read them with base R, keeping strings as strings.
read_dwd_semicolon <- function(f) {
  utils::read.delim(
    file = f,
    sep = ";",
    stringsAsFactors = FALSE,
    strip.white = TRUE
  )
}

# We collect:
# - wind_list: measurement tables from produkt_ff_stunde_* files
# - station_geo_list: station geo metadata from Metadaten_Geographie_* files
wind_list <- list()
station_geo_list <- list()

# Iterate all ZIPs and extract exactly one product file and one geo-meta file
# per ZIP (if present). The DWD packages are consistent in practice.
for (zf in zip_files) {
  # List ZIP contents without extracting everything (fast).
  zlist <- utils::unzip(zf, list = TRUE)
  
  # Identify the relevant file(s) within the ZIP by basename pattern.
  prod_name <- zlist$Name[grepl("^produkt_ff_stunde_", basename(zlist$Name))]
  geo_name  <- zlist$Name[grepl("^Metadaten_Geographie_", basename(zlist$Name))]
  
  # Create a deterministic unzip folder per ZIP to avoid collisions.
  outdir <- file.path(unzip_root, tools::file_path_sans_ext(basename(zf)))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract + read product file (hourly measurements)
  if (length(prod_name) >= 1) {
    utils::unzip(zf, files = prod_name[1], exdir = outdir, junkpaths = TRUE)
    prod_path <- file.path(outdir, basename(prod_name[1]))
    wind_list[[length(wind_list) + 1]] <- read_dwd_semicolon(prod_path)
  }
  
  # Extract + read geo metadata file (station coordinates etc.)
  if (length(geo_name) >= 1) {
    utils::unzip(zf, files = geo_name[1], exdir = outdir, junkpaths = TRUE)
    geo_path <- file.path(outdir, basename(geo_name[1]))
    station_geo_list[[length(station_geo_list) + 1]] <- read_dwd_semicolon(geo_path)
  }
}

# Bind all measurements into one table. For your dataset this is large (~millions rows).
wind_df <- bind_rows(wind_list)

# -------------------------------------------------------------------
# IMPORTANT: Harmonize DWD geo metadata column names
# -------------------------------------------------------------------
# Motivation:
# - The wind product tables contain station ID as "STATIONS_ID".
# - The geo metadata tables often contain station ID as "Stations_id"
#   and coordinate columns as "Geogr.Breite" / "Geogr.Laenge".
#
# Without the rename() below, the join between wind_df and station_geo will fail.
# That would break:
# - creation of stations_sf
# - creation of wind_raw sf table
# - AOI-buffer filtering (inside_idx)
station_geo <- bind_rows(station_geo_list) %>%
  rename(
    STATIONS_ID = Stations_id,
    geoBreite   = `Geogr.Breite`,
    geoLaenge   = `Geogr.Laenge`
  ) %>%
  distinct(STATIONS_ID, .keep_all = TRUE)

# Parse datetime from DWD field MESS_DATUM (hourly timestamp encoded as YYYYMMDDHH).
# tz="UTC" is used here to avoid implicit local timezone assumptions in summaries/plots.
wind_df <- wind_df %>%
  mutate(datetime = lubridate::ymd_h(as.character(MESS_DATUM), tz = "UTC"))

# Convert station geo table to sf points (WGS84 lon/lat).
# Note: geoBreite/geoLaenge are retained as numeric columns AND used for geometry.
stations_sf <- station_geo %>%
  transmute(
    STATIONS_ID  = STATIONS_ID,
    Stationsname = Stationsname,
    geoBreite    = as.numeric(geoBreite),
    geoLaenge    = as.numeric(geoLaenge)
  ) %>%
  st_as_sf(coords = c("geoLaenge", "geoBreite"), crs = 4326, remove = FALSE)

# Join measurements with station attributes and geometry.
# Result: an sf object (EPSG:4326) where each measurement row carries the station geometry.
# This is intentionally denormalized (many rows per station) because later filtering
# uses st_within on the measurement points (effectively station membership replicated).
wind_raw <- wind_df %>%
  left_join(st_drop_geometry(stations_sf), by = "STATIONS_ID") %>%
  left_join(stations_sf %>% select(STATIONS_ID, geometry), by = "STATIONS_ID") %>%
  st_as_sf(crs = 4326)

# Quick reporting: station IDs and names observed in the ZIP collection.
message("Stations in RAW wind ZIPs (STATIONS_ID):")
print(unique(wind_raw$STATIONS_ID))
message("Stations in RAW wind ZIPs (Stationsname):")
print(unique(wind_raw$Stationsname))

## -------------------------------------------------------------------
## 2) Helper: circular mean for wind direction (degrees)
## -------------------------------------------------------------------
# Wind direction is circular (0° == 360°). A linear mean is wrong.
# This function returns a mean direction in [0, 360).
mean_wdir_deg <- function(wd_deg) {
  theta <- wd_deg * pi / 180
  sin_mean <- mean(sin(theta), na.rm = TRUE)
  cos_mean <- mean(cos(theta), na.rm = TRUE)
  ang <- atan2(sin_mean, cos_mean) * 180 / pi
  if (ang < 0) ang <- ang + 360
  ang
}

# Direction bins used for frequency tables (30-degree sectors).
# Used in A5/B6 for direction-frequency tables.
dir_breaks <- seq(0, 360, by = 30)

## ===================================================================
## PART A: ALL STATIONS
## ===================================================================

## -------------------------------------------------------------------
## A1) Prepare openair input table (all stations)
## -------------------------------------------------------------------
# openair::windRose expects:
# - date: datetime
# - ws: wind speed
# - wd: wind direction
# - optional: station grouping
#
# Mapping from DWD product columns:
# - F (wind speed) -> ws
# - D (wind direction) -> wd
wind_for_openair_all <- wind_raw %>%
  st_drop_geometry() %>%
  transmute(
    date    = datetime,
    ws      = as.numeric(F),
    wd      = as.numeric(D),
    station = Stationsname
  ) %>%
  filter(!is.na(ws), !is.na(wd))

# Explicit iconv: ensure station labels are valid UTF-8 for plotting and stable RDS.
wind_for_openair_all <- wind_for_openair_all %>%
  mutate(station = iconv(station, from = "", to = "UTF-8", sub = ""))

# Persist the prepared table as an artefact (S2).
# This is the canonical "openair-ready" input (all stations).
saveRDS(wind_for_openair_all, paths[["wind_for_openair_all"]])

## -------------------------------------------------------------------
## A2–A4) Wind rose figures for ALL stations (PNG artefacts)
## -------------------------------------------------------------------

# (A2) single wind rose combining all stations
png(paths[["windrose_all_combined"]], width = 1800, height = 1400, res = 150)
windRose(
  mydata = wind_for_openair_all,
  ws     = "ws",
  wd     = "wd",
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)
dev.off()

# (A3) one panel per station
png(paths[["windrose_all_by_station"]], width = 2200, height = 1800, res = 150)
windRose(
  mydata = wind_for_openair_all,
  ws     = "ws",
  wd     = "wd",
  type   = "station",
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)
dev.off()

# (A4) combined panel set: "All stations" + each station as its own panel
wind_for_grid_all <- bind_rows(
  wind_for_openair_all %>% mutate(panel = "All stations"),
  wind_for_openair_all %>% mutate(panel = station)
)
wind_for_grid_all$panel <- factor(
  wind_for_grid_all$panel,
  levels = c("All stations", sort(unique(wind_for_openair_all$station)))
)

png(paths[["windrose_all_panel"]], width = 2200, height = 1800, res = 150)
windRose(
  mydata = wind_for_grid_all,
  ws     = "ws",
  wd     = "wd",
  type   = "panel",
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)
dev.off()

## -------------------------------------------------------------------
## A5) Tabular summaries for ALL stations (RDS artefacts)
## -------------------------------------------------------------------
# This block derives:
# - overall mean ws and circular mean wd (across all stations)
# - per-station mean ws and circular mean wd
# - direction frequency (overall and per station) in 30° bins
#
# These tables are typically used to:
# - document prevailing wind regimes (direction + speed)
# - check consistency between stations
# - support later decisions for windwardness descriptors / stratification

overall_summary_all <- wind_for_openair_all %>%
  summarise(
    n       = n(),
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd)
  )

station_summary_all <- wind_for_openair_all %>%
  group_by(station) %>%
  summarise(
    n       = n(),
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd),
    .groups = "drop"
  ) %>%
  arrange(station)

dir_table_overall_all <- wind_for_openair_all %>%
  mutate(dir_class = cut(wd, breaks = dir_breaks, include.lowest = TRUE, right = FALSE)) %>%
  count(dir_class, name = "n") %>%
  mutate(rel_freq = n / sum(n))

dir_table_by_station_all <- wind_for_openair_all %>%
  mutate(dir_class = cut(wd, breaks = dir_breaks, include.lowest = TRUE, right = FALSE)) %>%
  count(station, dir_class, name = "n") %>%
  group_by(station) %>%
  mutate(rel_freq = n / sum(n)) %>%
  ungroup() %>%
  arrange(station, dir_class)

# Persist all ALL-stations tables
saveRDS(overall_summary_all, paths[["overall_summary_all"]])
saveRDS(station_summary_all, paths[["station_summary_all"]])
saveRDS(dir_table_overall_all, paths[["dir_table_overall_all"]])
saveRDS(dir_table_by_station_all, paths[["dir_table_by_station_all"]])

## ===================================================================
## PART B: STATIONS INSIDE BURGWALD BUFFER (20 km)
## ===================================================================

## -------------------------------------------------------------------
## B1) Compute 20 km buffer in a metric CRS and spatially filter stations
## -------------------------------------------------------------------
# Rationale:
# - Buffer distances must be computed in meters => switch from EPSG:4326 to EPSG:25832.
# - wind_raw contains station geometry replicated per measurement row.
# - inside_idx flags those measurement rows whose station point lies within the buffer.
burgwald_buf_m <- st_transform(aoi_burgwald_wgs, 25832) %>%
  st_buffer(dist = 20000)

wind_raw_m <- st_transform(wind_raw, 25832)

inside_idx <- st_within(wind_raw_m, burgwald_buf_m, sparse = FALSE)[, 1]

wind_burgwald_sf <- wind_raw[inside_idx, ]

message("Stations inside Burgwald buffer (STATIONS_ID):")
print(unique(wind_burgwald_sf$STATIONS_ID))
message("Stations inside Burgwald buffer (Stationsname):")
print(unique(wind_burgwald_sf$Stationsname))

## -------------------------------------------------------------------
## B2) Prepare openair input table (buffer-only)
## -------------------------------------------------------------------
# Same mapping as A1, but restricted to the buffer subset.
wind_for_openair_buf <- wind_burgwald_sf %>%
  st_drop_geometry() %>%
  transmute(
    date    = datetime,
    ws      = as.numeric(F),
    wd      = as.numeric(D),
    station = Stationsname
  ) %>%
  filter(!is.na(ws), !is.na(wd))

wind_for_openair_buf <- wind_for_openair_buf %>%
  mutate(station = iconv(station, from = "", to = "UTF-8", sub = ""))

saveRDS(wind_for_openair_buf, paths[["wind_for_openair_buf"]])

## -------------------------------------------------------------------
## B3–B5) Wind rose figures for BUFFER-only stations (PNG artefacts)
## -------------------------------------------------------------------

png(paths[["windrose_buf_combined"]], width = 1800, height = 1400, res = 150)
windRose(
  mydata = wind_for_openair_buf,
  ws     = "ws",
  wd     = "wd",
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)
dev.off()

png(paths[["windrose_buf_by_station"]], width = 2200, height = 1800, res = 150)
windRose(
  mydata = wind_for_openair_buf,
  ws     = "ws",
  wd     = "wd",
  type   = "station",
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)
dev.off()

wind_for_grid_buf <- bind_rows(
  wind_for_openair_buf %>% mutate(panel = "All buffer stations"),
  wind_for_openair_buf %>% mutate(panel = station)
)
wind_for_grid_buf$panel <- factor(
  wind_for_grid_buf$panel,
  levels = c("All buffer stations", sort(unique(wind_for_openair_buf$station)))
)

png(paths[["windrose_buf_panel"]], width = 2200, height = 1800, res = 150)
windRose(
  mydata = wind_for_grid_buf,
  ws     = "ws",
  wd     = "wd",
  type   = "panel",
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)
dev.off()

## -------------------------------------------------------------------
## B6) Tabular summaries for BUFFER-only stations (RDS artefacts)
## -------------------------------------------------------------------

overall_summary_buf <- wind_for_openair_buf %>%
  summarise(
    n       = n(),
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd)
  )

station_summary_buf <- wind_for_openair_buf %>%
  group_by(station) %>%
  summarise(
    n       = n(),
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd),
    .groups = "drop"
  ) %>%
  arrange(station)

dir_table_overall_buf <- wind_for_openair_buf %>%
  mutate(dir_class = cut(wd, breaks = dir_breaks, include.lowest = TRUE, right = FALSE)) %>%
  count(dir_class, name = "n") %>%
  mutate(rel_freq = n / sum(n))

dir_table_by_station_buf <- wind_for_openair_buf %>%
  mutate(dir_class = cut(wd, breaks = dir_breaks, include.lowest = TRUE, right = FALSE)) %>%
  count(station, dir_class, name = "n") %>%
  group_by(station) %>%
  mutate(rel_freq = n / sum(n)) %>%
  ungroup() %>%
  arrange(station, dir_class)

saveRDS(overall_summary_buf, paths[["overall_summary_buf"]])
saveRDS(station_summary_buf, paths[["station_summary_buf"]])
saveRDS(dir_table_overall_buf, paths[["dir_table_overall_buf"]])
saveRDS(dir_table_by_station_buf, paths[["dir_table_by_station_buf"]])

## -------------------------------------------------------------------
## B7) Buffer-mean time series and seasonal split
## -------------------------------------------------------------------
# This creates a "buffer-mean" series:
# - For each timestamp, average ws across stations (linear mean)
# - For each timestamp, average wd across stations (circular mean)
#
# season2 rule (fixed):
# - Summer: Apr–Sep
# - Winter: Oct–Mar
wind_buf_season <- wind_for_openair_buf %>%
  mutate(
    month = lubridate::month(date),
    season2 = if_else(month %in% 4:9, "Summer (Apr–Sep)", "Winter (Oct–Mar)")
  )

wind_buf_mean_ts <- wind_buf_season %>%
  group_by(date, season2) %>%
  summarise(
    ws = mean(ws, na.rm = TRUE),
    wd = mean_wdir_deg(wd),
    .groups = "drop"
  ) %>%
  mutate(station = "Buffer mean")

saveRDS(wind_buf_mean_ts, paths[["wind_buf_mean_ts"]])

## -------------------------------------------------------------------
## B8) Wind rose for buffer-mean split by season (PNG artefact)
## -------------------------------------------------------------------

png(paths[["windrose_buf_mean_season"]], width = 1800, height = 1400, res = 150)
windRose(
  mydata = wind_buf_mean_ts,
  ws     = "ws",
  wd     = "wd",
  type   = "season2",
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)
dev.off()

## -------------------------------------------------------------------
## B9) Seasonal + annual mean wind table from buffer-mean series (RDS)
## -------------------------------------------------------------------
# These means are computed on the *buffer-mean* series (not on raw station rows),
# i.e. each timestamp has equal weight after the within-timestamp averaging in B7.
summer_mean_wind <- wind_buf_mean_ts %>%
  filter(season2 == "Summer (Apr–Sep)") %>%
  summarise(
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd),
    n       = n()
  )

winter_mean_wind <- wind_buf_mean_ts %>%
  filter(season2 == "Winter (Oct–Mar)") %>%
  summarise(
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd),
    n       = n()
  )

annual_mean_wind <- wind_buf_mean_ts %>%
  summarise(
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd),
    n       = n()
  )

wind_means_summary <- tibble::tibble(
  period  = c("Summer", "Winter", "Annual"),
  mean_ws = c(summer_mean_wind$mean_ws, winter_mean_wind$mean_ws, annual_mean_wind$mean_ws),
  mean_wd = c(summer_mean_wind$mean_wd, winter_mean_wind$mean_wd, annual_mean_wind$mean_wd)
)

saveRDS(wind_means_summary, paths[["wind_means_summary"]])
