#!/usr/bin/env Rscript

#######################################################################
# Script:   03_dwd_wind_rose_all_and_burgwald_buffer.R
# Author:   [Your Name]
# Project:  [Your Project Name]
#
# Purpose:
# --------
# - Use ALL stations contained in the DWD hourly wind RDS file.
# - Additionally, select only those stations whose locations lie
#   inside a buffered Burgwald AOI.
#
# For BOTH sets (all stations, buffer-only):
#   * Create wind roses:
#       (1) One wind rose for all stations combined.
#       (2) One wind rose per station in a lattice/grid plot.
#       (3) A combined panel with:
#             - "All stations" in one panel, and
#             - each station in its own panel.
#   * Produce tabular summaries:
#       - Overall mean wind speed and circular mean direction.
#       - Mean wind speed and direction per station.
#       - Direction-frequency tables (overall and per station).
#
# Requirements:
# -------------
# - RDS file with DWD hourly wind data:
#     data/processed/dwd-stations/burgwald_hourly_wind_F-D_20231212_20251211.rds
#   containing at least:
#     * datetime (POSIXct)
#     * F (wind speed)
#     * D (wind direction, degrees)
#     * geometry (POINT, sf)
#     * Stationsname (station name)
# - An sf polygon aoi_burgwald_wgs (Burgwald AOI, EPSG:4326) must be
#   available in the current environment (e.g. from 00-setup-burgwald.R).
#######################################################################

## -------------------------------------------------------------------
## 0) Libraries
## -------------------------------------------------------------------

library(here)
library(dplyr)
library(sf)
library(openair)
library(lubridate)

## -------------------------------------------------------------------
## 1) Load hourly DWD wind data (ALL stations in file)
## -------------------------------------------------------------------
# sourcing of the setup and specific used functions
source("src/_core/00-setup-burgwald.R")
wind_raw <- readRDS(
  here::here("data", "processed", "dwd-stations",
             "burgwald_hourly_wind_F-D_20231212_20251211.rds")
)

# Optional: check that the AOI object exists
if (!exists("aoi_burgwald_wgs")) {
  stop("Object 'aoi_burgwald_wgs' is not found. Please load or source it before running this script.")
}

# CRS check: station data and AOI must share the same CRS
if (!identical(st_crs(wind_raw), st_crs(aoi_burgwald_wgs))) {
  stop("CRS mismatch: wind_raw and aoi_burgwald_wgs have different CRS.")
}

# Quick info: which stations are present in the file?
message("Stations in file (STATIONS_ID):")
print(unique(wind_raw$STATIONS_ID))
message("Stations in file (Stationsname):")
print(unique(wind_raw$Stationsname))


## -------------------------------------------------------------------
## 2) Helper: circular mean of wind direction in degrees
## -------------------------------------------------------------------

mean_wdir_deg <- function(wd_deg) {
  theta <- wd_deg * pi / 180
  sin_mean <- mean(sin(theta), na.rm = TRUE)
  cos_mean <- mean(cos(theta), na.rm = TRUE)
  ang <- atan2(sin_mean, cos_mean) * 180 / pi
  if (ang < 0) ang <- ang + 360
  return(ang)
}

# Direction classes for frequency tables (30° sectors)
dir_breaks <- seq(0, 360, by = 30)


## ===================================================================
## PART A: ALL STATIONS IN THE FILE
## ===================================================================

## -------------------------------------------------------------------
## A1) Prepare data for openair::windRose() (all stations)
## -------------------------------------------------------------------

wind_for_openair_all <- wind_raw %>%
  st_drop_geometry() %>%
  transmute(
    date    = datetime,          # POSIXct timestamp
    ws      = as.numeric(F),     # wind speed [m/s]
    wd      = as.numeric(D),     # wind direction [°]
    station = Stationsname       # station name for faceting
  ) %>%
  filter(!is.na(ws), !is.na(wd))

## -------------------------------------------------------------------
## A2) WINDROSE: all stations combined
## -------------------------------------------------------------------

windRose(
  mydata = wind_for_openair_all,
  ws     = "ws",
  wd     = "wd",
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)

## -------------------------------------------------------------------
## A3) WINDROSE: one panel per station (all stations)
## -------------------------------------------------------------------

windRose(
  mydata = wind_for_openair_all,
  ws     = "ws",
  wd     = "wd",
  type   = "station",   # facet by station name
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
  # you can add layout = c(nrow, ncol) if many stations
)

## -------------------------------------------------------------------
## A4) WINDROSE: combined panel ("All stations" + each station)
## -------------------------------------------------------------------

wind_for_grid_all <- bind_rows(
  wind_for_openair_all %>% mutate(panel = "All stations"),
  wind_for_openair_all %>% mutate(panel = station)
)

# Panel order: "All stations" first, then alphabetical by station
wind_for_grid_all$panel <- factor(
  wind_for_grid_all$panel,
  levels = c("All stations", sort(unique(wind_for_openair_all$station)))
)

windRose(
  mydata = wind_for_grid_all,
  ws     = "ws",
  wd     = "wd",
  type   = "panel",
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)

## -------------------------------------------------------------------
## A5) Tabular summaries for ALL stations
## -------------------------------------------------------------------

# A5.1 Overall summary (all stations combined)
overall_summary_all <- wind_for_openair_all %>%
  summarise(
    n       = n(),
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd)
  )

message("Overall wind summary (ALL stations in file):")
print(overall_summary_all)

# A5.2 Per-station summary
station_summary_all <- wind_for_openair_all %>%
  group_by(station) %>%
  summarise(
    n       = n(),
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd),
    .groups = "drop"
  ) %>%
  arrange(station)

message("Per-station wind summary (ALL stations in file):")
print(station_summary_all)

# A5.3 Direction-frequency table (overall)
dir_table_overall_all <- wind_for_openair_all %>%
  mutate(
    dir_class = cut(
      wd,
      breaks = dir_breaks,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  count(dir_class, name = "n") %>%
  mutate(
    rel_freq = n / sum(n)
  )

message("Overall direction-frequency table (ALL stations):")
print(dir_table_overall_all)

# A5.4 Direction-frequency table per station
dir_table_by_station_all <- wind_for_openair_all %>%
  mutate(
    dir_class = cut(
      wd,
      breaks = dir_breaks,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  count(station, dir_class, name = "n") %>%
  group_by(station) %>%
  mutate(
    rel_freq = n / sum(n)
  ) %>%
  ungroup() %>%
  arrange(station, dir_class)

message("Direction-frequency table per station (ALL stations):")
print(dir_table_by_station_all)


## ===================================================================
## PART B: ONLY STATIONS INSIDE BURGWALD BUFFER
## ===================================================================

## -------------------------------------------------------------------
## B1) Spatial filter: keep only records within Burgwald buffer
## -------------------------------------------------------------------

# 20 km buffer around Burgwald AOI (units in metres)
burgwald_buf <- st_buffer(aoi_burgwald_wgs, dist = 20000)

inside_idx <- st_within(wind_raw, burgwald_buf, sparse = FALSE)[, 1]

wind_burgwald_sf <- wind_raw[inside_idx, ]

message("Stations inside Burgwald buffer (STATIONS_ID):")
print(unique(wind_burgwald_sf$STATIONS_ID))
message("Stations inside Burgwald buffer (Stationsname):")
print(unique(wind_burgwald_sf$Stationsname))

## -------------------------------------------------------------------
## B2) Prepare data for openair::windRose() (buffer-only)
## -------------------------------------------------------------------

wind_for_openair_buf <- wind_burgwald_sf %>%
  st_drop_geometry() %>%
  transmute(
    date    = datetime,
    ws      = as.numeric(F),
    wd      = as.numeric(D),
    station = Stationsname
  ) %>%
  filter(!is.na(ws), !is.na(wd))

## -------------------------------------------------------------------
## B3) WINDROSE: all buffer stations combined
## -------------------------------------------------------------------

windRose(
  mydata = wind_for_openair_buf,
  ws     = "ws",
  wd     = "wd",
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)

## -------------------------------------------------------------------
## B4) WINDROSE: one panel per station (buffer-only)
## -------------------------------------------------------------------

windRose(
  mydata = wind_for_openair_buf,
  ws     = "ws",
  wd     = "wd",
  type   = "station",
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)

## -------------------------------------------------------------------
## B5) WINDROSE: combined panel ("All buffer stations" + each station)
## -------------------------------------------------------------------

wind_for_grid_buf <- bind_rows(
  wind_for_openair_buf %>% mutate(panel = "All buffer stations"),
  wind_for_openair_buf %>% mutate(panel = station)
)

wind_for_grid_buf$panel <- factor(
  wind_for_grid_buf$panel,
  levels = c("All buffer stations", sort(unique(wind_for_openair_buf$station)))
)

windRose(
  mydata = wind_for_grid_buf,
  ws     = "ws",
  wd     = "wd",
  type   = "panel",
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)

## -------------------------------------------------------------------
## B6) Tabular summaries for BUFFER-ONLY stations
## -------------------------------------------------------------------

# B6.1 Overall summary (buffer-only)
overall_summary_buf <- wind_for_openair_buf %>%
  summarise(
    n       = n(),
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd)
  )

message("Overall wind summary (stations inside Burgwald buffer):")
print(overall_summary_buf)

# B6.2 Per-station summary (buffer-only)
station_summary_buf <- wind_for_openair_buf %>%
  group_by(station) %>%
  summarise(
    n       = n(),
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd),
    .groups = "drop"
  ) %>%
  arrange(station)

message("Per-station wind summary (stations inside Burgwald buffer):")
print(station_summary_buf)

# B6.3 Direction-frequency table (buffer-only, overall)
dir_table_overall_buf <- wind_for_openair_buf %>%
  mutate(
    dir_class = cut(
      wd,
      breaks = dir_breaks,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  count(dir_class, name = "n") %>%
  mutate(
    rel_freq = n / sum(n)
  )

message("Overall direction-frequency table (stations inside Burgwald buffer):")
print(dir_table_overall_buf)

# B6.4 Direction-frequency table per station (buffer-only)
dir_table_by_station_buf <- wind_for_openair_buf %>%
  mutate(
    dir_class = cut(
      wd,
      breaks = dir_breaks,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  count(station, dir_class, name = "n") %>%
  group_by(station) %>%
  mutate(
    rel_freq = n / sum(n)
  ) %>%
  ungroup() %>%
  arrange(station, dir_class)
## -------------------------------------------------------------------
## B7) Spatial mean over buffer stations + split into summer / winter
## -------------------------------------------------------------------
# Definition:
#   Here: Summer = April–September (months 4–9)
#          Winter = October–March (10–12, 1–3)
#   Adjust the month sets if you prefer another definition.

wind_buf_season <- wind_for_openair_buf %>%
  dplyr::mutate(
    month = lubridate::month(date),
    season2 = dplyr::if_else(
      month %in% 4:9,
      "Summer (Apr–Sep)",
      "Winter (Oct–Mar)"
    )
  )

# B7.1 Hourly time series: mean over all buffer stations
#      - ws: arithmetic mean
#      - wd: circular mean (mean_wdir_deg())
wind_buf_mean_ts <- wind_buf_season %>%
  dplyr::group_by(date, season2) %>%
  dplyr::summarise(
    ws = mean(ws, na.rm = TRUE),
    wd = mean_wdir_deg(wd),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    station = "Buffer mean"  # synthetic station name
  )

# Optional: quick check
message("First rows of buffer-mean time series (with summer/winter flag):")
print(head(wind_buf_mean_ts))

## -------------------------------------------------------------------
## B8) Wind roses for buffer-mean, split by summer / winter
## -------------------------------------------------------------------

# B8.1 One wind rose per season (separate plots)
windRose(
  mydata = wind_buf_mean_ts,
  ws     = "ws",
  wd     = "wd",
  type   = "season2",   # separate panel for Summer/Winter
  angle  = 10,
  paddle = FALSE,
  key.position = "right"
)

# If you prefer separate calls instead of faceting, you can also do:
# windRose(
#   mydata = dplyr::filter(wind_buf_mean_ts, season2 == "Summer (Apr–Sep)"),
#   ws     = "ws",
#   wd     = "wd",
#   angle  = 10,
#   paddle = FALSE,
#   key.position = "right",
#   main   = "Buffer mean – Summer"
# )
#
# windRose(
#   mydata = dplyr::filter(wind_buf_mean_ts, season2 == "Winter (Oct–Mar)"),
#   ws     = "ws",
#   wd     = "wd",
#   angle  = 10,
#   paddle = FALSE,
#   key.position = "right",
#   main   = "Buffer mean – Winter"
# )

## -------------------------------------------------------------------
## B9) Seasonal summary tables for buffer-mean series
## -------------------------------------------------------------------

# B9.1 Overall seasonal statistics (buffer mean)
buffer_season_summary <- wind_buf_mean_ts %>%
  dplyr::group_by(season2) %>%
  dplyr::summarise(
    n       = n(),
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd),
    .groups = "drop"
  )


message("Seasonal wind summary for buffer-mean series:")
print(buffer_season_summary)


#######################################################################
# Seasonal and annual mean wind for Burgwald-buffer stations
# (computed from wind_buf_mean_ts = spatially averaged time series)
#######################################################################

# Ensure seasons exist
if(!all(c("Summer (Apr–Sep)", "Winter (Oct–Mar)") %in% unique(wind_buf_mean_ts$season2))) {
  stop("Seasons missing in wind_buf_mean_ts")
}

# --- Helper for circular mean again (safe) ---
mean_wdir_deg <- function(wd_deg) {
  theta <- wd_deg * pi / 180
  sin_mean <- mean(sin(theta), na.rm = TRUE)
  cos_mean <- mean(cos(theta), na.rm = TRUE)
  ang <- atan2(sin_mean, cos_mean) * 180 / pi
  if (ang < 0) ang <- ang + 360
  return(ang)
}

# ------------------------------
# 1) SUMMER mean wind
# ------------------------------
summer_mean_wind <- wind_buf_mean_ts %>%
  filter(season2 == "Summer (Apr–Sep)") %>%
  summarise(
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd),
    n       = n()
  )

# ------------------------------
# 2) WINTER mean wind
# ------------------------------
winter_mean_wind <- wind_buf_mean_ts %>%
  filter(season2 == "Winter (Oct–Mar)") %>%
  summarise(
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd),
    n       = n()
  )

# ------------------------------
# 3) ANNUAL mean wind (all values)
# ------------------------------
annual_mean_wind <- wind_buf_mean_ts %>%
  summarise(
    mean_ws = mean(ws, na.rm = TRUE),
    mean_wd = mean_wdir_deg(wd),
    n       = n()
  )

# ------------------------------
# 4) Combine in one clean table
# ------------------------------
wind_means_summary <- tibble::tibble(
  period = c("Summer", "Winter", "Annual"),
  mean_ws = c(summer_mean_wind$mean_ws,
              winter_mean_wind$mean_ws,
              annual_mean_wind$mean_ws),
  mean_wd = c(summer_mean_wind$mean_wd,
              winter_mean_wind$mean_wd,
              annual_mean_wind$mean_wd)
)

message("=== Mean Wind Summary (Burgwald Buffer) ===")
print(wind_means_summary)
#######################################################################


message("Direction-frequency table per station (stations inside Burgwald buffer):")
print(dir_table_by_station_buf)
