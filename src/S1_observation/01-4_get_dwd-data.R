#!/usr/bin/env Rscript

############################################################
# Script:  01-4-get_DWD-data.R
# Author:  [Your Name]
# Project: [Your Project Name]
#
# Purpose
# -------
# Thin wrapper script that **controls the time window** for DWD data
# and **calls the Burgwald-specific helper functions** to download
# and preprocess:
#
#   1) Hourly wind (FF = speed, DD = direction)
#   2) Hourly precipitation (R1)
#   3) Sub-hourly precipitation (10-min and 5-min)
#
# All heavy lifting (station selection by AOI, downloads, unzipping,
# time filtering, joining metadata, writing CSV/RDS) is done inside:
#
#   - burgwald_get_hourly_dwd()
#   - burgwald_get_subhourly_precip()
#
# This script is therefore the **orchestrator** for DWD forcing data:
# you only set the period (`start_date`, `end_date`) and then call the
# high-level functions once.
#
# Output (handled by the helper functions)
# ---------------------------------------
# - Raw and processed DWD station files under:
#       data/raw/dwd-stations/
#       data/processed/dwd-stations/
#
# Prerequisites
# -------------
# - 00-setup-burgwald.R has been run and:
#     * loads all required packages
#     * defines aoi_burgwald_wgs and project paths
#     * defines burgwald_get_hourly_dwd()
#       and burgwald_get_subhourly_precip()
############################################################


# Source the global setup.
# This must define:
#   - aoi_burgwald_wgs
#   - path / directory settings
#   - the DWD helper functions used below.
source("src","_core/00-setup-burgwald.R")


# ----------------------------------------------------------
# 1) Define the DWD download period
# ----------------------------------------------------------
# Here: last 2 years, from 00:00:00 on start_date
#       to 23:59:59 on end_date (inclusive day range).
start_date <- Sys.Date() - 365 * 2
end_date   <- Sys.Date()

# POSIX versions (UTC) – useful if you need exact timestamps
# for further processing or logging.
start_dt <- as.POSIXct(paste0(start_date, " 00:00:00"), tz = "UTC")
end_dt   <- as.POSIXct(paste0(end_date,   " 23:59:59"), tz = "UTC")



# ----------------------------------------------------------
# 2) Hourly wind: FF (speed), DD (direction)
# ----------------------------------------------------------
FD_bw <- burgwald_get_hourly_dwd(
  var        = "wind",           # hourly wind product
  params     = c("FF", "DD"),    # FF -> F (speed), DD -> D (direction)
  start_date = start_date,
  end_date   = end_date
)
# → Returns an sf object with wind time series for all stations
#   inside Burgwald AOI + buffer; also writes CSV/RDS to disk
#   (depending on helper defaults).


# ----------------------------------------------------------
# 3) Hourly precipitation: R1
# ----------------------------------------------------------
rr_60min_bw <- burgwald_get_hourly_dwd(
  var        = "precipitation",  # hourly precipitation product
  params     = "R1",             # hourly rain sum (mm)
  start_date = start_date,
  end_date   = end_date
)
# → Again: sf object with hourly precipitation for AOI stations,
#   plus CSV/RDS written by the helper.


# ----------------------------------------------------------
# 4) Sub-hourly precipitation: 10-minute sums
# ----------------------------------------------------------
rr_10min <- burgwald_get_subhourly_precip(
  resolution = "10min",          # 10-minute precipitation product
  start_date = start_date,
  end_date   = end_date
)
# → sf object with 10-min rain sums (RWS_10) for AOI stations.


# ----------------------------------------------------------
# 5) Sub-hourly precipitation: 5-minute sums
# ----------------------------------------------------------
rr_5min <- burgwald_get_subhourly_precip(
  resolution = "5min",           # 5-minute precipitation product
  start_date = start_date,
  end_date   = end_date
)
# → sf object with 5-min rain sums (RS_05) for AOI stations.

# At this point, all relevant DWD forcing data for the chosen
# time period are available on disk and in memory (FD_bw,
# rr_60min_bw, rr_10min, rr_5min) for further modelling.
