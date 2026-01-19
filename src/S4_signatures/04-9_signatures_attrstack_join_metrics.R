#!/usr/bin/env Rscript
############################################################
# Script:  04-9_signatures_attrstack_join_metrics.R
# Project: Burgwald
#
# Purpose:
#   Build the full segment signature attrstack (S4_signatures) by
#   joining all S4 signature layers onto the stable segment geometry.
#
#   - NO computation of new metrics
#   - NO decisions
#   - NO clustering
#
# Inputs (productive, from outputs.tsv via paths[]):
#   - layer0_segments                    (S3_structure, gpkg)  [geometry backbone]
#   - layer0_attr_it_metrics             (S4_signatures, gpkg)
#   - layer0_attr_physio_metrics         (S4_signatures, gpkg)
#   - layer0_attr_hydro_metrics          (S4_signatures, gpkg)
#   - layer0_attr_cover_metrics          (S4_signatures, gpkg)
#   - layer0_attr_biostructure_metrics   (S4_signatures, gpkg)
#
# Output (productive, from outputs.tsv via paths[]):
#   - layer0_segments_attrstack_metrics  (S4_signatures, gpkg)
#
# Contract:
#   - All inputs must have 'segment_id'
#   - Geometry comes ONLY from layer0_segments
############################################################

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(dplyr)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))
source(here::here("src", "r-libs", "metrics-fun.R"))

## -------------------------------------------------------------------
## 1) Productive input (from outputs.tsv)
## -------------------------------------------------------------------

seg_file    <- paths[["layer0_segments"]]
it_file     <- paths[["layer0_attr_it_metrics"]]
physio_file <- paths[["layer0_attr_physio_metrics"]]
hydro_file  <- paths[["layer0_attr_hydro_metrics"]]
cover_file  <- paths[["layer0_attr_cover_metrics"]]
bio_file    <- paths[["layer0_attr_biostructure_metrics"]]

stopifnot(file.exists(seg_file))
stopifnot(file.exists(it_file))
stopifnot(file.exists(physio_file))
stopifnot(file.exists(hydro_file))
stopifnot(file.exists(cover_file))
stopifnot(file.exists(bio_file))

## -------------------------------------------------------------------
## 2) Productive outputs (from outputs.tsv)
## -------------------------------------------------------------------

out_file <- paths[["layer0_segments_attrstack_metrics"]]

## -------------------------------------------------------------------
## 2b) tmp output folder (NON-productive)
## -------------------------------------------------------------------

tmp_root <- here::here("data", "tmp", "layer0_segments", "tmp_attrstack_join")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Read inputs
## -------------------------------------------------------------------

segs_sf <- sf::read_sf(seg_file)
stopifnot("segment_id" %in% names(segs_sf))

it_sf     <- sf::read_sf(it_file)
physio_sf <- sf::read_sf(physio_file)
hydro_sf  <- sf::read_sf(hydro_file)
cover_sf  <- sf::read_sf(cover_file)
bio_sf    <- sf::read_sf(bio_file)

stopifnot("segment_id" %in% names(it_sf))
stopifnot("segment_id" %in% names(physio_sf))
stopifnot("segment_id" %in% names(hydro_sf))
stopifnot("segment_id" %in% names(cover_sf))
stopifnot("segment_id" %in% names(bio_sf))

## -------------------------------------------------------------------
## 4) Drop geometries from signature layers (keep only attributes)
## -------------------------------------------------------------------

it_df     <- sf::st_drop_geometry(it_sf)
physio_df <- sf::st_drop_geometry(physio_sf)
hydro_df  <- sf::st_drop_geometry(hydro_sf)
cover_df  <- sf::st_drop_geometry(cover_sf)
bio_df    <- sf::st_drop_geometry(bio_sf)

## -------------------------------------------------------------------
## 5) Join onto backbone geometry
## -------------------------------------------------------------------

out_sf <- segs_sf %>%
  dplyr::select(segment_id, geometry) %>%
  dplyr::left_join(it_df,     by = "segment_id") %>%
  dplyr::left_join(physio_df, by = "segment_id") %>%
  dplyr::left_join(hydro_df,  by = "segment_id") %>%
  dplyr::left_join(cover_df,  by = "segment_id") %>%
  dplyr::left_join(bio_df,    by = "segment_id")

## -------------------------------------------------------------------
## 6) Write productive output
## -------------------------------------------------------------------

out_dir <- dirname(out_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sf::st_write(out_sf, out_file, delete_dsn = TRUE, quiet = TRUE)
