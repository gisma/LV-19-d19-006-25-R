#!/usr/bin/env Rscript

# -------------------------------------------------------------------
# DEM Processing Pipeline — Technical Header
#
# Purpose
# -------
# This script derives terrain and hydrological predictor layers from a
# high-resolution DEM using SAGA GIS via Rsagacmd. All numerical processing
# is executed in SAGA; R is used only for orchestration and metadata handling.
#
# The output consists of two canonical multiband GeoTIFF stacks:
#   (1) relief_stack_10m   — terrain morphology and topographic context
#   (2) hydro_stack_10m    — hydrological structure and connectivity
#
#
# Spatial Assumptions
# -------------------
# • Coordinate reference system is inherited from the input DEM.
# • All processing is performed on a uniform 10 m grid.
# • Grid alignment and extent are defined exclusively by the resampled DEM.
# • No reprojection or snapping is applied inside this workflow.
#
#
# Resampling Semantics
# --------------------
# • DEM resampling is performed using SAGA grid_tools$resampling.
# • scale_down = 1 activates resampling in the Rsagacmd binding.
# • target_user_size = 10 defines the target cell size (meters).
# • keep_type = 1 preserves the original raster data type.
# • Resampling is deterministic and produces a single canonical grid
#   (dem10m.sgrd) used by all downstream steps.
#
#
# Hydrological Preprocessing
# --------------------------
# • Sink filling is performed using the Wang & Liu XXL algorithm.
# • A minimum slope threshold of 0.1 is enforced to avoid flat artefacts.
# • All hydrological derivatives operate exclusively on the filled DEM.
#
#
# Morphometric Derivatives
# ------------------------
# The following terrain metrics are computed:
# • Slope
# • Aspect
# • Curvatures:
#     - General
#     - Profile
#     - Plan
#     - Tangential
#     - Longitudinal
#     - Cross-sectional
#     - Minimal
#     - Maximal
#     - Total
#     - Flowline
#
# Units:
# • Slope and aspect units are explicitly controlled via SAGA parameters.
#
#
# Topographic Position Index (TPI)
# --------------------------------
# • TPI is computed independently at three spatial radii:
#     - 30 m
#     - 50 m
#     - 100 m
# • Each radius represents a distinct terrain context scale.
# • No aggregation or smoothing between scales is applied.
#
#
# Hydrological Derivatives
# ------------------------
# The following hydrological layers are produced:
# • Flow accumulation (top-down)
# • Strahler stream order
# • Drainage basin ID
# • Overland flow distance to channel network
#
# Channel Definition:
# • The Strahler order grid is used directly as the channel mask.
# • Non-zero cells define the channel network.
# • No vector rasterization is involved in distance computation.
#
#
# Output Semantics
# ----------------
# • All intermediate products are written as SAGA grids in a temporary
#   working directory and deleted after completion.
# • Final products are exported as multiband GeoTIFF using SAGA io_gdal.
# • Band names are assigned post-export for semantic clarity.
# • Pixel values are not modified during metadata renaming.
#
#
# Determinism and Reproducibility
# -------------------------------
# • The pipeline is deterministic given identical input DEM and parameters.
# • All intermediate filenames are fixed and externally inspectable.
# • No randomization, tiling heuristics, or adaptive rescaling is applied.
#
# -------------------------------------------------------------------


suppressPackageStartupMessages({
  library(here)
  library(Rsagacmd)
})

# Project setup: paths registry, helper functions, environment
source(here::here("src", "_core", "01-setup-burgwald.R"))

quiet <- FALSE
msg <- function(...) if (!quiet) message(...)

# -------------------------------------------------------------------
# Registry I/O (canonical keys, no new keys introduced)
# -------------------------------------------------------------------

# Input DEM (AOI extent)
dem_in           <- paths[["aoi_dgm"]]

# Output stacks (final GeoTIFF products)
out_relief_stack <- paths[["relief_stack_10m"]]
out_hydro_stack  <- paths[["hydro_stack_10m"]]

# -------------------------------------------------------------------
# Working directory for all intermediate SAGA grids
# -------------------------------------------------------------------

# All intermediate .sgrd/.sdat files are written here and deleted at the end
work_dir <- file.path(dirname(out_flowacc))
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# SAGA session
# -------------------------------------------------------------------

# SAGA runs all computations.
# The terra backend is used only as a GDAL adapter inside Rsagacmd,
# not for raster processing in the workflow itself.
saga <- Rsagacmd::saga_gis(
  raster_backend = "terra",
  vector_backend = "sf",
  all_outputs    = FALSE
)

# -------------------------------------------------------------------
# Canonical intermediate grid filenames (stable and inspectable)
# -------------------------------------------------------------------

# DEM chain
dem10m_sgrd      <- file.path(work_dir, "dem10m.sgrd")
dem10m_filled    <- file.path(work_dir, "dem10m_filled.sgrd")

# Morphometry outputs
slope_sgrd       <- file.path(work_dir, "Slope.sgrd")
aspect_sgrd      <- file.path(work_dir, "Aspect.sgrd")
c_gene_sgrd      <- file.path(work_dir, "GeneralCurvature.sgrd")
c_prof_sgrd      <- file.path(work_dir, "ProfileCurvature.sgrd")
c_plan_sgrd      <- file.path(work_dir, "PlanCurvature.sgrd")
c_tang_sgrd      <- file.path(work_dir, "TangentialCurvature.sgrd")
c_long_sgrd      <- file.path(work_dir, "LongitudinalCurvature.sgrd")
c_cros_sgrd      <- file.path(work_dir, "CrossSectionalCurvature.sgrd")
c_mini_sgrd      <- file.path(work_dir, "MinimalCurvature.sgrd")
c_maxi_sgrd      <- file.path(work_dir, "MaximalCurvature.sgrd")
c_tota_sgrd      <- file.path(work_dir, "TotalCurvature.sgrd")
c_roto_sgrd      <- file.path(work_dir, "FlowLineCurvature.sgrd")

# Multi-scale TPI outputs
tpi_r30_sgrd     <- file.path(work_dir, "tpi_r30.sgrd")
tpi_r50_sgrd     <- file.path(work_dir, "tpi_r50.sgrd")
tpi_r100_sgrd    <- file.path(work_dir, "tpi_r100.sgrd")

# Hydrology outputs
flowacc_sgrd     <- file.path(work_dir, "flowacc_10m.sgrd")
order_sgrd       <- file.path(work_dir, "strahler_order_10m.sgrd")
basin_sgrd       <- file.path(work_dir, "watershed_id_10m.sgrd")
dist_sgrd        <- file.path(work_dir, "dist_to_stream_10m.sgrd")

# -------------------------------------------------------------------
# 0) DEM resampling to 10 m grid
# -------------------------------------------------------------------
# - scale_down = 1 activates resampling in this Rsagacmd binding
# - target_user_size = 10 defines the target cell size in meters
# - keep_type = 1 preserves the original data type

msg("Import DEM to SAGA grid: ", dem10m_sgrd)
saga$grid_tools$resampling(
  input      = dem_in,
  output     = dem10m_sgrd,
  keep_type  = 1,
  scale_down = 1,
  target_user_size = 10,
  .verbose   = TRUE
)

# -------------------------------------------------------------------
# 1) Sink filling (hydrologic preprocessing)
# -------------------------------------------------------------------
# Produces a hydraulically corrected DEM for all downstream analyses.

msg("Fill sinks: ", dem10m_filled)
saga$ta_preprocessor$fill_sinks_xxl_wang_liu(
  elev     = dem10m_sgrd,
  filled   = dem10m_filled,
  minslope = 0.1,
  .verbose = TRUE
)

# -------------------------------------------------------------------
# 2) Morphometric derivatives
# -------------------------------------------------------------------
# Computes slope, aspect and multiple curvature metrics from the filled DEM.

msg("Morphometry (slope/aspect/curvature)")
saga$ta_morphometry$slope_aspect_curvature(
  elevation = dem10m_filled,
  slope     = slope_sgrd,
  aspect    = aspect_sgrd,
  c_gene    = c_gene_sgrd,
  c_prof    = c_prof_sgrd,
  c_plan    = c_plan_sgrd,
  c_tang    = c_tang_sgrd,
  c_long    = c_long_sgrd,
  c_cros    = c_cros_sgrd,
  c_mini    = c_mini_sgrd,
  c_maxi    = c_maxi_sgrd,
  c_tota    = c_tota_sgrd,
  c_roto    = c_roto_sgrd,
  method       = 6,
  unit_slope   = 1,
  unit_aspect  = 1,
  .verbose     = TRUE
)

# -------------------------------------------------------------------
# 3) Multi-scale Topographic Position Index (TPI)
# -------------------------------------------------------------------
# Three radii are evaluated independently to capture terrain context
# at different spatial scales.

msg("TPI r=30/50/100")
saga$ta_morphometry$topographic_position_index_tpi(
  dem        = dem10m_filled,
  tpi        = tpi_r30_sgrd,
  radius_min = 30,
  .verbose   = TRUE
)
saga$ta_morphometry$topographic_position_index_tpi(
  dem        = dem10m_filled,
  tpi        = tpi_r50_sgrd,
  radius_min = 50,
  .verbose   = TRUE
)
saga$ta_morphometry$topographic_position_index_tpi(
  dem        = dem10m_filled,
  tpi        = tpi_r100_sgrd,
  radius_min = 100,
  .verbose   = TRUE
)

# -------------------------------------------------------------------
# 4) Flow accumulation
# -------------------------------------------------------------------
# Computes total contributing area using the filled DEM.

msg("Flow accumulation")
saga$ta_hydrology$flow_accumulation_top_down(
  elevation  = dem10m_filled,
  accu_total = flowacc_sgrd,
  .verbose   = TRUE
)

# -------------------------------------------------------------------
# 5) Channel network and drainage basins
# -------------------------------------------------------------------
# Generates:
#   - Strahler order grid
#   - Basin ID grid
#   - Optional vector outputs (channels, drainage basins)

msg("Channels+Basins (GRIDS only)")
saga$ta_channels$channel_network_and_drainage_basins(
  dem       = dem10m_filled,
  threshold = 5,
  order     = order_sgrd,
  basin     = basin_sgrd,
  segments  = file.path(work_dir, "Channels.gpkg"),
  basins    = file.path(work_dir, "DrainageBasins.gpkg"),
  .verbose  = TRUE
)

# -------------------------------------------------------------------
# 6) Distance to channel network
# -------------------------------------------------------------------
# Uses the Strahler order grid directly as the channel mask.
# Non-zero cells are interpreted as channel network.

msg("Distance to channel network (channels = order grid)")
saga$ta_channels$overland_flow_distance_to_channel_network(
  elevation = dem10m_filled,
  channels  = order_sgrd,
  distance  = dist_sgrd,
  method    = 1,
  boundary  = FALSE,
  flow_k_default = 20,
  flow_r_default = 0.05,
  .verbose  = TRUE
)

# -------------------------------------------------------------------
# 7) Export final products to canonical GeoTIFF stacks
# -------------------------------------------------------------------
# Export is performed by SAGA (io_gdal), not by terra.

msg("Export Hydro Stack to: ", dirname(out_flowacc))
saga$io_gdal$export_raster(
  grids = list(flowacc_sgrd, order_sgrd, basin_sgrd, dist_sgrd),
  multiple   = 1,
  file       = out_hydro_stack,
  format     = 1,
  type       = 0,
  set_nodata = 1,
  nodata     = -99999,
  .verbose   = TRUE
)

msg("Export Topo Stack to: ", dirname(out_relief_stack))
saga$io_gdal$export_raster(
  grids = vars <- list(
    dem10m_sgrd,
    dem10m_filled,
    slope_sgrd,
    aspect_sgrd,
    c_gene_sgrd,
    c_prof_sgrd,
    c_plan_sgrd,
    c_tang_sgrd,
    c_long_sgrd,
    c_cros_sgrd,
    c_mini_sgrd,
    c_maxi_sgrd,
    c_tota_sgrd,
    c_roto_sgrd,
    tpi_r30_sgrd,
    tpi_r50_sgrd,
    tpi_r100_sgrd
  ),
  multiple   = 1,
  file       = out_relief_stack,
  format     = 1,
  type       = 0,
  set_nodata = 1,
  nodata     = -99999,
  .verbose   = TRUE
)

# -------------------------------------------------------------------
# 8) Band renaming (metadata only)
# -------------------------------------------------------------------
# GeoTIFF bands are renamed for semantic clarity.
# This does not alter pixel values.

r <- terra::rast(out_relief_stack)
names(r) <- c(
  "dem10m", "dem10m_filled", "slope", "aspect",
  "c_gene", "c_prof", "c_plan", "c_tang",
  "c_long", "c_cros", "c_mini", "c_maxi",
  "c_tota", "c_roto",
  "tpi_r30", "tpi_r50", "tpi_r100"
)
tmp_file <- paste0(out_relief_stack, ".tmp")
terra::writeRaster(r, tmp_file, overwrite = TRUE)
file.rename(tmp_file, out_relief_stack)

r <- terra::rast(out_hydro_stack)
names(r) <- c("flowacc", "strahler", "watershed", "dist_to_water")
tmp_file <- paste0(out_hydro_stack, ".tmp")
terra::writeRaster(r, tmp_file, overwrite = TRUE)
file.rename(tmp_file, out_hydro_stack)

# -------------------------------------------------------------------
# 9) Cleanup of intermediate SAGA artifacts
# -------------------------------------------------------------------
# Removes all temporary SAGA grids and auxiliary files.

delfiles = list.files(
  path    = dirname(out_hydro_stack),
  pattern = "\\.(sgrd|mgrd|sdat|xml|prj)$",
  full.names = TRUE
)
unlink(delfiles)

msg("Done.")
