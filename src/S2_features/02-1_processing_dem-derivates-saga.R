#!/usr/bin/env Rscript

# -------------------------------------------------------------------
# DEM Processing Pipeline — Technical Header
# Script:  02-1_processing_dem-derivates-saga.R
# Project:  Burgwald
#
# Purpose
# -------
# Derive terrain + hydrology predictor layers from a DEM using SAGA GIS
# (via Rsagacmd). All raster computations are executed in SAGA; R only
# orchestrates module calls and canonical I/O via `paths`.
#
# Final outputs (canonical GeoTIFF stacks)
# ----------------------------------------
# 1) out_relief_stack (relief_stack_10m):
#    Terrain morphology + topographic context layers on a 10 m grid:
#      - dem10m, dem10m_filled
#      - slope, aspect
#      - curvatures: c_gene, c_prof, c_plan, c_tang, c_long, c_cros,
#                   c_mini, c_maxi, c_tota, c_roto
#      - TPI: tpi_r30, tpi_r50, tpi_r100
#
# 2) out_hydro_stack (hydro_stack_10m):
#    Hydrological structure/connectivity layers on a 10 m grid:
#      - flowacc (flow accumulation)
#      - strahler (Strahler order)
#      - watershed (drainage basin id)
#      - dist_to_water (overland flow distance to channel network)
#
# -------------------------------------------
# All intermediate products are written as SAGA grid pairs in work_dir:
#   *.sgrd (header) + *.sdat (binary data) + aux files (*.prj, *.xml, ...)
# Filenames are fixed and inspectable for debugging (stable naming). 
# However by default they will be removed.
#
# Dependencies / requirements
# ---------------------------
# R:
#   - here     : project-relative paths
#   - Rsagacmd : R interface to SAGA GIS (calls saga_cmd under the hood)
# External:
#   - SAGA GIS installed and callable (saga_cmd available)
#   - GDAL stack for reading/writing GeoTIFF/NetCDF as used by SAGA io_gdal
#
# Spatial / grid semantics (as implemented here)
# ----------------------------------------------
# - The workflow defines a uniform 10 m grid by resampling the input DEM
#   with SAGA grid_tools$resampling.
# - All downstream steps use that resampled grid (dem10m.sgrd) and its
#   sink-filled version (dem10m_filled.sgrd).
# - No reprojection is performed here; CRS is inherited from dem_in.
#
# Resampling (important for reproducibility)
# ------------------------------------------
# The grid_tools$resampling call defines the canonical 10 m grid:
#   - scale_down = 1  : activates downscaling/resampling (binding semantics)
#   - target_user_size = 10 : target cell size in map units (meters in EPSG:25832)
#   - keep_type = 1   : preserve original data type
# This produces dem10m.sgrd which anchors extent, alignment, resolution.
#
# Hydrology preprocessing
# -----------------------
# Sink filling is performed on dem10m with:
#   ta_preprocessor$fill_sinks_xxl_wang_liu
# Parameter:
#   - minslope = 0.1 enforces a minimum slope to avoid flat artefacts.
# All hydrological derivatives operate on dem10m_filled (not raw dem10m).
#
# Channel definition for distance computation
# -------------------------------------------
# The Strahler order grid is used directly as channel mask:
#   - channels = order_sgrd
# Interpretation:
#   - non-zero cells define the channel network
# No vector rasterization is needed for the distance tool call.
#
# Output semantics
# ----------------
# - Export to GeoTIFF stacks is done using SAGA io_gdal$export_raster
#   (not terra). This keeps the processing chain SAGA-centric.
# - terra is used only at the end to assign band names (metadata) and
#   rewrite the file to persist those names.
# - Intermediate SAGA grids in work_dir are deleted at the end.
# -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(Rsagacmd)
})

# Project setup:
# - loads `paths` registry (from outputs.tsv)
# - provides helper functions used across the project
source(here::here("src", "_core", "01-setup-burgwald.R"))

quiet <- FALSE
msg <- function(...) if (!quiet) message(...)

# -------------------------------------------------------------------
# Registry I/O (canonical keys, no new keys)
# -------------------------------------------------------------------

# Input DEM: AOI DEM (usually DEM mosaic+clip from earlier pipeline steps)
dem_in           <- paths[["aoi_dgm"]]

# Final output stacks (multiband GeoTIFF)
out_relief_stack <- paths[["relief_stack_10m"]]
out_hydro_stack  <- paths[["hydro_stack_10m"]]

# -------------------------------------------------------------------
# Working directory for intermediate SAGA grids
# -------------------------------------------------------------------
# NOTE: All intermediate .sgrd/.sdat products are written here.
# They are removed at the end of the script (cleanup block).
#
# IMPORTANT: This line references `out_hydro_stack`, which must exist in the
# environment/registry in your project context. In your earlier working
# version, `work_dir` was defined as `dirname(out_hydro_stack)` to place the
# intermediate grids next to the canonical hydrology outputs.
work_dir <- file.path(dirname(out_hydro_stack))
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# SAGA session
# -------------------------------------------------------------------
# - raster_backend="terra" here is only the adapter layer inside Rsagacmd.
#   The actual computations are performed by SAGA modules via saga_cmd.
# - vector_backend="sf" enables vector output handling if modules write shapes.
saga <- Rsagacmd::saga_gis(
  raster_backend = "terra",
  vector_backend = "sf",
  all_outputs    = FALSE
)

# -------------------------------------------------------------------
# Canonical intermediate grid filenames (fixed naming)
# -------------------------------------------------------------------

# DEM chain
dem10m_sgrd      <- file.path(work_dir, "dem10m.sgrd")         # resampled 10 m DEM
dem10m_filled    <- file.path(work_dir, "dem10m_filled.sgrd")  # sink-filled 10 m DEM

# Morphometry outputs (names match your canonical SAGA naming convention)
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

# Multi-scale TPI outputs (independent scales, no aggregation)
tpi_r30_sgrd     <- file.path(work_dir, "tpi_r30.sgrd")
tpi_r50_sgrd     <- file.path(work_dir, "tpi_r50.sgrd")
tpi_r100_sgrd    <- file.path(work_dir, "tpi_r100.sgrd")

# Hydrology outputs
flowacc_sgrd     <- file.path(work_dir, "flowacc_10m.sgrd")
order_sgrd       <- file.path(work_dir, "strahler_order_10m.sgrd")
basin_sgrd       <- file.path(work_dir, "watershed_id_10m.sgrd")
dist_sgrd        <- file.path(work_dir, "dist_to_stream_10m.sgrd")

# -------------------------------------------------------------------
# 0) DEM resampling to 10 m grid (SAGA)
# -------------------------------------------------------------------
# Input:
#   - dem_in (from paths registry)
# Output:
#   - dem10m_sgrd (canonical 10 m SAGA grid)
# Parameters (binding-specific semantics you already validated):
#   - keep_type = 1          keep original type (avoid implicit conversion)
#   - scale_down = 1         enable downscale/resampling
#   - target_user_size = 10  target cell size
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
# 1) Sink filling (SAGA)
# -------------------------------------------------------------------
# Input:
#   - dem10m_sgrd
# Output:
#   - dem10m_filled
# Parameter:
#   - minslope=0.1 ensures hydraulic connectivity (avoids flat artefacts)
msg("Fill sinks: ", dem10m_filled)
saga$ta_preprocessor$fill_sinks_xxl_wang_liu(
  elev     = dem10m_sgrd,
  filled   = dem10m_filled,
  minslope = 0.1,
  .verbose = TRUE
)

# -------------------------------------------------------------------
# 2) Morphometry: slope/aspect/curvature (SAGA)
# -------------------------------------------------------------------
# Input:
#   - dem10m_filled (not raw dem10m)
# Outputs:
#   - slope_sgrd, aspect_sgrd
#   - all curvature grids listed above
# Parameters:
#   - method=6      SAGA morphometry method option used in your workflow
#   - unit_slope=1  slope unit control (SAGA parameter)
#   - unit_aspect=1 aspect unit control (SAGA parameter)
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
# 3) Multi-scale Topographic Position Index (TPI) (SAGA)
# -------------------------------------------------------------------
# Input:
#   - dem10m_filled
# Outputs:
#   - tpi_r30_sgrd, tpi_r50_sgrd, tpi_r100_sgrd
# Parameter:
#   - radius_min defines the spatial neighborhood radius (map units)
# Note:
#   - each radius is computed independently (no merging between scales)
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
# 4) Flow accumulation (SAGA)
# -------------------------------------------------------------------
# Input:
#   - dem10m_filled
# Output:
#   - flowacc_sgrd
# Tool:
#   ta_hydrology$flow_accumulation_top_down
msg("Flow accumulation")
saga$ta_hydrology$flow_accumulation_top_down(
  elevation  = dem10m_filled,
  accu_total = flowacc_sgrd,
  .verbose   = TRUE
)

# -------------------------------------------------------------------
# 5) Channel network and drainage basins (SAGA)
# -------------------------------------------------------------------
# Input:
#   - dem10m_filled
# Outputs:
#   - order_sgrd (Strahler order grid)
#   - basin_sgrd (basin id grid)
# Plus (optional) vector outputs written into work_dir:
#   - Channels.gpkg
#   - DrainageBasins.gpkg
# Parameter:
#   - threshold controls channel initiation (your workflow value here: 5)
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
# 6) Overland flow distance to channel network (SAGA)
# -------------------------------------------------------------------
# Inputs:
#   - elevation = dem10m_filled
#   - channels  = order_sgrd (Strahler order grid as channel mask)
# Output:
#   - dist_sgrd
# Parameters:
#   - method=1, boundary=FALSE, flow_k_default=20, flow_r_default=0.05
# Semantics:
#   - channels mask is derived from non-zero cells of `order_sgrd`
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
# 7) Export final products to canonical GeoTIFF stacks (SAGA io_gdal)
# -------------------------------------------------------------------
# Export is done by SAGA (io_gdal) to write multiband GeoTIFF stacks.
# - `multiple = 1` writes a multiband raster to `file`.
# - `format = 1` selects GeoTIFF in this binding.
# - nodata handling is forced to -99999 to keep missing values explicit.
msg("Export Hydro Stack to: ", dirname(out_hydro_stack))
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

r = rast(out_hydro_stack)
names(r) = c("flowacc", "strahler", "watershed", "dist_to_water")
tmp_file <- paste0(out_hydro_stack, ".tmp")
terra::writeRaster(r, tmp_file, filetype= "GTiff",overwrite = TRUE)
file.rename(tmp_file, out_hydro_stack)


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
# 8) Band renaming (metadata only; values unchanged)
# -------------------------------------------------------------------
# GeoTIFF written by SAGA gets band names set via terra to provide
# stable semantic layer naming downstream (signatures, modeling, plots).
r <- terra::rast(out_relief_stack)

# existing canonical names
nm <- c(
  "dem10m", "dem10m_filled", "slope", "aspect",
  "c_gene", "c_prof", "c_plan", "c_tang",
  "c_long", "c_cros", "c_mini", "c_maxi",
  "c_tota", "c_roto",
  "tpi_r30", "tpi_r50", "tpi_r100"
)
names(r) <- nm

# Southness (0..1): cos(aspect - 180°) scaled to [0,1]
# aspect is degrees
aspect_rad <- r[["aspect"]] * pi / 180
southness  <- cos(aspect_rad - pi)        # [-1..1]
southness01 <- (southness + 1) / 2        # [0..1]
names(southness01) <- "southness"

# append derived band
r2 <- c(r, southness01)

tmp_file <- paste0(out_relief_stack, ".tmp")
terra::writeRaster(r2, tmp_file, filetype= "GTiff",overwrite = TRUE)
file.rename(tmp_file, out_relief_stack)

# -------------------------------------------------------------------
# 9) Cleanup: remove intermediate SAGA artifacts
# -------------------------------------------------------------------
# Removes SAGA grid pairs and sidecar files from the output directory.
# This keeps the repo clean and prevents accidental Git tracking of
# intermediate grids (especially large .sdat files).
delfiles = list.files(
  path    = dirname(out_hydro_stack),
  pattern = "\\.(sgrd|mgrd|sdat|xml|prj)$",
  full.names = TRUE
)
unlink(delfiles)

msg("Done.")
