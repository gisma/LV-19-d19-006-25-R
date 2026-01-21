#!/usr/bin/env Rscript

# src/02-1_processing_dem-derivates-saga.R
#
# Terrain and hydrological derivatives using SAGA GIS (via Rsagacmd).
#
# Design principles:
# - All SAGA tools are called with EXPLICIT output file arguments.
# - Intermediate products are passed as FILE PATHS, not as terra objects.
# - SAGA writes *.sgrd/*.sdat pairs; terra always reads the *.sdat side.
# - Only final products are exported as GeoTIFF stacks.
#
# This avoids implicit wrapper behaviour, version-dependent naming,
# and terra <-> SAGA ping-pong artefacts.

source("src/_core/01-setup-burgwald.R")

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(Rsagacmd)
})

# ---------------------------------------------------------------------------
# Global SAGA instance
# ---------------------------------------------------------------------------

saga <- Rsagacmd::saga_gis(
  raster_backend = "terra",
  vector_backend = "sf",
  all_outputs    = FALSE
)

# ===========================================================================
# RELIEF DERIVATIVES
# ===========================================================================

#' Derive relief layers from a DEM and write a multi-band GeoTIFF stack.
#'
#' Processing chain:
#'   1) Optional resampling
#'   2) Sink filling (Wang & Liu)
#'   3) Slope / aspect / curvature (explicit outputs)
#'   4) Multi-scale TPI (explicit outputs)
#'   5) Stack assembly and export
#'
#' All SAGA outputs are written as *.sgrd/*.sdat and only loaded into terra
#' once they exist on disk.

saga_derive_relief_stack <- function(dem_1m,
                                     prefix,
                                     out_stack_file,
                                     outdir              = dirname(out_stack_file),
                                     target_res_m        = NULL,
                                     resample_method     = c("bilinear", "near"),
                                     minslope            = 0.1,
                                     tpi_scales_m        = c(30, 50, 250),
                                     saga_obj            = NULL) {
  
  resample_method <- match.arg(resample_method)
  
  if (is.null(saga_obj)) {
    saga_obj <- Rsagacmd::saga_gis(
      raster_backend = "terra",
      vector_backend = "sf",
      all_outputs    = FALSE
    )
  }
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # -------------------------------------------------------------------------
  # 0) Optional resampling of the base DEM
  # -------------------------------------------------------------------------
  dem <- dem_1m
  
  if (!is.null(target_res_m)) {
    tpl <- terra::rast(
      extent     = terra::ext(dem_1m),
      crs        = terra::crs(dem_1m),
      resolution = target_res_m
    )
    dem <- terra::resample(dem_1m, tpl, method = resample_method)
  }
  
  # -------------------------------------------------------------------------
  # 1) Fill sinks (Wang & Liu XXL)
  # -------------------------------------------------------------------------
  dem_filled_file <- file.path(outdir, paste0(prefix, "_dem_filled.sgrd"))
  
  saga_obj$ta_preprocessor$fill_sinks_xxl_wang_liu(
    elev     = dem,
    filled   = dem_filled_file,
    minslope = minslope
  )
  
  stopifnot(file.exists(dem_filled_file))
  
  dem_filled <- terra::rast(sub("\\.sgrd$", ".sdat", dem_filled_file))
  
  # -------------------------------------------------------------------------
  # 2) Slope / aspect / curvatures (all outputs explicit)
  # -------------------------------------------------------------------------
  lsps_files <- list(
    slope  = file.path(outdir, paste0(prefix, "_slope.sgrd")),
    aspect = file.path(outdir, paste0(prefix, "_aspect.sgrd")),
    c_gene = file.path(outdir, paste0(prefix, "_c_gene.sgrd")),
    c_prof = file.path(outdir, paste0(prefix, "_c_prof.sgrd")),
    c_plan = file.path(outdir, paste0(prefix, "_c_plan.sgrd")),
    c_tang = file.path(outdir, paste0(prefix, "_c_tang.sgrd")),
    c_long = file.path(outdir, paste0(prefix, "_c_long.sgrd")),
    c_cros = file.path(outdir, paste0(prefix, "_c_cros.sgrd")),
    c_mini = file.path(outdir, paste0(prefix, "_c_mini.sgrd")),
    c_maxi = file.path(outdir, paste0(prefix, "_c_maxi.sgrd")),
    c_tota = file.path(outdir, paste0(prefix, "_c_tota.sgrd")),
    c_roto = file.path(outdir, paste0(prefix, "_c_roto.sgrd"))
  )
  
  saga_obj$ta_morphometry$slope_aspect_curvature(
    elevation   = dem_filled_file,
    slope        = lsps_files$slope,
    aspect       = lsps_files$aspect,
    c_gene       = lsps_files$c_gene,
    c_prof       = lsps_files$c_prof,
    c_plan       = lsps_files$c_plan,
    c_tang       = lsps_files$c_tang,
    c_long       = lsps_files$c_long,
    c_cros       = lsps_files$c_cros,
    c_mini       = lsps_files$c_mini,
    c_maxi       = lsps_files$c_maxi,
    c_tota       = lsps_files$c_tota,
    c_roto       = lsps_files$c_roto,
    method       = 6,
    unit_slope   = 0,
    unit_aspect  = 0,
    .verbose     = TRUE
  )
  
  stopifnot(all(file.exists(unlist(lsps_files))))
  
  lsps_stack <- terra::rast(unlist(sub("\\.sgrd$", ".sdat", lsps_files)))
  
  # -------------------------------------------------------------------------
  # 3) Multi-scale TPI (explicit outputs)
  # -------------------------------------------------------------------------
  tpi_files <- setNames(
    file.path(outdir, paste0(prefix, "_tpi_r", tpi_scales_m, "m.sgrd")),
    paste0("tpi_r", tpi_scales_m, "m")
  )
  
  for (i in seq_along(tpi_scales_m)) {
    saga_obj$ta_morphometry$topographic_position_index_tpi(
      dem        = dem_filled_file,
      tpi        = tpi_files[i],
      radius_min = 0,
      radius_max = tpi_scales_m[i],
      standard   = TRUE,
      .verbose   = TRUE
    )
    stopifnot(file.exists(tpi_files[i]))
  }
  
  tpi_stack <- terra::rast(unlist(sub("\\.sgrd$", ".sdat", tpi_files)))
  
  # -------------------------------------------------------------------------
  # 4) Assemble and export relief stack
  # -------------------------------------------------------------------------
  out_stack <- c(dem_filled, lsps_stack, tpi_stack)
  names(out_stack) <- make.names(names(out_stack), unique = TRUE)
  
  terra::writeRaster(out_stack, out_stack_file, overwrite = TRUE)
  invisible(out_stack_file)
}

# ===========================================================================
# HYDROLOGY DERIVATIVES
# ===========================================================================

#' Derive hydrological layers from a DEM and export both single products
#' and a combined hydro stack.

saga_derive_hydro_rasters <- function(
    dem_1m,
    out_flowacc_file,
    out_strahler_file,
    out_watershed_file,
    out_dist2stream_file,
    outdir,
    target_res_m = 10,
    resample_method = "bilinear",
    minslope = 0.1,
    saga_obj,
    threshold = 1000,
    .verbose = TRUE
) {
  
  stopifnot(!is.null(saga_obj))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # -------------------------------------------------------------------------
  # 1) DEM preparation (written as SAGA grid)
  # -------------------------------------------------------------------------
  dem_10m_file <- file.path(outdir, "dem10m.sdat")
  terra::writeRaster(dem_1m, dem_10m_file, overwrite = TRUE)
  
  # -------------------------------------------------------------------------
  # 2) Fill sinks
  # -------------------------------------------------------------------------
  dem_filled_file <- file.path(outdir, "dem10m_filled.sgrd")
  
  saga_obj$ta_preprocessor$fill_sinks_xxl_wang_liu(
    elev     = dem_10m_file,
    filled   = dem_filled_file,
    minslope = minslope,
    .verbose = .verbose
  )
  stopifnot(file.exists(dem_filled_file))
  
  # -------------------------------------------------------------------------
  # 3) Flow accumulation
  # -------------------------------------------------------------------------
  flowacc_sgrd <- file.path(outdir, "flowacc_10m.sgrd")
  
  saga_obj$ta_hydrology$flow_accumulation_top_down(
    elevation  = dem_filled_file,
    accu_total = flowacc_sgrd,
    .verbose   = .verbose
  )
  stopifnot(file.exists(flowacc_sgrd))
  
  terra::writeRaster(
    terra::rast(sub("\\.sgrd$", ".sdat", flowacc_sgrd)),
    out_flowacc_file,
    overwrite = TRUE
  )
  
  # -------------------------------------------------------------------------
  # 4) Channel network and drainage basins
  # -------------------------------------------------------------------------
  channels_gpkg <- file.path(outdir, "channels.gpkg")
  order_sgrd    <- file.path(outdir, "strahler_order_10m.sgrd")
  basin_sgrd   <- file.path(outdir, "watershed_id_10m.sgrd")
  basins_gpkg   <- file.path(outdir, "watershed_id_10m.gpkg")

  saga_obj$ta_channels$channel_network_and_drainage_basins(
    dem       = dem_filled_file,
    threshold = threshold,
    segments  = channels_gpkg,
    order     = order_sgrd,
    basin     = basin_sgrd,
    basins    = basins_gpkg


  )
  
  terra::writeRaster(terra::rast(sub("\\.sgrd$", ".sdat", order_sgrd)),
                     out_strahler_file, overwrite = TRUE)
  terra::writeRaster(terra::rast(sub("\\.sgrd$", ".sdat", basin_sgrd)),
                     out_watershed_file, overwrite = TRUE)
  
  
  
  
  # -------------------------------------------------------------------------
  # 5) Overland flow distance to channel network (explicit outputs)
  # -------------------------------------------------------------------------
  # ------------------------------------------------------------
  # Rasterize channel network (vector → grid for SAGA input)
  # ------------------------------------------------------------
  
  channels_gpkg  <- file.path(outdir, "channels.gpkg")          # existing vector output
  channels_sdat <- file.path(outdir, "channels_10m.sdat")     # terra output
  channels_sgrd <- file.path(outdir, "channels_10m.sgrd")     # SAGA header
  

  
  # reference grid: use filled DEM grid geometry
  dem_ref <- terra::rast(sub("\\.sgrd$", ".sdat", dem_filled_file))
  
  channels_vect <- terra::vect(channels_gpkg)
  
  # rasterize: presence / binary mask (1 = channel, NA = no channel)
  channels_rast <- terra::rasterize(
    x     = channels_vect,
    y     = dem_ref,
    field = 1,
    background = NA_real_
  )
  
  # write as SAGA-compatible raster (terra writes .sdat)
  terra::writeRaster(channels_rast, channels_sdat, overwrite = TRUE)
  

  # SAGA expects the .sgrd header; Rsagacmd will create it implicitly
  # when the file is used as input, but some setups prefer it explicitly:
  if (!file.exists(channels_sgrd)) {
    file.copy(channels_sdat, channels_sgrd)
  }
  
  
  dist_files <- list(
    distance = file.path(outdir, "dist_to_stream_10m.sgrd"),
    distvert = file.path(outdir, "distvert_to_stream_10m.sgrd"),
    cahnnels = file.path(outdir, "channels_10m.sgrd"),
    disthorz = file.path(outdir, "disthorz_to_stream_10m.sgrd"),
    time     = file.path(outdir, "flow_time_10m.sgrd"),
    sdr      = file.path(outdir, "sdr_10m.sgrd"),
    passes   = file.path(outdir, "passes_10m.sgrd")
  )
  
  saga_obj$ta_channels$overland_flow_distance_to_channel_network(
    elevation = dem_filled_file,
    channels  = channels_sgrd,
    distance  = dist_files$distance,
    distvert  = dist_files$distvert,
    disthorz  = dist_files$disthorz,
    time      = dist_files$time,
    sdr       = dist_files$sdr,
    passes    = dist_files$passes,
    method           = 1,
    boundary         = FALSE,
    flow_k_default   = 20,
    flow_r_default   = 0.05,
    .verbose = TRUE
  )
  
  stopifnot(file.exists(dist_files$distance))
  
  terra::writeRaster(
    terra::rast(sub("\\.sgrd$", ".sdat", dist_files$distance)),
    out_dist2stream_file,
    overwrite = TRUE
  )
  
  # -------------------------------------------------------------------------
  # 6) Assemble hydro stack
  # -------------------------------------------------------------------------
  hydro_files <- c(
    flowacc   = flowacc_sgrd,
    strahler  = order_sgrd,
    watershed = basin_sgrd,
    dist      = dist_files$distance
  )
  
  hydro_stack <- terra::rast(sub("\\.sgrd$", ".sdat", hydro_files))
  names(hydro_stack) <- names(hydro_files)
  
  hydro_stack_file <- file.path(outdir, "hydro_stack_10m.tif")
  terra::writeRaster(hydro_stack, hydro_stack_file, overwrite = TRUE)
  
  invisible(list(
    dem_10m      = dem_10m_file,
    dem_filled   = dem_filled_file,
    flowacc_sgrd = flowacc_sgrd,
    channels     = channels_gpkg,
    order        = order_sgrd,
    basin        = basin_sgrd,
    basins       = basins_gpkg,
    dist2stream  = dist_files$distance,
    hydro_stack  = hydro_stack_file
  ))
}

# ===========================================================================
# EXECUTION
# ===========================================================================

dem_1m_file <- paths[["aoi_dgm"]]
stopifnot(file.exists(dem_1m_file))
dem_1m <- terra::rast(dem_1m_file)

# --- Relief stack -----------------------------------------------------------
out_10m <- paths[["relief_stack_10m"]]

saga_derive_relief_stack(
  dem_1m          = dem_1m,
  prefix          = "relief_10m_burgwald",
  out_stack_file  = out_10m,
  outdir          = dirname(out_10m),
  target_res_m    = 10,
  resample_method = "bilinear",
  tpi_scales_m    = c(30, 50, 250),
  minslope        = 0.1,
  saga_obj        = saga
)

message("Wrote relief stack: ", out_10m)

# --- Hydrology --------------------------------------------------------------
out_flowacc   <- paths[["flowacc_10m"]]
out_strahler  <- paths[["strahler_order_10m"]]
out_watershed <- paths[["watershed_id_10m"]]
out_dist2str  <- paths[["dist_to_stream_10m"]]

saga_derive_hydro_rasters(
  dem_1m               = dem_1m,
  out_flowacc_file     = out_flowacc,
  out_strahler_file    = out_strahler,
  out_watershed_file   = out_watershed,
  out_dist2stream_file = out_dist2str,
  outdir               = dirname(out_flowacc),
  target_res_m         = 10,
  resample_method      = "bilinear",
  minslope             = 0.1,
  saga_obj             = saga
)

message("Wrote hydrology: ", out_flowacc)
message("Wrote hydrology: ", out_strahler)
message("Wrote hydrology: ", out_watershed)
message("Wrote hydrology: ", out_dist2str)


# out_gpkg <- paths[["hydrology_vectors_10m"]]
# stopifnot(!is.null(out_gpkg))

# # channels (line network)
# stopifnot(file.exists(channels_gpkg))
# ch <- sf::st_read(channels_gpkg, quiet = TRUE)
# sf::st_write(ch, out_gpkg, layer = "channels", delete_layer = TRUE, quiet = TRUE)
# 
# # basins (polygons)
# stopifnot(file.exists(basins_gpkg))
# ba <- sf::st_read(basins_gpkg, quiet = TRUE)
# sf::st_write(ba, out_gpkg, layer = "basins", delete_layer = TRUE, quiet = TRUE)
# 
# # optional nodes
# if (exists("nodes_shp", inherits = FALSE) && file.exists(nodes_shp)) {
#   nd <- sf::st_read(nodes_gpkg, quiet = TRUE)
#   sf::st_write(nd, out_gpkg, layer = "nodes", delete_layer = TRUE, quiet = TRUE)
# }
# 
# stopifnot(file.exists(out_gpkg))