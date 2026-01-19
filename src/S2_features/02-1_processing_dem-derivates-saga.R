#!/usr/bin/env Rscript

# src/02-1_processing_dem-derivates-saga.R
#
# Terrain derivatives via SAGA (Rsagacmd), written as ONE multi-band GeoTIFF stack.
# Base input is the 1 m DEM; a target resolution can be selected via resampling.
#
# Output logic:
# - Canonical outputs are resolved via `paths[[<key>]]` (from metadata/outputs.tsv).
# - This script writes ONLY the final stack(s). No “work_*” artefacts.

source("src/_core/01-setup-burgwald.R")

3# --- SAGA wrapper (one global instance) --------------------------------------
saga <- Rsagacmd::saga_gis(
  raster_backend = "terra",
  vector_backend = "sf",
  all_outputs    = FALSE
)

# ensure_parent_dir <- function(path) {
#   dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
#   invisible(TRUE)
# }

#' Derive relief layers from a DEM (SAGA) and write them as a single GeoTIFF stack
#'
#' What this does (high level):
#' 1) Optionally resample the base DEM to a target resolution (e.g., 10 m).
#' 2) Fill sinks (Wang & Liu).
#' 3) Compute slope / aspect / curvatures (SAGA "Slope, Aspect, Curvature").
#' 4) Compute multi-scale TPI for given radii in meters.
#' 5) Write ONE multi-band GeoTIFF and (optionally) 
#'
#' Notes on “SAGA temp files”:
#' Rsagacmd usually manages SAGA I/O through temp grids, but depending on backend
#' and module behaviour you can still end up with .sdat/.sgrd sidecars in outdir.
#' If that happens, `cleanup_saga_sidecars = TRUE` removes them by prefix match.
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
  
  # Use provided SAGA instance or create one here
  if (is.null(saga_obj)) {
    saga_obj <- Rs.player <- Rsagacmd::saga_gis(
      raster_backend = "terra",
      vector_backend = "sf",
      all_outputs    = TRUE
    )
  }
  #ensure_parent_dir(out_stack_file)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # --- 0) Optional resampling from the 1 m base DEM --------------------------
  dem <- dem_1m
  
  if (!is.null(target_res_m)) {
    # Build a target grid that shares extent+CRS but has the requested resolution.
    # We use a template raster and then resample.
    tpl <- terra::rast(
      extent = terra::ext(dem_1m),
      crs    = terra::crs(dem_1m),
      resolution = target_res_m
    )
    dem <- terra::resample(dem_1m, tpl, method = resample_method)
  }
  
  # --- 1) Fill sinks (Wang & Liu) -------------------------------------------
  dem_filled <- saga_obj$ta_preprocessor$fill_sinks_xxl_wang_liu(
    elev     = dem,
    minslope = minslope
  )
  
  # --- 2) Slope / aspect / curvatures --------------------------------------
  lsps <- saga_obj$ta_morphometry$slope_aspect_curvature(
    elevation   = dem_filled,
    method      = 6,  # SAGA method id as in your current script
    unit_slope  = 0,  # radians
    unit_aspect = 0   # radians
  )
  lsps_stack <- terra::rast(lsps)
  
  # --- 3) Multi-scale TPI (radii in meters) ---------------------------------
  tpi_rasters <- vector("list", length(tpi_scales_m))
  for (i in seq_along(tpi_scales_m)) {
    Rm <- tpi_scales_m[i]
    tpi_rasters[[i]] <- saga_obj$ta_morphometry$topographic_position_index_tpi(
      dem        = dem_filled,
      radius_min = 0,
      radius_max = Rm,
      standard   = TRUE
    )
  }
  tpi_stack <- terra::rast(tpi_rasters)
  
  # --- 4) Assemble final stack (one file) -----------------------------------
  out_stack <- c(dem_filled, lsps_stack, tpi_stack)
  
  # Layer naming: stable + explicit
  n_lsps <- terra::nlyr(lsps_stack)
  lsps_names <- names(lsps_stack)
  if (is.null(lsps_names) || any(!nzchar(lsps_names))) {
    lsps_names <- paste0("lsps_", seq_len(n_lsps))
  }
  
  names(out_stack) <- c(
    "dem_filled",
    lsps_names,
    paste0("tpi_r", tpi_scales_m, "m")
  )
  
  terra::writeRaster(out_stack, out_stack_file, overwrite = TRUE)
  

  invisible(out_stack_file)
}

# ---------------------------------------------------------------------------
# Base DEM (1 m) comes from canonical S1 output: key = aoi_dgm
# ---------------------------------------------------------------------------
dem_1m_file <- paths[["aoi_dgm"]]
stopifnot(file.exists(dem_1m_file))
dem_1m <- terra::rast(dem_1m_file)

# ---------------------------------------------------------------------------
# YOUR CHOSEN SETUP:
# - Base model = 1 m DEM
# - Target resolution for derivatives = 10 m
# - TPI radii (meters): micro=30, meso=50, macro=250
# - Output = ONE stack GeoTIFF (canonical key from outputs.tsv)
# ---------------------------------------------------------------------------

out_10m <- paths[["relief_stack_10m"]]

saga_derive_relief_stack(
  dem_1m          = dem_1m,
  prefix          = "relief_10m_burgwald",
  out_stack_file  = out_10m,
  outdir          = dirname(out_10m),
  target_res_m    = 10,
  resample_method = "bilinear",
  tpi_scales_m    = c(30, 50, 250),
  minslope        = 0.1
)

message("Wrote relief stack: ", out_10m)
