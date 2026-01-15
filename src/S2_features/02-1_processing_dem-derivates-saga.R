#!/usr/bin/env Rscript

############################################################
# Script:  02-1_processing_DEM-derivates-SAGA.R
# Author:   [Your Name]
# Project:  [Your Project Name]
#
# Purpose:
# --------
## - Derive a set of SAGA-GIS based terrain / relief parameters
##   (filled DEM, slope/aspect/curvatures, multi-scale TPI)
##   from one or more DEMs using Rsagacmd.
##
## Workflow:
##  1) Initialise a global SAGA-GIS wrapper object (`saga`).
##  2) Define a function `saga_derive_relief()` that:
##       - fills sinks in the DEM,
##       - computes slope/aspect/curvatures,
##       - computes TPI for multiple radii,
##       - writes all outputs to disk and returns them invisibly.
##  3) Call this function for:
##       - an (intended) 1 m DEM,
##       - an aggregated 100 m DEM.
##
## Wrapper concept (API calls in mature GIS software):
## - `Rsagacmd::saga_gis()` creates an R **wrapper object** (`saga`)
##   that exposes SAGA tools as R functions:
##       `saga$ta_morphometry$slope_aspect_curvature(...)`
## - Internally, these wrappers:
##     * build the correct command-line / API call,
##     * manage parameters and temporary files,
##     * convert results back to R objects (here: `terra::SpatRaster`).
## - This pattern is similar to:
##     * `qgisprocess::qgis_run_algorithm()` (QGIS),
##     * `rgrass7::execGRASS()` (GRASS GIS),
##     * ArcPy (Python wrapper around ArcGIS).
## - Advantage: you keep full GIS functionality, but you orchestrate
##   everything reproducibly from R.
##
## NOTE:
## - The **executable code below is unchanged**; only comments and this
##   header were added for documentation.
## ============================================================

# -------------------------------------------------------
# SAGA-GIS Relief-Derivate mit Rsagacmd (sagacmd)
# -------------------------------------------------------

source(here::here("src", "00-setup-burgwald.R"))


# einmal global initialisieren
# -> this creates the SAGA wrapper object that provides access
#    to SAGA modules via R (see wrapper concept above).
saga <- saga_gis(
  raster_backend = "terra",
  vector_backend = "sf",
  all_outputs    = TRUE
)

saga_derive_relief <- function(dem,
                               prefix,
                               outdir     = here::here("data", "processed"),
                               minslope   = 0.1,
                               tpi_scales = c(100, 500, 1000),
                               saga       = saga) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  ## 1) Sinks füllen (Wang & Liu) – Rsagacmd gibt direkt das gefüllte DEM zurück
  dem_filled <- saga$ta_preprocessor$fill_sinks_xxl_wang_liu(
    elev     = dem,
    minslope = minslope
  )
  
  terra::writeRaster(
    dem_filled,
    file.path(outdir, paste0(prefix, "_filled_dem.tif")),
    overwrite = TRUE
  )
  
  ## 2) Standard-Ableitungen (Slope, Aspect, Curvatures)
  lsps <- saga$ta_morphometry$slope_aspect_curvature(
    elevation   = dem_filled,  # nicht dem = ...
    method      = 6,           # z.B. "maximum slope" (siehst du in available_opts)
    unit_slope  = 0,           # 0 = radians
    unit_aspect = 0            # 0 = radians
  )
  
  # lsps ist eine LISTE von SpatRaster-Objekten -> Stack bauen:
  lsps_stack <- terra::rast(lsps)
  
  terra::writeRaster(
    lsps_stack,
    here::here("data", "processed", paste0(prefix, "_slope_aspect_curvature.tif")),
    overwrite = TRUE
  )
  
  ## 3) Multi-scale TPI
  tpi_list <- list()
  for (R in tpi_scales) {
    tpi_r <- saga$ta_morphometry$topographic_position_index_tpi(
      dem        = dem_filled,
      radius_min = 0,
      radius_max = R,
      standard   = TRUE   # standardisierte TPI-Werte
    )
    out_tpi <- file.path(outdir, paste0(prefix, "_tpi_r", R, "m.tif"))
    terra::writeRaster(tpi_r, out_tpi, overwrite = TRUE)
    tpi_list[[paste0("tpi_r", R)]] <- tpi_r
  }
  
  invisible(list(
    dem_filled = dem_filled,
    lsps       = lsps,
    tpi        = tpi_list
  ))
}

# -------------------------------------------------------
# 1) Originalauflösung (DGM1 Burgwald, ~1 m)
# -------------------------------------------------------
dem_1m_file <- here::here("data", "raw", "AOI_Burgwald", "dem", "dem_dgm1_burgwald.tif")
stopifnot(file.exists(dem_1m_file))

dem_1m <- terra::rast(dem_1m_file)



# Originalauflösung
relief_1m <- saga_derive_relief(
  dem    = dem_1m,
  prefix = "dem1m_burgwald",
  outdir = here::here("data", "processed", "relief_1m")
)


# 100 m DEM in R aggregieren (Beispiel)
dem_100m <- terra::aggregate(dem_1m, fact = 100)



relief_100m <- saga_derive_relief(
  dem    = dem_100m,
  prefix = "dem100m_burgwald",
  outdir = here::here("data", "processed", "relief_100m")
)
