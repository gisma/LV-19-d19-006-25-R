############################################################
# 00_setup-burgwald.R
#
# Tasks:
# 1) Load all packages
# 2) Create project folders
# 3) Define AOI for Burgwald
# 4) Load helper functions
#
# IMPORTANT:
# - Must be run inside an RStudio Project
# - Paths are handled via here::here()
# - No setwd() anywhere
############################################################


# -----------------------------
# 0) Packages & project root


if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)

pacman::p_load(
  
  # --- project / filesystem / dates
  here,
  fs,
  lubridate,
  
  # --- spatial core (order matters)
  sf,
  terra,
  
  # --- cubes / STAC / remote sensing backends
  gdalcubes,
  stars,
  tidyterra,
  rstac,
  CDSE,
  
  # --- RS utilities
  exactextractr,
  RStoolbox,
  
  # --- tidy data / plotting
  dplyr,
  tidyr,
  ggplot2,
  
  # --- visualisation / maps
  mapview,
  mapedit,
  tmap,
  tmaptools,
  colorspace,
  
  # --- web / OGC / APIs
  httr,
  OpenStreetMap,
  ows4R,
  jsonlite,
  osmdata, 
  rvest, 
  data.table,
  
  # --- ML / classification (infrastructure only)
  randomForest,
  ranger,
  e1071,
  caret,
  Rsagacmd,
  bfast,
  sp
  
)

remotes::install_github("r-spatial/link2GI",ref = "master")
library(link2GI)

# # raster bewusst NICHT attachen
# if (!requireNamespace("raster", quietly = TRUE)) {
#   install.packages("raster")
# }


message("Project root: ", here::here())
options(timeout = 600)

# -------------------------------------------------------
# 1) Create mandatory folders
# -------------------------------------------------------
aoi_root <- here::here("data", "raw", "AOI_Burgwald")
fs::dir_create(aoi_root)
fs::dir_create(file.path(aoi_root, "dem"))
fs::dir_create(file.path(aoi_root, "clc"))
fs::dir_create(file.path(aoi_root, "osm_by_key"))
# -----------------------------
fs::dir_create(here("data", "raw"))
fs::dir_create(here("data", "processed"))
fs::dir_create(here("outputs", "figures"))
fs::dir_create(here("metadata"))
fs::dir_create(here("docs"))
fs::dir_create(here("src"))
fs::dir_create(here("data", "raw", "dwd-stations"))
fs::dir_create(here("data", "processed", "dwd-stations"))
# -----------------------------
# Global paths for DWD pipeline
# -----------------------------
path_dwd_raw      <- here::here("data", "raw", "dwd-stations")
path_dwd_processed <- here::here("data", "processed", "dwd-stations")
root_folder <- here::here()



# -----------------------------
# 2) Define Burgwald AOI (WGS84)
# -----------------------------

source(here::here("_core", "01-fun-data-retrieval.R"))
burgwald_bbox <- c(
  xmin = 8.70,
  xmax = 9.00,
  ymin = 50.85,
  ymax = 51.05
)

aoi_burgwald_50km <- aoi_with_buffer(burgwald_bbox, buffer_km = 50)

aoi_burgwald_wgs <- sf::st_as_sfc(
  sf::st_bbox(burgwald_bbox, crs = 4326)
)

# Save AOI for inspection / reuse
aoi_file <- here("data", "processed", "aoi_burgwald.gpkg")
sf::st_write(
  sf::st_sf(geometry = aoi_burgwald_wgs),
  aoi_file,
  delete_dsn = TRUE
)

cat("AOI saved to:", aoi_file, "\n")



clc_legend <- data.frame(
  class_id = 1:44,  # WICHTIG: 1–44 = Rasterwerte
  code = c(
    111,112,121,122,123,124,131,132,133,141,142,
    211,212,213,221,222,223,231,241,242,243,244,
    311,312,313,321,322,323,324,331,332,333,334,335,
    411,412,421,422,423,
    511,512,521,522,523
  ),
  name = c(
    "Continuous urban fabric",
    "Discontinuous urban fabric",
    "Industrial or commercial units",
    "Road and rail networks and associated land",
    "Port areas",
    "Airports",
    "Mineral extraction sites",
    "Dump sites",
    "Construction sites",
    "Green urban areas",
    "Sport and leisure facilities",
    
    "Non-irrigated arable land",
    "Permanently irrigated land",
    "Rice fields",
    "Vineyards",
    "Fruit trees and berry plantations",
    "Olive groves",
    "Pastures",
    "Annual crops associated with permanent crops",
    "Complex cultivation patterns",
    "Agriculture with significant areas of natural vegetation",
    "Agro-forestry areas",
    
    "Broad-leaved forest",
    "Coniferous forest",
    "Mixed forest",
    "Natural grasslands",
    "Moors and heathland",
    "Sclerophyllous vegetation",
    "Transitional woodland-shrub",
    "Beaches, dunes, sands",
    "Bare rocks",
    "Sparsely vegetated areas",
    "Burnt areas",
    "Glaciers and perpetual snow",
    
    "Inland marshes",
    "Peatbogs",
    "Salt marshes",
    "Salines",
    "Intertidal flats",
    
    "Water courses",
    "Water bodies",
    "Coastal lagoons",
    "Estuaries",
    "Sea and ocean"
  ),
  color = c(
    "#E6004D","#FF0000","#CC4DF2","#CC0000","#E68000","#E6E600",
    "#A6A6A6","#A6A6A6","#FFA6FF","#A6FFA6","#A6E6A6",
    "#FFFF00","#E6E600","#E6CC00","#E6CC7F","#E6B300","#E6A600",
    "#E68000","#FFFF7F","#FFFF00","#E6E64D","#E6E64D",
    "#00A600","#008000","#70A64D","#A6FF00","#A6A600","#A6A600",
    "#A6CC00","#E6E6E6","#DADADA","#BFBFBF","#F2F2F2","#FFFFFF",
    "#00CCF2","#00A6F2","#00CCE6","#66CCE6","#A6E6F2",
    "#0000CC","#3333FF","#6680FF","#00A6CC","#99CCE6"
  ),
  stringsAsFactors = FALSE
)

# Color palettes for vegetation indices:
ndvi.col <- function(n) rev(colorspace::sequential_hcl(n, "Green-Yellow"))
ano.col  <- colorspace::diverging_hcl(7, palette = "Red-Green", register = "rg")

# ---- 0) Project-local temp dirs (portable, user-space) -----------------

# Base temp dir inside the project (user-space, under version control only as .gitignore)
proj_tmp_root <- here::here("tmp")

# Allow override via environment variable (e.g. on cluster)
custom_tmp <- Sys.getenv("BURGWALD_TMP_DIR", unset = NA_character_)
if (!is.na(custom_tmp) && nzchar(custom_tmp)) {
  proj_tmp_root <- custom_tmp
}

# Subdirs for R / GDAL and terra
r_tmp_dir     <- file.path(proj_tmp_root, "Rtmp")
terra_tmp_dir <- file.path(proj_tmp_root, "terra")

# Create dirs if they don't exist (user-space only)
dir.create(r_tmp_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(terra_tmp_dir, showWarnings = FALSE, recursive = TRUE)

# Point R / GDAL temp to project-local dir
Sys.setenv(
  TMPDIR = r_tmp_dir,
  TEMP   = r_tmp_dir,
  TMP    = r_tmp_dir
)

# Tell terra to use project-local temp dir
terra::terraOptions(
  tempdir = terra_tmp_dir,
  memfrac = 0.2,   # fraction of RAM terra may use
  todisk  = TRUE   # write big objects to disk instead of RAM
)

# Optional: clean old terra tmp-files from previous runs
# Clean ALL old terra tmp-files (recommended)
#terra::tmpFiles(orphan = TRUE, current = FALSE, old = FALSE, remove = TRUE)

