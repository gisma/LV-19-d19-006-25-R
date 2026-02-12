############################################################
# 01-setup-burgwald.R
#
# Tasks:
# 1) Load packages
# 2) Create folder skeleton (single source: src/_core/00-folders.R)
# 3) Define AOI (WGS84)
# 4) Build canonical output paths from core/outputs.tsv
############################################################

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
# --- packages (Windows + Linux, fehlertolerant) ----------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load_tolerant <- function(pkgs, install = TRUE, repos = getOption("repos")) {
  ok <- character(0)
  fail <- list()
  
  for (p in pkgs) {
    res <- try({
      if (!requireNamespace(p, quietly = TRUE)) {
        if (!install) stop("not installed")
        install.packages(p, repos = repos, dependencies = TRUE)
      }
      suppressPackageStartupMessages(library(p, character.only = TRUE))
      TRUE
    }, silent = TRUE)
    
    if (isTRUE(res)) {
      ok <- c(ok, p)
    } else {
      fail[[p]] <- as.character(res)
      message("Package FAILED: ", p)
    }
  }
  
  message("\n--- Package load summary ---")
  message("OK:   ", paste(ok, collapse = ", "))
  message("FAIL: ", if (length(fail)) paste(names(fail), collapse = ", ") else "<none>")
  
  invisible(list(ok = ok, fail = fail))
}

pkgs <- c(
  "here","fs","lubridate",
  "sf","terra",
  "gdalcubes","stars","tidyterra","rstac","CDSE",
  "exactextractr","RStoolbox",
  "dplyr","tidyr","ggplot2","curl",
  "mapview","mapedit","tmap","tmaptools","colorspace",
  "httr","ows4R","jsonlite","osmdata","rvest","data.table",
  "randomForest","ranger","e1071","caret","Rsagacmd","bfast","sp","link2GI"
)

pkg_status <- p_load_tolerant(pkgs, install = TRUE)

# --- optional: OpenStreetMap (nicht kritisch) ------------------------------
# Unter Windows oft trouble; daher separat + tolerant.
if (!requireNamespace("OpenStreetMap", quietly = TRUE)) {
  try(install.packages("OpenStreetMap", dependencies = TRUE), silent = TRUE)
}
if (requireNamespace("OpenStreetMap", quietly = TRUE)) {
  suppressPackageStartupMessages(library(OpenStreetMap))
} else {
  message("OpenStreetMap konnte nicht installiert/geladen werden (optional).")
}

# --- access failed list ----------------------------------------------------
# names(pkg_status$fail)
# pkg_status$fail


message("Project root: ", here::here())
options(timeout = 600)
Sys.setenv(OSM_TIMEOUT = 600)


# --- folder skeleton (single source) --------------------------------------
source(file.path(root_folder, "src", "_core","00-folders.R"))
source(file.path(root_folder, "src", "_core","02-helper.R"))

root_folder <- here::here()
dirs <- createFolders_simple(root_folder = root_folder, folders = burgwald_folders(), create_folders = TRUE)

# --- canonical outputs (S,key -> absolute path) ---------------------------
outputs_df <- burgwald_outputs()              # metadata/outputs.tsv
paths <- compile_outputs(root_folder, outputs_df)
for (p in unname(paths)) {
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
}
# --- AOI (WGS84) ----------------------------------------------------------
source(here::here("src", "r-libs", "01-fun-data-retrieval.R"))

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

# --- raw AOI roots --------------------------------------------------------
aoi_root <- here::here("data", "raw", "AOI_Burgwald")

# --- DWD paths (raw/provider + processed/productive handled elsewhere) ----
path_dwd_raw <- here::here("data", "raw", "dwd-stations")
path_dwd_processed <- here::here("data", "processed", "dwd-stations")

# --- palettes / legend (unchanged from your setup) ------------------------
clc_legend <- data.frame(
  class_id = 1:44,
  code = c(
    111,112,121,122,123,124,131,132,133,141,142,
    211,212,213,221,222,223,231,241,242,243,244,
    311,312,313,321,322,323,324,331,332,333,334,335,
    411,412,421,422,423,
    511,512,521,522,523
  ),
  name = c(
    "Continuous urban fabric","Discontinuous urban fabric","Industrial or commercial units",
    "Road and rail networks and associated land","Port areas","Airports",
    "Mineral extraction sites","Dump sites","Construction sites",
    "Green urban areas","Sport and leisure facilities",
    "Non-irrigated arable land","Permanently irrigated land","Rice fields",
    "Vineyards","Fruit trees and berry plantations","Olive groves",
    "Pastures","Annual crops associated with permanent crops","Complex cultivation patterns",
    "Agriculture with significant areas of natural vegetation","Agro-forestry areas",
    "Broad-leaved forest","Coniferous forest","Mixed forest",
    "Natural grasslands","Moors and heathland","Sclerophyllous vegetation",
    "Transitional woodland-shrub","Beaches, dunes, sands","Bare rocks",
    "Sparsely vegetated areas","Burnt areas","Glaciers and perpetual snow",
    "Inland marshes","Peatbogs","Salt marshes","Salines","Intertidal flats",
    "Water courses","Water bodies","Coastal lagoons","Estuaries","Sea and ocean"
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

ndvi.col <- function(n) rev(colorspace::sequential_hcl(n, "Green-Yellow"))
ano.col  <- colorspace::diverging_hcl(7, palette = "Red-Green", register = "rg")

# --- temp dirs (unchanged logic) ------------------------------------------
proj_tmp_root <- here::here("tmp")
custom_tmp <- Sys.getenv("BURGWALD_TMP_DIR", unset = NA_character_)
if (!is.na(custom_tmp) && nzchar(custom_tmp)) proj_tmp_root <- custom_tmp

r_tmp_dir     <- file.path(proj_tmp_root, "Rtmp")
terra_tmp_dir <- file.path(proj_tmp_root, "terra")

dir.create(r_tmp_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(terra_tmp_dir, showWarnings = FALSE, recursive = TRUE)

Sys.setenv(TMPDIR = r_tmp_dir, TEMP = r_tmp_dir, TMP = r_tmp_dir)

terra::terraOptions(tempdir = terra_tmp_dir, memfrac = 0.2, todisk = TRUE)


# --- Climate interpolation paths (minimal envrmt compatibility) ------------
# These keys are required by fun_clim_* helpers (historically envimaR::createEnvi()).
# We keep this strictly as a path registry; no administrative extraction here.

envrmt <- list(
  # raw roots
  path_data_lev0 = here::here("data", "raw"),
  # derived products (interpolation outputs)
  path_data_lev1 = here::here("data", "productive", "S2_features"),
  # optional / not used in the simplified RSK pipeline but kept for compatibility
  path_data_lev2 = here::here("data", "productive", "S3_structure"),

  # providers
  path_CDC_KL    = here::here("data", "raw", "providers", "dwd"),
  path_GhcnDaily = here::here("data", "raw", "providers", "ghcn_daily"),
  path_GhcnMonthly = here::here("data", "raw", "providers", "ghcn_monthly")
)

# ensure provider dirs exist (no-ops if already there)
dir.create(envrmt$path_CDC_KL, recursive = TRUE, showWarnings = FALSE)
dir.create(envrmt$path_GhcnDaily, recursive = TRUE, showWarnings = FALSE)
dir.create(envrmt$path_GhcnMonthly, recursive = TRUE, showWarnings = FALSE)

