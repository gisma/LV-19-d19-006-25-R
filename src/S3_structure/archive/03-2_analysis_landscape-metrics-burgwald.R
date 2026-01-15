#!/usr/bin/env Rscript
############################################################
# 03-2_analysis_landscape-metrics-burgwald.R
#
# Purpose
# -------
# Apply information-theoretic and classical landscape metrics
# plus pattern-based spatial analysis (motif) to the
# Burgwald land-cover raster created in 01-1_get_base-geodata.R
#
# Conceptual background
# ---------------------
# 1) Information-theoretic (IT) landscape metrics
#    A categorical raster (land cover) is treated as a spatial signal.
#    From its co-occurrence matrix (how often class i is adjacent to j),
#    we derive:
#
#    - Entropy H(x): compositional diversity
#        * Measures uncertainty in the distribution of classes.
#        * Low H(x): one dominant class (monothematic).
#        * High H(x): many classes in similar proportions (multithematic).
#
#    - Relative Mutual Information U: configurational order / clumpiness
#        * Measures how strongly the spatial arrangement deviates from
#          random mixing.
#        * Low U: classes are intermixed in a near-random way
#                 → fragmented, noisy neighbourhoods.
#        * High U: classes form ordered blocks
#                 → aggregated, clumped landscape structure.
#
#    Together, H(x) (composition) and U (configuration) form a
#    2-dimensional pattern space:
#        - H low, U low   → few classes, randomly mixed
#        - H low, U high  → one/few classes, strongly clumped
#        - H high, U low  → many classes, highly fragmented mosaic
#        - H high, U high → many classes, but spatially ordered
#
#    This 2D description reduces redundancy compared to many classical
#    metrics (e.g. SHDI, AI, PD, ED) that are often strongly correlated.
#
# 2) Pattern-based spatial analysis with motif
#    Instead of single scalar metrics, motif builds spatial signatures:
#    - For each local window (e.g. 25×25 raster cells), it computes a
#      co-occurrence-based signature (COVE).
#    - Each signature is a vector describing adjacency probabilities
#      of land-cover classes inside that window.
#    - Distances between signatures (e.g. Jensen–Shannon) indicate how
#      similar or different local patterns are.
#    - Clustering these signatures yields landscape pattern types
#      (typologies) across the study area.
#
# Inputs
# ------
#   - AOI: data/processed/aoi_burgwald.gpkg
#   - CLC: data/raw/AOI_Burgwald/clc/clc5_2018_burgwald.tif
#
# Outputs (objects in memory)
# ---------------------------
#   it_metrics_global    : global H(x) and U for the entire Burgwald AOI
#   grid_all_metrics_sf  : 2×2 km grid with ent, relmutinf, shdi, ai
#   cor_metrics          : correlation matrix of ent, relmutinf, shdi, ai
#   lsp_cove_burgwald    : motif COVE signatures for ~2.5 km windows
#   lsp_clusters_sf      : sf polygons with "cluster" (pattern types)
############################################################

# ---------------------------------------------------------
# 0) Packages
# ---------------------------------------------------------
library(here)
library(terra)
library(sf)
library(landscapemetrics)
library(motif)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

message("Project root: ", here::here())

# ---------------------------------------------------------
# 1) Paths & input data
# ---------------------------------------------------------

# AOI root as defined in the setup script (01_setup-burgwald.R)
aoi_root <- here::here("data", "raw", "AOI_Burgwald")

# Land-cover raster (CLC5 2018, already clipped to a rough Burgwald extent)
# This is a categorical raster with CLC classes at ~100 m resolution.
clc_file <- file.path(aoi_root, "clc", "clc5_2018_burgwald.tif")

# AOI polygon (WGS84) from the setup script
aoi_file <- here::here("data", "processed", "aoi_burgwald.gpkg")

# Read data
lc_burgwald  <- terra::rast(clc_file)
aoi_burgwald <- sf::read_sf(aoi_file)

# Harmonise CRS:
# Ensure raster and AOI polygon are in the same coordinate reference system.
if (!sf::st_crs(aoi_burgwald) == sf::st_crs(lc_burgwald)) {
  aoi_burgwald <- sf::st_transform(aoi_burgwald, sf::st_crs(lc_burgwald))
}

# Crop and mask the raster to the AOI.
# This guarantees that all metrics refer exactly to the AOI extent.
lc_burgwald <- lc_burgwald |>
  terra::crop(terra::vect(aoi_burgwald)) |>
  terra::mask(terra::vect(aoi_burgwald))

# Ensure that land-cover classes are integers (important for landscapemetrics).
terra::values(lc_burgwald) <- round(terra::values(lc_burgwald))

# ---------------------------------------------------------
# 2) Global information-theoretic metrics (whole Burgwald)
# ---------------------------------------------------------

# Compute global IT metrics for the entire Burgwald AOI:
# - lsm_l_ent       → entropy H(x), compositional diversity
# - lsm_l_relmutinf → relative mutual information U, configurational order
#
# H(x): how many classes exist and how evenly they occur.
# U   : how strongly the spatial arrangement deviates from random mixing.
it_ent_global <- lsm_l_ent(lc_burgwald)
it_u_global   <- lsm_l_relmutinf(lc_burgwald)

it_metrics_global <- tibble(
  metric = c("entropy_Hx", "relmutinf_U"),
  value  = c(it_ent_global$value, it_u_global$value)
)

# ---------------------------------------------------------
# 3) Grid-based metrics: IT + classical (SHDI, AI)
# ---------------------------------------------------------

# 3.1 Create a regular grid over the AOI.
#     Use a metric CRS (ETRS89 / UTM 32N) so that cellsize is in metres.
aoi_utm <- sf::st_transform(aoi_burgwald, 25832)

# Grid resolution (metres) – landscape unit size for metrics.
# Each grid cell will be treated as its own 'landscape'.
grid_cellsize <- 2000  # 2 km × 2 km

grid_utm <- sf::st_make_grid(
  aoi_utm,
  cellsize = grid_cellsize,
  what     = "polygons",
  square   = TRUE
) |>
  sf::st_as_sf() |>
  dplyr::mutate(grid_id = dplyr::row_number())

# Transform the grid back to the raster CRS so that sample_lsm() aligns
# with the land-cover raster.
grid <- sf::st_transform(grid_utm, sf::st_crs(lc_burgwald))

# 3.2 Information-theoretic metrics per grid cell (landscape level).
# Each 2×2 km grid cell is a small landscape for which we compute:
# - entropy H(x)          → lsm_l_ent
# - relative mutual info U → lsm_l_relmutinf
lsm_grid_it <- sample_lsm(
  landscape = lc_burgwald,
  y         = grid,
  what      = c("lsm_l_ent", "lsm_l_relmutinf"),
  level     = "landscape"
)

lsm_grid_it_wide <- lsm_grid_it |>
  dplyr::select(grid_id = plot_id, metric, value) |>
  tidyr::pivot_wider(names_from = metric, values_from = value)

# 3.3 Classical metrics per grid cell:
# - SHDI (Shannon Diversity Index)      → compositional diversity
# - AI   (Aggregation Index)            → configurational aggregation
lsm_grid_classic <- sample_lsm(
  landscape = lc_burgwald,
  y         = grid,
  what      = c("lsm_l_shdi", "lsm_l_ai"),
  level     = "landscape"
)

lsm_grid_classic_wide <- lsm_grid_classic |>
  dplyr::select(grid_id = plot_id, metric, value) |>
  tidyr::pivot_wider(names_from = metric, values_from = value)

# 3.4 Merge all metrics into one sf object.
# This yields a 2×2 km grid with H(x), U, SHDI, and AI per cell.
grid_all_metrics_sf <- grid |>
  dplyr::left_join(lsm_grid_it_wide,      by = "grid_id") |>
  dplyr::left_join(lsm_grid_classic_wide, by = "grid_id") 


# ---------------------------------------------------------
# 4) Correlation between metrics (redundancy check)
# ---------------------------------------------------------

# Prepare a numeric table of the four metrics:
# - ent       = H(x), IT composition
# - relmutinf = U, IT configuration
# - shdi      = classical composition
# - ai        = classical configuration
metrics_num <- grid_all_metrics_sf |>
  sf::st_drop_geometry() |>
  dplyr::select(ent, relmutinf, shdi, ai)

# Correlation matrix shows how strongly the metrics are related:
# - ent vs shdi → compositional redundancy
# - relmutinf vs ai → configurational redundancy
cor_metrics <- stats::cor(metrics_num, use = "pairwise.complete.obs")

# ---------------------------------------------------------
# 5) Pattern-based spatial analysis with motif
# ---------------------------------------------------------

# 5.1 Compute COVE (co-occurrence) signatures over the Burgwald raster.
#     window = 25 means:
#       - 25 × 25 raster cells, NOT metres.
#       - For CLC at ~100 m, this is ~2500 m × 2500 m (2.5 km × 2.5 km).
#     Each signature is a probability distribution (PDF) describing
#     adjacency patterns of land-cover classes within that window.
lsp_cove_burgwald <- lsp_signature(
  x            = lc_burgwald,
  type         = "cove",
  window       = 20,
  normalization = "pdf"
)

# 5.2 Compute distances between all signatures.
#     dist_fun = "jensen-shannon" computes Jensen–Shannon divergence
#     between the co-occurrence PDFs, i.e. a symmetric measure of
#     difference between local pattern structures.
dist_mat <- lsp_to_dist(
  x        = lsp_cove_burgwald,
  dist_fun = "jensen-shannon"
)

# 5.3 Cluster signatures into k pattern types.
#     Here, hierarchical clustering (Ward) is used on the distance matrix.
#     Each ~2.5 km window is assigned to one of k landscape pattern types.
set.seed(123)
hc <- hclust(dist_mat, method = "ward.D2")
cluster_ids <- cutree(hc, k = 6)  # number of pattern types

# Attach cluster IDs as an attribute to each signature.
lsp_cove_burgwald$cluster <- cluster_ids

# 5.4 Convert signatures to sf polygons.
#     lsp_add_sf() reconstructs spatial footprints (tiles/windows) and
#     returns an sf object with geometry and attributes, including "cluster".
lsp_clusters_sf <- lsp_add_sf(lsp_cove_burgwald)

# Example interactive visualization:
# Each polygon represents a ~2.5 km COVE window, colored by pattern type.
# Pattern types can be interpreted as:
#   - homogeneous forest blocks,
#   - agricultural mosaics,
#   - mixed-use transition zones, etc.
mapview::mapview(lsp_clusters_sf["cluster"])

# ---------------------------------------------------------
# 6) Object overview
# ---------------------------------------------------------

message("Global IT metrics (Burgwald):")
print(it_metrics_global)
# -> A single H(x) and U value describing composition and configuration
#    of the entire Burgwald AOI.

message("\nGrid-based metrics (first rows, 2×2 km cells):")
print(head(sf::st_drop_geometry(grid_all_metrics_sf)))
# -> For each 2×2 km cell: ent, relmutinf, shdi, ai
#    This shows how landscape structure varies across the AOI.

message("\nCorrelation matrix of metrics (ent, relmutinf, shdi, ai):")
print(cor_metrics)
# -> Indicates redundancy or independence between IT and classical metrics.

# Additional suggested interactive checks (uncomment when needed):
mapview::mapview(grid_all_metrics_sf, zcol = "ent")       # composition H(x)
mapview::mapview(grid_all_metrics_sf, zcol = "relmutinf") # configuration U
mapview::mapview(grid_all_metrics_sf, zcol = "shdi")      # classical diversity
mapview::mapview(grid_all_metrics_sf, zcol = "ai")        # classical aggregation
#mapview::mapview(lsp_clusters_sf, zcol = "cluster")       # pattern typologies

############################################################
# Objects created in the workspace:
#   it_metrics_global    - tibble with global H(x) and U
#   grid_all_metrics_sf  - sf grid with ent, relmutinf, shdi, ai per 2×2 km cell
#   cor_metrics          - 4×4 correlation matrix of ent, relmutinf, shdi, ai
#   lsp_cove_burgwald    - motif COVE signatures for ~2.5 km windows
#   lsp_clusters_sf      - sf polygons with "cluster" field (pattern types)
############################################################
