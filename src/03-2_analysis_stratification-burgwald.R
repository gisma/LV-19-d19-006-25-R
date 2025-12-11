#!/usr/bin/env Rscript
############################################################
# 04_station-stratification-burgwald.R
#
# Purpose
# -------
# Derive objective spatial strata for rain gauge placement in
# the Burgwald catchment, based on:
#
#   1) Information-theoretic landscape metrics
#      (entropy H(x), relative mutual information U)
#      computed on a 2×2 km grid.
#   2) Pattern-based landscape signatures (motif COVE)
#      on ~2.5 km windows, clustered into pattern types.
#
# Resulting objects (IN MEMORY ONLY, no files written):
#   - grid_it_strata_sf:
#       2×2 km grid with IT metrics and stratum_id
#   - station_candidates_sf:
#       several centroid points per IT stratum (candidate gauges)
#   - lsp_cove_burgwald:
#       motif COVE signatures (~2.5 km windows)
#   - lsp_pattern_sf:
#       sf polygons with pattern_cluster (pattern types)
#   - station_with_pattern_sf:
#       station candidates annotated with pattern_cluster
#
############################################################
# Conceptual link to network design (DETAILED EXPLANATION)
# ----------------------------------------------------------
# This script is part of a three-layer network design philosophy
# used in modern hydrological observatories. The logic below
# explains *why* we compute IT metrics, create strata, and later
# refine gauge placement. All three layers are required to design
# a scientifically defensible rain-gauge network.
#
# ==========================================================
# 1) PHYSIOGRAPHIC STRATIFICATION  (STRUCTURAL REPRESENTATIVENESS)
# ==========================================================
# Goal:
#   Ensure that rainfall stations represent the dominant
#   structural units of the landscape. These units differ in:
#     - land-cover composition
#     - forest–openland transitions
#     - fragmentation and clumpiness
#     - canopy interception potential
#     - micro-topographic exposure
#
# Traditional stratification used classes such as:
#   elevation bands, slope/aspect, geology, or forest types.
#
# In this script, we use *information-theoretic landscape metrics*:
#
#   - ent (entropy H(x)):
#       Measures thematic diversity of land cover.
#       Low → homogeneous landscape (e.g., closed forest).
#       High → heterogeneous, mixed, patchy.
#
#   - relmutinf (relative mutual information U):
#       Measures configurational order (clumpiness).
#       Low → fragmented patterns.
#       High → ordered, aggregated patterns.
#
# The combination (H(x), U) gives each 2×2 km grid cell a
# *structural fingerprint*. Clustering these fingerprints produces
# **IT strata**, i.e., groups of areas that share similar spatial
# structure. Each stratum represents a structurally distinct
# landscape region that may support different rainfall regimes.
#
# This mirrors the logic of process-based observatories such as
# CAOS (Attert), Reynolds Creek, and Walnut Gulch, but is more
# objective, reproducible, and data-driven.
#
# ==========================================================
# 2) INFORMATION-GAIN OPTIMISATION  (ANALYTICAL EFFICIENCY)
# ==========================================================
# Goal:
#   Once physiographic strata exist, we determine *where inside*
#   these strata additional stations produce the highest
#   information gain.
#
# This step typically uses:
#   - kriging variance of interpolated rainfall,
#   - radar–gauge residuals,
#   - representativeness error,
#   - redundancy analysis between stations.
#
# Why this is not done first:
#   Pure statistical optimisation is blind to structural
#   landscape drivers, especially in forested terrain.
#
#   → If we optimise without stratification,
#     the model tends to place gauges in open fields or in
#     "statistically convenient" areas, not where rainfall
#     processes actually differ.
#
# Therefore, *information-gain optimisation must operate
# within and across the physiographic strata*, NOT on the full
# catchment at once.
#
# This mirrors the HYREX and Henriksen (2024) approach, where
# initial stratification or geometry is refined with explicit
# uncertainty reduction metrics.
#
# NOTE:
#   This script does **not** perform the optimisation step yet.
#   It prepares the required structural layer.
#
# ==========================================================
# 3) HYDROLOGICAL COUPLING  (FUNCTIONAL ADEQUACY)
# ==========================================================
# Goal:
#   Validate and refine the network using hydrological response.
#   A station is only useful if it helps explain or predict:
#     - streamflow peaks (P–Q coherence),
#     - water balance components,
#     - spatial storm organisation relevant for runoff.
#
# Even an optimally placed station (statistically) may be
# useless hydrologically if it does not capture the rainfall
# regime controlling the hydrograph at key outlets.
#
# Examples from observatories:
#   - Walnut Gulch: one gauge per ephemeral channel.
#   - Reynolds Creek: stations align with snow/rain transitions.
#   - CAOS: gauges positioned to represent "functional units."
#
# Hydrological coupling is the third layer because:
#   - Landscape structure defines where rainfall differs.
#   - Statistical optimisation defines where data are most useful.
#   - Hydrology validates whether the network captures the
#     processes that matter for runoff and storage.
#
# This script does **not** yet integrate discharge or water
# balance data, but the resulting strata and station candidates
# provide the required structure for such coupling.
#
# ==========================================================
# SUMMARY OF THE THREE-LAYER LOGIC
# ==========================================================
#
#   (1) STRUCTURE (IT-METRICS):
#         Identify objective landscape strata using entropy and
#         mutual information → “Where is the landscape different?”
#
#   (2) OPTIMISATION:
#         Refine station density and locations within these
#         strata using kriging variance or radar uncertainty
#         → “Where do we gain the most information?”
#
#   (3) HYDROLOGY:
#         Validate the network using streamflow, storage, or
#         water-balance signals → “Where does rainfall matter
#         for catchment function?”
#
# The script implements Layer 1 entirely:
#   - Builds 2×2 km structural grid
#   - Computes H(x) and U
#   - Clusters cells into structural strata
#   - Identifies representative grid cells
#   - Produces multiple station candidates per stratum
#
# Layers 2 and 3 can be added directly on top of these strata.
#
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

# Root folder for AOI-specific raw data (created in setup script)
aoi_root <- here::here("data", "raw", "AOI_Burgwald")

# Land-cover raster (CLC-like, clipped to Burgwald region)
# Assumed to be a categorical raster with land-cover classes
# at ~100 m resolution (e.g., CLC 2018).
clc_file <- file.path(aoi_root, "clc", "clc5_2018_burgwald.tif")

# AOI polygon (Burgwald outline)
aoi_file <- here::here("data", "processed", "aoi_burgwald.gpkg")

# Read data
lc_burgwald  <- terra::rast(clc_file)
aoi_burgwald <- sf::read_sf(aoi_file)

# Harmonise CRS of AOI and raster so that operations are coherent
if (!sf::st_crs(aoi_burgwald) == sf::st_crs(lc_burgwald)) {
  aoi_burgwald <- sf::st_transform(aoi_burgwald, sf::st_crs(lc_burgwald))
}

# Crop & mask raster to AOI to ensure all metrics refer exactly
# to the Burgwald catchment, not a larger CLC extent
lc_burgwald <- lc_burgwald |>
  terra::crop(terra::vect(aoi_burgwald)) |>
  terra::mask(terra::vect(aoi_burgwald))

# Force integer land-cover classes (landscapemetrics expects
# categorical data as integers)
terra::values(lc_burgwald) <- round(terra::values(lc_burgwald))

# ---------------------------------------------------------
# 2) Macro-level IT metrics on a 2×2 km grid
# ---------------------------------------------------------
# Concept:
#   - Create a regular 2×2 km grid across the Burgwald.
#   - Treat each grid cell as a small "landscape".
#   - Compute two information-theoretic metrics:
#       * ent       = entropy H(x)
#            → compositional diversity (how many classes, how even)
#       * relmutinf = relative mutual information U
#            → configurational order (clumpiness vs randomness)
#   - Each cell gets a coordinate in (H, U) space, which will be
#     used to define strata for rain gauge placement.

# Transform AOI to a metric CRS (ETRS89 / UTM 32N) so that
# grid cell size is specified in metres
aoi_utm <- sf::st_transform(aoi_burgwald, 25832)  # EPSG:25832

# Grid resolution (metres)
grid_cellsize <- 2000  # 2 km × 2 km

# Create 2×2 km grid covering the AOI.
# Each polygon is a candidate "unit landscape" for IT metrics.
grid_utm <- sf::st_make_grid(
  aoi_utm,
  cellsize = grid_cellsize,
  what     = "polygons",
  square   = TRUE
) |>
  sf::st_as_sf() |>
  dplyr::mutate(grid_id = dplyr::row_number())

# Transform the grid back to the raster CRS to align with lc_burgwald
grid <- sf::st_transform(grid_utm, sf::st_crs(lc_burgwald))

# Compute IT metrics per grid cell (landscape-level metrics).
# 'sample_lsm()' cuts out each grid cell as a mini-landscape and
# calculates:
#   - lsm_l_ent       → metric name "ent"
#   - lsm_l_relmutinf → metric name "relmutinf"
lsm_grid_it <- sample_lsm(
  landscape = lc_burgwald,
  y         = grid,
  what      = c("lsm_l_ent", "lsm_l_relmutinf"),
  level     = "landscape"
)

# Convert metrics to wide format:
#   grid_id | ent | relmutinf
lsm_grid_it_wide <- lsm_grid_it |>
  dplyr::select(grid_id = plot_id, metric, value) |>
  tidyr::pivot_wider(names_from = metric, values_from = value)

# Merge IT metrics with grid geometry to obtain an sf object
# where each grid cell has ent and relmutinf values
grid_it_sf <- grid |>
  dplyr::left_join(lsm_grid_it_wide, by = "grid_id")

# ---------------------------------------------------------
# 3) Cluster grid cells into IT-based strata
# ---------------------------------------------------------
# Concept:
#   - Use (ent, relmutinf) to cluster 2×2 km grid cells into a
#     small number of structural strata.
#   - Each stratum is a region with similar land-cover composition
#     and configuration (H(x), U).
#   - Later, we will select several representative cells per stratum
#     as candidate gauge locations.

# Extract numeric IT metrics (ent, relmutinf) without geometry
it_clust_df <- grid_it_sf |>
  sf::st_drop_geometry() |>
  dplyr::select(grid_id, ent, relmutinf) |>
  dplyr::filter(!is.na(ent), !is.na(relmutinf))

# Standardise ent and relmutinf (z-scores) so they have equal weight
it_clust_scaled <- it_clust_df |>
  dplyr::mutate(
    ent_z       = as.numeric(scale(ent)),
    relmutinf_z = as.numeric(scale(relmutinf))
  )

# Number of strata (you can change this if needed)
k_strata <- 5

set.seed(123)
km_it <- stats::kmeans(
  it_clust_scaled[, c("ent_z", "relmutinf_z")],
  centers = k_strata,
  nstart  = 50
)

# Attach stratum IDs to the scaled table
it_clust_scaled$stratum_id <- km_it$cluster

# Attach stratum IDs back to the grid sf (geometry + IT metrics)
grid_it_strata_sf <- grid_it_sf |>
  dplyr::left_join(
    it_clust_scaled |>
      dplyr::select(grid_id, stratum_id),
    by = "grid_id"
  )

# ---------------------------------------------------------
# 4) Derive multiple candidate station locations per IT stratum
# ---------------------------------------------------------
# Concept:
#   - For each stratum, compute the centroid in (ent_z, relmutinf_z) space.
#   - For each cell, compute its distance to the stratum centroid.
#   - Select several nearest cells per stratum → "most typical" cells.
#   - Use the centroid of each selected polygon as station candidate.
#
# This allows, for example, 20 candidates by choosing 4 cells
# per stratum if k_strata = 5.

# 4.1 Compute stratum centroids in z-score space
strata_centroids <- it_clust_scaled |>
  dplyr::group_by(stratum_id) |>
  dplyr::summarise(
    ent_z_mean       = mean(ent_z, na.rm = TRUE),
    relmutinf_z_mean = mean(relmutinf_z, na.rm = TRUE),
    .groups = "drop"
  )

# 4.2 Distance of each cell to its stratum centroid
it_clust_scaled <- it_clust_scaled |>
  dplyr::left_join(strata_centroids, by = "stratum_id") |>
  dplyr::mutate(
    dist_to_center = sqrt(
      (ent_z       - ent_z_mean)^2 +
        (relmutinf_z - relmutinf_z_mean)^2
    )
  )

# 4.3 Select several representative cells per stratum
#     Here: aim for ~20 stations in total
target_n_stations <- 20
n_per_stratum     <- ceiling(target_n_stations / k_strata)

rep_cells_basic <- it_clust_scaled |>
  dplyr::group_by(stratum_id) |>
  dplyr::slice_min(
    order_by  = dist_to_center,
    n         = n_per_stratum,
    with_ties = FALSE
  ) |>
  dplyr::ungroup() |>
  dplyr::select(
    grid_id,
    stratum_id,
    ent_z,
    relmutinf_z
  )

# 4.4 Bring back geometries and raw IT metrics for the representative cells
rep_cells_sf <- grid_it_sf |>
  dplyr::inner_join(
    rep_cells_basic,
    by = "grid_id"
  )
# rep_cells_sf now contains:
#   grid_id, ent, relmutinf, stratum_id, ent_z, relmutinf_z, geometry

# 4.5 Compute centroids of representative cells as station candidates
# 4.5 Compute centroids of representative cells as station candidates
# NOTE:
#   - In rep_cells_sf the geometry column is called 'x' (see printout).
#   - sf keeps that as the active geometry column.
#   - We therefore DO NOT mention 'geometry' explicitly in select();
#     sf will keep the geometry column automatically.

station_candidates_sf <- rep_cells_sf |>
  sf::st_centroid() |>
  dplyr::mutate(
    station_id = dplyr::row_number()
  ) |>
  dplyr::select(
    station_id,
    stratum_id,
    grid_id,
    ent,
    relmutinf,
    ent_z,
    relmutinf_z
    # geometry column is kept implicitly by sf
  )


# station_candidates_sf now contains:
#   - station_id: running ID of candidate station
#   - stratum_id: IT-based stratum
#   - grid_id   : ID of the representative 2×2 km cell
#   - ent, relmutinf: raw IT metrics of that cell
#   - ent_z, relmutinf_z: z-scores (for interpretation)
#   - geometry: point (cell centroid) as gauge candidate location

# ---------------------------------------------------------
# 5) Pattern-based analysis with motif (optional refinement)
# ---------------------------------------------------------
# Concept:
#   - motif computes co-occurrence (COVE) signatures on moving windows.
#   - window = 25 → 25×25 cells; at ~100 m resolution this is
#     ≈ 2.5 km × 2.5 km.
#   - Each window yields a probability distribution of adjacency
#     patterns between classes.
#   - Distances between signatures (e.g. Jensen–Shannon) allow
#     clustering into pattern types ("texture regimes").
#   - These pattern types can be used to interpret or further
#     refine the IT strata and station placement.

# Compute COVE signatures on the land-cover raster
lsp_cove_burgwald <- lsp_signature(
  x             = lc_burgwald,
  type          = "cove",
  window        = 25,        # 25 cells ≈ 2.5 km
  normalization = "pdf"      # probability distribution
)

# Compute distances between signatures using Jensen–Shannon divergence,
# a symmetric distance between probability distributions
dist_mat <- lsp_to_dist(
  x        = lsp_cove_burgwald,
  dist_fun = "jensen-shannon"
)

# Choose number of pattern clusters (pattern types)
k_pattern <- 6
set.seed(123)
hc <- hclust(dist_mat, method = "ward.D2")
pattern_ids <- cutree(hc, k = k_pattern)

# Attach pattern_cluster to the lsp object
lsp_cove_burgwald$pattern_cluster <- pattern_ids

# Convert signatures to sf polygons
lsp_pattern_sf <- lsp_add_sf(lsp_cove_burgwald)
# lsp_pattern_sf contains:
#   - geometry: polygon (window footprint)
#   - pattern_cluster: ID of the pattern type
#   - plus additional list-columns with signatures

# ---------------------------------------------------------
# 6) Attach pattern types to station candidates (optional)
# ---------------------------------------------------------
# Concept:
#   - For each station candidate, determine which pattern_cluster
#     it falls into (point-in-polygon).
#   - This provides additional context: which structural pattern
#     regime the candidate station is representing.

station_with_pattern_sf <- sf::st_join(
  station_candidates_sf,
  lsp_pattern_sf |> dplyr::select(pattern_cluster),
  join = sf::st_intersects,
  left = TRUE
)

# station_with_pattern_sf:
#   - one row per station candidate
#   - station_id, stratum_id (IT stratum)
#   - ent, relmutinf, ent_z, relmutinf_z
#   - pattern_cluster (motif-based pattern type, if any intersection)

# ---------------------------------------------------------
# 7) Object overview (console only)
# ---------------------------------------------------------

message("IT-based grid strata (first rows, no geometry):")
print(
  grid_it_strata_sf |>
    sf::st_drop_geometry() |>
    dplyr::select(grid_id, ent, relmutinf, stratum_id) |>
    head()
)

message("\nStation candidates per IT stratum (with pattern types if available):")
print(
  station_with_pattern_sf |>
    sf::st_drop_geometry() |>
    dplyr::select(
      station_id,
      stratum_id,
      grid_id,
      ent,
      relmutinf,
      ent_z,
      relmutinf_z,
      pattern_cluster
    )
)

############################################################
# Objects created in the workspace:
#
#   lc_burgwald
#       Land-cover raster (Burgwald AOI, categorical classes).
#
#   aoi_burgwald
#       AOI polygon.
#
#   grid_it_strata_sf
#       2×2 km sf grid with:
#         - ent (entropy H(x))
#         - relmutinf (relative mutual information U)
#         - stratum_id (IT-based cluster ID)
#
#   station_candidates_sf
#       Several centroids per IT stratum, with:
#         - station_id
#         - ent, relmutinf (raw IT metrics)
#         - ent_z, relmutinf_z (z-scores)
#         - geometry (point coordinates)
#
#   lsp_cove_burgwald
#       motif COVE signatures with:
#         - pattern_cluster (pattern-type ID)
#
#   lsp_pattern_sf
#       sf polygons (signature windows) with:
#         - pattern_cluster
#         - additional list-columns for signatures
#
#   station_with_pattern_sf
#       station_candidates_sf with pattern_cluster attached
#       via spatial join (point-in-polygon).
#
# Suggested interactive plotting (run manually in RStudio):
#
#   library(mapview)
#
#   # 1) IT strata as 2×2 km grid
#   mapview::mapview(grid_it_strata_sf, zcol = "stratum_id")
#
#   # 2) Station candidates on top of strata
   mapview::mapview(grid_it_strata_sf, zcol = "stratum_id") +
     mapview::mapview(station_with_pattern_sf, zcol = "stratum_id", cex = 4)
#
#   # 3) Pattern clusters (COVE) – remove list-columns for mapview
   geom_col <- attr(lsp_pattern_sf, "sf_column")
   lsp_pattern_sf_plot <- lsp_pattern_sf[, c("pattern_cluster", geom_col)]
   lsp_pattern_sf_plot$pattern_cluster <- as.factor(lsp_pattern_sf_plot$pattern_cluster)
   mapview::mapview(lsp_pattern_sf_plot, zcol = "pattern_cluster")
#
############################################################
