#!/usr/bin/env Rscript
############################################################
# Script:  05-9_decisions_selector_settings.R
# Project: Burgwald
#
# Purpose:
#   S5 selector: produce spatially distributed candidate sets for
#   different settings from existing S5 domain strata outputs.
#
# Core logic:
#   A) Filter (hard constraints, e.g. cov_forest >= 0.8, watershed_major)
#   B) Coverage targets (e.g. physio_stratum, coverage_stratum, hydro_stratum)
#   C) Quality ranking (dist_to_center per domain)
#   D) Spatial thinning (min distance, greedy) -> points-in-segments
#
# Inputs (productive, from outputs.tsv via paths[]):
#   - layer0_it_strata_segments
#   - layer0_physio_strata_segments
#   - layer0_hydro_strata_segments
#   - layer0_biostructure_strata_segments
#   - layer0_coverage_strata_segments
#
# Outputs:
#   - written to tmp_root (NON-productive) as separate gpkg per setting:
#       segments_<setting>.gpkg
#       points_<setting>.gpkg
############################################################

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))
source(here::here("src", "r-libs", "metrics-fun.R"))

## -------------------------------------------------------------------
## 0) Parameters
## -------------------------------------------------------------------

d_min_m <- 500  # minimum spacing in meters

# default per-setting quotas (can be edited in-place)
n_per_physio   <- 3L
n_per_coverage <- 3L
n_per_hydro    <- 3L

set.seed(42)

## -------------------------------------------------------------------
## 1) Productive inputs
## -------------------------------------------------------------------

it_seg_file     <- paths[["layer0_it_strata_segments"]]
physio_seg_file <- paths[["layer0_physio_strata_segments"]]
hydro_seg_file  <- paths[["layer0_hydro_strata_segments"]]
bio_seg_file    <- paths[["layer0_biostructure_strata_segments"]]
cov_seg_file    <- paths[["layer0_coverage_strata_segments"]]

stopifnot(file.exists(it_seg_file))
stopifnot(file.exists(physio_seg_file))
stopifnot(file.exists(hydro_seg_file))
stopifnot(file.exists(bio_seg_file))
stopifnot(file.exists(cov_seg_file))

## -------------------------------------------------------------------
## 2) tmp output folder (NON-productive)
## -------------------------------------------------------------------

tmp_root <- here::here("data", "tmp", "layer0_decisions", "tmp_selector_settings")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Read + build a single segment table (geometry backbone = IT layer)
## -------------------------------------------------------------------

it_sf     <- sf::read_sf(it_seg_file)
physio_sf <- sf::read_sf(physio_seg_file)
hydro_sf  <- sf::read_sf(hydro_seg_file)
bio_sf    <- sf::read_sf(bio_seg_file)
cov_sf    <- sf::read_sf(cov_seg_file)

stopifnot("segment_id" %in% names(it_sf))
stopifnot("it_stratum" %in% names(it_sf))
stopifnot("it_dist_to_center" %in% names(it_sf))

stopifnot("physio_stratum" %in% names(physio_sf))
stopifnot("physio_dist_to_center" %in% names(physio_sf))

stopifnot("hydro_stratum" %in% names(hydro_sf))
stopifnot("hydro_dist_to_center" %in% names(hydro_sf))

stopifnot("bio_stratum" %in% names(bio_sf))
stopifnot("bio_dist_to_center" %in% names(bio_sf))

stopifnot("cov_stratum" %in% names(cov_sf))
stopifnot("cov_dist_to_center" %in% names(cov_sf))

# Keep all attributes from IT geometry layer, then join other strata/dists
seg_sf <- it_sf %>%
  dplyr::select(segment_id, geometry, dplyr::everything())

seg_sf <- seg_sf %>%
  dplyr::left_join(
    physio_sf %>% sf::st_drop_geometry() %>%
      dplyr::select(segment_id, physio_stratum, physio_dist_to_center),
    by = "segment_id"
  ) %>%
  dplyr::left_join(
    hydro_sf %>% sf::st_drop_geometry() %>%
      dplyr::select(segment_id, hydro_stratum, hydro_dist_to_center),
    by = "segment_id"
  ) %>%
  dplyr::left_join(
    bio_sf %>% sf::st_drop_geometry() %>%
      dplyr::select(segment_id, bio_stratum, bio_dist_to_center),
    by = "segment_id"
  ) %>%
  dplyr::left_join(
    cov_sf %>% sf::st_drop_geometry() %>%
      dplyr::select(segment_id, cov_stratum, cov_dist_to_center),
    by = "segment_id"
  )

# Contract from you:
stopifnot("cov_forest" %in% names(seg_sf))

## -------------------------------------------------------------------
## 4) Helper: greedy spatial thinning on points
## -------------------------------------------------------------------

thin_points_greedy <- function(pts_sf, d_min_m) {
  if (nrow(pts_sf) <= 1) return(pts_sf)
  
  keep <- rep(FALSE, nrow(pts_sf))
  for (i in seq_len(nrow(pts_sf))) {
    if (!any(keep)) {
      keep[i] <- TRUE
      next
    }
    # check distance to already kept points
    ok <- TRUE
    for (j in which(keep)) {
      if (sf::st_is_within_distance(pts_sf[i, ], pts_sf[j, ], dist = d_min_m)[[1]] |> length() > 0) {
        ok <- FALSE
        break
      }
    }
    if (ok) keep[i] <- TRUE
  }
  pts_sf[keep, ]
}

## -------------------------------------------------------------------
## 5) Helper: select top-n per target, then thin
## -------------------------------------------------------------------

select_setting <- function(seg_sf, setting_id, filter_expr, target_col, rank_cols, n_per_target, d_min_m) {
  # 1) filter
  x <- seg_sf
  if (!is.null(filter_expr)) {
    x <- x %>% dplyr::filter(!!rlang::parse_expr(filter_expr))
  }
  
  # drop segments without target or rank columns
  needed <- c(target_col, rank_cols)
  stopifnot(all(needed %in% names(x)))
  
  x_tbl <- x %>% sf::st_drop_geometry() %>%
    dplyr::select(segment_id, dplyr::all_of(target_col), dplyr::all_of(rank_cols))
  
  x_tbl <- x_tbl %>% tidyr::drop_na()
  
  # 2) ranking (lexicographic: rank_cols in given order)
  #    smaller dist = better
  ord_expr <- lapply(rank_cols, function(cc) rlang::sym(cc))
  names(ord_expr) <- NULL
  
  pick_ids <- x_tbl %>%
    dplyr::group_by(.data[[target_col]]) %>%
    dplyr::arrange(!!!ord_expr) %>%
    dplyr::slice_head(n = n_per_target) %>%
    dplyr::ungroup() %>%
    dplyr::select(segment_id, target = .data[[target_col]], dplyr::all_of(rank_cols))
  
  # 3) segments
  sel_seg <- x %>% dplyr::inner_join(pick_ids, by = "segment_id") %>%
    dplyr::mutate(setting = setting_id)
  
  # 4) points in segments + thinning
  sel_pts <- sf::st_point_on_surface(sel_seg) %>%
    dplyr::mutate(setting = setting_id) %>%
    dplyr::select(segment_id, setting, target, dplyr::all_of(rank_cols), geometry)
  
  sel_pts_thin <- thin_points_greedy(sel_pts, d_min_m = d_min_m)
  
  # keep only segments that survived thinning
  sel_seg_thin <- sel_seg %>%
    dplyr::semi_join(sel_pts_thin %>% sf::st_drop_geometry() %>% dplyr::select(segment_id), by = "segment_id")
  
  list(segments = sel_seg_thin, points = sel_pts_thin)
}

## -------------------------------------------------------------------
## 6) Settings (edit in-place)
## -------------------------------------------------------------------
# Setting A: forest-only, cover physio strata, quality = physio representative
res_A <- select_setting(
  seg_sf      = seg_sf,
  setting_id  = "A_forest_physio",
  filter_expr = "cov_forest >= 0.8",
  target_col  = "physio_stratum",
  rank_cols   = c("physio_dist_to_center"),
  n_per_target= n_per_physio,
  d_min_m     = d_min_m
)

# Setting B: across LUCC types (kmeans coverage_stratum), quality = LUCC representative,
#            tie-break by physio representative (spatial distribution + process nuance)
res_B <- select_setting(
  seg_sf      = seg_sf,
  setting_id  = "B_lucc_coverage",
  filter_expr = NULL,
  target_col  = "cov_stratum",
  rank_cols   = c("cov_dist_to_center", "physio_dist_to_center"),
  n_per_target= n_per_coverage,
  d_min_m     = d_min_m
)

# Setting C: within each watershed_major (if available), cover hydro strata
#            (writes one gpkg per watershed into tmp_root)
do_C <- "watershed_major" %in% names(seg_sf)

## -------------------------------------------------------------------
## 7) Write outputs (NON-productive)
## -------------------------------------------------------------------

write_setting <- function(res, tmp_root) {
  seg_out <- file.path(tmp_root, paste0("segments_", unique(res$segments$setting), ".gpkg"))
  pts_out <- file.path(tmp_root, paste0("points_",   unique(res$points$setting), ".gpkg"))
  
  sf::st_write(res$segments, seg_out, delete_dsn = TRUE, quiet = TRUE)
  sf::st_write(res$points,   pts_out, delete_dsn = TRUE, quiet = TRUE)
  
  message("Wrote: ", seg_out)
  message("Wrote: ", pts_out)
}

write_setting(res_A, tmp_root)
write_setting(res_B, tmp_root)

if (do_C) {
  ws_ids <- sort(unique(seg_sf$watershed_major))
  ws_ids <- ws_ids[!is.na(ws_ids)]
  
  for (w in ws_ids) {
    res_Cw <- select_setting(
      seg_sf      = seg_sf,
      setting_id  = paste0("C_ws_", w, "_hydro"),
      filter_expr = paste0("watershed_major == ", w),
      target_col  = "hydro_stratum",
      rank_cols   = c("hydro_dist_to_center"),
      n_per_target= n_per_hydro,
      d_min_m     = d_min_m
    )
    write_setting(res_Cw, tmp_root)
  }
} else {
  message("Skip setting C (watershed_major not present in segment table).")
}
