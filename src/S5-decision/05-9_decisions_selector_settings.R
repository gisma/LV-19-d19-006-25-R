#!/usr/bin/env Rscript
############################################################
# Script:  05-9_decisions_selector_settings.R
# Project: Burgwald
#
# Purpose:
#   Constraint-based candidate selection ("selector") for decision support.
#
# What this script does:
#   - Reads a single segment attribute truth source:
#       layer0_segments_attrstack_metrics  (S4_signatures, gpkg)
#
#   - Executes multiple selection "settings" (scenarios). Each setting is:
#       1) Hard filter:        filter_expr  (forest-only, watershed-only, etc.)
#       2) Coverage target:    target_col + n_per_target (e.g., physio strata)
#       3) Ranking model:      rank_mode = "lex" OR "mcda"
#          - lex:  rank_cols in explicit order (primary, tie-breaker, ...)
#          - mcda: compute mcda_score from multiple dist columns + weights
#       4) Spatial spacing:    greedy thinning with min distance d_min_m
#       5) Points:             point-on-surface inside selected polygons
#
# Outputs:
#   This script is typically used to create non-productive scenario outputs
#   in tmp/ for rapid iteration. You can later promote selected outputs to
#   productive keys if you want.
#
# Notes:
#   - No changes to segmentation geometry.
#   - Selector is the place for constraints + coverage + spacing.
#   - MCDA (if used) is only a scoring model to rank candidates.
############################################################

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(rlang)
  library(units)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))
source(here::here("src", "r-libs", "metrics-fun.R"))   # must contain thin_points_greedy() and mcda_add_score()

## -------------------------------------------------------------------
## 0) Parameters (global knobs for this run)
## -------------------------------------------------------------------

set.seed(42)

# default spacing for scenarios (override per setting if needed)
d_min_default_m <- 500

# default quota per target (override per setting if needed)
n_per_target_default <- 3L

## -------------------------------------------------------------------
## 1) Productive input (from outputs.tsv)
## -------------------------------------------------------------------

attr_file <- paths[["layer0_segments_attrstack_metrics"]]
stopifnot(file.exists(attr_file))

## -------------------------------------------------------------------
## 2) NON-productive output folder (scenario results)
## -------------------------------------------------------------------

tmp_root <- here::here("data", "tmp", "layer0_decisions", "selector_settings")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## 3) Read inputs
## -------------------------------------------------------------------

seg_sf <- sf::read_sf(attr_file)
stopifnot("segment_id" %in% names(seg_sf))

## -------------------------------------------------------------------
## 4) Helper: evaluate filter safely on attribute table
## -------------------------------------------------------------------

apply_filter_expr <- function(x_sf, filter_expr) {
  if (is.null(filter_expr) || !nzchar(filter_expr)) return(x_sf)
  
  # Evaluate expression in the data context (no geometry ops)
  x_tbl <- x_sf %>% sf::st_drop_geometry()
  
  keep <- tryCatch(
    eval(parse(text = filter_expr), envir = x_tbl, enclos = baseenv()),
    error = function(e) stop("filter_expr evaluation failed: ", filter_expr, "\n", e$message)
  )
  
  if (!is.logical(keep) || length(keep) != nrow(x_tbl)) {
    stop("filter_expr must return a logical vector of length nrow(data). Got: ", typeof(keep))
  }
  
  x_sf[which(keep), ]
}

## -------------------------------------------------------------------
## 5) Core selector function (one setting -> segments + points)
## -------------------------------------------------------------------
# Inputs:
#   seg_sf       : sf polygons with segment_id + signature columns
#   setting_id   : string identifier
#   filter_expr  : string expression evaluated in attribute context
#   target_col   : column to cover (strata / classes)
#   n_per_target : quota per target value
#   d_min_m      : minimum spacing distance in meters
#
# Ranking options:
#   rank_mode="lex":
#     rank_cols  : ordered vector of columns (smaller is better)
#
#   rank_mode="mcda":
#     mcda_dist_cols: vector of dist columns to fuse
#     mcda_weights  : named numeric vector with weights for mcda_dist_cols
#     mcda_norm_scope: "filtered" (recommended) or "global"
#
select_setting <- function(seg_sf,
                           setting_id,
                           filter_expr,
                           target_col,
                           n_per_target,
                           d_min_m,
                           rank_mode = c("lex", "mcda"),
                           rank_cols = NULL,
                           mcda_dist_cols = NULL,
                           mcda_weights = NULL,
                           mcda_norm_scope = c("filtered", "global")) {
  
  rank_mode <- match.arg(rank_mode)
  mcda_norm_scope <- match.arg(mcda_norm_scope)
  
  stopifnot("segment_id" %in% names(seg_sf))
  stopifnot(target_col %in% names(seg_sf))
  
  # -----------------------------------------------------------------
  # 1) Hard filter (allowed design space)
  # -----------------------------------------------------------------
  x <- apply_filter_expr(seg_sf, filter_expr)
  
  if (nrow(x) == 0) {
    stop("Setting '", setting_id, "': filter_expr leaves no segments.")
  }
  
  # -----------------------------------------------------------------
  # 2) Ranking preselection: pick n_per_target best per target value
  # -----------------------------------------------------------------
  if (rank_mode == "mcda") {
    stopifnot(!is.null(mcda_dist_cols), length(mcda_dist_cols) >= 1)
    stopifnot(!is.null(mcda_weights))
    stopifnot(all(mcda_dist_cols %in% names(x)))
    
    x <- mcda_add_score(
      seg_sf         = x,
      dist_cols      = mcda_dist_cols,
      weights        = mcda_weights,
      norm_scope     = mcda_norm_scope,
      keep_norm_cols = FALSE,
      score_col      = "mcda_score"
    )
    
    x_tbl <- x %>%
      sf::st_drop_geometry() %>%
      dplyr::select(segment_id, dplyr::all_of(target_col), mcda_score) %>%
      tidyr::drop_na()
    
    pick_ids <- x_tbl %>%
      dplyr::group_by(.data[[target_col]]) %>%
      dplyr::arrange(mcda_score) %>%
      dplyr::slice_head(n = n_per_target) %>%
      dplyr::ungroup() %>%
      dplyr::rename(target = .data[[target_col]])
    
    sel_seg <- x %>%
      dplyr::inner_join(pick_ids, by = "segment_id") %>%
      dplyr::mutate(setting = setting_id)
    
    sel_pts <- sf::st_point_on_surface(sel_seg) %>%
      dplyr::mutate(setting = setting_id) %>%
      dplyr::select(segment_id, setting, target, mcda_score, geometry)
    
  } else {
    stopifnot(!is.null(rank_cols), length(rank_cols) >= 1)
    stopifnot(all(rank_cols %in% names(x)))
    
    x_tbl <- x %>%
      sf::st_drop_geometry() %>%
      dplyr::select(segment_id, dplyr::all_of(target_col), dplyr::all_of(rank_cols)) %>%
      tidyr::drop_na()
    
    ord_expr <- lapply(rank_cols, rlang::sym)
    names(ord_expr) <- NULL
    
    pick_ids <- x_tbl %>%
      dplyr::group_by(.data[[target_col]]) %>%
      dplyr::arrange(!!!ord_expr) %>%
      dplyr::slice_head(n = n_per_target) %>%
      dplyr::ungroup() %>%
      dplyr::rename(target = .data[[target_col]])
    
    sel_seg <- x %>%
      dplyr::inner_join(pick_ids, by = "segment_id") %>%
      dplyr::mutate(setting = setting_id)
    
    sel_pts <- sf::st_point_on_surface(sel_seg) %>%
      dplyr::mutate(setting = setting_id) %>%
      dplyr::select(segment_id, setting, target, dplyr::all_of(rank_cols), geometry)
  }
  
  # -----------------------------------------------------------------
  # 3) Spatial thinning (enforce spacing)
  # -----------------------------------------------------------------
  if (!is.null(d_min_m) && is.finite(d_min_m) && d_min_m > 0) {
    sel_pts_thin <- thin_points_greedy(sel_pts, d_min_m = d_min_m)
    sel_seg_thin <- sel_seg %>% dplyr::semi_join(sel_pts_thin %>% sf::st_drop_geometry(),
                                                 by = "segment_id")
  } else {
    sel_pts_thin <- sel_pts
    sel_seg_thin <- sel_seg
  }
  
  list(
    segments = sel_seg_thin,
    points   = sel_pts_thin
  )
}

## -------------------------------------------------------------------
## 6) Settings (scenarios)
## -------------------------------------------------------------------
# You can keep this as a short declarative list; each run writes two outputs.
#
# Minimal recipe (the 4 knobs):
#   filter_expr  : allowed space
#   target_col   : what you want to cover
#   n_per_target : quota per target
#   d_min_m      : spacing
#
# Ranking:
#   - lex: rank_cols in priority order
#   - mcda: mcda_dist_cols + mcda_weights (rank by mcda_score)
# -------------------------------------------------------------------

settings <- list(
  
  # Example A: forest-only, cover physio strata, lex ranking (physio distance only)
  list(
    id            = "A_forest_physio_lex",
    filter_expr   = "cov_forest >= 0.8",
    target_col    = "physio_stratum",
    n_per_target  = n_per_target_default,
    d_min_m       = d_min_default_m,
    rank_mode     = "lex",
    rank_cols     = c("physio_dist_to_center")
  ),
  
  # Example B: LUCC coverage, cover cov strata, lex ranking (coverage distance then physio tie-break)
  list(
    id            = "B_lucc_cov_lex",
    filter_expr   = "",
    target_col    = "cov_stratum",
    n_per_target  = n_per_target_default,
    d_min_m       = d_min_default_m,
    rank_mode     = "lex",
    rank_cols     = c("cov_dist_to_center", "physio_dist_to_center")
  ),
  
  # Example C: forest-only, cover physio strata, MCDA ranking across multiple domains
  list(
    id              = "C_forest_physio_mcda",
    filter_expr     = "cov_forest >= 0.8",
    target_col      = "physio_stratum",
    n_per_target    = n_per_target_default,
    d_min_m         = d_min_default_m,
    rank_mode       = "mcda",
    mcda_dist_cols  = c("physio_dist_to_center", "it_dist_to_center",
                        "hydro_dist_to_center", "bio_dist_to_center",
                        "cov_dist_to_center"),
    mcda_weights    = c(
      physio_dist_to_center = 1.0,
      it_dist_to_center     = 0.5,
      hydro_dist_to_center  = 1.0,
      bio_dist_to_center    = 0.5,
      cov_dist_to_center    = 1.0
    ),
    mcda_norm_scope = "filtered"
  )
)

## -------------------------------------------------------------------
## 7) Run settings + write outputs (tmp)
## -------------------------------------------------------------------

for (s in settings) {
  message("Running setting: ", s$id)
  
  res <- select_setting(
    seg_sf          = seg_sf,
    setting_id      = s$id,
    filter_expr     = s$filter_expr,
    target_col      = s$target_col,
    n_per_target    = s$n_per_target,
    d_min_m         = s$d_min_m,
    rank_mode       = s$rank_mode,
    rank_cols       = s$rank_cols,
    mcda_dist_cols  = s$mcda_dist_cols,
    mcda_weights    = s$mcda_weights,
    mcda_norm_scope = if (!is.null(s$mcda_norm_scope)) s$mcda_norm_scope else "filtered"
  )
  
  out_seg <- file.path(tmp_root, paste0("sel_", s$id, "_segments.gpkg"))
  out_pts <- file.path(tmp_root, paste0("sel_", s$id, "_points.gpkg"))
  
  sf::st_write(res$segments, out_seg, delete_dsn = TRUE, quiet = TRUE)
  sf::st_write(res$points,   out_pts, delete_dsn = TRUE, quiet = TRUE)
  
  message("  wrote: ", out_seg)
  message("  wrote: ", out_pts)
}
