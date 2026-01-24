#!/usr/bin/env Rscript

############################################################
# Script:   02-5_processing_seg_clc_trainDF.R
# Project:  Burgwald
#
# What this script does (high level)
# ----------------------------------
# This script creates a *segment-based* training table by transferring land-cover
# supervision from a CLC raster to your segment polygons, and then attaching
# segment-wise predictor signatures (means from a reduced S2 predictor stack).
#
# Core outputs:
#   (1) out_seg_labeled (vector): all segments + (class, purity, area_m2, keep_area)
#   (2) out_train_rds   (RDS):    *hard-labeled* segments + predictor means
#
# Why "purity" exists
# -------------------
# For each segment we compute class shares (from CLC area fractions within the segment).
# The dominant class share is called "purity".
#
# purity = max_k( share_k )   where share_k are aggregated class shares within the segment.
#
# Purity is used as a policy gate to reduce label noise:
#   - segments with low purity are mixed/transition segments and are excluded from training
#     (but they remain in out_seg_labeled for diagnostics / downstream mapping).
#
# Downstream implications (important)
# -----------------------------------
# - out_train_rds is the *only* object used for supervised model fitting downstream.
# - Therefore, purity_min directly controls:
#     * training set size (n rows)
#     * class balance (which classes survive the hard filter)
#     * label noise (how many mixed segments are allowed)
# - Downstream scripts will train on cleaner labels if purity_min is higher,
#   at the cost of fewer training samples and potentially reduced spatial coverage.
#
# FIXED VERSION note
# ------------------
# Header claims a "normalise CLC raster to TRUE 3-digit codes" step exists.
# The intention is: avoid category-index semantics and ensure codes are real CLC3
# values (111,112,...,523) before mapping.
# In the current code body, the extraction/mapping step uses clc_rast directly.
# The mapping policy assumes the raster values are *already* CLC3 codes.
############################################################

suppressPackageStartupMessages({
  library(here)
  library(terra)
  library(sf)
  library(dplyr)
  library(tibble)
  library(exactextractr)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))

quiet <- FALSE
msg <- function(...) if (!quiet) message(...)

# --------------------------------------------------------------------
# 1) Canonical IO (NO invented keys; derive filenames from existing ones)
# --------------------------------------------------------------------
# pred_file: reduced predictor stack (already feature-selected upstream)
# seg_file : polygon segments (must contain segment_id)
# clc_file : CLC raster used for supervision transfer
# out_train_rds     : RDS with the *training table* (hard-labeled)
# out_seg_labeled   : vector layer with per-segment label diagnostics (soft + hard info)
pred_file <- paths[["s2_pred_stack_2021_rf"]]
seg_file <- paths[["layer0_segments"]]
out_train_rds = paths[["train_seg_rds"]]
clc_file <- paths[["aoi_clc"]]
out_seg_labeled = paths[["layer0_segments_clc_supervised"]]

# --------------------------------------------------------------------
# 2) Policy parameters
# --------------------------------------------------------------------
# crs_m: metric CRS for area computation (m²)
# purity_min:
#   - minimum dominant-class share required to accept a segment as "hard-labeled"
#   - applies only to training set creation (train_df), not to labeled segments export
# area_quantile_keep:
#   - defines a minimum segment area threshold via quantile cut
#   - keep_area = segment_area >= A_min, where A_min = quantile(area, 1 - area_quantile_keep)
crs_m <- 25832
purity_min <- 0.66
area_quantile_keep <- 0.95

# --------------------------------------------------------------------
# 3) Mapping policy (CLC3 -> model classes)
# --------------------------------------------------------------------
# map_clc_to_class():
#   maps *CLC Level-3* codes (integers like 311/312/313/...) to your model classes.
#
# IMPORTANT:
#   This mapping assumes raster values are TRUE CLC3 codes.
#   If the raster stores category indices (1..N) you must normalise to codes first,
#   otherwise this mapping will silently produce NA/incorrect classes.
map_clc_to_class <- function(code) {
  code <- as.integer(code)
  
  dplyr::case_when(
    # water
    code >= 500 & code < 600 ~ "water",
    
    # forest split
    code == 311 ~ "forest_broadleaf",
    code == 312 ~ "forest_coniferous",
    code == 313 ~ "forest_mixed",
    
    # cropland
    code >= 200 & code < 300 ~ "cropland",
    
    # grassland
    code %in% c(231, 321) ~ "grassland",
    
    # mixed green
    code %in% c(141, 142) ~ "mixed_green",
    
    # urban sealed
    code >= 100 & code < 200 ~ "urban_sealed",
    
    TRUE ~ NA_character_
  )
}

# --------------------------------------------------------------------
# 4) Load data
# --------------------------------------------------------------------
# pred_stack: predictor stack (raster) from which we compute segment means
# segments  : segment polygons; required: segment_id
# clc_rast  : CLC raster used for class shares / purity
msg("Loading predictor stack ...")
pred_stack <- terra::rast(pred_file)

msg("Loading segments ...")
segments <- sf::st_read(seg_file, quiet = TRUE)

if (!("segment_id" %in% names(segments))) {
  stop("Segments must contain 'segment_id' column.")
}

msg("Loading raw CLC raster (category indices) ...")
clc_rast <- terra::rast(clc_file)

# --------------------------------------------------------------------
# 5) Prepare geometries
# --------------------------------------------------------------------
# - compute segment area in a metric CRS (m²)
# - derive minimum area threshold from quantile (area_quantile_keep)
# - keep_area is an additional policy filter (independent of purity)
segments_m <- sf::st_transform(segments, crs_m)
seg_area_m2 <- as.numeric(sf::st_area(segments_m))
A_min <- as.numeric(stats::quantile(seg_area_m2, probs = 1 - area_quantile_keep, na.rm = TRUE, names = FALSE))

msg(sprintf("Segment area filter: A_min = %.2f m² (cut at Q%.0f)", A_min, (1 - area_quantile_keep) * 100))

keep_area <- seg_area_m2 >= A_min

# Reproject segments to match rasters for exact extraction:
# - segments_clc  for extracting class shares from CLC raster
# - segments_pred for extracting predictor means from predictor stack
segments_clc  <- sf::st_transform(segments_m, terra::crs(clc_rast))
segments_pred <- sf::st_transform(segments_m, terra::crs(pred_stack))

# --------------------------------------------------------------------
# 6) CLC area shares per segment (exactextractr)
# --------------------------------------------------------------------
# exact_extract() is used with a custom fun() that receives:
#   df: values under the polygon coverage
#   coverage_fraction: per-cell coverage weights for the polygon
#
# For each segment_id:
#   - group by CLC code value v
#   - sum coverage weights w per code (area share proxy)
#   - normalise to shares that sum to 1
#
# Output clc_shares_df is a long table:
#   segment_id | code | share
msg("Computing CLC3 shares per segment (exactextractr) ...")

clc_shares_df <- exactextractr::exact_extract(
  clc_rast,
  segments_clc,
  include_cols = "segment_id",
  fun = function(df, coverage_fraction) {
    
    v   <- as.integer(df$value)
    sid <- df$segment_id[1]
    w   <- as.numeric(coverage_fraction)
    
    ok <- is.finite(v) & is.finite(w) & w > 0
    if (!any(ok)) {
      return(data.frame(
        segment_id = sid,
        code  = NA_integer_,
        share = NA_real_
      ))
    }
    
    tab <- tapply(w[ok], v[ok], sum, default = 0)
    out <- data.frame(
      segment_id = sid,
      code  = as.integer(names(tab)),
      share = as.numeric(tab)
    )
    
    s <- sum(out$share)
    if (is.finite(s) && s > 0) out$share <- out$share / s
    out
  }
)

# Flatten list-of-dfs returned by exact_extract() into one data frame
clc_shares_df <- dplyr::bind_rows(clc_shares_df)

# Drop NA codes (segments with no valid overlap)
clc_shares_df <- clc_shares_df %>% dplyr::filter(!is.na(code))

# Map CLC3 -> model class and keep only valid mapped classes
clc_long <- clc_shares_df %>%
  mutate(
    code  = as.integer(code),
    share = as.numeric(share),
    class = map_clc_to_class(code)
  ) %>%
  filter(!is.na(class), is.finite(share), share > 0)

# Aggregate shares by model class within each segment (collapse multiple codes to one class)
clc_class_shares <- clc_long %>%
  group_by(segment_id, class) %>%
  summarise(share = sum(share), .groups = "drop")

# Dominant class per segment + purity diagnostics:
# - class     : class with largest share
# - purity    : share of that dominant class
# - n_classes : number of model classes present (after aggregation)
seg_label <- clc_class_shares %>%
  group_by(segment_id) %>%
  arrange(desc(share), .by_group = TRUE) %>%
  summarise(
    class     = first(class),
    purity    = first(share),
    n_classes = n(),
    .groups   = "drop"
  )

# Attach geometry-derived diagnostics (area + keep_area policy flag)
seg_label <- seg_label %>%
  left_join(
    tibble(segment_id = segments$segment_id, area_m2 = seg_area_m2, keep_area = keep_area),
    by = "segment_id"
  )

# Hard-label selection (THIS is where purity impacts training):
# Keep only segments that are:
#   (a) large enough (keep_area == TRUE)
#   (b) sufficiently "pure" (purity >= purity_min)
seg_label_hard <- seg_label %>%
  filter(keep_area, purity >= purity_min)

msg(sprintf("Labeled segments (hard): %d / %d", nrow(seg_label_hard), nrow(segments)))

# --------------------------------------------------------------------
# 7) Segment predictor signatures (means)
# --------------------------------------------------------------------
# Compute per-segment predictor means from pred_stack via exactextractr.
# fun="mean" returns one row per polygon with columns mean.<layername>.
msg("Computing segment predictor means (exactextractr) ...")

ext_mean <- exactextractr::exact_extract(
  pred_stack,
  segments_pred,
  fun = "mean"
)

# Strip exactextractr "mean." prefix to match predictor layer names
names(ext_mean) <- sub("^mean\\.", "", names(ext_mean))

# Attach segment_id (must align with the segments order used in exact_extract)
seg_feat <- tibble::as_tibble(ext_mean) %>%
  mutate(segment_id = segments$segment_id)

pred_names <- setdiff(names(seg_feat), "segment_id")

# --------------------------------------------------------------------
# 8) Build training table
# --------------------------------------------------------------------
# Training table is the join of:
#   seg_label_hard (segment_id + class/purity diagnostics)
#   seg_feat       (segment_id + predictor means)
# plus: filtering finite predictor values (no NA/Inf/NaN)
#
# Downstream: this is the table used for supervised learning.
out_seg_labeled
train_df <- seg_label_hard %>%
  left_join(seg_feat, by = "segment_id") %>%
  filter(if_all(all_of(pred_names), ~ is.finite(.x)))

train_df$class <- factor(train_df$class)

msg(sprintf("Training table rows: %d", nrow(train_df)))

# --------------------------------------------------------------------
# 9) Write outputs
# --------------------------------------------------------------------
# (A) Labeled segments (ALL segments): this includes:
#     - segments with low purity
#     - segments failing keep_area
#     These are kept for diagnostics / mapping, but they are NOT used in training.
msg("Writing labeled segments ...")
seg_out <- segments %>%
  left_join(
    seg_label %>% dplyr::select(segment_id, class, purity, area_m2, keep_area),
    by = "segment_id"
  )

dir.create(dirname(out_seg_labeled), recursive = TRUE, showWarnings = FALSE)
sf::st_write(seg_out, out_seg_labeled, delete_dsn = TRUE, quiet = TRUE)

# (B) Training table (ONLY hard-labeled segments):
#     This RDS is the training contract consumed by downstream model scripts.
#     purity_min and keep_area are "baked in" here.
msg("Writing training table (RDS) ...")
saveRDS(train_df, out_train_rds)

msg("DONE.")
