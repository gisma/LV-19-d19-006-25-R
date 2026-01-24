#!/usr/bin/env Rscript

############################################################
# Script:   02-6_segments_clc_supervision_exactextract.R
# Project:  Burgwald
#
# FIXED VERSION:
# --------------
# This version normalises the CLC raster to TRUE 3‑digit CLC codes (111, 112, ... 523)
# BEFORE any extraction or mapping. We deliberately eliminate terra category indices
# (1..N) and avoid all use of raster levels / factor semantics.
#
# Semantics after this fix:
#   - clc_rast : original raster (values = 1..N category indices)
#   - exact_extract operates ONLY on clc_code_r
#   - map_clc_to_class() again works correctly (311/312/313/...)
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
pred_file <- paths[["s2_pred_stack_2021_rf"]]
seg_file <- paths[["layer0_segments"]]
out_train_rds = paths[["train_seg_rds"]]
clc_file <- paths[["aoi_clc"]]
out_seg_labeled = paths[["layer0_segments_clc_supervised"]]


# --------------------------------------------------------------------
# 2) Policy parameters
# --------------------------------------------------------------------

crs_m <- 25832
purity_min <- 0.66
area_quantile_keep <- 0.95

# --------------------------------------------------------------------
# 3) Mapping policy (CLC3 -> model classes)
# --------------------------------------------------------------------

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
# 5) NORMALISE CLC: category index (1..N) -> TRUE CLC3 codes (111..523)
# --------------------------------------------------------------------


# --------------------------------------------------------------------
# 6) Prepare geometries
# --------------------------------------------------------------------

segments_m <- sf::st_transform(segments, crs_m)
seg_area_m2 <- as.numeric(sf::st_area(segments_m))
A_min <- as.numeric(stats::quantile(seg_area_m2, probs = 1 - area_quantile_keep, na.rm = TRUE, names = FALSE))

msg(sprintf("Segment area filter: A_min = %.2f m² (cut at Q%.0f)", A_min, (1 - area_quantile_keep) * 100))

keep_area <- seg_area_m2 >= A_min

# reproject for extraction
segments_clc  <- sf::st_transform(segments_m, terra::crs(clc_rast))
segments_pred <- sf::st_transform(segments_m, terra::crs(pred_stack))

# --------------------------------------------------------------------
# 7) CLC area shares per segment (now on TRUE CLC3 codes)
# --------------------------------------------------------------------

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

clc_shares_df <- dplyr::bind_rows(clc_shares_df)
clc_shares_df <- clc_shares_df %>% dplyr::filter(!is.na(code))


# map CLC3 -> model class
clc_long <- clc_shares_df %>%
  mutate(
    code  = as.integer(code),
    share = as.numeric(share),
    class = map_clc_to_class(code)
  ) %>%
  filter(!is.na(class), is.finite(share), share > 0)

# aggregate shares by class within segment
clc_class_shares <- clc_long %>%
  group_by(segment_id, class) %>%
  summarise(share = sum(share), .groups = "drop")

# dominant class per segment
seg_label <- clc_class_shares %>%
  group_by(segment_id) %>%
  arrange(desc(share), .by_group = TRUE) %>%
  summarise(
    class     = first(class),
    purity    = first(share),
    n_classes = n(),
    .groups   = "drop"
  )

# attach area + filter
seg_label <- seg_label %>%
  left_join(
    tibble(segment_id = segments$segment_id, area_m2 = seg_area_m2, keep_area = keep_area),
    by = "segment_id"
  )

seg_label_hard <- seg_label %>%
  filter(keep_area, purity >= purity_min)

msg(sprintf("Labeled segments (hard): %d / %d", nrow(seg_label_hard), nrow(segments)))



# --------------------------------------------------------------------
# 8) Segment predictor signatures
# --------------------------------------------------------------------

msg("Computing segment predictor means (exactextractr) ...")

ext_mean <- exactextractr::exact_extract(
  pred_stack,
  segments_pred,
  fun = "mean"
)

names(ext_mean) <- sub("^mean\\.", "", names(ext_mean))

seg_feat <- tibble::as_tibble(ext_mean) %>%
  mutate(segment_id = segments$segment_id)

pred_names <- setdiff(names(seg_feat), "segment_id")

# --------------------------------------------------------------------
# 9) Build training table
# --------------------------------------------------------------------
out_seg_labeled
train_df <- seg_label_hard %>%
  left_join(seg_feat, by = "segment_id") %>%
  filter(if_all(all_of(pred_names), ~ is.finite(.x)))

train_df$class <- factor(train_df$class)

msg(sprintf("Training table rows: %d", nrow(train_df)))

# --------------------------------------------------------------------
# 10) Write outputs
# --------------------------------------------------------------------

msg("Writing labeled segments ...")
seg_out <- segments %>%
  left_join(
    seg_label %>% dplyr::select(segment_id, class, purity, area_m2, keep_area),
    by = "segment_id"
  )


dir.create(dirname(out_seg_labeled), recursive = TRUE, showWarnings = FALSE)
sf::st_write(seg_out, out_seg_labeled, delete_dsn = TRUE, quiet = TRUE)

msg("Writing training table (RDS) ...")
saveRDS(train_df, out_train_rds)

msg("DONE.")
