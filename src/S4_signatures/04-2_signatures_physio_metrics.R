#!/usr/bin/env Rscript
############################################################
# Script:  04-2_signatures_physio_metrics.R
# Project: Burgwald
#
# Purpose:
#   Build S4_signatures per segment by extracting a COMPLETE signature
#   from canonical S2 stacks (no re-derivation):
#     - relief_stack_10m
#     - hydro_stack_10m
#
# Stats per continuous band:
#   mean, stdev, min, max, coefficient_of_variation, q25, q50, q75
#
# Categorical / ordinal hydro bands:
#   - watershed: mode + mode_share (area-weighted)
#   - strahler:  mode + mode_share (area-weighted) + max
#
# Inputs (paths registry):
#   - layer0_segments
#   - relief_stack_10m
#   - hydro_stack_10m
#   - wind_means_summary (optional)
#
# Output:
#   - layer0_attr_physio_metrics
############################################################

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(terra)
  library(dplyr)
  library(tibble)
  library(exactextractr)
})

source(here::here("src","_core","01-setup-burgwald.R"))

seg_file    <- paths[["layer0_segments"]]
relief_file <- paths[["relief_stack_10m"]]
hydro_file  <- paths[["hydro_stack_10m"]]
out_file    <- paths[["layer0_attr_physio_metrics"]]

message("Input segments:     ", seg_file)
message("Input relief stack: ", relief_file)
message("Input hydro stack:  ", hydro_file)
message("Output:             ", out_file)

stopifnot(file.exists(seg_file), file.exists(relief_file), file.exists(hydro_file))

segments_sf <- sf::read_sf(seg_file)
if (!"segment_id" %in% names(segments_sf)) stop("Missing segment_id in: ", seg_file)

relief <- terra::rast(relief_file)
hydro  <- terra::rast(hydro_file)

# CRS: segments -> raster CRS (relief anchors)
crs_relief <- terra::crs(relief, proj = TRUE)
if (is.na(crs_relief) || crs_relief == "") stop("relief_stack_10m has no CRS: ", relief_file)
if (is.na(sf::st_crs(segments_sf))) stop("Segments have no CRS: ", seg_file)

if (!identical(sf::st_crs(segments_sf)$wkt, crs_relief)) {
  segments_sf <- sf::st_transform(segments_sf, crs_relief)
}
if (any(!sf::st_is_valid(segments_sf))) {
  message("Fixing invalid segment geometries with st_make_valid() ...")
  segments_sf <- sf::st_make_valid(segments_sf)
}

# Hydro CRS must match (pipeline invariant)
crs_hydro <- terra::crs(hydro, proj = TRUE)
if (!identical(crs_hydro, crs_relief)) {
  stop("CRS mismatch: relief_stack_10m vs hydro_stack_10m. Fix in S2 export.")
}

# Relief must contain at least these (for windwardness + sanity); southness is expected from S2
req_relief <- c("dem10m", "slope", "aspect", "southness")
missing_req <- setdiff(req_relief, names(relief))
if (length(missing_req) > 0) {
  stop(
    "relief_stack_10m missing required bands: ",
    paste(missing_req, collapse = ", "),
    "\nCurrent bands: ", paste(names(relief), collapse = ", ")
  )
}

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

# Robust quantile column finder across exactextractr versions
find_quantile_col <- function(colnames_vec, q) {
  # accept: quantile_25, quantile_0.25, quantile.25, quantile0.25, etc.
  pct <- as.integer(round(q * 100))
  # candidates containing 0.25 or 25
  cand <- grep(paste0("^quantile.*(", q, "|", pct, ")$"), colnames_vec, value = TRUE)
  if (length(cand) == 1) return(cand[[1]])
  # fallback: anything starting with quantile that ends with pct
  cand <- grep(paste0("^quantile.*", pct, "$"), colnames_vec, value = TRUE)
  if (length(cand) == 1) return(cand[[1]])
  NA_character_
}

extract_signature_continuous_allbands <- function(r, polys_sf, prefix) {
  bn <- names(r)
  if (length(bn) == 0) stop("Raster has no band names for prefix: ", prefix)
  
  stats_fun <- c("mean", "stdev", "min", "max", "coefficient_of_variation", "quantile")
  qs <- c(0.25, 0.50, 0.75)
  
  out <- tibble::tibble(segment_id = polys_sf$segment_id)
  
  for (b in bn) {
    x <- exactextractr::exact_extract(
      r[[b]],
      polys_sf,
      fun = stats_fun,
      quantiles = qs
    )
    
    nms <- names(x)
    need <- c("mean", "stdev", "min", "max", "coefficient_of_variation")
    if (!all(need %in% nms)) {
      stop("Missing expected stats from exact_extract for band '", b, "': ",
           paste(setdiff(need, nms), collapse = ", "))
    }
    
    q25 <- find_quantile_col(nms, 0.25)
    q50 <- find_quantile_col(nms, 0.50)
    q75 <- find_quantile_col(nms, 0.75)
    
    out[[paste0(prefix, "__", b, "__mean")]] <- as.numeric(x[["mean"]])
    out[[paste0(prefix, "__", b, "__sd")]]   <- as.numeric(x[["stdev"]])
    out[[paste0(prefix, "__", b, "__min")]]  <- as.numeric(x[["min"]])
    out[[paste0(prefix, "__", b, "__max")]]  <- as.numeric(x[["max"]])
    out[[paste0(prefix, "__", b, "__cv")]]   <- as.numeric(x[["coefficient_of_variation"]])
    
    if (!is.na(q25)) out[[paste0(prefix, "__", b, "__q25")]] <- as.numeric(x[[q25]])
    if (!is.na(q50)) out[[paste0(prefix, "__", b, "__q50")]] <- as.numeric(x[[q50]])
    if (!is.na(q75)) out[[paste0(prefix, "__", b, "__q75")]] <- as.numeric(x[[q75]])
  }
  
  out
}

extract_mode_share <- function(values, coverage_fraction) {
  ok <- is.finite(values) & is.finite(coverage_fraction)
  v <- values[ok]
  w <- coverage_fraction[ok]
  
  v <- v[!is.na(v)]
  w <- w[!is.na(w)]
  
  if (length(v) == 0) return(c(mode = NA_real_, mode_share = NA_real_))
  
  tab <- tapply(w, v, sum)  # area-weighted vote
  mode_val <- as.numeric(names(tab)[which.max(tab)])
  mode_w   <- as.numeric(max(tab))
  total_w  <- as.numeric(sum(tab))
  
  c(mode = mode_val, mode_share = if (total_w > 0) mode_w / total_w else NA_real_)
}

# ------------------------------------------------------------------
# 1) Relief: treat as continuous stack (full signature for all bands)
# ------------------------------------------------------------------
relief_sig <- extract_signature_continuous_allbands(relief, segments_sf, prefix = "relief")

# ------------------------------------------------------------------
# 2) Hydro: mixed handling
#    - watershed: mode + mode_share
#    - strahler:  mode + mode_share + max
#    - everything else: continuous signature
# ------------------------------------------------------------------
bn_h <- names(hydro)
if (length(bn_h) == 0) stop("hydro_stack_10m has no band names")

# build out with segment_id
hydro_sig <- tibble::tibble(segment_id = segments_sf$segment_id)

for (b in bn_h) {
  
  if (b %in% c("watershed", "watershed_id", "basin", "basin_id")) {
    
    mode_df <- exactextractr::exact_extract(
      hydro[[b]],
      segments_sf,
      fun = function(values, coverage_fraction) {
        extract_mode_share(values, coverage_fraction)
      }
    )
    
    hydro_sig[[paste0("hydro__", b, "__mode")]]       <- as.numeric(mode_df[["mode"]])
    hydro_sig[[paste0("hydro__", b, "__mode_share")]] <- as.numeric(mode_df[["mode_share"]])
    
    next
  }
  
  if (b %in% c("strahler", "order", "strahler_order")) {
    
    mode_df <- exactextractr::exact_extract(
      hydro[[b]],
      segments_sf,
      fun = function(values, coverage_fraction) {
        ms <- extract_mode_share(values, coverage_fraction)
        ok <- is.finite(values)
        vmax <- if (any(ok)) suppressWarnings(max(values[ok], na.rm = TRUE)) else NA_real_
        c(mode = ms[["mode"]], mode_share = ms[["mode_share"]], max = vmax)
      }
    )
    
    hydro_sig[[paste0("hydro__", b, "__mode")]]       <- as.numeric(mode_df[["mode"]])
    hydro_sig[[paste0("hydro__", b, "__mode_share")]] <- as.numeric(mode_df[["mode_share"]])
    hydro_sig[[paste0("hydro__", b, "__max")]]        <- as.numeric(mode_df[["max"]])
    
    next
  }
  
  # default: treat as continuous
  tmp <- extract_signature_continuous_allbands(hydro[[b]], segments_sf, prefix = "hydro")
  # tmp has: segment_id + hydro__<b>__* columns (because names(hydro[[b]]) == b)
  # merge into hydro_sig
  keep_cols <- setdiff(names(tmp), "segment_id")
  for (cc in keep_cols) hydro_sig[[cc]] <- tmp[[cc]]
}

# ------------------------------------------------------------------
# 3) Optional: windwardness (uses relief band 'aspect' only)
# ------------------------------------------------------------------
has_wind <- ("wind_means_summary" %in% names(paths)) && file.exists(paths[["wind_means_summary"]])
wind_df <- NULL

if (has_wind) {
  wind_means_summary <- readRDS(paths[["wind_means_summary"]])
  
  if (all(c("period","mean_wd") %in% names(wind_means_summary)) && "aspect" %in% names(relief)) {
    
    wd_summer <- wind_means_summary$mean_wd[wind_means_summary$period == "Summer"]
    wd_winter <- wind_means_summary$mean_wd[wind_means_summary$period == "Winter"]
    
    if (length(wd_summer) == 1L && length(wd_winter) == 1L) {
      
      aspect_deg <- relief[["aspect"]]
      windward_summer_r <- cos((aspect_deg - wd_summer) * pi / 180)  # [-1..1]
      windward_winter_r <- cos((aspect_deg - wd_winter) * pi / 180)  # [-1..1]
      
      w_s <- exactextractr::exact_extract(windward_summer_r, segments_sf, "mean")
      w_w <- exactextractr::exact_extract(windward_winter_r, segments_sf, "mean")
      
      wind_df <- tibble::tibble(
        segment_id = segments_sf$segment_id,
        windward__summer__mean = as.numeric(w_s),
        windward__winter__mean = as.numeric(w_w)
      )
      
      message("Windwardness enabled (uses relief band 'aspect').")
    } else {
      message("wind_means_summary missing exactly one Summer and one Winter mean_wd. Skipping windwardness.")
    }
  } else {
    message("wind_means_summary present but missing required columns or relief lacks 'aspect'. Skipping windwardness.")
  }
} else {
  message("wind_means_summary not available. Windwardness omitted.")
}

# ------------------------------------------------------------------
# 4) Join + write
# ------------------------------------------------------------------
geom_col <- attr(segments_sf, "sf_column")
stopifnot(is.character(geom_col), geom_col %in% names(segments_sf))

out_sf <- segments_sf %>%
  dplyr::select(segment_id, dplyr::all_of(geom_col)) %>%
  dplyr::left_join(relief_sig, by = "segment_id") %>%
  dplyr::left_join(hydro_sig,  by = "segment_id")

if (!is.null(wind_df)) out_sf <- out_sf %>% dplyr::left_join(wind_df, by = "segment_id")

out_dir <- dirname(out_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sf::st_write(out_sf, out_file, delete_dsn = TRUE, quiet = TRUE)

message("Wrote: ", out_file)
message("n_segments: ", nrow(out_sf))
message("n_cols: ", ncol(out_sf))
message("Example cols: ", paste(head(names(out_sf), 30), collapse = ", "))
