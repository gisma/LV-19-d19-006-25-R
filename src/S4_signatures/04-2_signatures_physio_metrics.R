#!/usr/bin/env Rscript
############################################################
# Script:  04-2_signatures_physio_metrics.R
# Project: Burgwald
#
# Purpose
# -------
# Build S4 (Signatures) from canonical S2 raster stacks by aggregating values
# over S3 segments. No physical feature derivation is allowed here.
#
# Output concept
# --------------
# A segment is a polygonal spatial unit (S3). A signature is a vector of
# statistics per segment (S4). The output is a GeoPackage containing:
#   - geometry (segment polygons)
#   - a wide attribute table: one column per (band × statistic)
#
# Continuous vs. categorical/ordinal layers
# -----------------------------------------
# For continuous raster bands, a "complete" signature is extracted:
#   mean, stdev, min, max, coefficient_of_variation, quantiles (25/50/75)
#
# Some hydro layers are not continuous:
#   - watershed: categorical basin ID
#   - strahler : ordinal stream order
# For these, mean/quantiles are not meaningful; instead:
#   - watershed: majority (mode) + purity (mode_share)
#   - strahler : majority (mode) + purity (mode_share) + max order present
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

# ---------------------------------------------------------
# 0) Inputs / Outputs (canonical path registry)
# ---------------------------------------------------------
seg_file    <- paths[["layer0_segments"]]
relief_file <- paths[["relief_stack_10m"]]
hydro_file  <- paths[["hydro_stack_10m"]]
out_file    <- paths[["layer0_attr_physio_metrics"]]

message("Input segments:     ", seg_file)
message("Input relief stack: ", relief_file)
message("Input hydro stack:  ", hydro_file)
message("Output:             ", out_file)

stopifnot(file.exists(seg_file), file.exists(relief_file), file.exists(hydro_file))

# ---------------------------------------------------------
# 1) Read inputs
# ---------------------------------------------------------
segments_sf <- sf::read_sf(seg_file)
if (!"segment_id" %in% names(segments_sf)) stop("Missing segment_id in: ", seg_file)

relief <- terra::rast(relief_file)
hydro  <- terra::rast(hydro_file)

# ---------------------------------------------------------
# 2) CRS alignment (mandatory for exact extraction)
# ---------------------------------------------------------
# exactextractr expects polygons and rasters in the same CRS. Segments are
# transformed to the raster CRS (relief anchors the reference).
crs_relief <- terra::crs(relief, proj = TRUE)
if (is.na(crs_relief) || crs_relief == "") stop("Relief stack has no CRS: ", relief_file)
if (is.na(sf::st_crs(segments_sf))) stop("Segments have no CRS: ", seg_file)

if (!identical(sf::st_crs(segments_sf)$wkt, crs_relief)) {
  segments_sf <- sf::st_transform(segments_sf, crs_relief)
}

# Invalid geometries can break extraction or produce unstable results.
# st_make_valid() repairs topology without changing the conceptual unit.
if (any(!sf::st_is_valid(segments_sf))) {
  message("Fixing invalid segment geometries with st_make_valid() ...")
  segments_sf <- sf::st_make_valid(segments_sf)
}

# Hydro stack must match the same CRS; this is a pipeline invariant.
crs_hydro <- terra::crs(hydro, proj = TRUE)
if (!identical(crs_hydro, crs_relief)) {
  stop("CRS mismatch: relief stack vs hydro stack. Fix in S2 export.")
}

# ---------------------------------------------------------
# 3) Helper functions (why they exist and what they do)
# ---------------------------------------------------------

# 3.1 Quantile column resolver
# exactextractr returns quantiles with version-dependent column names
# (e.g., 'quantile_25' vs 'quantile_0.25'). This helper locates the correct
# column for a requested quantile in a robust way.
find_quantile_col <- function(colnames_vec, q) {
  pct <- as.integer(round(q * 100))
  cand <- grep(paste0("^quantile.*(", q, "|", pct, ")$"), colnames_vec, value = TRUE)
  if (length(cand) == 1) return(cand[[1]])
  cand <- grep(paste0("^quantile.*", pct, "$"), colnames_vec, value = TRUE)
  if (length(cand) == 1) return(cand[[1]])
  NA_character_
}

# 3.2 Continuous-band signature extraction for a full stack
# Iterates over all bands of a raster stack. For each band, exactextractr
# computes summary statistics within each polygon (area-weighted at boundaries).
#
# The output is a segment-level table (one row per segment) with stable column
# naming: <prefix>__<band>__<stat>.
extract_signature_continuous_allbands <- function(r, polys_sf, prefix) {
  bn <- names(r)
  if (length(bn) == 0) stop("Raster has no band names for prefix: ", prefix)
  
  stats_fun <- c("mean", "stdev", "min", "max", "coefficient_of_variation", "quantile")
  qs <- c(0.25, 0.50, 0.75)
  
  out <- tibble::tibble(segment_id = polys_sf$segment_id)
  
  for (b in bn) {
    x <- exactextractr::exact_extract(r[[b]], polys_sf, fun = stats_fun, quantiles = qs)
    
    # These columns are expected for continuous extraction.
    nms <- names(x)
    need <- c("mean", "stdev", "min", "max", "coefficient_of_variation")
    if (!all(need %in% nms)) {
      stop("Missing expected stats from exact_extract for band '", b, "': ",
           paste(setdiff(need, nms), collapse = ", "))
    }
    
    # Quantiles are optional if the version does not provide named columns.
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

# 3.3 Mode + purity for categorical / ordinal rasters
# For categorical IDs (watershed) and ordinal classes (strahler), the mean is
# not meaningful. Instead, majority voting is used, weighted by the fraction of
# each raster cell covered by the polygon (coverage_fraction).
#
# mode       = category value with the largest area-weighted support
# mode_share = that support divided by total covered area -> a "purity" measure
extract_mode_share <- function(values, coverage_fraction) {
  ok <- is.finite(values) & is.finite(coverage_fraction)
  v <- values[ok]; w <- coverage_fraction[ok]
  v <- v[!is.na(v)]; w <- w[!is.na(w)]
  
  if (length(v) == 0) return(c(mode = NA_real_, mode_share = NA_real_))
  
  tab <- tapply(w, v, sum)
  mode_val <- as.numeric(names(tab)[which.max(tab)])
  mode_w   <- as.numeric(max(tab))
  total_w  <- as.numeric(sum(tab))
  
  c(mode = mode_val, mode_share = if (total_w > 0) mode_w / total_w else NA_real_)
}

# 3.4 Robust normalization for custom exact_extract outputs
# exactextractr custom functions can return different container types across
# versions/configuration (data.frame, matrix, list-of-vectors). Downstream code
# expects a data.frame with named columns (e.g., 'mode', 'mode_share').
#
# This helper converts the returned object into a predictable data.frame. If
# that is not possible, it stops with a diagnostic error rather than producing
# silent wrong indexing.
normalize_custom_extract <- function(x, expected) {
  
  # already a data.frame with correct names
  if (is.data.frame(x) && all(expected %in% names(x))) return(x)
  
  # matrix / array
  if (is.matrix(x) || (is.array(x) && length(dim(x)) == 2)) {
    
    d <- dim(x)
    dn <- dimnames(x)
    
    # Case 1: fields are rows (your case) -> transpose
    if (!is.null(dn) && length(dn) >= 1 && !is.null(dn[[1]]) && all(expected %in% dn[[1]])) {
      xt <- t(x)
      colnames(xt) <- dn[[1]]  # use existing rownames
      # ensure expected column order
      xt <- xt[, expected, drop = FALSE]
      return(as.data.frame(xt))
    }
    
    # Case 2: fields are columns -> just set colnames (or reorder)
    if (!is.null(dn) && length(dn) >= 2 && !is.null(dn[[2]]) && all(expected %in% dn[[2]])) {
      colnames(x) <- dn[[2]]
      x <- x[, expected, drop = FALSE]
      return(as.data.frame(x))
    }
    
    # Case 3: no dimnames, but shape matches expected in columns
    if (d[2] == length(expected)) {
      colnames(x) <- expected
      return(as.data.frame(x))
    }
    
    # Case 4: no dimnames, but shape matches expected in rows -> transpose
    if (d[1] == length(expected)) {
      xt <- t(x)
      colnames(xt) <- expected
      return(as.data.frame(xt))
    }
    
    stop("Custom exact_extract output matrix has unexpected shape: ", paste(d, collapse = "x"))
  }
  
  # list-of-vectors fallback
  if (is.list(x) && !is.data.frame(x)) {
    df <- dplyr::bind_rows(lapply(x, function(v) {
      v <- as.numeric(v)
      names(v) <- expected
      as.data.frame(as.list(v))
    }))
    return(df)
  }
  
  stop("Custom exact_extract output could not be normalized. Got type: ", paste(class(x), collapse = "/"))
}



# ---------------------------------------------------------
# 4) Relief signatures: all bands treated as continuous
# ---------------------------------------------------------
relief_sig <- extract_signature_continuous_allbands(relief, segments_sf, prefix = "relief")

# ---------------------------------------------------------
# 5) Hydro signatures: mixed handling by band semantics
# ---------------------------------------------------------
bn_h <- names(hydro)
if (length(bn_h) == 0) stop("hydro_stack_10m has no band names")

hydro_sig <- tibble::tibble(segment_id = segments_sf$segment_id)

for (b in bn_h) {
  
  # watershed: categorical basin ID -> mode + mode_share
  if (b == "watershed") {
    
    raw <- exactextractr::exact_extract(
      hydro[[b]],
      segments_sf,
      fun = function(values, coverage_fraction) {
        extract_mode_share(values, coverage_fraction)
      }
    )
    
    mode_df <- normalize_custom_extract(raw, expected = c("mode", "mode_share"))
    
    hydro_sig[[paste0("hydro__", b, "__mode")]]       <- as.numeric(mode_df[["mode"]])
    hydro_sig[[paste0("hydro__", b, "__mode_share")]] <- as.numeric(mode_df[["mode_share"]])
    next
  }
  
  # strahler: ordinal stream order -> mode + mode_share + max order present
  if (b == "strahler") {
    
    raw <- exactextractr::exact_extract(
      hydro[[b]],
      segments_sf,
      fun = function(values, coverage_fraction) {
        ms <- extract_mode_share(values, coverage_fraction)
        ok <- is.finite(values)
        vmax <- if (any(ok)) suppressWarnings(max(values[ok], na.rm = TRUE)) else NA_real_
        c(mode = ms[["mode"]], mode_share = ms[["mode_share"]], max = vmax)
      }
    )
    
    mode_df <- normalize_custom_extract(raw, expected = c("mode", "mode_share", "max"))
    
    hydro_sig[[paste0("hydro__", b, "__mode")]]       <- as.numeric(mode_df[["mode"]])
    hydro_sig[[paste0("hydro__", b, "__mode_share")]] <- as.numeric(mode_df[["mode_share"]])
    hydro_sig[[paste0("hydro__", b, "__max")]]        <- as.numeric(mode_df[["max"]])
    next
  }
  
  # default: continuous hydro layers -> full signature extraction
  tmp <- extract_signature_continuous_allbands(hydro[[b]], segments_sf, prefix = "hydro")
  keep_cols <- setdiff(names(tmp), "segment_id")
  for (cc in keep_cols) hydro_sig[[cc]] <- tmp[[cc]]
}

# ---------------------------------------------------------
# 6) Optional windwardness (requires 'aspect' band)
# ---------------------------------------------------------
# This is computed only if both a seasonal wind-direction summary exists and
# the relief stack contains an 'aspect' band. The computation is a deterministic
# transformation of an existing S2 variable (aspect) with an external parameter
# (mean wind direction) and produces a segment-level mean.
has_wind <- ("wind_means_summary" %in% names(paths)) && file.exists(paths[["wind_means_summary"]])
wind_df <- NULL

if (has_wind) {
  wind_means_summary <- readRDS(paths[["wind_means_summary"]])
  
  if (all(c("period","mean_wd") %in% names(wind_means_summary)) && "aspect" %in% names(relief)) {
    
    wd_summer <- wind_means_summary$mean_wd[wind_means_summary$period == "Summer"]
    wd_winter <- wind_means_summary$mean_wd[wind_means_summary$period == "Winter"]
    
    if (length(wd_summer) == 1L && length(wd_winter) == 1L) {
      
      aspect_deg <- relief[["aspect"]]
      windward_summer_r <- cos((aspect_deg - wd_summer) * pi / 180)
      windward_winter_r <- cos((aspect_deg - wd_winter) * pi / 180)
      
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

# ---------------------------------------------------------
# 7) Join signatures to geometry and write output
# ---------------------------------------------------------
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
