# =====================================================================
# Helper functions for 03-0_analysis_base-segmentation_otb.R
# Place into: src/r-libs/metrics-fun.R
# =====================================================================

# ---- Build a safe BandMathX z-score expression for an N-band raster ----
# Returns a shell-quoted expression string that OTB BandMathX can execute.
build_bandmathx_zscore_expr <- function(nb) {
  stopifnot(length(nb) == 1, nb >= 1)
  
  expr_vec <- character(nb)
  for (i in seq_len(nb)) {
    expr_vec[i] <- paste0(
      "(im1b", i, " - im1b", i, "Mean) / sqrt(im1b", i, "Var)"
    )
  }
  
  expr <- paste(expr_vec, collapse = ";")
  
  # Critical: must be shell-quoted because it contains ';' and '()'
  shQuote(expr)
}


# ---- Run OTB BandMathX to create a z-score standardized multi-band stack ----
# Self-describing API, validated keys. Writes to out_file.
run_otb_bandmathx_zscore <- function(otb, in_file, out_file, ram = 256, quiet = FALSE) {
  # Determine band count (metadata only; terra does not load full raster)
  r  <- terra::rast(in_file)
  nb <- terra::nlyr(r)
  
  expr <- build_bandmathx_zscore_expr(nb)
  
  cmd <- link2GI::otb_build_cmd(
    "BandMathX",
    otb,
    include_optional = "defaults",
    require_output   = TRUE
  )
  
  cmd[["il"]]  <- in_file
  cmd[["exp"]] <- expr
  
  if (!is.null(ram)) cmd[["ram"]] <- as.character(ram)
  
  cmd <- link2GI::otb_set_out(cmd, otb, key = "out", path = out_file)
  
  link2GI::runOTB(cmd, otb, quiet = quiet)
  invisible(out_file)
}


# ---- Run OTB DimensionalityReduction (PCA) to create a feature stack ----
run_otb_pca <- function(otb, in_file, out_file, ncomp = 6L, ram = 256, quiet = FALSE) {
  stopifnot(ncomp >= 1)
  
  cmd <- link2GI::otb_build_cmd(
    "DimensionalityReduction",
    otb,
    include_optional = "defaults",
    require_output   = TRUE
  )
  
  cmd[["in"]]     <- in_file
  cmd[["method"]] <- "pca"
  cmd[["nbcomp"]] <- as.character(ncomp)
  
  if (!is.null(ram)) cmd[["ram"]] <- as.character(ram)
  
  cmd <- link2GI::otb_set_out(cmd, otb, key = "out", path = out_file)
  
  link2GI::runOTB(cmd, otb, quiet = quiet)
  invisible(out_file)
}


# ---- Run OTB LargeScaleMeanShift in raster mode (label raster) ----
# Returns the written label raster path.
run_otb_meanshift_labels <- function(otb, in_raster, out_label_raster,
                                     spatialr, ranger, minsize,
                                     tilesize = 256L,
                                     ram = 256,
                                     quiet = TRUE) {
  
  cmd <- link2GI::otb_build_cmd(
    "LargeScaleMeanShift",
    otb,
    include_optional = "defaults",
    require_output   = TRUE
  )
  
  cmd[["in"]]       <- in_raster
  cmd[["mode"]]     <- "raster"
  cmd[["spatialr"]] <- as.character(spatialr)
  cmd[["ranger"]]   <- as.character(ranger)
  cmd[["minsize"]]  <- as.character(minsize)
  
  # raster-output key in OTB app
  cmd <- link2GI::otb_set_out(cmd, otb, key = "mode.raster.out", path = out_label_raster)
  
  cmd[["tilesizex"]] <- as.character(tilesize)
  cmd[["tilesizey"]] <- as.character(tilesize)
  
  if (!is.null(ram)) cmd[["ram"]] <- as.character(ram)
  
  link2GI::runOTB(cmd, otb, quiet = quiet)
  
  stopifnot(file.exists(out_label_raster))
  invisible(out_label_raster)
}


# ---- Run OTB LargeScaleMeanShift in vector mode (final polygons) ----
run_otb_meanshift_vector <- function(otb, in_raster, out_vector,
                                     spatialr, ranger, minsize,
                                     tilesize = 256L,
                                     ram = 256,
                                     quiet = FALSE) {
  
  cmd <- link2GI::otb_build_cmd(
    "LargeScaleMeanShift",
    otb,
    include_optional = "defaults",
    require_output   = TRUE
  )
  
  cmd[["in"]]       <- in_raster
  cmd[["mode"]]     <- "vector"
  cmd[["spatialr"]] <- as.character(spatialr)
  cmd[["ranger"]]   <- as.character(ranger)
  cmd[["minsize"]]  <- as.character(minsize)
  
  # vector-output key in OTB app
  cmd <- link2GI::otb_set_out(cmd, otb, key = "mode.vector.out", path = out_vector)
  
  cmd[["tilesizex"]] <- as.character(tilesize)
  cmd[["tilesizey"]] <- as.character(tilesize)
  
  if (!is.null(ram)) cmd[["ram"]] <- as.character(ram)
  
  link2GI::runOTB(cmd, otb, quiet = quiet)
  
  stopifnot(file.exists(out_vector))
  invisible(out_vector)
}


# ---- Compute label-size statistics for a label raster (pixels per label) ----
label_size_stats <- function(label_raster) {
  lab <- if (inherits(label_raster, "SpatRaster")) label_raster else terra::rast(label_raster)
  
  fr <- terra::freq(lab)
  fr <- fr[!is.na(fr[, 1]), , drop = FALSE]
  if (is.null(fr) || nrow(fr) == 0) {
    return(list(freq = data.frame(label = integer(0), n_px = integer(0)),
                quantile = NA, shares = NA))
  }
  
  names(fr) <- c("label", "n_px")
  n_px <- fr$n_px
  
  list(
    freq = fr,
    summary = summary(n_px),
    quantile = quantile(n_px, probs = c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1), na.rm = TRUE),
    shares = c(
      lt2  = mean(n_px < 2,  na.rm = TRUE),
      lt5  = mean(n_px < 5,  na.rm = TRUE),
      lt10 = mean(n_px < 10, na.rm = TRUE),
      lt20 = mean(n_px < 20, na.rm = TRUE)
    )
  )
}


# ---- Choose min_px threshold from stats (explicit, audit-friendly) ----
choose_min_px <- function(stats,
                          onepx_threshold = 0.05,
                          min_px_if_triggered = 5L,
                          min_px_else = 0L) {
  
  if (is.null(stats$freq) || nrow(stats$freq) == 0) return(min_px_else)
  
  onepx_share <- mean(stats$freq$n_px <= 1.5, na.rm = TRUE)
  
  if (onepx_share >= onepx_threshold) {
    min_px_if_triggered
  } else {
    min_px_else
  }
}


# ---- Sieve / absorb small labels (<min_px) into dominant neighboring label ----
# Pure terra-based cleanup. Works on integer label rasters.
sieve_labels <- function(label_raster, min_px = 5L, max_iter = 50L) {
  lab <- if (inherits(label_raster, "SpatRaster")) label_raster else terra::rast(label_raster)
  lab <- terra::as.int(lab)
  
  for (iter in seq_len(max_iter)) {
    fr <- terra::freq(lab)
    fr <- fr[!is.na(fr[, 1]), , drop = FALSE]
    if (is.null(fr) || nrow(fr) == 0) break
    names(fr) <- c("label", "n_px")
    
    small <- fr$label[fr$n_px < min_px]
    if (length(small) == 0) break
    
    for (lbl in small) {
      m <- (lab == lbl)
      
      # fast skip if label disappeared
      s <- terra::global(m, "sum", na.rm = TRUE)[1, 1]
      if (is.na(s) || s == 0) next
      
      # dilate mask by 1 pixel (8-neighborhood)
      dil <- terra::focal(m, w = matrix(1, 3, 3), fun = max, na.rm = TRUE, fillvalue = 0)
      edge <- (dil == 1) & (m == 0)
      
      neigh_vals <- terra::values(lab, mat = FALSE)[terra::values(edge, mat = FALSE) == 1]
      neigh_vals <- neigh_vals[!is.na(neigh_vals) & neigh_vals != lbl]
      
      if (length(neigh_vals) == 0) next
      
      to_lbl <- as.integer(names(which.max(table(neigh_vals))))
      
      v <- terra::values(lab, mat = FALSE)
      v[v == lbl] <- to_lbl
      lab <- terra::setValues(lab, v)
    }
  }
  
  lab
}

# ---- Adjusted Rand Index (memory-safe implementation) ------------------------
# Computes ARI without constructing a full contingency table (table(x,y)),
# which would explode for many unique labels.
adjusted_rand_index <- function(x, y) {
  
  # remove NA pairs
  ok <- !is.na(x) & !is.na(y)
  x  <- x[ok]
  y  <- y[ok]
  
  n <- length(x)
  if (n < 2) return(NA_real_)
  
  # relabel to compact integer levels
  xf <- factor(x)
  yf <- factor(y)
  
  comb2 <- function(v) v * (v - 1) / 2
  
  # joint counts for observed (x,y) pairs only
  ij  <- interaction(xf, yf, drop = TRUE)
  nij <- tabulate(as.integer(ij))
  
  # marginal counts
  ai <- tabulate(as.integer(xf))
  bj <- tabulate(as.integer(yf))
  
  sum_nij <- sum(comb2(nij))
  sum_ai  <- sum(comb2(ai))
  sum_bj  <- sum(comb2(bj))
  n2      <- comb2(n)
  
  if (n2 == 0) return(NA_real_)
  
  expected  <- (sum_ai * sum_bj) / n2
  max_index <- 0.5 * (sum_ai + sum_bj)
  denom     <- max_index - expected
  
  if (denom == 0) return(NA_real_)
  
  (sum_nij - expected) / denom
}
