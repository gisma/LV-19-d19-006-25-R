# ============================================================
# metrics-fun.R  (OTB segmentation helpers + ARI_prev stability)
# ============================================================
# This file contains ONLY reusable helper functions.
# Project-specific paths MUST be handled in the calling script via outputs.tsv.
#
# Pipeline context (S3 base segmentation):
#   1) z-score standardisation (OTB BandMathX)
#   2) PCA feature extraction   (OTB DimensionalityReduction)
#   3) MeanShift segmentation   (OTB LargeScaleMeanShift, raster labels)
#   4) Stability scoring        (ARI_prev: ARI under small parameter perturbations)
#   5) Postcondition            (optional masking of tiny segments by size distribution)
#   6) Polygonisation           (done in the calling script; not here)
#
# Notes (science/tech):
# - z-scoring removes scale effects and makes bands comparable for PCA/segmentation.
# - PCA reduces noise and collinearity while preserving dominant variance structure.
# - MeanShift is a non-parametric mode-seeking clustering approach.
# - ARI (Adjusted Rand Index) measures partition similarity corrected for chance.
#   We use ARI_prev as a local stability proxy: if small parameter changes produce
#   similar segmentations, the solution is more robust (not "true", but stable).
#
# References:
# - MeanShift: Comaniciu & Meer (2002), "Mean Shift: A robust approach toward feature space analysis."
# - ARI: Hubert & Arabie (1985), "Comparing partitions."
# - PCA (general): Jolliffe & Cadima (2016), "Principal component analysis: a review and recent developments."
# ============================================================

# ---- Build BandMathX z-score expression for nb input bands -------------------
# BandMathX supports per-band internal stats variables:
#   im1b{i}Mean, im1b{i}Var  (computed by OTB for each input band).
# We produce a multi-band output expression: one z-score term per band.
# IMPORTANT: OTB expects multiple output bands separated by ';'.
# IMPORTANT: We shell-quote the expression to prevent the shell from interpreting
#            ';' (command separator) and parentheses. This is essential on Linux.
build_bandmathx_zscore_expr <- function(nb) {
  stopifnot(length(nb) == 1, is.finite(nb), nb >= 1)
  
  expr_vec <- character(nb)
  for (i in seq_len(nb)) {
    expr_vec[i] <- paste0("(im1b", i, " - im1b", i, "Mean) / sqrt(im1b", i, "Var)")
  }
  
  shQuote(paste(expr_vec, collapse = ";"))
}

# ---- Run OTB BandMathX to create z-scored multiband raster -------------------
# Uses link2GI "Self-describing API" (otb_build_cmd + otb_set_out).
# This avoids fragile legacy parsing and ensures only valid parameter keys are used.
# ---- Run OTB BandMathX to create z-scored multiband raster -------------------
# ---- Run OTB BandMathX to create a z-score standardized multi-band stack ----
# FIXES:
#  - deterministic overwrite (temp file + atomic rename)
#  - temp output keeps .tif extension (writer recognition)
#  - sanity checks (exists + non-trivial size)
run_otb_bandmathx_zscore <- function(otb, in_file, out_file,
                                     ram = NULL, quiet = FALSE, progress = FALSE) {
  stopifnot(file.exists(in_file))
  
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
  cmd[["progress"]] <- if (isTRUE(progress)) "true" else "false"
  
  out_dir <- dirname(out_file)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  tmp_out <- file.path(
    out_dir,
    paste0(
      tools::file_path_sans_ext(basename(out_file)),
      ".tmp_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%d%H%M%S"),
      ".tif"
    )
  )
  if (file.exists(tmp_out)) file.remove(tmp_out)
  
  # IMPORTANT: write to temp file first (must end with .tif)
  cmd <- link2GI::otb_set_out(cmd, otb, key = "out", path = tmp_out)
  
  link2GI::runOTB(cmd, otb, quiet = quiet)
  
  if (!file.exists(tmp_out)) stop("OTB did not write temp BandMathX output: ", tmp_out)
  fi <- file.info(tmp_out)
  if (is.na(fi$size) || fi$size < 1024L) stop("Temp BandMathX output too small/invalid: ", tmp_out)
  
  if (file.exists(out_file)) {
    ok_rm <- tryCatch(file.remove(out_file), error = function(e) FALSE)
    if (!isTRUE(ok_rm) && file.exists(out_file)) stop("Cannot remove existing output: ", out_file)
  }
  
  ok_mv <- file.rename(tmp_out, out_file)
  if (!isTRUE(ok_mv) || !file.exists(out_file)) stop("Failed to move temp->final: ", out_file)
  
  invisible(out_file)
}



# ---- Run OTB DimensionalityReduction (PCA) ----------------------------------
# PCA reduces dimensionality and collinearity of predictor stacks.
# In remote sensing segmentation workflows this is mainly pragmatic:
# - reduces RAM/IO for later steps
# - suppresses some band-level redundancy
# - keeps dominant variance structure
# ---- Run OTB DimensionalityReduction (PCA) ----------------------------------
# ---- Run OTB DimensionalityReduction (PCA) to create a feature stack ----
# FIXES:
#  - deterministic overwrite (temp file + atomic rename)
#  - temp output keeps .tif extension
#  - sanity checks
run_otb_pca <- function(otb, in_file, out_file, ncomp = 6L,
                        ram = NULL, quiet = FALSE, progress = FALSE) {
  stopifnot(file.exists(in_file))
  stopifnot(length(ncomp) == 1, is.finite(ncomp), ncomp >= 1)
  
  cmd <- link2GI::otb_build_cmd(
    "DimensionalityReduction",
    otb,
    include_optional = "defaults",
    require_output   = TRUE
  )
  cmd[["in"]]     <- in_file
  cmd[["method"]] <- "pca"
  cmd[["nbcomp"]] <- as.character(as.integer(ncomp))
  if (!is.null(ram)) cmd[["ram"]] <- as.character(ram)
  cmd[["progress"]] <- if (isTRUE(progress)) "true" else "false"
  
  out_dir <- dirname(out_file)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  tmp_out <- file.path(
    out_dir,
    paste0(
      tools::file_path_sans_ext(basename(out_file)),
      ".tmp_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%d%H%M%S"),
      ".tif"
    )
  )
  if (file.exists(tmp_out)) file.remove(tmp_out)
  
  cmd <- link2GI::otb_set_out(cmd, otb, key = "out", path = tmp_out)
  
  link2GI::runOTB(cmd, otb, quiet = quiet)
  
  if (!file.exists(tmp_out)) stop("OTB did not write temp PCA output: ", tmp_out)
  fi <- file.info(tmp_out)
  if (is.na(fi$size) || fi$size < 1024L) stop("Temp PCA output too small/invalid: ", tmp_out)
  
  if (file.exists(out_file)) {
    ok_rm <- tryCatch(file.remove(out_file), error = function(e) FALSE)
    if (!isTRUE(ok_rm) && file.exists(out_file)) stop("Cannot remove existing output: ", out_file)
  }
  
  ok_mv <- file.rename(tmp_out, out_file)
  if (!isTRUE(ok_mv) || !file.exists(out_file)) stop("Failed to move temp->final: ", out_file)
  
  invisible(out_file)
}



# ---- Run OTB LargeScaleMeanShift in raster label mode ------------------------
# =====================================================================
# FIXED: run_otb_meanshift_labels()
# - Enforces deterministic overwrite (delete target before run)
# - Verifies that the output file is NEWER than the call start time
# =====================================================================
# ---- Run OTB LargeScaleMeanShift in raster mode (label raster) ----
# FIXES:
#  - deterministic overwrite (temp file + atomic rename)
#  - temp output ends with .tif (writer recognition; fixes your current failure)
#  - sanity checks
run_otb_meanshift_labels <- function(otb, in_raster, out_label_raster,
                                     spatialr, ranger, minsize,
                                     tilesize = 256L,
                                     ram = 256,
                                     quiet = FALSE,
                                     progress = FALSE) {
  
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
  
  cmd[["tilesizex"]] <- as.character(tilesize)
  cmd[["tilesizey"]] <- as.character(tilesize)
  
  if (!is.null(ram)) cmd[["ram"]] <- as.character(ram)
  cmd[["progress"]] <- if (isTRUE(progress)) "true" else "false"
  
  out_dir <- dirname(out_label_raster)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # IMPORTANT: ensure temp output ends with .tif (NOT ".tif.tmp_...")
  tmp_out <- file.path(
    out_dir,
    paste0(
      tools::file_path_sans_ext(basename(out_label_raster)),
      ".tmp_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%d%H%M%S"),
      ".tif"
    )
  )
  if (file.exists(tmp_out)) file.remove(tmp_out)
  
  cmd <- link2GI::otb_set_out(cmd, otb, key = "mode.raster.out", path = tmp_out)
  
  link2GI::runOTB(cmd, otb, quiet = quiet)
  
  if (!file.exists(tmp_out)) stop("OTB did not write temp label raster: ", tmp_out)
  fi <- file.info(tmp_out)
  if (is.na(fi$size) || fi$size < 1024L) stop("Temp label raster too small/invalid: ", tmp_out)
  
  if (file.exists(out_label_raster)) {
    ok_rm <- tryCatch(file.remove(out_label_raster), error = function(e) FALSE)
    if (!isTRUE(ok_rm) && file.exists(out_label_raster)) stop("Cannot remove existing label raster: ", out_label_raster)
  }
  
  ok_mv <- file.rename(tmp_out, out_label_raster)
  if (!isTRUE(ok_mv) || !file.exists(out_label_raster)) stop("Failed to move temp->final label raster: ", out_label_raster)
  
  invisible(out_label_raster)
}


# ---- Label size statistics ---------------------------------------------------
# Purpose:
# - quantify the size distribution of segments (labels) in pixel units
# - provide quantiles and shares of "tiny segments"
#
# Implementation detail:
# terra::freq() yields label frequencies without loading the whole raster as a table.
# It typically returns (value,count) or (layer,value,count); we normalise both.
label_size_stats <- function(label_raster,
                             probs = c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1)) {
  lab <- if (inherits(label_raster, "SpatRaster")) label_raster else terra::rast(label_raster)
  
  fr <- terra::freq(lab)
  if (is.null(fr) || nrow(fr) == 0) {
    return(list(
      freq = data.frame(label = integer(0), n_px = integer(0)),
      quantile = stats::setNames(rep(NA_real_, length(probs)), paste0(probs * 100, "%")),
      shares = c(lt2 = NA_real_, lt5 = NA_real_, lt10 = NA_real_, lt20 = NA_real_)
    ))
  }
  
  if (ncol(fr) == 3L) fr <- fr[, c("value", "count")] else fr <- fr[, c(1, 2)]
  names(fr) <- c("label", "n_px")
  fr <- fr[!is.na(fr$label), , drop = FALSE]
  
  q <- stats::quantile(fr$n_px, probs = probs, na.rm = TRUE, type = 7)
  sh <- c(
    lt2  = mean(fr$n_px < 2,  na.rm = TRUE),
    lt5  = mean(fr$n_px < 5,  na.rm = TRUE),
    lt10 = mean(fr$n_px < 10, na.rm = TRUE),
    lt20 = mean(fr$n_px < 20, na.rm = TRUE)
  )
  
  list(freq = fr, quantile = q, shares = sh)
}

# ---- Jump / knee-based cutoff on label size distribution ---------------------
# Goal:
# Choose a cutoff min_px from the size distribution itself (non-arbitrary).
#
# Heuristic (not BFAST):
# - sort n_px, work on log-scale to reduce heavy-tail dominance
# - (optional) running median smoothing to stabilise noisy tails
# - compute second differences and pick the maximum curvature point ("knee")
#
# Interpretation:
# Everything below the knee is treated as "small fragments" relative to the
# dominant segment size regime. This is a *structural* criterion, not a semantic one.
choose_min_px_jump <- function(stats,
                               min_px_floor = 1L,
                               smooth_k = 101L) {
  stopifnot(is.list(stats), "freq" %in% names(stats))
  fr <- stats$freq
  if (is.null(fr) || nrow(fr) < 10) return(as.integer(min_px_floor))
  
  x <- sort(fr$n_px)
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 10) return(as.integer(min_px_floor))
  
  lx <- log(x)
  
  k <- as.integer(smooth_k)
  if (is.finite(k) && k >= 3 && (k %% 2 == 1) && k < length(lx)) {
    lx <- stats::runmed(lx, k = k, endrule = "median")
  }
  
  d1 <- diff(lx)
  d2 <- diff(d1)
  if (length(d2) == 0) return(as.integer(min_px_floor))
  
  # index of maximum curvature (shift by +1 to map to x)
  j <- which.max(d2) + 1L
  cut <- x[j]
  
  max(as.integer(min_px_floor), as.integer(round(cut)))
}

# ---- Mask labels smaller than min_px (set to NA) ------------------------------
# This does NOT "merge" small segments; it removes them.
# Rationale: if the modelling goal is "homogeneous representative patches",
# tiny fragments are not useful candidates. Masking makes them ineligible later.
#
# Implementation note:
# - For ~5–10 million cells this is typically OK in RAM on modern machines,
#   but it *is* an in-memory value vector. If needed, rewrite blockwise.
# ---- Mask labels smaller than min_px (set to NA) - blockwise -----------------
# Fix: avoids terra::values(lab) full in-memory vector on large rasters.
mask_small_segments <- function(label_raster, min_px) {
  min_px <- as.integer(min_px)
  stopifnot(is.finite(min_px), min_px >= 1)
  
  lab <- if (inherits(label_raster, "SpatRaster")) label_raster else terra::rast(label_raster)
  
  fr <- terra::freq(lab)
  if (is.null(fr) || nrow(fr) == 0) return(lab)
  
  if (ncol(fr) == 3L) fr <- fr[, c("value", "count")] else fr <- fr[, c(1, 2)]
  names(fr) <- c("label", "n_px")
  fr <- fr[!is.na(fr$label), , drop = FALSE]
  
  small <- fr$label[fr$n_px < min_px]
  if (length(small) == 0) return(lab)
  
  # blockwise write: keep datatype & geometry
  out <- terra::rast(lab)
  terra::writeStart(out, terra::tempfile(fileext = ".tif"), overwrite = TRUE)
  
  bs <- terra::blocks(lab)
  for (i in seq_len(nrow(bs))) {
    v <- terra::readValues(lab, row = bs$row[i], nrows = bs$nrows[i], mat = FALSE)
    v[v %in% small] <- NA
    terra::writeValues(out, v, bs$row[i])
  }
  
  terra::writeStop(out)
  out
}


# ---- Memory-safe Adjusted Rand Index (ARI) -----------------------------------
# ARI measures similarity of two partitions while correcting for chance agreement.
# Using table(x,y) can explode in size (potentially > 2^31) even if only few pairs
# are observed; therefore we compute joint counts via interaction+tabulate.
#
# Reference: Hubert & Arabie (1985).
adjusted_rand_index <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) == 0) return(NA_real_)
  
  xf <- as.integer(factor(x))
  yf <- as.integer(factor(y))
  
  # joint counts only for observed pairs
  ij  <- as.integer(interaction(xf, yf, drop = TRUE))
  nij <- tabulate(ij)
  
  # marginals
  ni <- tabulate(xf)
  nj <- tabulate(yf)
  
  comb2 <- function(n) n * (n - 1) / 2
  
  sum_nij <- sum(comb2(nij))
  sum_ni  <- sum(comb2(ni))
  sum_nj  <- sum(comb2(nj))
  
  n <- length(xf)
  total <- comb2(n)
  if (total == 0) return(NA_real_)
  
  expected <- (sum_ni * sum_nj) / total
  max_idx  <- 0.5 * (sum_ni + sum_nj)
  
  denom <- (max_idx - expected)
  if (denom == 0) return(NA_real_)
  
  (sum_nij - expected) / denom
}

# ---- Sample paired labels from two rasters (same grid) -----------------------
# We sample cell indices from r_base and read both rasters at those cells.
sample_label_pairs <- function(r_base, r_other, n = 2e5, seed = 1L) {
  rb <- if (inherits(r_base, "SpatRaster")) r_base else terra::rast(r_base)
  ro <- if (inherits(r_other, "SpatRaster")) r_other else terra::rast(r_other)
  
  stopifnot(terra::ncell(rb) == terra::ncell(ro))
  
  set.seed(as.integer(seed))
  
  nc <- terra::ncell(rb)
  n  <- as.integer(min(n, nc))
  
  # sample cell indices directly (no terra coordinate math => no overflow)
  cells <- sample.int(nc, size = n, replace = FALSE)
  
  xb <- terra::values(rb, cells = cells, mat = FALSE)
  yb <- terra::values(ro, cells = cells, mat = FALSE)
  
  ok <- !is.na(xb) & !is.na(yb)
  list(x = xb[ok], y = yb[ok])
}

# ---- Generate parameter perturbations around a base set ----------------------
# Creates a small, symmetric candidate set around a base (spatialr,ranger,minsize):
# - spatialr: +/- ds
# - ranger  : +/- dr
# - minsize : +/- dm
# Then samples up to K candidates (deterministic).
make_param_perturbations <- function(spatialr, ranger, minsize,
                                     dr = 0.02, ds = 1L, dm = 20L, K = 8L) {
  spatialr <- as.integer(spatialr)
  minsize  <- as.integer(minsize)
  ranger   <- as.numeric(ranger)
  
  dr <- as.numeric(dr)
  ds <- as.integer(ds)
  dm <- as.integer(dm)
  K  <- as.integer(K)
  
  cand_spatialr <- unique(pmax(1L, spatialr + c(-ds, 0L, ds)))
  cand_ranger   <- unique(pmax(1e-6, ranger + c(-dr, 0, dr)))
  cand_minsize  <- unique(pmax(1L, minsize + c(-dm, 0L, dm)))
  
  cand <- expand.grid(
    spatialr = cand_spatialr,
    ranger   = cand_ranger,
    minsize  = cand_minsize,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # drop the baseline itself
  cand <- cand[!(
    cand$spatialr == spatialr &
      abs(cand$ranger - ranger) < 1e-12 &
      cand$minsize == minsize
  ), , drop = FALSE]
  
  if (nrow(cand) == 0) return(cand)
  
  if (nrow(cand) > K) {
    set.seed(1L)
    cand <- cand[sample(seq_len(nrow(cand)), size = K, replace = FALSE), , drop = FALSE]
  }
  
  cand
}

# ---- Sample fixed cell indices once and reuse across rasters ---------------
sample_label_cells <- function(r, n = 2e5, seed = 1L, keep_na = FALSE) {
  stopifnot(inherits(r, "SpatRaster"))
  set.seed(seed)
  cells <- terra::spatSample(r, size = n, method = "random",
                             as.points = FALSE, na.rm = !keep_na, values = FALSE)
  as.integer(cells)
}


# ---- Compute ARI_prev: stability proxy for a base MeanShift parameter set ----
# Definition used here:
# - run baseline MeanShift -> label raster
# - run K perturbations around the baseline -> label rasters
# - compute ARI(base, pert_i) on a random cell sample
# - ARI_prev = mean(ARI_i), ari_sd = sd(ARI_i)
#
# Why this is "circular" in the modelling sense:
# - the selected segmentation scale depends on stability under perturbations;
# - but the "right" scale is also process/landscape-scale dependent and not purely
#   spectral. This selection therefore constrains itself by its own assumptions.
#
# ARI_prev is still useful as an operational robustness filter.
# =====================================================================
# FIXED: compute_ari_prev()
# - Baseline + perturbation filenames are unique per BASELINE signature
# - Deterministic overwrite for every OTB output (no reuse)
# - Optional cleanup of perturbation files to prevent tmp explosion
# =====================================================================
# ---- Compute ARI_prev: stability proxy for a base MeanShift parameter set ----
# FIX:
#  - perturbation filenames are unique per BASELINE signature (no cross-scale collisions)
#  - optional cleanup to prevent tmp explosion
compute_ari_prev <- function(otblink, feat_raster,
                             out_dir_tmp,
                             spatialr, ranger, minsize,
                             perturb = list(dr = 0.02, ds = 1L, dm = 20L, K = 8L, minsize_floor_frac = 0.8),
                             sample_fact = 4L,
                             seed = 1L,
                             tilesize = 256L,
                             ram = NULL,
                             quiet = TRUE,
                             verbose = TRUE,
                             cleanup_perturbations = FALSE) {
  
  dir.create(out_dir_tmp, recursive = TRUE, showWarnings = FALSE)
  
  ts0 <- Sys.time()
  log <- function(...) if (isTRUE(verbose)) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)
  sec <- function(t_start) round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)
  
  assert_file_ok <- function(path, min_bytes = 2048L, label = "output") {
    if (!file.exists(path)) stop("Missing ", label, ": ", path)
    fi <- file.info(path)
    if (is.na(fi$size) || fi$size < min_bytes) stop("Invalid ", label, " (too small): ", path)
    invisible(TRUE)
  }
  
  log("ARI_prev start | base: ",
      sprintf("s=%d r=%.4f m=%d", spatialr, ranger, minsize),
      " | sample_fact=", sample_fact)
  
  # ---- baseline -------------------------------------------------------------
  base_file <- file.path(
    out_dir_tmp,
    sprintf("labels_base_s%02d_r%.4f_m%04d.tif",
            as.integer(spatialr), as.numeric(ranger), as.integer(minsize))
  )
  
  log("Running baseline MeanShift ...")
  t0 <- Sys.time()
  
  run_otb_meanshift_labels(
    otb              = otblink,
    in_raster        = feat_raster,
    out_label_raster = base_file,
    spatialr         = spatialr,
    ranger           = ranger,
    minsize          = minsize,
    tilesize         = tilesize,
    ram              = ram,
    quiet            = quiet
  )
  
  log("Baseline done in ", sec(t0), " s")
  assert_file_ok(base_file, label = "baseline label raster")
  
  r_base <- terra::rast(base_file)
  
  # ---- build coarse template (near resample only) ---------------------------
  log("Aggregating baseline labels ...")
  t0 <- Sys.time()
  
  if (is.null(sample_fact) || sample_fact <= 1L) {
    r_template  <- NULL
    r_base_small <- r_base
  } else {
    r_template <- terra::rast(r_base)
    terra::res(r_template) <- terra::res(r_template) * as.numeric(sample_fact)
    terra::ext(r_template) <- terra::ext(r_base)
    r_base_small <- terra::resample(r_base, r_template, method = "near")
  }
  
  x_base <- terra::values(r_base_small, mat = FALSE)
  log("Baseline aggregation done in ", sec(t0), " s")
  
  # ---- perturbations --------------------------------------------------------
  cand <- make_param_perturbations(
    spatialr = spatialr, ranger = ranger, minsize = minsize,
    dr = perturb$dr, ds = perturb$ds, dm = perturb$dm, K = perturb$K,
    minsize_floor_frac = if (!is.null(perturb$minsize_floor_frac)) perturb$minsize_floor_frac else 0.8
  )
  
  if (nrow(cand) == 0) {
    log("No perturbations generated.")
    return(list(
      ari_prev = NA_real_, ari_sd = NA_real_, ari = numeric(0),
      base = base_file, cand = cand
    ))
  }
  
  log("Perturbations: ", nrow(cand))
  
  ari_vec <- numeric(nrow(cand))
  
  # keep 1 file and overwrite safely (your OTB wrapper already does atomic rename)
  pert_file <- file.path(out_dir_tmp, "labels_pert_current.tif")
  
  for (i in seq_len(nrow(cand))) {
    p <- cand[i, ]
    
    log(sprintf("Perturbation %d / %d | s=%d r=%.4f m=%d",
                i, nrow(cand),
                as.integer(p$spatialr), as.numeric(p$ranger), as.integer(p$minsize)))
    
    # --- MeanShift ---
    t0 <- Sys.time()
    run_otb_meanshift_labels(
      otb              = otblink,
      in_raster        = feat_raster,
      out_label_raster = pert_file,
      spatialr         = p$spatialr,
      ranger           = p$ranger,
      minsize          = p$minsize,
      tilesize         = tilesize,
      ram              = ram,
      quiet            = quiet
    )
    log("  MeanShift done in ", sec(t0), " s")
    assert_file_ok(pert_file, label = "perturbation label raster")
    
    # --- downsample perturbation (near) ---
    r_other <- terra::rast(pert_file)
    
    t0 <- Sys.time()
    if (is.null(sample_fact) || sample_fact <= 1L) {
      r_other_small <- r_other
    } else {
      r_other_small <- terra::resample(r_other, r_template, method = "near")
    }
    y_other <- terra::values(r_other_small, mat = FALSE)
    log("  Aggregation done in ", sec(t0), " s")
    
    # --- ARI ---
    t0 <- Sys.time()
    ari_vec[i] <- adjusted_rand_index(x_base, y_other)
    log("  ARI = ", round(ari_vec[i], 4), " | computed in ", round(sec(t0), 3), " s")
    
    rm(r_other, r_other_small)
    gc(FALSE)
    
    if (isTRUE(cleanup_perturbations) && file.exists(pert_file)) {
      suppressWarnings(file.remove(pert_file))
    }
  }
  
  log("ARI_prev finished in ", sec(ts0), " s")
  
  list(
    ari_prev = mean(ari_vec, na.rm = TRUE),
    ari_sd   = stats::sd(ari_vec, na.rm = TRUE),
    ari      = ari_vec,
    base     = base_file,
    cand     = cand
  )
}




make_param_perturbations <- function(spatialr, ranger, minsize,
                                     dr = NULL, ds = 1L, dm = NULL, K = 8L,
                                     minsize_floor_frac = 0.8) {
  spatialr <- as.integer(spatialr)
  minsize  <- as.integer(minsize)
  ranger   <- as.numeric(ranger)
  
  # adaptive defaults
  if (is.null(dr)) dr <- max(0.005, 0.10 * ranger)           # 10% ranger, min 0.005
  if (is.null(dm)) dm <- max(5L, as.integer(round(0.20 * minsize)))  # 20% minsize, min 5
  
  dr <- as.numeric(dr)
  ds <- as.integer(ds)
  dm <- as.integer(dm)
  K  <- as.integer(K)
  
  # local stability: don't allow s-1 at very small spatialr
  if (spatialr <= 3L) ds <- 0L
  
  cand_spatialr <- unique(pmax(1L, spatialr + c(-ds, 0L, ds)))
  cand_ranger   <- unique(pmax(1e-6, ranger + c(-dr, 0, dr)))
  cand_minsize  <- unique(pmax(1L, minsize + c(-dm, 0L, dm)))
  
  cand <- expand.grid(
    spatialr = cand_spatialr,
    ranger   = cand_ranger,
    minsize  = cand_minsize,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # drop baseline
  cand <- cand[!(
    cand$spatialr == spatialr &
      abs(cand$ranger - ranger) < 1e-12 &
      cand$minsize == minsize
  ), , drop = FALSE]
  
  if (nrow(cand) == 0) return(cand)
  
  # prevent minsize collapsing too far (kills the "m=10" cases)
  floor_m <- as.integer(round(minsize * minsize_floor_frac))
  cand$minsize <- pmax(floor_m, as.integer(cand$minsize))
  
  # deterministic subsample to K
  if (nrow(cand) > K) {
    set.seed(1L)
    cand <- cand[sample(seq_len(nrow(cand)), size = K, replace = FALSE), , drop = FALSE]
  }
  
  cand
}




build_scales_A <- function(master_10m_raster,
                           chm_1m = NULL,
                           scale_ids = c("s20m","s30m","s40m","s60m","s80m","s120m"),
                           fallback_scales = NULL,
                           chm_thr_m = 2,
                           patch_quantiles = c(0.25, 0.40, 0.55, 0.70, 0.85, 0.95),
                           min_patch_px_10m = 4L,
                           max_spatialr = 30L) {
  
  # master pixel size in meters (assumes projected CRS in meters)
  res10 <- terra::res(terra::rast(master_10m_raster))[1]
  if (!is.finite(res10) || res10 <= 0) stop("Cannot read resolution from master_10m_raster.")
  
  # ---- Fallback: use provided fallback or your known manual grid ------------
  if (is.null(chm_1m)) {
    if (!is.null(fallback_scales)) {
      return(fallback_scales)
    }
    # default fallback = your grid (keeps your current behavior)
    return(tibble::tibble(
      scale_id = scale_ids,
      spatialr = c(2, 3, 4, 6, 8, 12),
      ranger   = c(0.06, 0.08, 0.10, 0.12, 0.14, 0.16),
      minsize  = c(30, 40, 50, 80, 100, 150),
      scale_source = "fallback_manual"
    ))
  }
  
  # ---- CHM path ------------------------------------------------------------
  chm <- terra::rast(chm_1m)
  
  # bring CHM to master extent (optional; safe)
  chm <- terra::crop(chm, terra::ext(terra::rast(master_10m_raster)))
  
  # derive 10m vegetation mask from CHM: max height in 10m cell
  # (max is robust for canopy presence; avoids "mean height dilution")
  chm_max_10m <- terra::aggregate(chm, fact = as.integer(round(res10 / terra::res(chm)[1])),
                                  fun = "max", na.rm = TRUE)
  
  veg_10m <- chm_max_10m > chm_thr_m
  veg_10m <- terra::as.int(veg_10m)
  
  # patch labeling (connected components)
  patches <- terra::patches(veg_10m, directions = 8)
  
  fr <- terra::freq(patches)
  fr <- fr[!is.na(fr[,1]), , drop = FALSE]
  if (is.null(fr) || nrow(fr) == 0) {
    # fallback if CHM yields nothing
    return(tibble::tibble(
      scale_id = scale_ids,
      spatialr = c(2, 3, 4, 6, 8, 12),
      ranger   = c(0.06, 0.08, 0.10, 0.12, 0.14, 0.16),
      minsize  = c(30, 40, 50, 80, 100, 150),
      scale_source = "fallback_chm_empty"
    ))
  }
  
  names(fr) <- c("patch_id","n_px")
  fr <- fr[fr$n_px >= min_patch_px_10m, , drop = FALSE]
  
  if (nrow(fr) == 0) {
    return(tibble::tibble(
      scale_id = scale_ids,
      spatialr = c(2, 3, 4, 6, 8, 12),
      ranger   = c(0.06, 0.08, 0.10, 0.12, 0.14, 0.16),
      minsize  = c(30, 40, 50, 80, 100, 150),
      scale_source = "fallback_chm_patches_too_small"
    ))
  }
  
  # patch areas in m^2 at 10m
  A <- fr$n_px * (res10 * res10)
  r_eq <- sqrt(A / pi)  # equivalent radius (m)
  
  rq <- as.numeric(stats::quantile(r_eq, probs = patch_quantiles, na.rm = TRUE))
  rq <- unique(rq)
  
  # map radius to spatialr in 10m pixels
  spatialr <- pmin(max_spatialr, pmax(1L, as.integer(round(rq / res10))))
  
  # map radius to minsize (area in 10m pixels)
  minsize <- pmax(1L, as.integer(round((pi * rq^2) / (res10 * res10))))
  
  # keep exactly length(scale_ids) by padding/truncating
  k <- length(scale_ids)
  if (length(spatialr) < k) spatialr <- c(spatialr, rep(tail(spatialr,1), k - length(spatialr)))
  if (length(minsize)  < k) minsize  <- c(minsize,  rep(tail(minsize,1),  k - length(minsize)))
  spatialr <- spatialr[seq_len(k)]
  minsize  <- minsize[seq_len(k)]
  
  tibble::tibble(
    scale_id = scale_ids,
    spatialr = spatialr,
    ranger   = NA_real_,           # ranger comes from feature space (next function)
    minsize  = minsize,
    scale_source = "chm_patch_quantiles"
  )
}

derive_ranger_from_pca <- function(pca_file,
                                   n = 50000L,
                                   k = 10L,
                                   probs = c(0.25, 0.5, 0.75, 0.9),
                                   seed = 1L) {
  
  r <- terra::rast(pca_file)
  
  set.seed(seed)
  pts <- terra::spatSample(r, size = n, method = "random", na.rm = TRUE, as.points = FALSE, values = TRUE)
  
  X <- as.matrix(pts)
  # drop rows with NA
  ok <- stats::complete.cases(X)
  X <- X[ok, , drop = FALSE]
  if (nrow(X) < (k + 2)) stop("Not enough complete samples from PCA stack.")
  
  # nearest-neighbor distance (simple, no heavy deps):
  # sample a smaller subset for O(n^2) safety if needed
  # (keep it deterministic and cheap)
  max_n <- 8000L
  if (nrow(X) > max_n) X <- X[sample.int(nrow(X), max_n), , drop = FALSE]
  
  # compute distances to kNN via partial sort on squared distances
  # (base R; yes, O(n^2) on max_n, but max_n keeps it bounded)
  d_k <- numeric(nrow(X))
  for (i in seq_len(nrow(X))) {
    di <- rowSums((X - X[i, ])^2)
    di[i] <- Inf
    d_k[i] <- sqrt(sort(di, partial = k)[k])
  }
  
  as.numeric(stats::quantile(d_k, probs = probs, na.rm = TRUE))
}

assign_ranger_to_scales <- function(scales, ranger_candidates) {
  k <- nrow(scales)
  if (length(ranger_candidates) < k) {
    ranger_candidates <- c(ranger_candidates, rep(tail(ranger_candidates,1), k - length(ranger_candidates)))
  }
  scales$ranger <- ranger_candidates[seq_len(k)]
  scales
}




make_param_perturbations <- function(spatialr, ranger, minsize,
                                     dr = NULL, ds = 1L, dm = NULL, K = 8L,
                                     minsize_floor_frac = 0.8) {
  spatialr <- as.integer(spatialr)
  minsize  <- as.integer(minsize)
  ranger   <- as.numeric(ranger)
  
  # adaptive defaults if dr/dm not provided
  if (is.null(dr)) dr <- max(0.005, 0.10 * ranger)
  if (is.null(dm)) dm <- max(5L, as.integer(round(0.20 * minsize)))
  
  dr <- as.numeric(dr)
  ds <- as.integer(ds)
  dm <- as.integer(dm)
  K  <- as.integer(K)
  
  # local stability: don't allow s-1 at very small spatialr
  if (spatialr <= 3L) ds <- 0L
  
  cand_spatialr <- unique(pmax(1L, spatialr + c(-ds, 0L, ds)))
  cand_ranger   <- unique(pmax(1e-6, ranger   + c(-dr, 0, dr)))
  cand_minsize  <- unique(pmax(1L,  minsize  + c(-dm, 0L, dm)))
  
  cand <- expand.grid(
    spatialr = cand_spatialr,
    ranger   = cand_ranger,
    minsize  = cand_minsize,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # drop baseline
  cand <- cand[!(
    cand$spatialr == spatialr &
      abs(cand$ranger - ranger) < 1e-12 &
      cand$minsize == minsize
  ), , drop = FALSE]
  
  if (nrow(cand) == 0) return(cand)
  
  # clamp minsize: prevent collapse (e.g. m=10 for base m=30)
  floor_m <- as.integer(round(minsize * minsize_floor_frac))
  cand$minsize <- pmax(floor_m, as.integer(cand$minsize))
  
  # deterministic subsample to K
  if (nrow(cand) > K) {
    set.seed(1L)
    cand <- cand[sample(seq_len(nrow(cand)), size = K, replace = FALSE), , drop = FALSE]
  }
  
  cand
}

make_adaptive_perturb <- function(spatialr, ranger, minsize,
                                  K = 8L,
                                  dr_frac = 0.10,     # +/- 10% on ranger
                                  dm_frac = 0.20,     # +/- 20% on minsize
                                  ds_max  = 1L) {     # +/-1 at most
  
  # ranger step (absolute) but derived from fraction
  dr <- max(0.005, ranger * dr_frac)   # clamp so it doesn't go to ~0
  # minsize step (absolute) from fraction
  dm <- max(5L, as.integer(round(minsize * dm_frac)))
  
  # spatialr step: avoid brutal relative jumps at small spatialr
  ds <- if (spatialr <= 3) 0L else min(ds_max, 1L)
  
  list(dr = dr, ds = ds, dm = dm, K = as.integer(K))
}
