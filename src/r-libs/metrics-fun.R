# src/metrics-fun.R
# ------------------------------------------------------------------------------
# Metrics + stability utilities for OTB LargeScaleMeanShift scale selection
#
# Core idea:
#   ARI_prev(p) = mean Adjusted Rand Index between baseline segmentation S(p)
#                and K perturbed segmentations S(p_k).
#
# References (conceptual):
#   - Hubert, L. & Arabie, P. (1985). Comparing partitions. J. Classification.
#   - OTB Cookbook: LargeScaleMeanShift parameters and output modes.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(link2GI)
})

# ---- Adjusted Rand Index (Hubert & Arabie, 1985) -----------------------------
# x, y: integer/factor vectors of same length, labels per pixel.
adjusted_rand_index <- function(x, y) {
  stopifnot(length(x) == length(y))
  ok <- !(is.na(x) | is.na(y))
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2L) return(NA_real_)
  
  x <- as.integer(as.factor(x))
  y <- as.integer(as.factor(y))
  
  # contingency table
  tab <- table(x, y)
  
  # helper: nC2
  comb2 <- function(n) n * (n - 1) / 2
  
  sum_nij <- sum(comb2(tab))
  sum_ai  <- sum(comb2(rowSums(tab)))
  sum_bj  <- sum(comb2(colSums(tab)))
  n       <- sum(tab)
  total   <- comb2(n)
  
  if (total == 0) return(NA_real_)
  
  expected <- (sum_ai * sum_bj) / total
  maxind   <- 0.5 * (sum_ai + sum_bj)
  
  denom <- (maxind - expected)
  if (denom == 0) return(NA_real_)
  
  (sum_nij - expected) / denom
}

# ---- Create deterministic perturbations around a base parameter set ----------
# This is "true perturbation": it changes parameters and re-runs segmentation.
#
# Strategy:
#   - spatialr: +/- 1 (>=1)
#   - ranger:   +/- dr (>= small positive)
#   - minsize:  +/- dm (>= 1)
#
# K perturbations are built as a small factorial-ish neighborhood but capped.
make_param_perturbations <- function(spatialr, ranger, minsize,
                                     dr = 0.01, ds = 1L, dm = 20L,
                                     K = 8L) {
  base <- list(spatialr = as.integer(spatialr),
               ranger   = as.numeric(ranger),
               minsize  = as.integer(minsize))
  
  cand <- expand.grid(
    spatialr = c(base$spatialr - ds, base$spatialr + ds, base$spatialr),
    ranger   = c(base$ranger - dr,   base$ranger + dr,   base$ranger),
    minsize  = c(base$minsize - dm,  base$minsize + dm,  base$minsize),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # clean bounds
  cand$spatialr <- pmax(1L, as.integer(round(cand$spatialr)))
  cand$ranger   <- pmax(1e-6, as.numeric(cand$ranger))
  cand$minsize  <- pmax(1L, as.integer(round(cand$minsize)))
  
  # remove exact duplicates and the exact baseline
  cand <- unique(cand)
  cand <- cand[!(cand$spatialr == base$spatialr &
                   cand$ranger   == base$ranger &
                   cand$minsize  == base$minsize), , drop = FALSE]
  
  # cap to K (stable ordering for reproducibility)
  if (nrow(cand) > K) cand <- cand[seq_len(K), , drop = FALSE]
  cand
}

# ---- Run OTB LargeScaleMeanShift in raster mode (labels) ---------------------
# Output is the labeled raster of the "small region merging" step (OTB doc).
run_otb_meanshift_raster <- function(otblink, in_raster, out_raster,
                                     spatialr, ranger, minsize,
                                     tilesize = 500L,
                                     ram = NULL,
                                     quiet = FALSE) {
  algo <- "LargeScaleMeanShift"
  cmd  <- parseOTBFunction(algo = algo, gili = otblink)
  
  cmd[["in"]] <- in_raster
  cmd[["spatialr"]] <- as.character(spatialr)
  cmd[["ranger"]]   <- as.character(ranger)
  cmd[["minsize"]]  <- as.character(minsize)
  
  cmd[["mode"]] <- "raster"
  cmd[["mode.raster.out"]] <- out_raster
  
  cmd[["tilesizex"]] <- as.character(tilesize)
  cmd[["tilesizey"]] <- as.character(tilesize)
  
  if (!is.null(ram)) cmd[["ram"]] <- as.character(ram)
  
  ret <- runOTB(cmd, gili = otblink, quiet = quiet, retRaster = FALSE)
  stopifnot(file.exists(out_raster))
  invisible(out_raster)
}

# ---- Sample aligned pixels from two label rasters ----------------------------
# Using sampling avoids holding full rasters in RAM for ARI.
sample_label_pairs <- function(r_base, r_other, n = 2e6, seed = 1L) {
  stopifnot(terra::nlyr(r_base) == 1L, terra::nlyr(r_other) == 1L)
  
  # Ensure same geometry; if not, resample other to base
  if (!terra::compareGeom(r_base, r_other, stopOnError = FALSE)) {
    r_other <- terra::resample(r_other, r_base, method = "near")
  }
  
  set.seed(seed)
  # sample indices from base raster cells
  N <- terra::ncell(r_base)
  n <- min(as.integer(n), N)
  idx <- sample.int(N, size = n)
  
  x <- terra::values(r_base, cells = idx, mat = FALSE)
  y <- terra::values(r_other, cells = idx, mat = FALSE)
  list(x = x, y = y)
}

# ---- Compute ARI_prev for one base parameter set -----------------------------
compute_ari_prev <- function(otblink, feat_raster,
                             out_dir_tmp,
                             spatialr, ranger, minsize,
                             perturb = list(dr = 0.01, ds = 1L, dm = 20L, K = 8L),
                             sample_n = 2e6,
                             seed = 1L,
                             tilesize = 500L,
                             ram = NULL,
                             quiet = TRUE) {
  dir.create(out_dir_tmp, recursive = TRUE, showWarnings = FALSE)
  
  # 1) baseline
  base_file <- file.path(out_dir_tmp,
                         sprintf("ms_base_s%02d_r%.4f_m%04d.tif", spatialr, ranger, minsize))
  run_otb_meanshift_raster(
    otblink = otblink, in_raster = feat_raster, out_raster = base_file,
    spatialr = spatialr, ranger = ranger, minsize = minsize,
    tilesize = tilesize, ram = ram, quiet = quiet
  )
  
  r_base <- terra::rast(base_file)
  
  # 2) perturbed runs
  cand <- make_param_perturbations(spatialr, ranger, minsize,
                                   dr = perturb$dr, ds = perturb$ds, dm = perturb$dm,
                                   K = perturb$K)
  
  if (nrow(cand) == 0) {
    return(list(ari_prev = NA_real_, ari_sd = NA_real_, ari = numeric(0), base = base_file, cand = cand))
  }
  
  ari_vec <- numeric(nrow(cand))
  
  for (i in seq_len(nrow(cand))) {
    p <- cand[i, ]
    f <- file.path(out_dir_tmp,
                   sprintf("ms_pert%02d_s%02d_r%.4f_m%04d.tif",
                           i, p$spatialr, p$ranger, p$minsize))
    
    run_otb_meanshift_raster(
      otblink = otblink, in_raster = feat_raster, out_raster = f,
      spatialr = p$spatialr, ranger = p$ranger, minsize = p$minsize,
      tilesize = tilesize, ram = ram, quiet = quiet
    )
    
    r_other <- terra::rast(f)
    pairs <- sample_label_pairs(r_base, r_other, n = sample_n, seed = seed + i)
    ari_vec[i] <- adjusted_rand_index(pairs$x, pairs$y)
  }
  
  list(
    ari_prev = mean(ari_vec, na.rm = TRUE),
    ari_sd   = stats::sd(ari_vec, na.rm = TRUE),
    ari      = ari_vec,
    base     = base_file,
    cand     = cand
  )
}
