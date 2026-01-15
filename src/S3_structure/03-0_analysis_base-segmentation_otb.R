#!/usr/bin/env Rscript

# sourcing of the setup and specific used functions
source("_core/00-setup-burgwald.R")

# ---- helpers FIRST (create this file exactly as above) -----------------------
source("r-libs/metrics-fun.R")

# ---------------------------------------------------------
# 0) Inputs / outputs
# ---------------------------------------------------------
year <- 2018
pred_stack_file <- here::here("data", paste0("pred_stack_", year, ".tif"))
stopifnot(file.exists(pred_stack_file))

out_root <- here::here("data", "processed", "layer0_segments")
feat_dir <- file.path(out_root, "features")
seg_dir  <- file.path(out_root, "segments")
met_dir  <- file.path(out_root, "metrics")         # NEW
tmp_dir  <- file.path(out_root, "tmp_meanshift")   # NEW (raster outputs for ARI)
fs::dir_create(c(out_root, feat_dir, seg_dir, met_dir, tmp_dir))

# ---------------------------------------------------------
# 1) OTB environment + runner
# ---------------------------------------------------------
otblink <- link2GI::linkOTB(searchLocation = "~/apps/otb911")

# ---------------------------------------------------------
# 2) Standardize (z-score) with BandMathX
# ---------------------------------------------------------
r  <- terra::rast(pred_stack_file)
nb <- terra::nlyr(r)
rm(r); gc()

expr <- paste(
  vapply(seq_len(nb), function(i) {
    sprintf("(im1b%d - im1b%dMean) / sqrt(im1b%dVar)", i, i, i)
  }, character(1)),
  collapse = ";"
)
expr <- paste0("'", expr, "'")

algo <- "BandMathX"
std_file <- file.path(out_root, paste0(algo, "_zscore_stack.tif"))
cmd <- parseOTBFunction(algo = algo, gili = otblink)
cmd$il  <- pred_stack_file
cmd$out <- std_file
cmd$exp <- expr
cmd[["progress"]] <- 1
runOTB(cmd, gili = otblink, quiet = FALSE, retRaster = FALSE)
stopifnot(file.exists(std_file))
message("Standardized stack written: ", std_file)

# ---------------------------------------------------------
# 3) PCA (OTB) on standardized stack
# ---------------------------------------------------------
ncomp <- 6
pca_file <- file.path(feat_dir, paste0("features_otb_pca_y", year, "_pc", ncomp, ".tif"))

algo <- "DimensionalityReduction"
cmd <- parseOTBFunction(algo = algo, gili = otblink)
cmd[["in"]]     <- std_file
cmd[["method"]] <- "pca"
cmd[["nbcomp"]] <- as.character(ncomp)
cmd$out <- pca_file
runOTB(cmd, gili = otblink, quiet = FALSE, retRaster = FALSE)
stopifnot(file.exists(pca_file))

# ---------------------------------------------------------
# 4) LargeScaleMeanShift multiscale: METRICS FIRST (ARI_prev)
# ---------------------------------------------------------
param_grid <- data.frame(
  scale_id  = c("s30m", "s60m", "s120m"),
  spatialr  = c(3, 6, 12),
  ranger    = c(0.10, 0.12, 0.15),
  minsize   = c(80, 120, 200),
  stringsAsFactors = FALSE
)

# --- tuning knobs for ARI_prev (explicit) ------------------------------------
ari_cfg <- list(
  dr = 0.01,  # ranger perturbation
  ds = 1L,    # spatialr perturbation
  dm = 20L,   # minsize perturbation
  K  = 8L     # number of perturbed runs
)
sample_n <- 2e6   # pixels sampled for ARI (reduce if RAM tight)
tilesiz  <- 500L
ram_mb   <- NULL  # set e.g. 8192 if you want explicit
seed0    <- 42L

metrics <- vector("list", nrow(param_grid))

for (i in seq_len(nrow(param_grid))) {
  p <- param_grid[i, ]
  message(sprintf("ARI_prev: %s (spatialr=%d, ranger=%.3f, minsize=%d)",
                  p$scale_id, p$spatialr, p$ranger, p$minsize))
  
  out_tmp <- file.path(tmp_dir, paste0("ari_", p$scale_id, "_y", year))
  res <- compute_ari_prev(
    otblink    = otblink,
    feat_raster= pca_file,
    out_dir_tmp= out_tmp,
    spatialr   = p$spatialr,
    ranger     = p$ranger,
    minsize    = p$minsize,
    perturb    = ari_cfg,
    sample_n   = sample_n,
    seed       = seed0 + i * 100L,
    tilesize   = tilesiz,
    ram        = ram_mb,
    quiet      = TRUE
  )
  
  metrics[[i]] <- data.frame(
    scale_id = p$scale_id,
    spatialr = p$spatialr,
    ranger   = p$ranger,
    minsize  = p$minsize,
    ARI_prev = res$ari_prev,
    ARI_sd   = res$ari_sd
  )
}

metrics_df <- do.call(rbind, metrics)
metrics_file <- file.path(met_dir, paste0("meanshift_ari_prev_y", year, ".csv"))
write.csv(metrics_df, metrics_file, row.names = FALSE)
message("ARI_prev metrics written: ", metrics_file)
print(metrics_df)

# ---------------------------------------------------------
# 5) (optional) Final segmentation export (vector) ONLY for chosen scales
# ---------------------------------------------------------
# Example rule: take the best ARI_prev (you can replace by your own logic)
best_id <- metrics_df$scale_id[which.max(metrics_df$ARI_prev)]
message("Best (by ARI_prev): ", best_id)

for (i in seq_len(nrow(param_grid))) {
  p <- param_grid[i, ]
  if (p$scale_id != best_id) next
  
  out_vec <- file.path(seg_dir, paste0("segments_meanshift_", p$scale_id, "_y", year, ".gpkg"))
  
  algo <- "LargeScaleMeanShift"
  cmd <- parseOTBFunction(algo = algo, gili = otblink)
  
  cmd[["in"]] <- pca_file
  cmd[["spatialr"]] <- as.character(p$spatialr)
  cmd[["ranger"]]   <- as.character(p$ranger)
  cmd[["minsize"]]  <- as.character(p$minsize)
  
  cmd[["mode"]]            <- "vector"
  cmd[["mode.vector.out"]] <- out_vec
  
  cmd[["tilesizex"]] <- "500"
  cmd[["tilesizey"]] <- "500"
  
  runOTB(cmd, gili = otblink, quiet = FALSE, retRaster = FALSE)
  stopifnot(file.exists(out_vec))
  message("Final vector segmentation written: ", out_vec)
}

message("Layer 0 done.")
