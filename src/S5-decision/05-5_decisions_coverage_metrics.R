#!/usr/bin/env Rscript
#######################################################################
# Script: 05-5_decisions_coverage_metrics.R
# Project: Burgwald
#
# Purpose
# -------
# S5 decisions (Coverage/LUCC): build a "coverage strata" typification over
# S3 segments using landcover fractions already present in the S4 attrstack,
# then select representative candidates per stratum.
#
# Conceptual role
# --------------
# This is not a landcover classifier. It consumes S4 outputs (fractions) and
# builds a discrete typification (strata) for later decision steps.
#
# "Stratum" here means:
#   a partition of the coverage-fraction space that should be represented
#   by candidates (for network design evaluation, trade-off analysis, etc.).
#
# IMPORTANT: domain restriction BEFORE selection
# ----------------------------------------------
# Candidates are searched ONLY within the gauge-station hull polygon.
# The spatial predicate is applied BEFORE scaling, clustering, and selection.
#
# Outputs
# -------
# - layer0_coverage_strata_segments: ALL segments with coverage stratum fields
#   attached (outside hull -> NA, because the coverage decision space is not
#   evaluated there).
# - layer0_coverage_candidates_pts: candidate points (only inside hull).
#######################################################################

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stats)
  
  library(clustertend)  # hopkins()
  library(fpc)          # pamk(), kmeansruns(), clusterboot()
  library(cluster)      # silhouette()
  library(NbClust)      # NbClust()
})

source(here::here("src", "_core", "01-setup-burgwald.R"))

## -------------------------------------------------------------------
## Parameters (same structure as physio/hydro)
## -------------------------------------------------------------------
n_per_stratum <- 3L

k_min <- 2L
k_max <- 20L
k_scan <- k_min:k_max

cor_cutoff <- 0.9

kmeans_nstart_final <- 50L
kmeans_nstart_eval  <- 20L
clusterboot_B       <- 30L

k_eval_frac <- 0.2
crs_proj <- 25832

set.seed(42)

## -------------------------------------------------------------------
## Productive input
## -------------------------------------------------------------------
attr_file <- paths[["layer0_segments_attrstack_metrics"]]
stopifnot(file.exists(attr_file))

## -------------------------------------------------------------------
## Domain polygon (gauge-station hull)
## -------------------------------------------------------------------
hull_file <- paths[["gauges_hull"]]
stopifnot(file.exists(hull_file))

## -------------------------------------------------------------------
## Productive outputs
## -------------------------------------------------------------------
cov_seg_out <- paths[["layer0_coverage_strata_segments"]]
cov_pts_out <- paths[["layer0_coverage_candidates_pts"]]

dir.create(dirname(cov_seg_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cov_pts_out), recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## Helpers (copied from physio/hydro; unchanged)
## -------------------------------------------------------------------

drop_correlated_greedy <- function(Xz, cutoff = 0.95) {
  stopifnot(is.matrix(Xz))
  cm <- abs(stats::cor(Xz, use = "pairwise.complete.obs"))
  diag(cm) <- 0
  
  cols <- colnames(cm)
  keep <- rep(TRUE, length(cols))
  names(keep) <- cols
  
  for (j in seq_along(cols)) {
    if (!keep[j]) next
    too_corr <- which(cm[j, ] > cutoff & keep)
    if (length(too_corr) > 0) keep[too_corr] <- FALSE
  }
  
  names(keep)[keep]
}

calinski_harabasz <- function(km, n) {
  k <- as.integer(length(km$size))
  if (k <= 1L) return(NA_real_)
  bss <- as.numeric(km$betweenss)
  wss <- as.numeric(km$tot.withinss)
  (bss / (k - 1)) / (wss / (n - k))
}

davies_bouldin <- function(X, cl, centers) {
  X <- as.matrix(X)
  cl <- as.integer(cl)
  k <- nrow(centers)
  if (k <= 1L) return(NA_real_)
  
  S <- numeric(k)
  for (i in seq_len(k)) {
    idx <- which(cl == i)
    if (length(idx) == 0) { S[i] <- NA_real_; next }
    dif <- X[idx, , drop = FALSE] -
      matrix(centers[i, ], nrow = length(idx), ncol = ncol(X), byrow = TRUE)
    S[i] <- mean(sqrt(rowSums(dif^2)))
  }
  
  M <- as.matrix(dist(centers))
  diag(M) <- NA_real_
  
  R <- numeric(k)
  for (i in seq_len(k)) {
    rij <- (S[i] + S) / M[i, ]
    R[i] <- suppressWarnings(max(rij, na.rm = TRUE))
  }
  
  mean(R, na.rm = TRUE)
}

choose_k_internal_stability_clusterboot <- function(Xz, k_scan,
                                                    nstart = 20L,
                                                    boot_B = 30L,
                                                    bootmethod = "boot",
                                                    seed = 1L) {
  stopifnot(is.matrix(Xz), nrow(Xz) > 10, ncol(Xz) >= 2)
  
  dX <- dist(Xz)
  n  <- nrow(Xz)
  
  tab <- lapply(k_scan, function(k) {
    message("k-eval: ", k)
    
    km <- kmeans(Xz, centers = k, nstart = nstart)
    
    sil <- mean(cluster::silhouette(km$cluster, dX)[, 3])
    ch  <- calinski_harabasz(km, n = n)
    db  <- davies_bouldin(Xz, km$cluster, km$centers)
    wss <- as.numeric(km$tot.withinss)
    
    stab <- NA_real_
    cb <- try(
      fpc::clusterboot(
        Xz,
        B = boot_B,
        bootmethod = bootmethod,
        clustermethod = fpc::kmeansCBI,
        k = k,
        seed = seed
      ),
      silent = TRUE
    )
    if (!inherits(cb, "try-error") && !is.null(cb$result$bootmean)) {
      stab <- mean(cb$result$bootmean, na.rm = TRUE)
    }
    
    tibble::tibble(
      k = as.integer(k),
      sil = as.numeric(sil),
      ch  = as.numeric(ch),
      db  = as.numeric(db),
      wss = as.numeric(wss),
      stab = as.numeric(stab)
    )
  }) %>% dplyr::bind_rows()
  
  tab_n <- tab %>%
    dplyr::mutate(
      sil_n  = as.numeric(scale(sil)),
      ch_n   = as.numeric(scale(ch)),
      db_n   = -as.numeric(scale(db)),
      wss_n  = -as.numeric(scale(wss)),
      stab_n = as.numeric(scale(stab))
    ) %>%
    dplyr::mutate(
      score = 1.0*sil_n + 1.0*ch_n + 1.0*db_n + 0.5*wss_n + 2.0*stab_n
    )
  
  tab_ok <- tab_n %>% dplyr::filter(is.finite(stab_n))
  if (nrow(tab_ok) > 0) {
    k_best <- tab_ok$k[which.max(tab_ok$score)]
  } else {
    k_best <- tab_n$k[which.max(tab_n$score)]
  }
  
  list(k = as.integer(k_best), table = tab_n)
}

choose_k_aggregate <- function(k_vec, k_min, k_max) {
  kv <- as.integer(k_vec)
  kv <- kv[is.finite(kv)]
  kv <- kv[kv >= k_min & kv <= k_max]
  if (length(kv) == 0) stop("No valid k suggestions available after filtering.")
  tab <- sort(table(kv), decreasing = TRUE)
  if (length(tab) >= 1 && tab[1] >= 2) return(as.integer(names(tab)[1]))
  as.integer(stats::median(kv))
}

get_k_from_nbclust <- function(nb) {
  if (is.null(nb) || is.null(nb$Best.nc)) return(NA_integer_)
  bn <- nb$Best.nc
  k <- suppressWarnings(as.integer(bn[1, 1]))
  if (!is.finite(k)) return(NA_integer_)
  k
}

## -------------------------------------------------------------------
## 1) Read input + hull; apply domain restriction (WITHIN) EARLY
## -------------------------------------------------------------------
segs_sf_all <- sf::read_sf(attr_file)
stopifnot(inherits(segs_sf_all, "sf"))
stopifnot("segment_id" %in% names(segs_sf_all))

# Project convention: geometry column is "geom"
stopifnot("geom" %in% names(segs_sf_all))

hull <- sf::read_sf(hull_file)
stopifnot(inherits(hull, "sf"))
hull <- sf::st_make_valid(hull)
hull <- sf::st_union(hull)

if (sf::st_crs(segs_sf_all) != sf::st_crs(hull)) {
  hull <- sf::st_transform(hull, sf::st_crs(segs_sf_all))
}

# Strict domain restriction: only segments inside hull are eligible.
segs_sf <- segs_sf_all[sf::st_within(segs_sf_all, hull, sparse = FALSE), ]
message("Segments within hull: ", nrow(segs_sf), " / ", nrow(segs_sf_all))
stopifnot(nrow(segs_sf) > 0)

## -------------------------------------------------------------------
## 2) Coverage decision space (ONLY difference to physio/hydro)
## -------------------------------------------------------------------
# In your attrstack, landcover fractions are stored as p_* columns.
# These are interpreted as compositional "coverage" descriptors.
cover_cols <- c(
  "p_cropland",
  "p_grassland",
  "p_urban_sealed",
  "p_water",
  "p_forest_broadleaf",
  "p_forest_coniferous",
  "p_forest_mixed"
)
stopifnot(all(cover_cols %in% names(segs_sf)))

cov_df <- segs_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(segment_id, dplyr::all_of(cover_cols)) %>%
  tidyr::drop_na()

stopifnot(nrow(cov_df) > 0)

X  <- as.matrix(cov_df[, cover_cols, drop = FALSE])
Xz <- scale(X)

# If coverage fractions are highly collinear, drop redundant dimensions.
keep_cols <- drop_correlated_greedy(Xz, cutoff = cor_cutoff)
stopifnot(length(keep_cols) >= 2)

data_for_k <- as.matrix(Xz[, keep_cols, drop = FALSE])
message("Coverage vars: ", length(cover_cols), " total; ",
        length(keep_cols), " kept after |cor| <= ", cor_cutoff)

## -------------------------------------------------------------------
## 3) k-evaluation on a subsample (same logic as physio/hydro)
## -------------------------------------------------------------------
n_all  <- nrow(data_for_k)
n_eval <- max(2000L, as.integer(round(k_eval_frac * n_all)))
n_eval <- min(n_eval, n_all)

set.seed(42)
idx_eval <- sample.int(n_all, size = n_eval, replace = FALSE)
data_for_k_eval <- data_for_k[idx_eval, , drop = FALSE]

message("k-evaluation sample: ", n_eval, " / ", n_all, " (frac=", k_eval_frac, ")")

h_n <- min(1000L, nrow(data_for_k_eval))
hop <- clustertend::hopkins(as.data.frame(data_for_k_eval), n = h_n)
message("Hopkins statistic (n=", h_n, "): ", format(hop, digits = 5))

pamk_best <- fpc::pamk(
  as.data.frame(data_for_k_eval),
  criterion = "asw",
  krange = k_scan,
  usepam = TRUE,
  ns = 2,
  critout = TRUE
)
k_pamk <- as.integer(pamk_best$nc)
message("pamk (ASW) suggested k: ", k_pamk)

kmeans_best <- fpc::kmeansruns(
  as.data.frame(data_for_k_eval),
  krange = k_scan,
  runs = 10,
  criterion = c("asw"),
  critout = TRUE
)
k_kmeansruns <- as.integer(kmeans_best$bestk)
message("kmeansruns (ASW) suggested k: ", k_kmeansruns)

nb_wardD2 <- NbClust::NbClust(
  as.data.frame(data_for_k_eval),
  distance = "euclidean",
  min.nc = k_min, max.nc = k_max,
  method = "ward.D2",
  index = "dindex",
  alphaBeale = 0.05
)
k_nb_wardD2 <- get_k_from_nbclust(nb_wardD2)
message("NbClust (ward.D2, dindex) suggested k: ", k_nb_wardD2)

cb_res <- choose_k_internal_stability_clusterboot(
  Xz      = data_for_k_eval,
  k_scan  = k_scan,
  nstart  = kmeans_nstart_eval,
  boot_B  = clusterboot_B,
  seed    = 1L
)
message("clusterboot internal+stability table (sample):")
print(cb_res$table)
k_cb <- cb_res$k
message("clusterboot suggested k (sample): ", k_cb)

k_suggestions <- c(k_pamk, k_kmeansruns, k_nb_wardD2, k_cb)
message("k suggestions (raw): ", paste(k_suggestions, collapse = ", "))

k_cov <- choose_k_aggregate(k_suggestions, k_min = k_min, k_max = k_max)
message("Aggregated operational k_cov: ", k_cov)

## -------------------------------------------------------------------
## 4) Final k-means stratification (FULL, within hull, complete cases)
## -------------------------------------------------------------------
cov_km <- kmeans(data_for_k, centers = k_cov, nstart = kmeans_nstart_final)

centers_assigned <- cov_km$centers[cov_km$cluster, , drop = FALSE]
cov_dist <- sqrt(rowSums((data_for_k - centers_assigned)^2))

cov_assign <- tibble::tibble(
  segment_id          = cov_df$segment_id,
  cov_stratum         = as.integer(cov_km$cluster),
  cov_dist_to_center  = as.numeric(cov_dist),
  cov_k               = as.integer(k_cov)
)

## -------------------------------------------------------------------
## 5) Attach strata back to ALL segments (product)
## -------------------------------------------------------------------
cov_seg_sf_all <- segs_sf_all %>%
  dplyr::left_join(cov_assign, by = "segment_id")

cov_seg_sf_in <- segs_sf %>%
  dplyr::left_join(cov_assign, by = "segment_id")

## -------------------------------------------------------------------
## 6) Candidate selection per stratum (ONLY within hull + assigned)
## -------------------------------------------------------------------
cov_candidates_ids <- cov_assign %>%
  dplyr::group_by(cov_stratum) %>%
  dplyr::slice_min(cov_dist_to_center, n = n_per_stratum, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(segment_id)

cov_candidates_seg <- cov_seg_sf_in %>%
  dplyr::semi_join(cov_candidates_ids, by = "segment_id") %>%
  dplyr::select(segment_id, cov_stratum, cov_dist_to_center, geom)

## -------------------------------------------------------------------
## 7) Candidate points (planar point_on_surface)
## -------------------------------------------------------------------
cov_candidates_pts <- cov_candidates_seg %>%
  sf::st_make_valid() %>%
  sf::st_transform(crs_proj) %>%
  sf::st_point_on_surface() %>%
  sf::st_transform(sf::st_crs(cov_candidates_seg)) %>%
  dplyr::mutate(space = "coverage") %>%
  dplyr::select(
    segment_id,
    space,
    stratum = cov_stratum,
    dist_to_center = cov_dist_to_center,
    geom
  )

## -------------------------------------------------------------------
## 8) Write outputs
## -------------------------------------------------------------------
sf::st_write(
  cov_seg_sf_all,
  dsn = cov_seg_out,
  layer = "coverage_strata_segments",
  driver = "GPKG",
  delete_dsn = TRUE,
  quiet = TRUE
)

sf::st_write(
  cov_candidates_pts,
  dsn = cov_pts_out,
  layer = "coverage_candidates_pts",
  driver = "GPKG",
  delete_dsn = TRUE,
  quiet = TRUE
)

message("Wrote segments:   ", cov_seg_out, " (rows=", nrow(cov_seg_sf_all), ")")
message("Wrote candidates: ", cov_pts_out, " (rows=", nrow(cov_candidates_pts), ")")
