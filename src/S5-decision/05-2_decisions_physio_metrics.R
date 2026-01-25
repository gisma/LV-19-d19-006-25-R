#!/usr/bin/env Rscript
#######################################################################
# Script: 05-2_decisions_physio_metrics.R
# Project: Burgwald
#
# Role in the pipeline
# --------------------
# This script implements an S5-style "decision support" step for the
# physio space: it builds a physio stratification (clusters/strata)
# over the already-segmented spatial units (S3 segments with S4 metrics),
# then selects representative candidate segments/points per stratum.
#
# Key idea: "stratum" ≠ truth
# --------------------------
# A "physio_stratum" is a discrete typification of continuous relief
# signatures. It is an operational partition of the candidate space
# used to:
#   (a) enforce coverage across different relief regimes,
#   (b) prevent selection collapsing into one dominant regime.
#
# IMPORTANT CHANGE (requested)
# ----------------------------
# Candidate search is constrained to a polygon derived from gauge stations
# (the hull). The "within" predicate is applied BEFORE candidate selection.
# This ensures candidates are only sought in the intended subdomain.
#
# NOTE on consistency filters
# ---------------------------
# Clustering uses complete cases of the physio variables (drop_na()).
# Therefore not every segment inside the hull necessarily receives a
# cluster assignment. Candidate selection must operate only on segments
# that actually have an assignment (physio_assign).
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
# metrics-fun.R assumed to exist in your repo; keep as in your project
# source(here::here("src", "r-libs", "metrics-fun.R"))

## -------------------------------------------------------------------
## Parameters (keep minimal, project-style)
## -------------------------------------------------------------------
n_per_stratum <- 5L

k_min <- 4L
k_max <- 25L
k_scan <- k_min:k_max

cor_cutoff <- 0.9

kmeans_nstart_final <- 50L
kmeans_nstart_eval  <- 20L
clusterboot_B       <- 30L

k_eval_frac <- 0.2

# Use projected CRS for geometric operations if lon/lat is present.
crs_proj <- 25832
set.seed(42)

## -------------------------------------------------------------------
## Input: segments with S4 metrics
## -------------------------------------------------------------------
attr_file <- paths[["layer0_segments_attrstack_metrics"]]
stopifnot(file.exists(attr_file))

## -------------------------------------------------------------------
## Input: hull polygon (gauge station hull)
## -------------------------------------------------------------------
# Use the key you actually have in your project.
# If the key differs: replace paths[["gauges_hull_gpkg"]] accordingly.
hull_file <- paths[["gauges_hull"]]
stopifnot(file.exists(hull_file))

## -------------------------------------------------------------------
## Outputs
## -------------------------------------------------------------------
# Output 1: all segments, with physio_* fields attached (outside hull -> NA)
physio_seg_out <- paths[["layer0_physio_strata_segments"]]

# Output 2: candidates points (only within hull by construction)
physio_pts_out <- paths[["layer0_physio_candidates_pts"]]

dir.create(dirname(physio_seg_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(physio_pts_out), recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## Helpers (kept compact; deterministic behaviour)
## -------------------------------------------------------------------

# Greedy removal of highly correlated predictors after standardisation.
# This is not "model selection"; it is a guard against redundant dimensions
# dominating Euclidean distance in k-means.
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

# Internal evaluation table across k: silhouette + CH + DB + WSS + bootstrap stability.
# This is a pragmatic "multiple weak signals" approach: no single index is decisive.
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
      # Weighted aggregate: stability gets extra weight because we don't want
      # "nice indices" that are unstable under resampling.
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
## 1) Read segments and hull; apply domain restriction (WITHIN) EARLY
## -------------------------------------------------------------------
segs_sf_all <- sf::read_sf(attr_file)
stopifnot(inherits(segs_sf_all, "sf"))
stopifnot("segment_id" %in% names(segs_sf_all))

# In your project the geometry column is called "geom".
# Keep that, do not invent other names.
stopifnot("geom" %in% names(segs_sf_all))

hull <- sf::read_sf(hull_file)
stopifnot(inherits(hull, "sf"))
hull <- sf::st_make_valid(hull)
hull <- sf::st_union(hull)

# Ensure same CRS for spatial predicate.
if (sf::st_crs(segs_sf_all) != sf::st_crs(hull)) {
  hull <- sf::st_transform(hull, sf::st_crs(segs_sf_all))
}

# STRICT domain restriction: keep only segments fully inside the hull.
# This is the requested behaviour: candidates are searched only inside.
segs_sf <- segs_sf_all[sf::st_within(segs_sf_all, hull, sparse = FALSE), ]

message("Segments within hull: ", nrow(segs_sf), " / ", nrow(segs_sf_all))
stopifnot(nrow(segs_sf) > 0)

## -------------------------------------------------------------------
## 2) Build physio feature matrix from within-hull segments
## -------------------------------------------------------------------
physio_cols <- grep("^relief__", names(segs_sf), value = TRUE)
stopifnot(length(physio_cols) > 0)

# Clustering space is defined over complete cases (no NA).
# This defines the effective segment set that can receive a stratum.
physio_df <- segs_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::select(segment_id, dplyr::all_of(physio_cols)) %>%
  tidyr::drop_na()

stopifnot(nrow(physio_df) > 0)

X  <- as.matrix(physio_df[, physio_cols, drop = FALSE])
Xz <- scale(X)

keep_cols <- drop_correlated_greedy(Xz, cutoff = cor_cutoff)
stopifnot(length(keep_cols) >= 2)

data_for_k <- as.matrix(Xz[, keep_cols, drop = FALSE])

message("Physio vars: ", length(physio_cols), " total; ",
        length(keep_cols), " kept after |cor| <= ", cor_cutoff)

## -------------------------------------------------------------------
## 3) k-evaluation on a subsample (expensive step)
## -------------------------------------------------------------------
n_all  <- nrow(data_for_k)
n_eval <- max(2000L, as.integer(round(k_eval_frac * n_all)))
n_eval <- min(n_eval, n_all)

set.seed(42)
idx_eval <- sample.int(n_all, size = n_eval, replace = FALSE)
data_for_k_eval <- data_for_k[idx_eval, , drop = FALSE]

message("k-evaluation sample: ", n_eval, " / ", n_all, " (frac=", k_eval_frac, ")")

# Cluster tendency check: Hopkins near 0.5 ~ random, near 0 ~ strong clustering.
h_n <- min(1000L, nrow(data_for_k_eval))
hop <- clustertend::hopkins(as.data.frame(data_for_k_eval), n = h_n)
message("Hopkins statistic (n=", h_n, "): ", format(hop, digits = 5))

# Multiple k-suggestions (no single oracle).
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

k_physio <- choose_k_aggregate(k_suggestions, k_min = k_min, k_max = k_max)
message("Aggregated operational k_physio: ", k_physio)

## -------------------------------------------------------------------
## 4) Final k-means stratification (FULL, within hull, complete cases)
## -------------------------------------------------------------------
physio_km <- kmeans(data_for_k, centers = k_physio, nstart = kmeans_nstart_final)

# Distance to assigned centroid in standardised variable space.
centers_assigned <- physio_km$centers[physio_km$cluster, , drop = FALSE]
physio_dist <- sqrt(rowSums((data_for_k - centers_assigned)^2))

physio_assign <- tibble::tibble(
  segment_id            = physio_df$segment_id,
  physio_stratum        = as.integer(physio_km$cluster),
  physio_dist_to_center = as.numeric(physio_dist),
  physio_k              = as.integer(k_physio)
)

## -------------------------------------------------------------------
## 5) Attach physio strata back to segments
## -------------------------------------------------------------------
# Two outputs are useful:
# - segs_sf_all: full AOI (outside hull -> NA because no assignment there)
# - segs_sf (within hull): should receive assignments where complete cases exist
#
# This keeps the "all segments" product while respecting the requested
# domain restriction for candidate search.
physio_seg_sf_all <- segs_sf_all %>%
  dplyr::left_join(physio_assign, by = "segment_id")

physio_seg_sf_in <- segs_sf %>%
  dplyr::left_join(physio_assign, by = "segment_id")

## -------------------------------------------------------------------
## 6) Candidate selection per stratum (ONLY within hull + assigned)
## -------------------------------------------------------------------
# Candidate selection operates on physio_assign, which already:
#   - is within hull (because physio_df came from segs_sf)
#   - is complete-case (because physio_df was drop_na())
#
# This is the minimal correct domain for selecting "closest to centroid".
physio_candidates_ids <- physio_assign %>%
  dplyr::group_by(physio_stratum) %>%
  dplyr::slice_min(physio_dist_to_center, n = n_per_stratum, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(segment_id)

physio_candidates_seg <- physio_seg_sf_in %>%
  dplyr::semi_join(physio_candidates_ids, by = "segment_id") %>%
  dplyr::select(segment_id, physio_stratum, physio_dist_to_center, geom)

## -------------------------------------------------------------------
## 7) Candidate points (use planar CRS for point_on_surface)
## -------------------------------------------------------------------
physio_candidates_pts <- physio_candidates_seg %>%
  sf::st_make_valid() %>%
  sf::st_transform(crs_proj) %>%
  sf::st_point_on_surface() %>%
  sf::st_transform(sf::st_crs(physio_candidates_seg)) %>%
  dplyr::mutate(space = "physio") %>%
  dplyr::select(
    segment_id,
    space,
    stratum = physio_stratum,
    dist_to_center = physio_dist_to_center,
    geom
  )

## -------------------------------------------------------------------
## 8) Write outputs
## -------------------------------------------------------------------
# Write: full segments with physio assignment (outside hull -> NA)
sf::st_write(
  physio_seg_sf_all,
  dsn = physio_seg_out,
  layer = "physio_strata_segments",
  driver = "GPKG",
  delete_dsn = TRUE,
  quiet = TRUE
)

sf::st_write(
  physio_candidates_pts,
  dsn = physio_pts_out,
  layer = "physio_candidates_pts",
  driver = "GPKG",
  delete_dsn = TRUE,
  quiet = TRUE
)

message("Wrote segments:  ", physio_seg_out, " (rows=", nrow(physio_seg_sf_all), ")")
message("Wrote candidates:", physio_pts_out, " (rows=", nrow(physio_candidates_pts), ")")
