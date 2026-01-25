#!/usr/bin/env Rscript

# -------------------------------------------------------------------
# Script: 04-0_run_S4_signatures.R
# Project: Burgwald
#
# Purpose
# -------
# Orchestrate S4 signature generation scripts and build the downstream
# "truth source" attribute stack:
#   layer0_segments_attrstack_metrics.gpkg
#
# Execution rule (contract)
# -------------------------
# Each S4 metrics script runs only if:
#   - overwrite == TRUE, OR
#   - its expected output file does not exist yet.
#
# Notes
# -----
# - The hydro part is intentionally excluded here. Hydro-related metrics are
#   evaluated within the physio block in this project state.
# - The attrstack builder joins only the metric files that exist.
# -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(dplyr)
})

source(here::here("src", "_core", "01-setup-burgwald.R"))

# ------------------------------------------------------------------
# Global overwrite switch
#   overwrite = TRUE  -> force recomputation even if outputs exist
#   overwrite = FALSE -> run only if expected output is missing
# ------------------------------------------------------------------
overwrite <- FALSE

# ------------------------------------------------------------------
# Block switches (TRUE = enabled)
# ------------------------------------------------------------------
run_it     <- TRUE
run_physio <- TRUE
run_cover  <- TRUE
run_bio    <- TRUE
run_attrstack <- TRUE

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
must_have_path <- function(key) {
  if (!(key %in% names(paths))) stop("Missing paths[['", key, "']] in setup.")
  paths[[key]]
}

should_run <- function(out_key) {
  out_file <- must_have_path(out_key)
  if (isTRUE(overwrite)) return(TRUE)
  !file.exists(out_file)
}

run_one <- function(name, file) {
  message("\n=== RUN: ", name, " ===\n", file)
  if (!file.exists(file)) stop("Missing script: ", file)
  source(file, local = new.env(parent = globalenv()))
}

# ------------------------------------------------------------------
# Scripts (hydro intentionally omitted in this project state)
# ------------------------------------------------------------------
scripts <- list(
  it     = here::here("src", "04-1_signatures_it_metrics.R"),
  physio = here::here("src", "04-2_signatures_physio_metrics.R"),
  cover  = here::here("src", "04-3_signatures_RF-segs_classification.R"),
  bio    = here::here("src", "04-4_signatures_biostructure_metrics.R")
)

# ------------------------------------------------------------------
# Run blocks under overwrite/missing-output rule
# ------------------------------------------------------------------
if (run_it && should_run("layer0_attr_it_metrics")) {
  tryCatch(run_one("S4_it", scripts$it),
           error = function(e) message("FAIL S4_it: ", e$message))
} else if (run_it) {
  message("SKIP S4_it (output exists and overwrite=FALSE): ", paths[["layer0_attr_it_metrics"]])
}

if (run_physio && should_run("layer0_attr_physio_metrics")) {
  tryCatch(run_one("S4_physio", scripts$physio),
           error = function(e) message("FAIL S4_physio: ", e$message))
} else if (run_physio) {
  message("SKIP S4_physio (output exists and overwrite=FALSE): ", paths[["layer0_attr_physio_metrics"]])
}

if (run_cover && should_run("layer0_attr_rf_metrics")) {
  tryCatch(run_one("S4_cover", scripts$cover),
           error = function(e) message("FAIL S4_cover: ", e$message))
} else if (run_cover) {
  message("SKIP S4_cover (output exists and overwrite=FALSE): ", paths[["layer0_attr_rf_metrics"]])
}

if (run_bio && should_run("layer0_attr_biostructure_metrics")) {
  tryCatch(run_one("S4_bio", scripts$bio),
           error = function(e) message("FAIL S4_bio: ", e$message))
} else if (run_bio) {
  message("SKIP S4_bio (output exists and overwrite=FALSE): ", paths[["layer0_attr_biostructure_metrics"]])
}

# ------------------------------------------------------------------
# FINAL: build the segment-level attribute "truth source" for S5
#        layer0_segments_attrstack_metrics.gpkg
# ------------------------------------------------------------------
build_attrstack <- function() {
  
  base_file <- must_have_path("layer0_segments")
  out_file  <- must_have_path("layer0_segments_attrstack_metrics")
  
  if (!file.exists(base_file)) stop("Missing base segments: ", base_file)
  
  base <- sf::read_sf(base_file)
  if (!("segment_id" %in% names(base))) stop("Base segments missing segment_id: ", base_file)
  
  # Candidate metric products produced by S4 scripts (join only what exists).
  # Hydro is intentionally not listed (handled inside physio in this project state).
  metric_keys <- c(
    "layer0_attr_it_metrics",
    "layer0_attr_physio_metrics",
    "layer0_attr_biostructure_metrics",
    "layer0_attr_rf_metrics"
  )
  
  metric_files <- vapply(metric_keys, function(k) {
    if (k %in% names(paths)) paths[[k]] else NA_character_
  }, character(1))
  
  metric_files <- metric_files[is.finite(nchar(metric_files)) & file.exists(metric_files)]
  
  if (length(metric_files) == 0) {
    stop("No metric GPKGs found to build attrstack. Checked keys: ",
         paste(metric_keys, collapse = ", "))
  }
  
  out <- base %>% dplyr::select(segment_id, geom)
  
  for (f in metric_files) {
    message("ATTRSTACK JOIN: ", f)
    
    m <- sf::read_sf(f)
    if (!("segment_id" %in% names(m))) stop("Metrics file missing segment_id: ", f)
    
    # keep base geometry stable
    m <- sf::st_drop_geometry(m)
    
    # avoid accidental duplicate column names (except segment_id)
    dup <- intersect(setdiff(names(m), "segment_id"), names(out))
    if (length(dup) > 0) {
      names(m)[match(dup, names(m))] <- paste0(dup, "__dup")
    }
    
    out <- out %>% dplyr::left_join(m, by = "segment_id")
  }
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(out, out_file, delete_dsn = TRUE, quiet = TRUE)
  message("WROTE ATTRSTACK: ", out_file, "  (n=", nrow(out), ", p=", ncol(out), ")")
}

if (run_attrstack && should_run("layer0_segments_attrstack_metrics")) {
  tryCatch(build_attrstack(), error = function(e) message("FAIL ATTRSTACK: ", e$message))
} else if (run_attrstack) {
  message("SKIP ATTRSTACK (output exists and overwrite=FALSE): ", paths[["layer0_segments_attrstack_metrics"]])
}
