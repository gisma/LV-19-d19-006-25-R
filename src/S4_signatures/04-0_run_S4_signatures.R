#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
})


source(here::here("src","_core","01-setup-burgwald.R"))

# global switches (TRUE = run)
run_it     <- TRUE   # 04-1 (IT / landscape metrics on segments)
run_physio <- TRUE   # 04-2
run_hydro  <- TRUE   # 04-3
run_cover  <- TRUE   # 04-4
run_bio    <- TRUE   # 04-5

scripts <- list(
  it     = here::here("src", "04-1_signatures_landscape-metrics.R"),
  physio = here::here("src", "04-2_signatures_physiographic-metrics.R"),
  hydro  = here::here("src", "04-5_signatures_hydrology_metrics.R"),
  cover  = here::here("src", "04-3_signatures_cover_metrics.R"),
  bio    = here::here("src", "04-4_signatures_biostructure_metrics.R")
)

run_one <- function(name, file) {
  message("\n=== RUN: ", name, " ===\n", file)
  if (!file.exists(file)) stop("Missing script: ", file)
  source(file, local = new.env(parent = globalenv()))
}

if (run_it)     tryCatch(run_one("S4_it",     scripts$it),     error = function(e) message("FAIL S4_it: ", e$message))
if (run_physio) tryCatch(run_one("S4_physio", scripts$physio), error = function(e) message("FAIL S4_physio: ", e$message))
if (run_hydro)  tryCatch(run_one("S4_hydro",  scripts$hydro),  error = function(e) message("FAIL S4_hydro: ", e$message))
if (run_cover)  tryCatch(run_one("S4_cover",  scripts$cover),  error = function(e) message("FAIL S4_cover: ", e$message))
if (run_bio)    tryCatch(run_one("S4_bio",    scripts$bio),    error = function(e) message("FAIL S4_bio: ", e$message))
