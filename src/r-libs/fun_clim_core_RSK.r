# fun_clim_core_RSK.R
# Minimal function bundle for DWD CDC KL daily -> RSK workflow (newpipe=TRUE).
# Extracted from fun_clim_data.R without refactoring.

fix_enc <- function(v) iconv(v, to = "UTF-8", sub = "")
  st_meta[, `:=`(
    Stationsname = fix_enc(Stationsname),
    Bundesland   = fix_enc(Bundesland),
    Abgabe       = fix_enc(Abgabe)
  )]
  
  st_meta_u <- st_meta[, .(
    Stationsname  = Stationsname[1],
    Bundesland    = Bundesland[1],
    Abgabe        = Abgabe[1],
    Stationshoehe = Stationshoehe[1],
    geoBreite     = geoBreite[1],
    geoLaenge     = geoLaenge[1]
  ), by = STATIONS_ID]
  
  st_span_u <- st_meta[, .(
    von_datum = as.integer(format(min(von_datum, na.rm = TRUE), "%Y%m%d")),
    bis_datum = as.integer(format(max(bis_datum, na.rm = TRUE), "%Y%m%d"))
  ), by = STATIONS_ID]
  
  setkey(st_span_u, STATIONS_ID)
  setkey(dt,        STATIONS_ID)
  dt <- st_span_u[dt]
  
  if (anyNA(dt$von_datum) || anyNA(dt$bis_datum)) {
    bad <- unique(dt[is.na(von_datum) | is.na(bis_datum), STATIONS_ID])
    stop("Missing von_datum/bis_datum for STATIONS_ID: ", paste(bad, collapse = ", "))
  }


list_remote_files <- function(url, pattern = NULL) {
    doc <- rvest::read_html(url)
    files <- doc |> rvest::html_elements("a") |> rvest::html_attr("href")
    files <- files[!is.na(files) & files != "../"]
    if (!is.null(pattern)) files <- files[grepl(pattern, files)]
    files <- unique(files)
    data.table(url = paste0(url, files), file = files)
  }


download_if_missing <- function(url, destfile) {
    if (file.exists(destfile)) return(invisible(FALSE))
    message("download zip: ", destfile)
    tryCatch({
      utils::download.file(url, destfile = destfile, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) {
      message("DOWNLOAD FAILED: ", url)
      message(conditionMessage(e))
      FALSE
    })
  }


unzip_zip <- function(zipfile, exdir) {
    suppressWarnings(utils::unzip(zipfile, exdir = exdir))
    invisible(TRUE)
  }


parse_prod_span <- function(files) {
    bn <- basename(files)
    m <- regexec("^produkt_klima_tag_(\\d{8})_(\\d{8})_(\\d{5})\\.txt$", bn)
    mm <- regmatches(bn, m)
    ok <- lengths(mm) == 4
    if (!all(ok)) stop("Unexpected product filenames encountered.")
    data.table(
      STATIONS_ID = vapply(mm, `[[`, character(1), 4),
      file_start  = as.Date(vapply(mm, `[[`, character(1), 2), "%Y%m%d"),
      file_end    = as.Date(vapply(mm, `[[`, character(1), 3), "%Y%m%d"),
      file        = bn
    )
  }


read_one_product_file <- function(f, source_tag) {
    bn <- basename(f)
    m <- regexec("^produkt_klima_tag_(\\d{8})_(\\d{8})_(\\d{5})\\.txt$", bn)
    mm <- regmatches(bn, m)[[1]]
    if (length(mm) != 4) stop("Unexpected product filename: ", bn)
    
    file_sid <- mm[4]
    
    dt <- data.table::fread(
      f,
      sep = ";",
      na.strings = "-999",
      strip.white = TRUE,
      showProgress = FALSE
    )
    
    req <- c("STATIONS_ID", "MESS_DATUM", "QN_4")
    miss <- setdiff(req, names(dt))
    if (length(miss) > 0) {
      message("SKIP (missing cols): ", bn, " -> ", paste(miss, collapse = ", "))
      return(NULL)
    }
    
    dt[, STATIONS_ID := sprintf("%05d", as.integer(STATIONS_ID))]
    dt[, MESS_DATUM := as.Date(as.character(MESS_DATUM), format = "%Y%m%d")]
    dt[, QN_4 := as.integer(QN_4)]
    
    # HARD station-id validation: filename id must match content (wie Script)
    if (any(dt$STATIONS_ID != file_sid, na.rm = TRUE)) {
      bad <- dt[STATIONS_ID != file_sid, .N]
      stop("STATION-ID MISMATCH: filename SID=", file_sid,
           " but file contains other STATIONS_ID values. File=", bn,
           " mismatching rows=", bad)
    }
    
    # date window filter only (wie Script)
    dt <- dt[!is.na(MESS_DATUM) & MESS_DATUM >= start_date & MESS_DATUM <= end_date]
    if (nrow(dt) == 0) return(NULL)
    
    present <- intersect(keep_cols, names(dt))
    dt <- dt[, ..present]
    for (cc in setdiff(keep_cols, names(dt))) dt[, (cc) := NA]
    setcolorder(dt, keep_cols)
    
    # provenance for duplicate diagnostics (wie Script)
    dt[, source := source_tag]
    dt[, file := bn]
    
    dt[]
  }


sanitize_climate_param <- function(sf_day, cVar, date = NULL) {
  x <- as.numeric(sf_day[[cVar]])
  
  # Fehldatenwert -999 zu NA
  x[x == -999] <- NA
  
  if (cVar == "SDK" && !is.null(date)) {
    # Sonnenauf- und untergang → maximale Tageslänge (in Stunden)
    dt <- suncalc::getSunlightTimes(date = as.Date(date), lat = 51.0, lon = 9.0, tz = "UTC")
    td <- dt$sunset - dt$sunrise
    maxDaylight <- ceiling(as.numeric(td))
    x[x > maxDaylight] <- maxDaylight
  }
  
  if (cVar == "PM") {
    x[x > 1060.6] <- 1060.6
    #x[x < 954.9]  <- 954.9
    x[x < 500.0]  <- 500.0
  }
  
  if (cVar == "UPM") {
    x[x > 100] <- 100
    x[x < 0]   <- 0
  }
  
  if (cVar %in% c("TXK", "TNK", "TMK")) {
    x[x > 42]    <- 42
    x[x < -46.0] <- -46.0
  }
  
  if (cVar == "NM") {
    x[x > 8.0] <- 8.0
    x[x < 0]   <- 0
  }
  
  if (cVar == "RSK") {
    x[x > 312] <- 312
    x[x < 0]   <- 0
  }
  
  sf_day$tmp <- x
  return(sf_day)
}


kl_daily_core_prepare <- function(
    # Zeitfenster
  start_date,
  end_date,
  
  # HARTE Regeln wie im direktdownload-Script
  max_missing_frac = 0.25,
  delta_t = 7L,
  
  # DWD URLs (wie direktdownload)
  base_hist   = "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/",
  base_recent = "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/recent/",
  
  # Spaltenschema (wie direktdownload)
  keep_cols = c(
    "STATIONS_ID","MESS_DATUM","QN_4",
    "RSK","RSKF","SDK","SHK_TAG","NM","VPM","PM","TMK","UPM","TXK","TNK","TGK","eor"
  ),
  
  # Pfade (wie direktdownload; nur parametrisierbar, kein anderer Kram)
  root_raw = file.path("data", "data_lev0", "kl_daily", "raw"),
  out_dir  = file.path("data", "data_lev0", "kl_daily", "extracted"),
  
  # Flags (wie direktdownload)
  do_download = FALSE,
  force_unzip = FALSE,
  force_merge = TRUE,
  auto_resolve_duplicates = TRUE,
  
  # Station-Meta-Datei (lokal; exakt diese Datei liest der Script-Teil 9.5)
  meta_file = file.path("data", "data_lev0", "kl_daily", "raw", "KL_Tageswerte_Beschreibung_Stationen.txt"),
  
  # Output merged csv (wie direktdownload)
  write_merged_csv = TRUE
) {
  
  suppressPackageStartupMessages({
    library(data.table)
    library(rvest)
    library(xml2)
  })
  
  # --- Date coercion (minimal, aber zwingend) ---
  if (is.character(start_date)) start_date <- as.Date(start_date)
  if (is.character(end_date))   end_date   <- as.Date(end_date)
  stopifnot(inherits(start_date, "Date"), inherits(end_date, "Date"))
  if (end_date < start_date) stop("end_date < start_date")
  
  # --- Folder structure (1:1 direktdownload) ---
  raw_hist <- file.path(root_raw, "historical")
  raw_rec  <- file.path(root_raw, "recent")
  ex_hist  <- file.path(root_raw, "extracted_historical")
  ex_rec   <- file.path(root_raw, "extracted_recent")
  
  dir.create(raw_hist, showWarnings = FALSE, recursive = TRUE)
  dir.create(raw_rec,  showWarnings = FALSE, recursive = TRUE)
  dir.create(ex_hist,  showWarnings = FALSE, recursive = TRUE)
  dir.create(ex_rec,   showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
  
  n_days_total <- as.integer(end_date - start_date) + 1L
  full_dates   <- seq.Date(start_date, end_date, by = "day") # (wie Script)
  
  pat_prod <- "^produkt_klima_tag_\\d{8}_\\d{8}_\\d{5}\\.txt$"
  pat_zip  <- "\\.zip$"
  
  out_merged_csv_gz <- file.path(out_dir, sprintf(
    "kl_daily_merged_%s_%s_stationfilter_miss%02d_gap%d.csv",
    format(start_date, "%Y%m%d"),
    format(end_date, "%Y%m%d"),
    as.integer(max_missing_frac * 100),
    as.integer(delta_t)
  ))
  
  # diagnostics (wie Script)
  prod_files_csv <- file.path(out_dir, "diag_product_files_found.csv")
  qn4_global_csv <- file.path(out_dir, "diag_qn4_global_counts.csv")
  qn4_station_csv<- file.path(out_dir, "diag_qn4_by_station.csv")
  dup_keys_csv   <- file.path(out_dir, "diag_duplicates_station_date.csv")
  dup_rows_csv   <- file.path(out_dir, "diag_duplicates_rows_detail.csv")
  st_diag_csv    <- file.path(out_dir, "diag_station_window_missing_gap.csv")
  
  # PRE-FLIGHT (wie Script)
  if (!isTRUE(force_merge) && file.exists(out_merged_csv_gz)) {
    message("FINAL OUTPUT EXISTS → SKIP MERGE (set force_merge=TRUE): ", out_merged_csv_gz)
    stop("force_merge=FALSE and final output exists (script would quit here).")
  }
  
  # ============================================================
  # HELPERS (1:1 aus direktdownload)
  # ============================================================
  list_remote_files <- function(url, pattern = NULL) {
    doc <- rvest::read_html(url)
    files <- doc |> rvest::html_elements("a") |> rvest::html_attr("href")
    files <- files[!is.na(files) & files != "../"]
    if (!is.null(pattern)) files <- files[grepl(pattern, files)]
    files <- unique(files)
    data.table(url = paste0(url, files), file = files)
  }
  
  download_if_missing <- function(url, destfile) {
    if (file.exists(destfile)) return(invisible(FALSE))
    message("download zip: ", destfile)
    tryCatch({
      utils::download.file(url, destfile = destfile, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) {
      message("DOWNLOAD FAILED: ", url)
      message(conditionMessage(e))
      FALSE
    })
  }
  
  unzip_zip <- function(zipfile, exdir) {
    suppressWarnings(utils::unzip(zipfile, exdir = exdir))
    invisible(TRUE)
  }
  
  parse_prod_span <- function(files) {
    bn <- basename(files)
    m <- regexec("^produkt_klima_tag_(\\d{8})_(\\d{8})_(\\d{5})\\.txt$", bn)
    mm <- regmatches(bn, m)
    ok <- lengths(mm) == 4
    if (!all(ok)) stop("Unexpected product filenames encountered.")
    data.table(
      STATIONS_ID = vapply(mm, `[[`, character(1), 4),
      file_start  = as.Date(vapply(mm, `[[`, character(1), 2), "%Y%m%d"),
      file_end    = as.Date(vapply(mm, `[[`, character(1), 3), "%Y%m%d"),
      file        = bn
    )
  }
  
  read_one_product_file <- function(f, source_tag) {
    bn <- basename(f)
    m <- regexec("^produkt_klima_tag_(\\d{8})_(\\d{8})_(\\d{5})\\.txt$", bn)
    mm <- regmatches(bn, m)[[1]]
    if (length(mm) != 4) stop("Unexpected product filename: ", bn)
    
    file_sid <- mm[4]
    
    dt <- data.table::fread(
      f,
      sep = ";",
      na.strings = "-999",
      strip.white = TRUE,
      showProgress = FALSE
    )
    
    req <- c("STATIONS_ID", "MESS_DATUM", "QN_4")
    miss <- setdiff(req, names(dt))
    if (length(miss) > 0) {
      message("SKIP (missing cols): ", bn, " -> ", paste(miss, collapse = ", "))
      return(NULL)
    }
    
    dt[, STATIONS_ID := sprintf("%05d", as.integer(STATIONS_ID))]
    dt[, MESS_DATUM := as.Date(as.character(MESS_DATUM), format = "%Y%m%d")]
    dt[, QN_4 := as.integer(QN_4)]
    
    # HARD station-id validation: filename id must match content (wie Script)
    if (any(dt$STATIONS_ID != file_sid, na.rm = TRUE)) {
      bad <- dt[STATIONS_ID != file_sid, .N]
      stop("STATION-ID MISMATCH: filename SID=", file_sid,
           " but file contains other STATIONS_ID values. File=", bn,
           " mismatching rows=", bad)
    }
    
    # date window filter only (wie Script)
    dt <- dt[!is.na(MESS_DATUM) & MESS_DATUM >= start_date & MESS_DATUM <= end_date]
    if (nrow(dt) == 0) return(NULL)
    
    present <- intersect(keep_cols, names(dt))
    dt <- dt[, ..present]
    for (cc in setdiff(keep_cols, names(dt))) dt[, (cc) := NA]
    setcolorder(dt, keep_cols)
    
    # provenance for duplicate diagnostics (wie Script)
    dt[, source := source_tag]
    dt[, file := bn]
    
    dt[]
  }
  
  # ============================================================
  # 1) Remote list zips (wie Script)
  # ============================================================
  z_hist <- list_remote_files(base_hist, pat_zip)
  z_rec  <- list_remote_files(base_recent, pat_zip)
  
  message("zip count historical: ", nrow(z_hist))
  message("zip count recent:     ", nrow(z_rec))
  
  z_hist[, dest := file.path(raw_hist, file)]
  z_rec[,  dest := file.path(raw_rec,  file)]
  
  # ============================================================
  # 2) Download (wie Script)
  # ============================================================
  if (isTRUE(do_download)) {
    miss_hist <- z_hist[!file.exists(dest)]
    miss_rec  <- z_rec [!file.exists(dest)]
    
    message("ZIPs missing (historical): ", nrow(miss_hist))
    message("ZIPs missing (recent):     ", nrow(miss_rec))
    
    for (i in seq_len(nrow(miss_hist))) download_if_missing(miss_hist$url[i], miss_hist$dest[i])
    for (i in seq_len(nrow(miss_rec)))  download_if_missing(miss_rec$url[i],  miss_rec$dest[i])
  } else {
    message("DOWNLOAD DISABLED (do_download=FALSE) → using local ZIPs only")
  }
  
  missing_hist <- z_hist[!file.exists(dest)]
  missing_rec  <- z_rec [!file.exists(dest)]
  if (nrow(missing_hist) + nrow(missing_rec) > 0) {
    stop(
      "Missing required ZIPs locally. Set do_download=TRUE or place zips in:\n",
      "  ", raw_hist, "  (historical)\n",
      "  ", raw_rec,  "  (recent)\n",
      "Missing counts: hist=", nrow(missing_hist), " recent=", nrow(missing_rec)
    )
  }
  
  # ============================================================
  # 3) Unzip (wie Script)
  # ============================================================
  have_prod_hist <- length(list.files(ex_hist, pattern = pat_prod, recursive = TRUE)) > 0
  have_prod_rec  <- length(list.files(ex_rec,  pattern = pat_prod, recursive = TRUE)) > 0
  
  if (isTRUE(force_unzip) || !have_prod_hist) {
    message("UNZIP historical -> ", ex_hist)
    for (zf in z_hist$dest) unzip_zip(zf, exdir = ex_hist)
  } else {
    message("UNZIP historical skipped")
  }
  
  if (isTRUE(force_unzip) || !have_prod_rec) {
    message("UNZIP recent -> ", ex_rec)
    for (zf in z_rec$dest) unzip_zip(zf, exdir = ex_rec)
  } else {
    message("UNZIP recent skipped")
  }
  
  # ============================================================
  # 4) Find product files (+ span table from filenames) (wie Script)
  # ============================================================
  prod_hist <- list.files(ex_hist, pattern = pat_prod, full.names = TRUE, recursive = TRUE)
  prod_rec  <- list.files(ex_rec,  pattern = pat_prod, full.names = TRUE, recursive = TRUE)
  
  message("product txt count hist:   ", length(prod_hist))
  message("product txt count recent: ", length(prod_rec))
  
  if (length(prod_hist) + length(prod_rec) == 0) {
    stop("No produkt_klima_tag_*.txt found under: ", ex_hist, " or ", ex_rec)
  }
  
  data.table::fwrite(
    rbind(
      data.table(path = prod_hist, source = "historical"),
      data.table(path = prod_rec,  source = "recent")
    ),
    prod_files_csv, sep = ";"
  )
  
  span_all <- data.table::rbindlist(list(parse_prod_span(prod_hist), parse_prod_span(prod_rec)), use.names = TRUE)
  st_span_files <- span_all[, .(
    min_file_start = min(file_start, na.rm = TRUE),
    max_file_end   = max(file_end,   na.rm = TRUE),
    n_files        = .N
  ), by = STATIONS_ID]
  
  # ============================================================
  # 5) Read all products (wie Script)
  # ============================================================
  res_hist <- vector("list", length(prod_hist))
  for (i in seq_along(prod_hist)) {
    if (i %% 250 == 0) message("read(hist) ", i, " / ", length(prod_hist))
    res_hist[[i]] <- read_one_product_file(prod_hist[i], "historical")
  }
  
  res_rec <- vector("list", length(prod_rec))
  for (i in seq_along(prod_rec)) {
    if (i %% 250 == 0) message("read(recent) ", i, " / ", length(prod_rec))
    res_rec[[i]] <- read_one_product_file(prod_rec[i], "recent")
  }
  
  dt <- data.table::rbindlist(c(res_hist, res_rec), use.names = TRUE, fill = TRUE)
  rm(res_hist, res_rec); gc()
  
  if (nrow(dt) == 0) stop("No rows after reading+date filter.")
  data.table::setDT(dt)  # **wichtig**: garantiert data.table
  
  # ============================================================
  # 6) QN_4 diagnostics (wie Script)
  # ============================================================
  data.table::fwrite(dt[, .N, by = QN_4][order(QN_4)], qn4_global_csv, sep=";")
  data.table::fwrite(dt[, .N, by = .(STATIONS_ID, QN_4)][order(STATIONS_ID, QN_4)], qn4_station_csv, sep=";")
  
  # ============================================================
  # 7) Duplicate check (wie Script)
  # ============================================================
  dup_keys <- dt[, .N, by = .(STATIONS_ID, MESS_DATUM)][N > 1L]
  if (nrow(dup_keys) > 0) {
    data.table::fwrite(dup_keys, dup_keys_csv, sep = ";")
    
    setkey(dt, STATIONS_ID, MESS_DATUM)
    dup_detail <- dt[dup_keys, on = .(STATIONS_ID, MESS_DATUM)]
    data.table::fwrite(dup_detail, dup_rows_csv, sep = ";")
    
    if (!isTRUE(auto_resolve_duplicates)) {
      stop("Stop: duplicates exist. See diagnostics in extracted/. Set auto_resolve_duplicates=TRUE to resolve.")
    }
    
    val_cols <- setdiff(keep_cols, c("STATIONS_ID","MESS_DATUM","QN_4","eor"))
    dt[, nnz := rowSums(!is.na(.SD)), .SDcols = val_cols]
    setorder(dt, STATIONS_ID, MESS_DATUM, -nnz, -QN_4)
    dt <- dt[, .SD[1], by = .(STATIONS_ID, MESS_DATUM)]
    dt[, nnz := NULL]
  }
  
  # ============================================================
  # 8) Station filter (wie Script)
  # ============================================================
  dt[, MESS_DATUM := as.Date(MESS_DATUM)]
  dt_dates <- unique(dt[, .(STATIONS_ID, MESS_DATUM)])
  setorder(dt_dates, STATIONS_ID, MESS_DATUM)
  
  st_cov <- dt_dates[, .(n_dates = .N, min_date = min(MESS_DATUM), max_date = max(MESS_DATUM)), by = STATIONS_ID]
  st_cov[, expected := n_days_total]
  st_cov[, coverage_frac := n_dates / expected]
  st_cov[, missing_frac  := 1 - coverage_frac]
  
  st_gap <- dt_dates[, {
    d <- MESS_DATUM
    gap_before <- as.integer(min(d) - start_date)
    gap_after  <- as.integer(end_date - max(d))
    dif <- diff(d)
    internal <- if (length(dif)) pmax(0L, as.integer(dif) - 1L) else integer()
    list(max_gap = max(c(gap_before, internal, gap_after), na.rm = TRUE))
  }, by = STATIONS_ID]
  
  st_diag <- merge(st_cov, st_span_files, by = "STATIONS_ID", all.x = TRUE)
  st_diag <- merge(st_diag, st_gap,       by = "STATIONS_ID", all.x = TRUE)
  
  data.table::fwrite(st_diag[order(-missing_frac, -max_gap)], st_diag_csv, sep=";")
  
  st_keep <- st_diag[
    min_file_start <= start_date &
      max_file_end   >= end_date &
      missing_frac   <= max_missing_frac &
      max_gap        <= delta_t,
    STATIONS_ID
  ]
  
  dt <- dt[STATIONS_ID %in% st_keep]
  if (nrow(dt) == 0) stop("After station filter no stations remain. See diag_station_window_missing_gap.csv")
  
  # ============================================================
  # 9) HARD POST-CHECKS (wie Script)
  # ============================================================
  stopifnot(all(nchar(dt$STATIONS_ID) == 5))
  stopifnot(all(grepl("^[0-9]{5}$", dt$STATIONS_ID)))
  stopifnot(dt[, anyDuplicated(paste(STATIONS_ID, MESS_DATUM))] == 0L)
  
  dt_dates2 <- unique(dt[, .(STATIONS_ID, MESS_DATUM)])
  st_gap2 <- dt_dates2[, {
    d <- sort(MESS_DATUM)
    gap_before <- as.integer(min(d) - start_date)
    gap_after  <- as.integer(end_date - max(d))
    dif <- diff(d)
    internal <- if (length(dif)) pmax(0L, as.integer(dif) - 1L) else integer()
    list(max_gap = max(c(gap_before, internal, gap_after), na.rm = TRUE),
         n_dates  = .N)
  }, by = STATIONS_ID]
  
  stopifnot(all(st_gap2$max_gap <= delta_t))
  stopifnot(all((1 - st_gap2$n_dates / n_days_total) <= max_missing_frac))
  
  # ============================================================
  # 9.5) Merge station metadata (wie Script-Block in direktdownload)
  # ============================================================
  suppressPackageStartupMessages({
    library(sf)
  })
  
  meta_file <- path.expand(meta_file)
  stopifnot(file.exists(meta_file))
  
  ln <- readLines(meta_file, warn = FALSE)
  x  <- ln[-c(1,2)]
  x  <- x[nzchar(trimws(x))]
  
  st_meta <- data.table::rbindlist(lapply(x, function(s){
    m <- regexec(
      "^\\s*([0-9]{5})\\s+([0-9]{8})\\s+([0-9]{8})\\s+(-?[0-9]+)\\s+([0-9]+\\.[0-9]+)\\s+([0-9]+\\.[0-9]+)\\s*(.*)$",
      s
    )
    z <- regmatches(s, m)[[1]]
    if (length(z) != 8) return(NULL)
    
    rest <- trimws(z[8])
    p <- trimws(unlist(strsplit(rest, "\\s{2,}", perl = TRUE)))
    p <- p[nzchar(p)]
    
    list(
      STATIONS_ID   = sprintf("%05d", as.integer(z[2])),
      von_datum     = as.Date(z[3], "%Y%m%d"),
      bis_datum     = as.Date(z[4], "%Y%m%d"),
      Stationshoehe = as.numeric(z[5]),
      geoBreite     = as.numeric(z[6]),
      geoLaenge     = as.numeric(z[7]),
      Stationsname  = if (length(p) >= 1) p[1] else NA_character_,
      Bundesland    = if (length(p) >= 2) p[2] else NA_character_,
      Abgabe        = if (length(p) >= 3) p[3] else NA_character_
    )
  }), fill = TRUE)
  
  fix_enc <- function(v) iconv(v, to = "UTF-8", sub = "")
  st_meta[, `:=`(
    Stationsname = fix_enc(Stationsname),
    Bundesland   = fix_enc(Bundesland),
    Abgabe       = fix_enc(Abgabe)
  )]
  
  st_meta_u <- st_meta[, .(
    Stationsname  = Stationsname[1],
    Bundesland    = Bundesland[1],
    Abgabe        = Abgabe[1],
    Stationshoehe = Stationshoehe[1],
    geoBreite     = geoBreite[1],
    geoLaenge     = geoLaenge[1]
  ), by = STATIONS_ID]
  
  st_span_u <- st_meta[, .(
    von_datum = as.integer(format(min(von_datum, na.rm = TRUE), "%Y%m%d")),
    bis_datum = as.integer(format(max(bis_datum, na.rm = TRUE), "%Y%m%d"))
  ), by = STATIONS_ID]
  
  setkey(st_span_u, STATIONS_ID)
  setkey(dt,        STATIONS_ID)
  dt <- st_span_u[dt]
  
  if (anyNA(dt$von_datum) || anyNA(dt$bis_datum)) {
    bad <- unique(dt[is.na(von_datum) | is.na(bis_datum), STATIONS_ID])
    stop("Missing von_datum/bis_datum for STATIONS_ID: ", paste(bad, collapse = ", "))
  }
  
  setkey(st_meta_u, STATIONS_ID)
  dt[, STATIONS_ID := sprintf("%05d", as.integer(STATIONS_ID))]
  setkey(dt, STATIONS_ID)
  dt <- st_meta_u[dt]
  
  # ============================================================
  # 10) Output dt_legacy (wie Script)
  # ============================================================
  dt[, c("source","file") := NULL]
  setorder(dt, STATIONS_ID, MESS_DATUM)
  setcolorder(dt, keep_cols)
  
  if (isTRUE(write_merged_csv)) {
    data.table::fwrite(dt, out_merged_csv_gz, sep = ";", na = "")
    message("merged saved: ", out_merged_csv_gz)
  }
  
  dt
}


kl_daily_extract_sf <- function(
    dt = NULL,
    
    # chain switch: TRUE => core wird mit do_download=TRUE ausgeführt
    force_download_chain = FALSE,
    
    # welche Param-Spalte als Wert in sf
    param = "PM",
    
    # core args (nur relevant wenn dt NULL oder force_download_chain=TRUE)
    start_date = as.Date("2003-01-01"),
    end_date   = as.Date("2024-12-31"),
    max_missing_frac = 0.25,
    delta_t = 7L,
    base_hist   = "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/",
    base_recent = "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/recent/",
    keep_cols = c(
      "STATIONS_ID","MESS_DATUM","QN_4",
      "RSK","RSKF","SDK","SHK_TAG","NM","VPM","PM","TMK","UPM","TXK","TNK","TGK","eor"
    ),
    root_raw = file.path("data", "data_lev0", "kl_daily", "raw"),
    out_dir  = file.path("data", "data_lev0", "kl_daily", "extracted"),
    meta_file = file.path("data", "data_lev0", "kl_daily", "raw", "KL_Tageswerte_Beschreibung_Stationen.txt"),
    do_download = FALSE,
    force_unzip = FALSE,
    force_merge = TRUE,
    auto_resolve_duplicates = TRUE,
    write_merged_csv = TRUE
) {
  
  suppressPackageStartupMessages({
    library(data.table)
    library(sf)
  })
  
  if (isTRUE(force_download_chain) || is.null(dt)) {
    dt <- kl_daily_core_prepare(
      start_date = start_date,
      end_date   = end_date,
      max_missing_frac = max_missing_frac,
      delta_t = delta_t,
      base_hist = base_hist,
      base_recent = base_recent,
      keep_cols = keep_cols,
      root_raw = root_raw,
      out_dir  = out_dir,
      do_download = if (isTRUE(force_download_chain)) TRUE else isTRUE(do_download),
      force_unzip = force_unzip,
      force_merge = force_merge,
      auto_resolve_duplicates = auto_resolve_duplicates,
      meta_file = meta_file,
      write_merged_csv = write_merged_csv
    )
  }
  
  data.table::setDT(dt) # garantiert data.table
  
  need <- c("STATIONS_ID","MESS_DATUM", param,
            "Stationshoehe","Stationsname","Bundesland","Abgabe",
            "geoBreite","geoLaenge")
  stopifnot(all(need %chin% names(dt)))
  
  dt_legacy <- data.table::copy(dt)
  
  dt_legacy[, STATIONS_ID   := as.numeric(STATIONS_ID)]
  dt_legacy[, Stationshoehe := as.integer(Stationshoehe)]
  dt_legacy[, MESS_DATUM    := as.Date(MESS_DATUM)]
  
  if (!("von_datum" %chin% names(dt_legacy))) dt_legacy[, von_datum := NA_integer_]
  if (!("bis_datum" %chin% names(dt_legacy))) dt_legacy[, bis_datum := NA_integer_]
  dt_legacy[, von_datum := as.integer(von_datum)]
  dt_legacy[, bis_datum := as.integer(bis_datum)]
  
  keep <- c("STATIONS_ID","von_datum","bis_datum","Stationshoehe",
            "Stationsname","Bundesland","Abgabe",
            "MESS_DATUM", param, "geoLaenge","geoBreite")
  dt_legacy <- dt_legacy[, ..keep]
  
  sf_legacy <- sf::st_as_sf(
    dt_legacy,
    coords = c("geoLaenge","geoBreite"),
    crs = 4326,
    remove = TRUE
  )
  sf_legacy$STATIONS_ID <- as.numeric(sf_legacy$STATIONS_ID)
  sf_legacy$geom <- sf::st_geometry(sf_legacy)
  sf_legacy <- sf::st_set_geometry(sf_legacy, "geom")
  if ("geometry" %in% names(sf_legacy)) sf_legacy$geometry <- NULL
  
  sf_legacy
}
