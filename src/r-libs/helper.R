#' Sample training points from a categorical CLC raster (optionally extract predictors)
#'
#' @description
#' Creates a point-based training set by drawing (stratified) random samples from
#' a categorical land-cover raster (e.g., CLC). The class label is taken directly
#' from the raster cell value at each sampled location. Optionally, the function
#' extracts predictor values (multi-band SpatRaster) at the same locations so the
#' result can be used directly for Random Forest training.
#'
#' This is a pragmatic bootstrapping method: it produces internally consistent
#' labels (from CLC) but inherits CLC's scale/accuracy limitations. Use it as
#' baseline training, not as final ground truth.
#'
#' @param clc A categorical `terra::SpatRaster` (one layer) with integer class codes.
#' @param n_per_class Integer. Number of points to sample per class.
#' @param classes Optional integer vector of class codes to sample. Default: all non-NA classes.
#' @param aoi Optional AOI polygon (`sf` or `SpatVector`). If provided, sampling is restricted to AOI.
#' @param predictors Optional `terra::SpatRaster` (multi-layer). Predictor values are extracted to points.
#' @param min_dist_m Optional numeric. If set, enforces an approximate minimum distance between points (meters).
#'                   Implemented via greedy thinning; may reduce the final sample size.
#' @param seed Integer seed for reproducibility.
#' @param drop_na_predictors Logical. If TRUE and predictors are provided, drop points with any NA predictor.
#' @return An `sf` POINT object with columns:
#'         - `class` (integer class code from CLC)
#'         - optionally predictor columns (if `predictors` is given)
#'
#' @examples
#' # pts <- sample_training_points_from_clc(clc, n_per_class = 200, aoi = aoi_sf)
#' # pts <- sample_training_points_from_clc(clc, 200, aoi_sf, predictors = s2_stack_10m)
sample_training_points_from_clc <- function(clc,
                                            n_per_class,
                                            classes = NULL,
                                            aoi = NULL,
                                            predictors = NULL,
                                            min_dist_m = NULL,
                                            seed = 1L,
                                            drop_na_predictors = TRUE) {
  stopifnot(inherits(clc, "SpatRaster"))
  stopifnot(terra::nlyr(clc) == 1)
  stopifnot(is.numeric(n_per_class) && length(n_per_class) == 1 && n_per_class > 0)
  
  if (!is.null(aoi)) {
    aoi_v <- if (inherits(aoi, "SpatVector")) aoi else terra::vect(aoi)
    # align AOI CRS to raster CRS
    if (!terra::same.crs(aoi_v, clc)) aoi_v <- terra::project(aoi_v, terra::crs(clc))
    clc <- terra::mask(terra::crop(clc, aoi_v), aoi_v)
  }
  
  # determine classes
  if (is.null(classes)) {
    vals <- terra::values(clc, mat = FALSE)
    classes <- sort(unique(vals[is.finite(vals)]))
  } else {
    classes <- sort(unique(as.integer(classes)))
  }
  if (length(classes) == 0) stop("No classes found to sample (CLC is empty or all NA).")
  
  set.seed(as.integer(seed))
  
  pts_list <- vector("list", length(classes))
  names(pts_list) <- as.character(classes)
  
  for (k in seq_along(classes)) {
    cls <- classes[k]
    m <- clc == cls
    
    # sample n_per_class points from mask; may return fewer if class is small
    p <- terra::spatSample(
      x = m,
      size = as.integer(n_per_class),
      method = "random",
      as.points = TRUE,
      na.rm = TRUE,
      values = FALSE
    )
    
    if (is.null(p) || terra::nrow(p) == 0) {
      pts_list[[k]] <- NULL
      next
    }
    
    # attach class label
    p$class <- cls
    pts_list[[k]] <- p
  }
  
  pts_v <- do.call(rbind, pts_list)
  if (is.null(pts_v) || terra::nrow(pts_v) == 0) stop("Sampling returned zero points.")
  
  # optional thinning by minimum distance (greedy; approximate)
  if (!is.null(min_dist_m)) {
    stopifnot(is.numeric(min_dist_m) && min_dist_m > 0)
    pts_sf <- sf::st_as_sf(pts_v)
    # ensure projected CRS for meter distance
    if (sf::st_is_longlat(pts_sf)) {
      stop("min_dist_m requires a projected CRS (meters). Reproject clc/aoi to a metric CRS first.")
    }
    keep <- rep(TRUE, nrow(pts_sf))
    for (i in seq_len(nrow(pts_sf))) {
      if (!keep[i]) next
      # mark all later points within distance as FALSE
      d <- sf::st_is_within_distance(pts_sf[i, ], pts_sf, dist = min_dist_m, sparse = TRUE)[[1]]
      d <- d[d > i]
      keep[d] <- FALSE
    }
    pts_sf <- pts_sf[keep, , drop = FALSE]
    pts_v <- terra::vect(pts_sf)
  }
  
  # optional predictor extraction
  if (!is.null(predictors)) {
    stopifnot(inherits(predictors, "SpatRaster"))
    # align predictors to clc CRS if needed
    if (!terra::same.crs(predictors, clc)) predictors <- terra::project(predictors, terra::crs(clc))
    
    pred_df <- terra::extract(predictors, pts_v)
    # terra::extract returns an ID column; drop it
    if ("ID" %in% names(pred_df)) pred_df$ID <- NULL
    
    # bind to points
    pts_df <- as.data.frame(pts_v)
    pts_df <- cbind(pts_df, pred_df)
    
    if (drop_na_predictors) {
      ok <- stats::complete.cases(pts_df[, names(pred_df), drop = FALSE])
      pts_df <- pts_df[ok, , drop = FALSE]
      pts_v  <- pts_v[ok]
    }
    
    # re-attach attributes
    pts_v <- terra::vect(pts_v)
    pts_v <- terra::setValues(pts_v, pts_df[, setdiff(names(pts_df), c("x","y")), drop = FALSE])
  }
  
  # return sf
  pts_sf <- sf::st_as_sf(pts_v)
  # explicit class column
  if (!("class" %in% names(pts_sf))) stop("Internal error: class column missing.")
  pts_sf
}
