split_mz_data = function(poly= NULL,mz_data=NULL,state_id = NULL,out_id =NULL){
for (s_id in state_id){
  # crop to state
  if (!file.exists(paste0(envrmt$path_data_lev1,"/",s_id,"_state.rds")))
  {state = st_crop( poly %>% filter(substr(AGS,1,2)==s_id),poly)} else {
    state= readRDS(paste0(envrmt$path_data_lev1,"/",s_id,"_state.rds"))
  }
  if (!file.exists(paste0(envrmt$path_data_lev1,"/",s_id,"_",out_id,"_state.rds")))
  mz_state = st_crop(mz_data,state)
  if (!file.exists(paste0(envrmt$path_data_lev1,"/",s_id,"_",out_id,"_state.rds")))
  saveRDS(mz_state,paste0(envrmt$path_data_lev1,"/",s_id,"_",out_id,"_state.rds"))
  if (!file.exists(paste0(envrmt$path_data_lev1,"/",s_id,"_state.rds")))
  saveRDS(state,paste0(envrmt$path_data_lev1,"/",s_id,"_state.rds"))
}}

cv.IDW  <- function(spatialDF, stat.formula = NULL,
                    seqNeighbors = NULL, seqBeta = NULL,
                    evalGridSize = NULL,
                    evalRaster = NULL,
                    verbose = TRUE){

  ### LOAD LIBRARIES ###
  library(sp)
  library(sf)
  library(gstat)
  library(raster)


  ### PROVIDE DEFAULT VALUES FOR FUNCTION ARGUMENTS ###
  if (is.null(seqNeighbors)){
    seqNeighbors <- round(seq(3, length(spatialDF), length.out = 10))
  }
  if (is.null(seqBeta)){
    seqBeta <- c(0.1, seq(0.5, 3, 0.5))
  }
  # if (is.null(evalGridSize)){
  #   x.interval <- extent.DE["xmax"] - extent.DE["xmin"]
  #   y.interval <- extent.DE["ymax"] - extent.DE["ymin"]
  #   evalGridSize <- round(min(x.interval, y.interval) *0.05)
  # }
  if (is.null(stat.formula)){
    print('Please provide a formula!!')
    return()
  }
  # if (is.null(evalRaster)){
  #   extent.evalGrid <- raster::extent(spatialDF)
  # }else{
  #   extent.evalGrid <- raster::extent(evalRaster)
  # }


  ### BUILD A GRID FOR PARAMETER COMBINATIONS ###
  cv.Grid <- expand.grid(Beta = seqBeta,
                         Neighbors = seqNeighbors)
  cv.Grid$RMSE <- NA
  spatialDF = spTransform(spatialDF,crs)
  ### LOOP THROUGH ALL PARAMETER COMBINATIONS ###
  for (i in 1:nrow(cv.Grid)){
    ### BUILD IDW MODEL ###
    idw <- gstat(formula = stat.formula,
                 data = spatialDF,
                 nmax = cv.Grid[i, 'Neighbors'],
                 set = list(idp = cv.Grid[i, 'Beta']))
    ### PERFORM LOOCV ###
    crossval <- gstat.cv(idw,
                         nmax = cv.Grid[i, 'Neighbors'],
                         beta = v.Grid[i, 'Beta'],
                         debug.level = 0)
    cv.Grid[i, 'RMSE'] <- RMSE(crossval$residual)
    if (verbose){
      print(paste('Function call', i, 'out of',  nrow(cv.Grid)))
      print(paste('Evaluating beta =',
                  cv.Grid[i, 'Beta'],
                  'and neighbors =',
                  cv.Grid[i, 'Neighbors']))
      print(paste('RMSE=', RMSE(crossval$residual)))
    }
  }

  ### GET BEST PARAMTER VALUES ###
  idx.min <- which.min(cv.Grid$RMSE)
  best.Beta <- cv.Grid$Beta[idx.min]
  best.Neighbors <- cv.Grid$Neighbors[idx.min]
  min.RMSE <- cv.Grid$RMSE[idx.min]

  ### BUILD IDW MODEL BASED ON BEST PARAMTER VALUES ###
  idw.best <- gstat(formula = stat.formula,
                    data = spatialDF,
                    nmax = best.Neighbors,
                    set = list(idp = best.Beta))

  ### PREPARE EVALUATION GRID ###
  # grid.evalGrid  <- expand.grid(x = seq(from = round(extent.evalGrid@xmin),
  #                                       to = round(extent.evalGrid@xmax),
  #                                       by = evalGridSize),
  #                               y = seq(from = round(extent.evalGrid@ymin),
  #                                       to = round(extent.evalGrid@ymax),
  #                                       by = evalGridSize))
  #
  #
  # coordinates(grid.evalGrid) <- ~x + y
  # proj4string(grid.evalGrid) <- proj4string(spatialDF)
  grid.evalGrid = grid.DE
  gridded(grid.evalGrid) <- TRUE
  idw.best$data$var1$data = spTransform(idw.best$data$var1$data,crs)

  ### INTERPOLATE VALUES FOR EVALUATION GRID USING THE BEST MODEL ###
  idw.best.predict <- predict(object = idw.best$data,
                              newdata = grid.evalGrid,
                              debug.level = 0)

  ### RETURN RESULTS AND OBJECTS ###
  return(list('idwBestModel' = idw.best,
              'idwBestRaster' = idw.best.predict,
              'bestBeta' = best.Beta,
              'bestNeighbors' = best.Neighbors,
              'bestRMSE' = min.RMSE,
              'gridCV' = cv.Grid))
}

RMSE <- function(residuals){
  sqrt(sum((residuals)^2)/length(residuals))
}





