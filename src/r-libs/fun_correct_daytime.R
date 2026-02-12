correct_daytime = function(fn){
  for (f in fn){
    current = terra::rast(f)
    dt=suncalc::getSunlightTimes(date = as.Date(substr(basename(f),1,10)), lat = 51.0, lon = 9.0, tz = "UTC")
    td=dt$sunset-dt$sunrise
    maxDaylight= ceiling(as.numeric(unlist(stringr::str_split(td,"Time difference of "))))
    m <- c(-999, 0,0, maxDaylight,99999,maxDaylight)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    current <- terra::classify(current, rclmat, include.lowest=TRUE)
    names(current) = xfun::sans_ext(basename(f))
    if (calc_bl) 
    {  
      # Calculate data frame of min and max precipitation for all months
      var <- cbind(bl_sf, exactextractr::exact_extract(terra::rast(f), bl_sf, c("min", "max","count","majority","median","quantile","minority","variance","stdev","coefficient_of_variation"),quantiles = c(0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.25,0.5,0.75,0.9)))
      var$date = substr(tools::file_path_sans_ext(basename(f)),1,10)
      #c("min",	"max",	"count",	"majority",	"median",	"q10",	"q20",	"q30",	"q40",	"q60",	"q70",	"q80",	"q25",	"q50",	"q75",	"q90",	"minority","variance","stdev","coefficient_of_variation")
      vr=st_drop_geometry(var[,c("min",	"max",	"count",	"majority",	"median",	"q10",	"q20",	"q30",	"q40",	"q60",	"q70",	"q80",	"q25",	"q50",	"q75",	"q90",	"minority")]) %>%
        mutate_if(is.numeric, round, digits=dig)
      var_fin=sjmisc::replace_columns(var,vr)
      #saveRDS(var,file.path(envrmt$path_data_lev2,cVar,paste0(tools::file_path_sans_ext(basename(clim_files[i])),".rds")))
      data.table::fwrite(st_drop_geometry(var_fin),file=file.path(envrmt$path_data_lev2,bl,cVar,paste0(tools::file_path_sans_ext(basename(clim_files[i])),".csv")),dec = ".")
      current = terra::mask(raster::raster(f), bl_sf)
      current = terra::crop(current,bl_sf)
      raster::writeRaster(current,file.path(envrmt$path_data_lev1,bl,cVar,paste0(bl,basename(f))),overwrite=TRUE)
    }  else {
      writeRaster(current, f,gdal=c("COMPRESS=NONE"),overwrite=TRUE)
    }
  }
}
