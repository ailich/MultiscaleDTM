#' Calculates surface area of a DTM
#'
#' Calculates surface area on a per cell basis of a DTM based on Jenness, 2004.
#' @param r DTM as a SpatRaster or RasterLayer in a projected coordinate system where map units match elevation/depth units
#' @param na.rm Logical indicating whether to remove NAs from calculations. When FALSE, the sum of the eight triangles is calculated. When TRUE, the mean of the created triangles is calculated and multiplied by 8 to scale it to the proper area.
#' @param filename character Output filename.
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @return a SpatRaster or RasterLayer
#' @examples 
#' r<- rast(volcano, extent= ext(2667400, 2667400 + 
#' ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
#' crs = "EPSG:27200")
#' sa<- SurfaceArea(r)
#' plot(sa)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @references
#' Jenness, J.S., 2004. Calculating landscape surface area from digital elevation models. Wildlife Society Bulletin 32, 829-839.
#' @export

SurfaceArea<- function(r, na.rm=FALSE, filename=NULL, overwrite=FALSE, wopt=list()){
  og_class<- class(r)[1]
  if(og_class=="RasterLayer"){
    r<- terra::rast(r) #Convert to SpatRaster
  }
  # Input checks
  if(!(og_class %in% c("RasterLayer", "SpatRaster"))){
    stop("Error: Input must be a 'SpatRaster' or 'RasterLayer'")
  }
  if(terra::nlyr(r)!=1){
    stop("Error: Input raster must be one layer.")
  }
  if(isTRUE(terra::is.lonlat(r, perhaps=FALSE))){
    stop("Error: Coordinate system is Lat/Lon. Coordinate system must be projected with elevation/depth units matching map units.")
  }
  if(terra::is.lonlat(r, perhaps=TRUE, warn=FALSE)){
    warning("Coordinate system may be Lat/Lon. Please ensure that the coordinate system is projected with elevation/depth units matching map units.")
  }
  SA<- terra::focalCpp(r, w=c(3,3), fun = C_SurfaceArea,  x_res = terra::res(r)[1], y_res = terra::res(r)[2], na_rm=na.rm, fillvalue=NA, wopt=wopt)
  names(SA)<- "SA"
  #Return
  if(og_class =="RasterLayer"){
    SA<- raster::raster(SA)
    if(!is.null(filename)){
      return(raster::writeRaster(SA, filename=filename, overwrite=overwrite))
    }
  }
  if(!is.null(filename)){
    return(terra::writeRaster(SA, filename=filename, overwrite=overwrite, wopt=wopt))
  }
  return(SA)
  }