#' Calculates surface area of a DEM
#'
#' Calculates surface area on a per cell basis of a DEM based on Jenness, 2004.
#' @param r DEM as a SpatRaster or RasterLayer in a projected coordinate system where map units match elevation/depth units
#' @param filename character Output filename.
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @return a SpatRaster or RasterLayer
#' @import terra
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @references
#' Jenness, J.S., 2004. Calculating landscape surface area from digital elevation models. Wildlife Society Bulletin 32, 829-839.
#' @export

SurfaceArea<- function(r, filename=NULL, overwrite=FALSE){
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
  SA<- terra::focalCpp(r, w=c(3,3), fun = C_SurfaceArea,  x_res = terra::res(r)[1], y_res = terra::res(r)[2])
  names(SA)<- "SA"
  #Return
  if(og_class =="RasterLayer"){
    SA<- raster::raster(SA)
    if(!is.null(filename)){
      return(raster::writeRaster(SA, filename=filename, overwrite=overwrite))
    }
  }
  if(!is.null(filename)){
    return(terra::writeRaster(SA, filename=filename, overwrite=overwrite))
  }
  return(SA)
  }