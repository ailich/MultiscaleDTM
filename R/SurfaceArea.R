#' Calculates surface area of a DEM
#'
#' Calculates surface area on a per cell basis of a DEM based on Jenness, 2004.
#' @param r DEM as a SpatRaster or RasterLayer in a projected coordinate system where map units match elevation/depth units
#' @return a SpatRaster or RasterLayer
#' @import terra
#' @importFrom raster raster
#' @references
#' Jenness, J.S., 2004. Calculating landscape surface area from digital elevation models. Wildlife Society Bulletin 32, 829-839. https://doi.org/10.2193/0091-7648(2004)032[0829:CLSAFD]2.0.CO;2
#' @export

SurfaceArea<- function(r, expand=FALSE){
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
  if(og_class=="RasterLayer"){SA<- raster::raster(SA)}
  return(SA)
  }