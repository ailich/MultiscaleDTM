#' Calculates Bathymetric Position Index
#'
#' Calculates Bathymetric Position Index (BPI). BPI is a measure of relative position that calculates the difference between the value of the focal cell and the mean of cells contained within an annulus shaped neighborhood. Positive values indicate local highs (i.e. peaks) and negative values indicate local lows (i.e. depressions). BPI can be expressed in units of the input DTM raster or can standardized relative to the local topography by dividing by the standard deviation or range of included elevation values in the focal window. BPI calls the function RelPos internally which serves as a general purpose and more flexible function for calculating relative position.
#' @param r r DTM as a SpatRaster or RasterLayer.
#' @param w Vector of length 2 specifying c(inner, outer) radii of the annulus in "cell" or "map" units or a focal weights matrix created by MultiscaleDTM::annulus_window. Inner radius must be less than or equal to outer radius. There is no default size.
#' @param stand Standardization method. Either "none" (the default), "range" or "sd" indicating whether the relative position should be standardized by dividing by the standard deviation or range of included values in the focal window. If stand is 'none' the layer name will be "bpi", otherwise it will be "sbpi" to indicate that the layer has been standardized.
#' @param unit Unit for w if it is a vector (default is unit="cell"). If w is a matrix, unit is ignored and extracted directly from w.
#' @param na.rm Logical indicating whether or not to remove NA values before calculations.
#' @param include_scale Logical indicating whether to append window size to the layer names (default = FALSE) or a character vector specifying the name you would like to append or a number specifying the number of significant digits. If include_scale = TRUE the appended scale will be the inner and outer radius. If unit="map" then window size will have "MU" after the number indicating that the number represents the scale in map units (note units can be extracted from w created with MultiscaleDTM::circle_window and MultiscaleDTM::annulus_window).
#' @param filename Character output filename.
#' @param overwrite Logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt List with named options for writing files as in writeRaster.
#' @return A SpatRaster or RasterLayer.
#' @examples
#' r<- rast(volcano, extent= ext(2667400, 2667400 + 
#' ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
#' crs = "EPSG:27200")
#' bpi<- BPI(r, w = c(2,4), stand= "none", unit = "cell", na.rm = TRUE)
#' plot(bpi)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @importFrom dplyr case_when
#' @references 
#' Lundblad, E.R., Wright, D.J., Miller, J., Larkin, E.M., Rinehart, R., Naar, D.F., Donahue, B.T., Anderson, S.M., Battista, T., 2006. A benthic terrain classification scheme for American Samoa. Marine Geodesy 29, 89–111. https://doi.org/10.1080/01490410600738021
#' @export

BPI<- function(r, w, stand="none", 
               unit="cell", 
               na.rm=FALSE, include_scale =FALSE, 
               filename=NULL, overwrite=FALSE, wopt=list()){
  og_class<- class(r)[1]
  if(og_class=="RasterLayer"){
    r<- terra::rast(r) #Convert to SpatRaster
  }
  unit<- tolower(unit)
  #Input checks
  if(!(og_class %in% c("RasterLayer", "SpatRaster"))){
    stop("Error: Input must be a 'SpatRaster' or 'RasterLayer'")
  }
  
  bpi<- MultiscaleDTM::RelPos(r, w=w, shape= "annulus", stand=stand, 
               exclude_center= FALSE, unit=unit, fun="mean", na.rm=na.rm, 
               include_scale =include_scale, filename=NULL, overwrite=FALSE, wopt=list())
  names(bpi)<- gsub(pattern = "rpos", replacement = "bpi", names(bpi))
  
  #Return
  if(og_class =="RasterLayer"){
    bpi<- raster::raster(bpi)
    if(!is.null(filename)){
      return(raster::writeRaster(bpi, filename=filename, overwrite=overwrite))
    }
  }
  if(!is.null(filename)){
    return(terra::writeRaster(bpi, filename=filename, overwrite=overwrite, wopt=wopt))
  }
  return(bpi)
}