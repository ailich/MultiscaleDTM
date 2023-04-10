#' Calculates Difference from Mean Value (DMV)
#'
#' Calculates Difference from Mean Value (DMV). DMV is a measure of relative position that calculates the difference between the value of the focal cell and the mean of all cells in a rectangular or circular neighborhood. Positive values indicate local highs (i.e. peaks) and negative values indicate local lows (i.e. depressions). DMV can be expressed in units of the input DTM raster or can standardized relative to the local topography by dividing by the standard deviation or range of elevation values in the focal window. DMV calls the function RelPos internally which serves as a general purpose and more flexible function for calculating relative position.
#' @param r DTM as a SpatRaster or RasterLayer.
#' @param w For a "rectangle" focal window, a vector of length 2 containing odd numbers specifying dimensions where the first number is the number of rows and the second is the number of columns (or a single number if the number of rows and columns is equal). For a "circle" shaped focal window, a single integer representing the radius in "cell" or "map" units or a focal weights matrix created by MultiscaleDTM::circle_window.
#' @param shape Character representing the shape of the focal window. Either "rectangle" (default) or "circle".
#' @param stand Standardization method. Either "none" (the default), "range" or "sd" indicating whether the TPI should be standardized by dividing by the standard deviation or range of included values in the focal window. If stand is 'none' the layer name will be "dmv", otherwise it will be "sdmv" to indicate that the layer has been standardized.
#' @param unit Unit for w if shape is 'circle' and it is a vector (default is unit="cell"). For circular windows specified with a matrix, unit is ignored and extracted directly from w. For rectangular and custom focal windows set unit='cell' or set unit to NA/NULL.
#' @param na.rm Logical indicating whether or not to remove NA values before calculations.
#' @param include_scale Logical indicating whether to append window size to the layer names (default = FALSE) or a character vector specifying the name you would like to append or a number specifying the number of significant digits. If include_scale = TRUE the number of rows and number of columns will be appended for rectangular windows. For circular windows it will be a single number representing the radius. If unit="map" then window size will have "MU" after the number indicating that the number represents the scale in map units (note units can be extracted from w created with MultiscaleDTM::circle_window).
#' @param filename Character output filename.
#' @param overwrite Logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt List with named options for writing files as in writeRaster.
#' @return a SpatRaster or RasterLayer.
#' @examples
#' r<- rast(volcano, extent= ext(2667400, 2667400 + 
#' ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
#' crs = "EPSG:27200")
#' dmv<- DMV(r, w=c(5,5), shape= "rectangle", stand="range", na.rm = TRUE)
#' plot(dmv)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @importFrom dplyr case_when
#' @references 
#' Lecours, V., Devillers, R., Simms, A.E., Lucieer, V.L., Brown, C.J., 2017. Towards a Framework for Terrain Attribute Selection in Environmental Studies. Environmental Modelling & Software 89, 19-30. https://doi.org/10.1016/j.envsoft.2016.11.027
#' Wilson, J.P., Gallant, J.C. (Eds.), 2000. Terrain Analysis: Principles and Applications. John Wiley & Sons, Inc.
#' @export

DMV<- function(r, w = dplyr::case_when(tolower(shape)=="rectangle" ~ 3,
                                       tolower(shape)=="circle" & isTRUE(tolower(unit)=="cell") ~ 1,
                                       tolower(shape)=="circle" & isTRUE(tolower(unit)=="map") ~ max(terra::res(r))), 
               shape= "rectangle", stand="none", unit="cell", na.rm=FALSE, 
               include_scale=FALSE, filename=NULL, overwrite=FALSE, wopt=list()){
  og_class<- class(r)[1]
  if(og_class=="RasterLayer"){
    r<- terra::rast(r) #Convert to SpatRaster
  }
  #Input checks
  if(!(og_class %in% c("RasterLayer", "SpatRaster"))){
    stop("Error: Input must be a 'SpatRaster' or 'RasterLayer'")
  }
  if (!(shape %in% c("rectangle", "circle"))){
    stop("Error: shape must be 'rectangle' or 'circle'")
  }
  dmv<- MultiscaleDTM::RelPos(r, w=w, shape= shape, stand=stand, exclude_center= FALSE, 
               unit=unit, fun="mean", na.rm=na.rm, include_scale =include_scale, filename=filename, overwrite=overwrite, wopt=wopt)
  names(dmv)<- gsub(pattern = "rpos", replacement = "dmv", names(dmv))
  
  #Return
  if(og_class =="RasterLayer"){
    dmv<- raster::raster(dmv)
    if(!is.null(filename)){
      return(raster::writeRaster(dmv, filename=filename, overwrite=overwrite))
    }
  }
  if(!is.null(filename)){
    return(terra::writeRaster(dmv, filename=filename, overwrite=overwrite, wopt=wopt))
  }
  return(dmv)
}
