#' Calculates Topographic Position Index
#'
#' Calculates Topographic Position Index (TPI). TPI is a measure of relative position that calculates the the difference between the value of the focal cell and the mean of mean of the surrounding cells (i.e. local mean but excluding the value of the focal cell).Positive values indicate local highs (i.e. peaks) and negative values indicate local lows (i.e. depressions). TPI can be expressed in units of the input DTM raster or can standardized relative to the local topography by dividing by the standard deviation or range of included elevation  values in the focal window. 
#' @param r DTM as a SpatRaster or RasterLayer. TPI calls the function RelPos internally which serves as a general purpose and more flexible function for calculating relative position.
#' @param w For a "rectangle" focal window, a vector of length 2 specifying dimensions where the first number is the number of rows and the second is the number of columns (or a single number if the number of rows and columns is equal). Window size must be an odd number, and the default is 3x3. For circle windows, w can be set to NA or NULL and radius can be used instead, or w can be specified using focal weights matrix created by MultiscaleDTM::circle_window.
#' @param shape Character representing the shape of the focal window. Either "rectangle" (default) or "circle".
#' @param radius For "circle" shaped focal windows, a single integer representing the radius in "cell" or "map" units. For a circle, the default radius is 1 cell if units= "cell" or the maximum of the x and y cell resolution if unit="map".
#' @param stand Standardization method. Either "none" (the default), "range" or "sd" indicating whether the TPI should be standardized by dividing by the standard deviation or range of included values in the focal window. If stand is 'none' the layer name will be "tpi", otherwise it will be "stpi" to indicate that the layer has been standardized.
#' @param unit Unit for radius. If a circle shaped window is specified with the radius parameter either "cell" (number of cells) or "map" for map units (e.g. meters). Otherwise it should be NA or NULL. If radius is specified the default is "cell" otherwise, if w is specified, the default is NA.
#' @param na.rm Logical indicating whether or not to remove NA values before calculations.
#' @param include_scale Logical indicating whether to append window size to the layer names (default = FALSE). If include_scale = TRUE the number of rows and number of columns will be appended for rectangular or custom windows. For circular windows it will be a single number representing the radius. If unit="map" then window size will have "MU" after the number indicating that the number represents the scale in map units.
#' @param filename Character output filename.
#' @param overwrite Logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt List with named options for writing files as in writeRaster.
#' @return SpatRaster or RasterLayer.
#' @examples 
#' r<- rast(volcano, extent= ext(2667400, 2667400 + 
#' ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
#' crs = "EPSG:27200")
#' tpi<- TPI(r, w=c(5,5), shape="rectangle", stand="none", na.rm = TRUE)
#' plot(tpi)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @importFrom dplyr case_when
#' @references 
#' Weiss, A., 2001. Topographic Position and Landforms Analysis. Presented at the ESRI user conference, San Diego, CA.
#' @export
#' 
TPI<- function(r, w=ifelse(tolower(shape)=="rectangle", c(3,3), NA_real_), shape= "rectangle",
               radius=dplyr::case_when(is.matrix(w) ~ NA_real_,
                                                    tolower(shape)=="circle" & tolower(unit)=="cell" ~ 1,
                                                    tolower(shape)=="circle" & tolower(unit)=="map" ~ max(terra::res(r)),
                                                    TRUE ~ NA_real_),
               stand="none", unit=ifelse((!is.matrix(w)) & all(is.null(w) | is.na(w)), "cell", NA_character_), 
               na.rm=FALSE, include_scale =FALSE, 
               filename=NULL, overwrite=FALSE, wopt=list()){
  
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
  tpi<- MultiscaleDTM::RelPos(r, w=w, shape= shape, radius=radius, stand=stand, exclude_center= TRUE, 
                        unit=unit, na.rm=na.rm, include_scale =include_scale, filename=NULL, 
                        overwrite=FALSE, wopt=list())
  names(tpi)<- gsub(pattern = "rpos", replacement = "tpi", names(tpi))
  
  #Return
  if(og_class =="RasterLayer"){
    tpi<- raster::raster(tpi)
    if(!is.null(filename)){
      return(raster::writeRaster(tpi, filename=filename, overwrite=overwrite))
    }
  }
  if(!is.null(filename)){
    return(terra::writeRaster(tpi, filename=filename, overwrite=overwrite, wopt=wopt))
  }
  return(tpi)
}