#' Calculates Difference from Mean Value (DMV)
#'
#' Calculates Difference from Mean Value (DMV) which is a measure of relative position. DMV calculates the difference between the focal cell and the mean of all cells in a rectangular or circular neighborhood. Positive values indicate local highs (i.e. peaks) and negative values indicate local lows (i.e. depressions). DMV can be expressed in units of the input DTM raster or can standardized relative to the local topography by dividing by the standard deviation or range of elevation values in the focal window. 
#' @param r DTM as a SpatRaster or RasterLayer
#' @param w For a "rectangle" focal window, a vector of length 2 specifying the dimensions where the first number is the number of rows and the second number is the number of columns (or a single number if the number of rows and columns is equal). Window size must be an odd number, and the default is 3x3. w can also be a focal weights matrix. Circular windows this can be created using MultiscaleDTM::circle_window or w can be set to NULL and radius can be used instead.
#' @param shape character representing the shape of the focal window. Either "rectangle" (default), or "circle".
#' @param radius A single integer representing the radius of the circle or a vector in "cell" or "map" units. This is ignored if w is not NULL.
#' @param stand standardization method. Either "none" (the default), "range" or "sd" indicating whether the relative position should be standardized by the standard deviation or range of values in the focal window. If stand is 'none' the layer name will be dmv, otherwise it will be sdmv to indicate that the layer has been standardized.
#' @param unit unit for radius. Either "cell" (number of cells, the default) or "map" for map units (e.g. meters).
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE). If unit="map" then window size will have "MU" after the number indicating that the number represents the window size in map units if a circular window is used.
#' @param filename character Output filename.
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @return a SpatRaster or RasterLayer
#' @examples
#' r<- rast(volcano, extent= ext(2667400, 2667400 + 
#' ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
#' crs = "EPSG:27200")
#' dmv<- DMV(r, w=c(5,5), shape= "rectangle", na.rm = TRUE, stand="range")
#' plot(dmv)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @references 
#' Lecours, V., Devillers, R., Simms, A.E., Lucieer, V.L., Brown, C.J., 2017. Towards a Framework for Terrain Attribute Selection in Environmental Studies. Environmental Modelling & Software 89, 19-30. https://doi.org/10.1016/j.envsoft.2016.11.027
#' Wilson, J.P., Gallant, J.C. (Eds.), 2000. Terrain Analysis: Principles and Applications. John Wiley & Sons, Inc.
#' @export

DMV<- function(r, w = ifelse(tolower(shape)=="rectangle", c(3,3), NA_real_), shape= "rectangle",
               radius=dplyr::case_when(is.matrix(w) ~ NA_real_,
                                       tolower(shape)=="circle" & tolower(unit)=="cell" ~ 1,
                                       tolower(shape)=="circle" & tolower(unit)=="map" ~ max(terra::res(r)),
                                       TRUE ~ NA_real_), 
               unit="cell", stand="none", na.rm=FALSE, include_scale=FALSE, filename=NULL, overwrite=FALSE, wopt=list()){
  ## shape check
  if (!(shape %in% c("rectangle", "circle"))){
    stop("Error: shape must be 'rectangle' or 'circle'")
  }
  dmv<- RelPos(r, w=w, shape= shape, radius=radius, stand=stand, exclude_center= FALSE, 
               unit=unit, na.rm=na.rm, include_scale =include_scale, filename=filename, overwrite=overwrite, wopt=wopt)
  return(dmv)
}
