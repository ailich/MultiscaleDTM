#' Directional Slope
#'
#' Calculates the slope along a specified direction. Upslope values are positive and downslope values are negative.
#' @param alpha Angle (in specified 'unit') at which you would like to calculate slope. 0 represents up in map direction (usually North) and it increases clockwise. This can be a single number or it can be a raster of cell values.
#' @param dz.dx The change in elevation per unit distance in the x direction as a SpatRaster, RasterLayer, or a single number. Positive is to the right. See details for more.
#' @param dz.dy The change in elevation per unit distance in the y direction as a SpatRaster or RasterLayer,or a single number Positive is up. See details for more.
#' @param unit "degrees" or "radians" (default is "degrees")
#' @param abs logical indicating whether or not to return the absolute value of slope (default is FALSE)
#' @param include_dir logical indicating whether to append direction to layer name (default is FALSE)
#' @param filename character Output filename. Can be a single filename, or as many filenames as there are layers to write a file for each layer
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @import terra
#' @importFrom raster raster
#' @importFrom raster stack
#' @importFrom raster writeRaster
#' @return a SpatRaster or RasterStack of slope and/or aspect (and components of aspect)
#' @examples
#' r<- erupt()
#' dz1<- SlpAsp(r, metrics = c("dz.dx", "dz.dy"))
#' dz2<- Qfit(r, metrics = c(), return_params = TRUE, as_derivs=TRUE)
#' dz3<- Pfit(r, metrics = c("dz.dx", "dz.dy"))
#' dirslp1<- DirSlp(alpha = 45, dz.dx= dz1$dz.dx, dz.dy= dz1$dz.dy)
#' dirslp2<- DirSlp(alpha = 45, dz.dx= dz2$zx, dz.dy= dz2$zy)
#' dirslp3<- DirSlp(alpha = 45, dz.dx= dz3$dz.dx, dz.dy= dz3$dz.dy)
#' @details dz.dx and dz.dy can be calculated at a specified scale via `SlpAsp`, `Pfit`, `Qfit` (zx and zy), or from an existing layer calculated by another program.
#' @references
#' Neteler, M., & Mitasova, H. (2008). Open source GIS: A GRASS GIS approach (3rd ed.). Springer.
#' @import terra
#' @importFrom raster raster
#' @importFrom raster stack
#' @export

DirSlp<- function(alpha, dz.dx, dz.dy, unit = "degrees", abs=FALSE, include_dir=FALSE, filename=NULL, overwrite=FALSE, wopt=list()){
  og_class<- c(class(alpha)[1], class(dz.dx)[1], class(dz.dy)[1])
  if(og_class[1]=="RasterLayer"){
    alpha<- terra::rast(alpha)
  }
  if(og_class[2]=="RasterLayer"){
    dz.dx<- terra::rast(dz.dx)
  }
  if(og_class[3]=="RasterLayer"){
    dz.dy<- terra::rast(dz.dy)
  }
  
  #Input checks
  if(any(!(og_class %in% c("RasterLayer", "SpatRaster", "integer", "numeric")))){
    stop("Error: Inputs 'alpha', 'dz.dx', and 'dz.dy' must be a 'SpatRaster', 'RasterLayer', 'integer', or 'numeric'")
  }
  
  unit<- tolower(unit) #make lowercase
  if(!(unit %in% c("degrees", "radians"))){
    stop("unit must be 'degrees' or 'radians'")
  }
  
  if(class(alpha)[1]=="SpatRaster"){
    direction<- "variable"
  } else{
    direction<- alpha
  }
  
  if(unit=="degrees"){
    alpha<- alpha*(pi/180) # convert to radians
  }
  
  if(class(alpha)[1]=="SpatRaster"){
    dir_slp<- terra::math((dz.dx * terra::math(alpha, fun="sin", wopt=wopt)) + (dz.dy * terra::math(alpha, fun="cos", wopt=wopt)), fun="atan", wopt=wopt)
  } else{
    dir_slp<- atan(dz.dx * sin(alpha) + (dz.dy) * cos(alpha))
  }
  if(unit=="degrees"){
    dir_slp<- dir_slp*(180/pi) # convert to degrees
    names(dir_slp)<- "dirslope"
  }
  
  if(abs){
    if(class(dir_slp)[1]=="SpatRaster"){
      dir_slp<- terra::math(dir_slp, fun="abs", wopt=wopt)
    } else{
      dir_slp<- abs(dir_slp)
    }
  }
  
  if(include_dir){names(dir_slp)<- paste0(names(dir_slp), "_", direction)}
  
  #Return
  if(og_class[2]=="RasterLayer"){
    dir_slp<- raster::raster(dir_slp)
    if(!is.null(filename)){
      return(raster::writeRaster(dir_slp, filename=filename, overwrite=overwrite))
      }
  }
  
  if(!is.null(filename)){
    return(terra::writeRaster(dir_slp, filename=filename, overwrite=overwrite, wopt=wopt))
  }
  return(dir_slp)
}
