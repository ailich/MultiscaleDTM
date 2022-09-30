#' Creates circular focal window
#'
#' Creates circular focal window around central pixel.
#' @param radius radius of circular window
#' @param resolution resolution of intended raster layer (one number or a vector of length 2). Only necessary if unit= "map"
#' @param unit unit for radius. Either "cell" (number of cells, the default) or "map" for map units (e.g. meters).
#' @param return_dismat logical, if TRUE return a matrix of distances from focal cell instead of a matrix to pass to terra::focal.
#' @return if a matrix of 1's and NA's showing which cells to include and exclude respectively in focal calculations, or if return_dismat=TRUE, a matrix indicating the distance from the focal cell.
#' @export
circle_window<- function(radius, unit= "cell", resolution, return_dismat = FALSE){
  unit<- tolower(unit)
  if (unit != "map" & unit != "cell"){
    stop("Error: unit must equal 'map' or 'cell'")
  }
  if(length(radius)!=1){
    stop("Error: radius must be a single integer")
  }
  if(unit=="cell"){
    resolution<- c(1,1)
  }
  if (length(resolution)==1){
    resolution<- rep(resolution,2)
  }
  nrows<- floor((radius/resolution[2])*2+1)
  ncols<- floor((radius/resolution[1])*2+1)
  if(nrows %% 2 ==0){
    nrows<- nrows+1
  }
  if(ncols %% 2 ==0){
    ncols<- ncols+1
  } #nrow and ncol must be odd to have central pixel
  x<- matrix(seq(1:ncols), nrow = nrows, ncol =ncols, byrow = TRUE)  - ((ncols+1)/2)
  y<- matrix(seq(1:nrows), nrow=nrows, ncol=ncols, byrow = FALSE) - ((nrows+1)/2)
  x<- x * resolution[1]
  y<- y * resolution[2]
  dis_mat<- sqrt((y^2)+(x^2)) #Distance from center of window
  if(return_dismat){
    return(dis_mat)
  }
  w<- matrix(NA, nrow = nrow(dis_mat), ncol= ncol(dis_mat))
  w[dis_mat <= radius]<- 1
  if(unit=="map"){
    comment(w)<- paste0(radius,"MUx", radius, "MU")
    } else{
      comment(w)<- paste0(radius,"x", radius)
    }
  return(w)
}

#' Creates annulus focal window
#'
#' Creates annulus focal window around central pixel.
#' @param radius radius of inner annulus c(inner,outer)
#' @param unit unit for radius. Either "cell" (number of cells, the default) or "map" for map units (e.g. meters).
#' @param resolution resolution of intended raster layer (one number or a vector of length 2). Only necessary if unit= "map"
#' @param return_dismat logical, if TRUE return a matrix of distances from focal cell instead of a matrix to pass to terra::focal (default FALSE)
#' @return if a matrix of 1's and NA's showing which cells to include and exclude respectively in focal calculations, or if return_dismat=TRUE, a matrix indicating the distance from the focal cell.
#' @export
annulus_window<- function(radius, unit= "cell", resolution, return_dismat=FALSE){
  unit<- tolower(unit)
  if(length(radius)==1){radius<- rep(radius, 2)}
  if(length(radius) > 2){
    stop("Specified radius exceeds 2 dimensions")
  }
  if(radius[1] > radius[2]){
    stop("Error: inner radius must be less than or equal to outer radius")
  }
  if((radius[1] < 1) & (unit == "cell")){
    stop("Error: inner radius must be at least 1")
  }
  if(unit=="cell"){
    resolution<- c(1,1)
  }
  if((radius[1] < max(resolution)) & (unit == "map")){
    stop("Error: inner radius must be >= resolution")
  }
  if(length(resolution) > 2){
    stop("Specified inner radius exceeds 2 dimensions")
  }
  
  dis_mat<- circle_window(radius = radius[2], unit = unit, resolution = resolution, return_dismat = TRUE)
  if(return_dismat){
    return(dis_mat)
  }
  w<- matrix(NA, nrow=nrow(dis_mat), ncol=ncol(dis_mat))
  w[(dis_mat >= radius[1]) & (dis_mat <= radius[2])]<- 1
  if(unit=="map"){
    comment(w)<- paste0(radius[1],"MUx", radius[2], "MU")
  } else{
    comment(w)<- paste0(radius[1],"x", radius[2])
  }
  return(w)
}

#' Calculates Bathymetric Position Index
#'
#' Calculates Bathymetric Position Index (BPI). This is the value of the focal pixel minus the mean of the surrounding pixels contained within an annulus shaped window.
#' @param r DTM as a SpatRaster or RasterLayer
#' @param radius a vector of length 2 specifying the inner and outer radii of the annulus c(inner,outer). This is ignored if w is provided.
#' @param unit unit for radius. Either "cell" (number of cells, the default) or "map" for map units (e.g. meters). This is ignored if w is provided.
#' @param w A focal weights matrix specifying which cells to include and exclude in the annulus focal window which can be created using MultiscaleDTM::annulus_window.
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE). If unit="map" then window size will have "MU" after the number indicating that the number represents the window size in map units.
#' @param filename character Output filename.
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @return a SpatRaster or RasterLayer
#' @examples
#' r<- rast(volcano, extent= ext(2667400, 2667400 + 
#' ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
#' crs = "EPSG:27200")
#' bpi<- BPI(r, radius = c(2,4), unit = "cell", na.rm = TRUE)
#' plot(bpi)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @references 
#' Lundblad, E.R., Wright, D.J., Miller, J., Larkin, E.M., Rinehart, R., Naar, D.F., Donahue, B.T., Anderson, S.M., Battista, T., 2006. A benthic terrain classification scheme for American Samoa. Marine Geodesy 29, 89–111. https://doi.org/10.1080/01490410600738021
#' @export

BPI<- function(r, radius=NULL, unit= "cell", w=NULL, na.rm=FALSE, include_scale=FALSE, filename=NULL, overwrite=FALSE, wopt=list()){
  og_class<- class(r)[1]
  if(og_class=="RasterLayer"){
    r<- terra::rast(r) #Convert to SpatRaster
  }
  unit<- tolower(unit)
  #Input checks
  if(!(og_class %in% c("RasterLayer", "SpatRaster"))){
    stop("Error: Input must be a 'SpatRaster' or 'RasterLayer'")
  }
  if(terra::nlyr(r)!=1){
    stop("Error: Input raster must be one layer.")
  }
  if(!is.null(w)){
    if(class(w)[1] != "matrix"){
      stop("Error: w must be a matrix")
    }
    if(any(0 == (dim(w) %% 2))){
      stop("Error: dimensions of w must be odd")
    }
  }
  if(length(radius)==1){radius<- rep(radius, 2)}
  
  resolution<- terra::res(r)
  
  if(is.null(w)){
    w<- annulus_window(radius = radius, unit = unit, resolution = resolution)
  }
  
  bpi<- r - terra::focal(x = r, w = w, fun = mean, na.rm = na.rm, wopt=wopt)
  names(bpi)<- "bpi"
  if(include_scale){names(bpi)<- paste0(names(bpi), "_", comment(w))} #Add scale to layer names
  
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