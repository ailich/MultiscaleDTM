#' Implementation of the Sappington et al., (2007) vector ruggedness measure
#'
#' Implementation of the Sappington et al., (2007) vector ruggedness measure, modified from Evans (2021). 
#' @param r DTM as a SpatRaster or RasterLayer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number. Default is 3x3.
#' @param na.rm A logical indicating whether or not to remove NA values before calculations. See details for more information.
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @param filename character Output filename.
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @return a RasterLayer
#' @examples
#' r<- rast(volcano, extent= ext(2667400, 2667400 + 
#' ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
#' crs = "EPSG:27200")
#' vrm<- VRM(r, w=c(5,5), na.rm = TRUE)
#' plot(vrm)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @importFrom utils packageVersion
#' @importFrom utils compareVersion
#' @details
#' If the crs is cartesian, when na.rm=TRUE, NA's will be removed from the slope/aspect calculations. When the crs is lat/lon, na.rm=TRUE will not affect the calculation of slope/aspect as terra::terrain will be used since it can calculate slope and aspect for spherical geometry but it does not support na.rm. In both cases when na.rm=TRUE, the x, y, and z components will be summed with na.rm=TRUE, and the N used in the denominator of the VRM equation will be the number of non-NA cells in the window rather than the total number of cells.
#' 
#' @references
#' Evans JS (2021). spatialEco. R package version 1.3-6, https://github.com/jeffreyevans/spatialEco.
#' 
#' Sappington, J.M., Longshore, K.M., Thompson, D.B., 2007. Quantifying Landscape Ruggedness for Animal Habitat Analysis: A Case Study Using Bighorn Sheep in the Mojave Desert. The Journal of Wildlife Management 71, 1419-1426. https://doi.org/10.2193/2005-723
#' @export

VRM<- function(r, w=c(3,3), na.rm = FALSE, include_scale=FALSE, filename=NULL, overwrite=FALSE, wopt=list()){
  og_class<- class(r)[1]
  if(og_class=="RasterLayer"){
    r<- terra::rast(r) #Convert to SpatRaster
  }
  #Input checks
  if(!(og_class %in% c("RasterLayer", "SpatRaster"))){
    stop("Error: Input must be a 'SpatRaster' or 'RasterLayer'")
  }
  if(terra::nlyr(r)!=1){
    stop("Error: Input raster must be one layer.")
  }
  if(length(w)==1){
    w<- rep(w,2)}
  if(length(w) > 2){
    stop("Specified window exceeds 2 dimensions")}
  if(any(0 == (w %% 2))){
    stop("Error: w must be odd")}
  if(all(w<3)){
    stop("Error: w must be greater or equal to 3 in at least one dimension")
  }
  
  if((compareVersion(as.character(packageVersion("terra")), "1.5.49")==-1) & isTRUE(terra::is.lonlat(r, perhaps=FALSE))){
    warning("Distance calculations conducted using Haversine (spheroid) rather geodesic formulas since terra version is < 1.5.49")
    }
  
  if(isTRUE(terra::is.lonlat(r, perhaps=FALSE)) | (!na.rm)){
    sa <- terra::terrain(r, v=c("slope", "aspect"), unit="radians", neighbors=8, wopt=wopt)
    } else{
      sa <- SlpAsp(r, w=c(3,3),unit="radians", method="queen", metrics= c("slope", "aspect"), na.rm=na.rm, include_scale=FALSE, mask_aspect=FALSE, wopt=wopt)
    } #SlpAsp can use na.rm in slope calculations, terra::terrain cannot but can handle spherical geometry.
  
  #Decompose vectors into components 
  sin.slp <- terra::math(sa$slope, fun="sin", wopt=wopt)
  Xrast <- terra::math(sa$aspect, fun="sin", wopt=wopt) * sin.slp # xRaster
  Yrast <- terra::math(sa$aspect, fun="cos", wopt=wopt) * sin.slp # yRaster
  Zrast <- terra::math(sa$slope, fun="cos", wopt=wopt) # zRaster 
  x.sum <- terra::focal(Xrast, w = w, fun=sum, na.rm=na.rm, wopt=wopt)
  y.sum <- terra::focal(Yrast, w = w, fun=sum, na.rm=na.rm, wopt=wopt) 
  z.sum <- terra::focal(Zrast, w = w, fun=sum, na.rm=na.rm, wopt=wopt)
  vrm.fun <- function(x, y, z) { 
    sqrt( (x^2) + (y^2) + (z^2) ) 
  }
  res_vect <- terra::lapp(c(x.sum, y.sum, z.sum), fun=vrm.fun, wopt=wopt) #resultant vector
  if(!na.rm){
    scale.factor <- round(w[1] * w[2], 0) #Constant scale factor
    } else{
      scale.factor<- terra::focalCpp(Zrast, w= w, fun = C_CountVals, wopt=wopt)
      } #If include na.rm option, adjust scale.factor to be number of non-NA cells in focal window
  out<- 1 - (res_vect / scale.factor) 
  names(out)<- "vrm"
  if(include_scale){names(out)<- paste0(names(out), "_", w[1],"x", w[2])} #Add scale to layer names
  
  #Return
  if(og_class =="RasterLayer"){
    out<- raster::raster(out)
    if(!is.null(filename)){
      return(raster::writeRaster(out, filename=filename, overwrite=overwrite))
    }
  }
  if(!is.null(filename)){
    return(terra::writeRaster(out, filename=filename, overwrite=overwrite, wopt=wopt))
  }
  return(out)
  }