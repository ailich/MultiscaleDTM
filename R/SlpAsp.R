#' Multiscale Slope and Aspect
#'
#' Calculates multiscale slope and aspect based on a modified version of the algorithm from Misiuk et al (2021) which extends classical formulations of slope restricted to a 3x3 window to multiple scales by using only cells on the edges of the focal window (see details for more information).
#' @param r DTM as a SpatRaster or RasterLayer in a projected coordinate system where map units match elevation/depth units
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param unit "degrees" or "radians"
#' @param method "rook", "queen" (default), or "boundary". The method indicates which cells to use to in computations. "rook" uses only the 4 edge cells directly up, down, left, and right; "queen" adds an additional four corner cells; "boundary" uses all edge cells (see details for more information).
#' @param metrics a character string or vector of character strings of which terrain attributes to return. Default is to return all available attributes which are c("slope", "aspect", "eastness", "northness").
#' @param na.rm Logical indicating whether or not to remove NA values before calculations. Not applicable for "rook" method.
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @param mask_aspect A logical. When mask_aspect is TRUE (the default), if slope evaluates to 0, aspect will be set to NA and both eastness and northness will be set to 0. When mask_aspect is FALSE, when slope is 0 aspect will be pi/2 radians or 90 degrees which is the behavior of raster::terrain, and northness and eastness will be calculated from that.
#' @param filename character Output filename. Can be a single filename, or as many filenames as there are layers to write a file for each layer
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @import terra
#' @importFrom raster raster
#' @importFrom raster stack
#' @importFrom raster writeRaster
#' @return a SpatRaster or RasterStack of slope and/or aspect (and components of aspect)
#' @examples
#' r<- rast(volcano, extent= ext(2667400, 2667400 + 
#' ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
#' crs = "EPSG:27200")
#' slp_asp<- SlpAsp(r = r, w = c(5,5), unit = "degrees", 
#' method = "queen", metrics = c("slope", "aspect", 
#' "eastness", "northness"))
#' plot(slp_asp)
#' @details Slope is calculated atan(sqrt(dz.dx^2 + dz.dy^2)) and aspect is calculated as (-pi/2)-atan_2(dz.dy, dz.dx) and then constrained from 0 to 2 pi/0 to 360 degrees. 
#' dz.dx is the difference in between the weighted mean of the right side of the focal window and weighted mean of the left side of the focal window divided by the x distance of the focal window in map units.
#' dz.dy is the difference in between the weighted mean of the top side of the focal window and weighted mean of the bottom side of the focal window divided by the y distance of the focal window in map units.
#' The cells used in these computations is dependent on the "method" chosen. For methods "queen" and "boundary", corner cells have half the weight of all other cells used in the computations.
#' @references
#' Fleming, M.D., Hoffer, R.M., 1979. Machine processing of landsat MSS data and DMA topographic data for forest cover type mapping (No. LARS Technical Report 062879). Laboratory for Applications of Remote Sensing, Purdue University, West Lafayette, Indiana.
#' 
#' Horn, B.K., 1981. Hill Shading and the Reflectance Map. Proceedings of the IEEE 69, 14-47.
#' 
#' Misiuk, B., Lecours, V., Dolan, M.F.J., Robert, K., 2021. Evaluating the Suitability of Multi-Scale Terrain Attribute Calculation Approaches for Seabed Mapping Applications. Marine Geodesy 44, 327-385. https://doi.org/10.1080/01490419.2021.1925789
#' 
#' Ritter, P., 1987. A vector-based slope and aspect generation algorithm. Photogrammetric Engineering and Remote Sensing 53, 1109-1111.
#' @import terra
#' @importFrom raster raster
#' @importFrom raster stack
#' @export

SlpAsp <- function(r, w=c(3,3), unit="degrees", method="queen", metrics= c("slope", "aspect", "eastness", "northness"), na.rm=FALSE, include_scale=FALSE, mask_aspect=TRUE, filename=NULL, overwrite=FALSE, wopt=list()){ 
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
  if(isTRUE(terra::is.lonlat(r, perhaps=FALSE))){
    stop("Error: Coordinate system is Lat/Lon. Coordinate system must be projected with elevation/depth units matching map units.")
  }
  if(terra::is.lonlat(r, perhaps=TRUE, warn=FALSE)){
    warning("Coordinate system may be Lat/Lon. Please ensure that the coordinate system is projected with elevation/depth units matching map units.")
  }
  if(length(w)==1){w<- rep(w,2)}
  if(length(w) > 2){
    stop("Specified window exceeds 2 dimensions")}
  if(any(0 == (w %% 2))){
    stop("w must be odd")
  }
  if(any(w<3)){
    stop("Error: w must be greater or equal to 3")
  }
  unit<- tolower(unit) #make lowercase
  if(!(unit %in% c("degrees", "radians"))){
    stop("unit must be 'degrees' or 'radians'")
  }
  
  if(method==4){method<- "rook"}
  if(method==8){method<- "queen"}
  
  if(!(method %in% c("queen", "rook", "boundary"))){
    stop("method must be 'queen', 'rook', 'boundary', '8', or '4'")
  }
  metrics<- tolower(metrics) #Make all lowercase
  if(any(!(metrics %in% c("slope", "aspect","northness", "eastness")))){
    stop("metrics must be 'slope', 'aspect', 'northness', and/or 'eastness'")
  }
  
  needed_metrics<- metrics
  if(any(c("eastness", "northness") %in% needed_metrics) & !("aspect" %in% needed_metrics)){
    needed_metrics<- c(needed_metrics, "aspect")
  }
  if(mask_aspect & ("aspect" %in% needed_metrics)){
    needed_metrics<- c(needed_metrics, "slope")
  }
  if(method=="rook" & na.rm){
    warning("na.rm=TRUE is only relevant if `method` is 'queen' or '8'")
  }
  
  #j is the number of cells on either side of the focal cell; l is used to generate the focal matrix
  jx <- floor((w[2]/2))
  jy <- floor(w[1]/2)
  
  if(method=="queen"){
    mat.l<- matrix(data = NA, nrow=w[1], ncol=w[2]) # left
    mat.l[1,1]<- 1
    mat.l[jy+1,1]<- 2
    mat.l[w[1],1]<- 1
    mat.l<- mat.l/sum(mat.l, na.rm=TRUE)
    
    mat.t<- matrix(data = NA, nrow=w[1], ncol=w[2]) # top
    mat.t[1,1]<- 1
    mat.t[1, jx+1]<- 2
    mat.t[1, w[2]]<- 1
    mat.t<- mat.t/sum(mat.t, na.rm=TRUE)
    mat.t
  } else if(method=="rook"){
    mat.l<- matrix(data = NA, nrow=w[1], ncol=w[2]) # left
    mat.l[jy+1,1]<- 1
    
    mat.t<- matrix(data = NA, nrow=w[1], ncol=w[2]) # top
    mat.t[1, jx+1]<- 1
  } else {
    mat.l<- matrix(data = NA, nrow=w[1], ncol=w[2]) # left
    mat.l[,1]<- 2
    mat.l[1,1]<- 1
    mat.l[w[1],1]<- 1
    mat.l<- mat.l/sum(mat.l, na.rm=TRUE)
    mat.l
    
    mat.t<- matrix(data = NA, nrow=w[1], ncol=w[2]) # top
    mat.t[1,]<- 2
    mat.t[1, 1]<- 1
    mat.t[1, w[2]]<- 1
    mat.t<- mat.t/sum(mat.t, na.rm=TRUE)
    mat.t
  }
  
  mat.r<-  mat.l[, ncol(mat.l):1] # right, flip matrix
  mat.b<-  mat.t[nrow(mat.t):1, ] # bottom, flip matrix
  
  if(!na.rm){
    mean.r<- focal(r, w=mat.r, fun=sum, na.rm=FALSE)
    mean.l<- focal(r, w=mat.l, fun=sum, na.rm=FALSE)
    mean.t<- focal(r, w=mat.t, fun=sum, na.rm=FALSE)
    mean.b<- focal(r, w=mat.b, fun=sum, na.rm=FALSE)
  } else{
    mean.r<- focal(r, w=mat.r, fun=mean, na.rm=TRUE)
    mean.l<- focal(r, w=mat.l, fun=mean, na.rm=TRUE)
    mean.t<- focal(r, w=mat.t, fun=mean, na.rm=TRUE)
    mean.b<- focal(r, w=mat.b, fun=mean, na.rm=TRUE)
  }
  
  dx<- xres(r) * jx * 2
  dy<- yres(r) * jy * 2
  
  dz.dx<- (mean.r-mean.l)/dx
  dz.dy<- (mean.t-mean.b)/dy
  
  out<- terra::rast() #initialize output
  
  if("slope" %in% needed_metrics){
    slope.k<- terra::math(terra::math(dz.dx^2 + dz.dy^2, fun= "sqrt", wopt=wopt), fun= "atan", wopt=wopt)
    if(unit=="degrees" & ("slope" %in% metrics)){
      slope.k<- slope.k * (180/pi)
    }
    names(slope.k)<- "slope"
    out<- c(out, slope.k, warn=FALSE)
  }
  
  if("aspect" %in% needed_metrics){
    aspect.k<- (-pi/2) - terra::atan_2(dz.dy, dz.dx, wopt=wopt) # aspect relative to North
    aspect.k<- ifel(aspect.k < 0, yes = aspect.k+(2*pi), no= aspect.k, wopt=wopt) # Constrain range so between 0 and 2pi
    
    if("eastness" %in% needed_metrics){
      eastness.k<- terra::math(aspect.k, fun="sin", wopt=wopt)
      if(mask_aspect){
        eastness.k<- terra::mask(eastness.k, mask= slope.k, maskvalues = 0, updatevalue = 0, wopt=wopt) #Set eastness to 0 where slope is zero
      }
      names(eastness.k)<- "eastness"
      out<- c(out, eastness.k, warn=FALSE)
    }
    
    if("northness" %in% needed_metrics){
      northness.k<- terra::math(aspect.k, fun="cos", wopt=wopt)
      if(mask_aspect){
        northness.k<- terra::mask(northness.k, mask= slope.k, maskvalues = 0, updatevalue = 0, wopt=wopt) #Set northenss to 0 where slope is zero
      }
      names(northness.k)<- "northness"
      out<- c(out, northness.k, warn=FALSE)
    }
    
    if(mask_aspect){
      aspect.k<- terra::mask(aspect.k, mask= slope.k, maskvalues = 0, updatevalue = NA, wopt=wopt) #Set aspect to undefined where slope is zero
    }
    if(unit=="degrees"){
      aspect.k<- aspect.k * (180/pi)
    }
    names(aspect.k)<- "aspect"
    out<- c(out, aspect.k, warn=FALSE)
  }
  
  if(!is.null(metrics)){out<- terra::subset(out, metrics, wopt=wopt)} #Subset needed metrics to requested metrics in proper order
  if(include_scale){names(out)<- paste0(names(out), "_", w[1], "x", w[2])}
  
  #Return
  if(og_class=="RasterLayer"){
    if(terra::nlyr(out) > 1){
      out<- raster::stack(out) #Convert to RasterStack
      if(!is.null(filename)){
        if(length(filename)==1){
          return(raster::writeRaster(out, filename=filename, overwrite=overwrite, bylayer=FALSE))
        } else{
          return(raster::writeRaster(out, filename=filename, overwrite=overwrite, bylayer=TRUE))
        }
      }
    } else{
      out<- raster::raster(out)
      if(!is.null(filename)){
        return(raster::writeRaster(out, filename=filename, overwrite=overwrite))
      }
    }
  }
  if(!is.null(filename)){
    return(terra::writeRaster(out, filename=filename, overwrite=overwrite, wopt=wopt))
  }
  return(out)
}

  
  
