#' Implementation of the Sappington et al., (2007) vector ruggedness measure
#'
#' Implementation of the Sappington et al., (2007) vector ruggedness measure, modified from Evans (2021). 
#' @param r DEM as a SpatRaster or RasterLayer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @return a RasterLayer
#' @import terra
#' @importFrom raster raster
#' @references
#' Evans JS (2021). spatialEco. R package version 1.3-6, https://github.com/jeffreyevans/spatialEco.
#' 
#' Sappington, J.M., Longshore, K.M., Thompson, D.B., 2007. Quantifying Landscape Ruggedness for Animal Habitat Analysis: A Case Study Using Bighorn Sheep in the Mojave Desert. The Journal of Wildlife Management 71, 1419-1426. https://doi.org/10.2193/2005-723
#' @export

VRM<- function(r, w, na.rm = FALSE, include_scale=FALSE){
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
  sa <- terra::terrain(r, v=c("slope", "aspect"), unit="radians", neighbors=8) 					
  sin.slp <- terra::app(sa[["slope"]], fun=sin)                 # xyRaster 
  cos.slp <- terra::app(sa[["slope"]], fun=cos)                 # zRaster 
  sin.asp <- terra::app(sa[["aspect"]], fun=sin) * sin.slp      # yRaster
  cos.asp <- terra::app(sa[["aspect"]], fun=cos) * sin.slp      # xRaster  
  x.sum <- terra::focal(sin.asp, w = w, fun=sum, na.rm=na.rm)
  y.sum <- terra::focal(cos.asp, w = w, fun=sum, na.rm=na.rm) 
  z.sum <- terra::focal(cos.slp, w = w, fun=sum, na.rm=na.rm)
  vrm.fun <- function(x, y, z) { 
    sqrt( (x^2) + (y^2) + (z^2) ) 
  }
  res_vect <- terra::lapp(c(x.sum, y.sum, z.sum), fun=vrm.fun) #resultant vector
  if(!na.rm){
    scale.factor <- round(w[1] * w[2], 0) #Constant scale factor
    } else{
      scale.factor<- terra::focalCpp(cos.slp, w= w, fun = C_CountVals)
      } #If include na.rm option, adjust scale.factor to be number of non-NA cells in focal window
  out<- 1 - (res_vect / scale.factor) 
  names(out)<- "vrm"
  if(include_scale){names(out)<- paste0(names(out), "_", w[1],"x", w[2])} #Add scale to layer names
  if(og_class=="RasterLayer"){out<- raster::raster(out)}
  return(out)
  }
