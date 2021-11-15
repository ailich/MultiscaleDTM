#' Implementation of the Sappington et al., (2007) vector ruggedness measure
#'
#' Implementation of the Sappington et al., (2007) vector ruggedness measure, modified from Evans (2021). 
#' @param r DEM as a raster layer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @return a RasterLayer
#' @import raster
#' @references
#' Evans JS (2021). spatialEco. R package version 1.3-6, https://github.com/jeffreyevans/spatialEco.
#' 
#' Sappington, J.M., Longshore, K.M., Thompson, D.B., 2007. Quantifying Landscape Ruggedness for Animal Habitat Analysis: A Case Study Using Bighorn Sheep in the Mojave Desert. The Journal of Wildlife Management 71, 1419-1426. https://doi.org/10.2193/2005-723
#' @export

VRM<- function(r, w, include_scale=FALSE){
  if (class(r)[1] != "RasterLayer") 
    stop("r must be a raster object")
  if(length(w)==1){
    w<- rep(w,2)}
  if(length(w) > 2){
    stop("Specified window exceeds 2 dimensions")}
  if(any(0 == (w %% 2))){
    stop("Error: w must be odd")}
  if(all(w<3)){
    stop("Error: w must be greater or equal to 3 in at least one dimension")
  }
  
  scale.factor <- round(w[1] * w[2], 0)
  sa <- raster::terrain(r, opt=c("slope", "aspect"), unit="radians", 
                        neighbors=8) 					
  sin.slp <- raster::calc(sa[["slope"]], fun=sin)                 # xyRaster 
  cos.slp <- raster::calc(sa[["slope"]], fun=cos)                 # zRaster 
  sin.asp <- raster::calc(sa[["aspect"]], fun=sin) * sin.slp      # yRaster
  cos.asp <- raster::calc(sa[["aspect"]], fun=cos) * sin.slp      # xRaster  
  f = matrix(1,w[1],w[2]) 
  x.sum <- raster::focal(sin.asp, w = f, fun=sum, na.rm=FALSE) #Maybe allow na.rm to be either TRUE or FALSE
  y.sum <- raster::focal(cos.asp, w = f, fun=sum, na.rm=FALSE) 
  z.sum <- raster::focal(cos.slp, w = f, fun=sum, na.rm=FALSE) 
  vrm.fun <- function(x, y, z) { 
    sqrt( (x^2) + (y^2) + (z^2) ) 
  }
  res_vect <- raster::overlay(x.sum, y.sum, z.sum, fun=vrm.fun) #resultant vector
  out<- 1 - (res_vect / scale.factor) #If include na.rm option, adjust scale.factor to be number of non-NA cells.
  names(out)<- "vrm"
  if(include_scale){names(out)<- paste0(names(out), "_", w[1],"x", w[2])} #Add scale to layer names
  return(out)
  }
