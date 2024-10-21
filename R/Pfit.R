#' Calculates multiscale slope using a local planar fit
#'
#' Calculates multiscale slope, aspect, curvature, and morphometric features of a DTM over a sliding rectangular window using a local quadratic fit to the surface (Evans, 1980; Wood, 1996).
#' @param r DTM as a SpatRaster (terra) or RasterLayer (raster) in a projected coordinate system where map units match elevation/depth units (up is assumed to be north for calculations of aspect, northness, and eastness).
#' @param w Vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number. Default is 3x3.
#' @param unit "degrees" or "radians".
#' @param na.rm Logical indicating whether or not to remove NA values before calculations.
#' @param filename character Output filename. Can be a single filename, or as many filenames as there are layers to write a file for each layer
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @return a SpatRaster (terra) or RasterStack/RasterLayer (raster)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster stack
#' @importFrom raster writeRaster
#' @references
#' Evans, I.S., 1980. An integrated system of terrain analysis and slope mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295.
#' 
#' Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
#' 
#' Wilson, M.F., O’Connell, B., Brown, C., Guinan, J.C., Grehan, A.J., 2007. Multiscale Terrain Analysis of Multibeam Bathymetry Data for Habitat Mapping on the Continental Slope. Marine Geodesy 30, 3-35. https://doi.org/10.1080/01490410701295962
#' 
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' @export

Pfit<- function(r, w=c(3,3), unit= "degrees",na.rm=FALSE, include_scale=FALSE, filename=NULL, overwrite=FALSE, wopt=list()){
  
  og_class<- class(r)[1]
  if(og_class=="RasterLayer"){
    r<- terra::rast(r) #Convert to SpatRaster
  }
  
  # Input checks
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
  if(length(w)==1){
    w<- rep(w,2)}
  if(any(0 == (w %% 2))){
    stop("Error: w must be odd")}
  if(length(w) > 2){
    stop("Specified window exceeds 2 dimensions")}
  if(all(w<3)){
    stop("Error: w must be greater or equal to 3 in at least one dimension")
  }
  
  unit<- tolower(unit) #make lowercase
  if (!any(unit==c("degrees", "radians"))){
    stop("unit must be 'degrees' or 'radians'")
  }
 
  #Define local coordinate system of window
 
  x_mat<- matrix(data = seq(from = (-xres(r) * floor(w[2]/2)), to = (xres(r) * floor(w[2]/2)), length.out = w[2]), nrow = w[1], ncol=w[2], byrow=TRUE)
  x_mat[-c(1, nrow(x_mat)), -c(1, ncol(x_mat))] <- NA
  x<- as.vector(t(x_mat)) #Transpose for focal
  
  y_mat<- matrix(data = seq(from = (yres(r) * floor(w[1]/2)), to = (-yres(r) * floor(w[1]/2)), length.out = w[1]), nrow = w[1], ncol=w[2])
  y_mat[-c(1, nrow(y_mat)), -c(1, ncol(y_mat))] <- NA
  y<- as.vector(t(y_mat)) #Transpose for focal
  
  idx<- !is.na(x)
  
  X<- cbind(x, y, 1) #Z = dx+ey+f
  X<- X[idx,]
  
  if(!na.rm){
    Xt<- t(X)
    XtX_inv<- solve(Xt %*% X)
  }
  
  # Calculate Regression Parameters

  params<- terra::focalCpp(r, w=w, fun = C_Pfit1_narmF, X= X, Xt= Xt, XtX_inv= XtX_inv, idx=idx, fillvalue=NA, wopt=wopt)
  slp<- terra::math(terra::math(params$d^2 + params$e^2, fun="sqrt", wopt=wopt), fun="atan", wopt=wopt)
  if(unit=="degrees"){
    slp<- slp*(180/pi)
  }
  return(slp)
}
