#' Calculates standard deviation of bathymetry (a measure of rugosity) adjusted for slope
#'
#' Calculates standard deviation of bathymetry (a measure of rugosity). Using a sliding rectangular window a plane is fit to the data and the standard deviation of the residuals is calculated.
#' @param r DEM as a raster layer in a projected coordinate system where map units match elevation/depth units
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @return a RasterLayer of adjusted rugosity
#' @import terra
#' @importFrom raster raster
#' @export

AdjSD2<- function(r, w=c(3,3), na.rm=FALSE, include_scale=FALSE){
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
  if(terra::is.lonlat(r, perhaps=FALSE)){
    stop("Error: Coordinate system is Lat/Lon. Coordinate system must be projected with elevation/depth units matching map units.")
  }
  if(terra::is.lonlat(r, perhaps=TRUE, warn=FALSE)){
    warning("Coordinate system may be Lat/Lon. Please ensure that the coordinate system is projected with elevation/depth units matching map units.")
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
  if(prod(w) < 4){
    stop("Error: Window size must have at least 4 cells to fit surface and have residuals")
  }
  
  #Define local coordinate system of window
  x_mat<- matrix(terra::res(r)[2], nrow = w[1], ncol=w[2])
  for (C in 1:w[2]) {
    x_mat[,C]<- x_mat[,C]*C
  }
  x_mat<- x_mat - mean(x_mat)
  x<- as.vector(t(x_mat)) #Transpose for focal
  
  y_mat<- matrix(terra::res(r)[1], nrow = w[1], ncol=w[2])
  for (R in 1:w[1]) {
    y_mat[R,]<- y_mat[R,]*(w[1]-(R-1)) #Have up be positive
  }
  y_mat<- y_mat - mean(y_mat)
  y<- as.vector(t(y_mat)) #Transpose for focal
  
  #Explanatory Variable matrix X for linear fit
  X<- cbind(x, y, 1) #Z = dx+ey+f
  
  #Fit Quadratic and Extract Residuals
  out<- terra::focalCpp(r, w=w, fun = C_AdjSD, X_full= X, na_rm=na.rm, fillvalue=NA, expand=FALSE)
  
  names(out)<- "adjSD"
  if(include_scale){names(out)<- paste0(names(out), "_", w[1],"x", w[2])} #Add scale to layer names
  
  if(og_class=="RasterLayer"){out<- raster::raster(out)}
  return(out)
}