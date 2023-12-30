#' Calculates standard deviation of bathymetry (a measure of rugosity) adjusted for slope
#'
#' Calculates standard deviation of bathymetry (a measure of rugosity). Using a sliding rectangular window a plane is fit to the data and the standard deviation of the residuals is calculated (Ilich et al., 2023)
#' @param r DTM as a SpatRaster or RasterLayer in a projected coordinate system where map units match elevation/depth units
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param na.rm A logical indicating whether or not to remove NA values before calculations
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @param filename character Output filename.
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @return a SpatRaster or RasterLayer of adjusted rugosity
#' @examples
#' r<- rast(volcano, extent= ext(2667400, 2667400 + 
#' ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
#' crs = "EPSG:27200")
#' adjsd<- AdjSD(r, w=c(5,5), na.rm = TRUE)
#' plot(adjsd)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @references
#' Ilich, A. R., Misiuk, B., Lecours, V., & Murawski, S. A. (2023). MultiscaleDTM: An open-source R package for multiscale geomorphometric analysis. Transactions in GIS, 27(4). https://doi.org/10.1111/tgis.13067
#' @export

AdjSD<- function(r, w=c(3,3), na.rm=FALSE, include_scale=FALSE, filename=NULL, overwrite=FALSE, wopt=list()){
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
  x_mat<- matrix(data = seq(from = (-xres(r) * floor(w[2]/2)), to = (xres(r) * floor(w[2]/2)), length.out = w[2]), nrow = w[1], ncol=w[2], byrow=TRUE)
  x<- as.vector(t(x_mat)) #Transpose for focal
  y_mat<- matrix(data = seq(from = (yres(r) * floor(w[1]/2)), to = (-yres(r) * floor(w[1]/2)), length.out = w[1]), nrow = w[1], ncol=w[2])
  y<- as.vector(t(y_mat)) #Transpose for focal
  
  #Explanatory Variable matrix X for linear fit
  X<- cbind(x, y, 1) #Z = dx+ey+f
  
  if(!na.rm){
    Xt<- t(X)
    XtX_inv<- solve(Xt %*% X)
  }
  
  #Fit Quadratic and Extract Residuals
  if(na.rm){
    out<- terra::focalCpp(r, w=w, fun = C_AdjSD_narmT, X_full= X, fillvalue=NA, wopt=wopt)
  } else{
    out<- terra::focalCpp(r, w=w, fun = C_AdjSD_narmF, X= X, Xt=Xt, XtX_inv=XtX_inv, fillvalue=NA, wopt=wopt)
    }
  
  names(out)<- "adjSD"
  if(include_scale){names(out)<- paste0(names(out), "_", w[1],"x", w[2])} #Add scale to layer names
  
  # Return
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
