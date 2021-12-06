#' Calculates Topographic Position Index
#'
#' Calculates Topographic Position Index (TPI). This is the value of the focal pixel minus the mean of the surrounding pixels (i.e. local mean but excluding the value of the focal pixel).
#' @param r DEM as a SpatRaster or RasterLayer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number. Default is 3x3.
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @return a SpatRaster or RasterLayer
#' @import terra
#' @importFrom raster raster
#' @references 
#' Weiss, A., 2001. Topographic Position and Landforms Analysis. Presented at the ESRI user conference, San Diego, CA.
#' @export
#' 

TPI2<- function(r, w=c(3,3), na.rm=FALSE, include_scale=FALSE){
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
  if(length(w)==1){w<- rep(w, 2)}
  if(length(w) > 2){
    stop("Specified window exceeds 2 dimensions")}
  if(any(0 == (w %% 2))){
    stop("Error: w must be odd")}
  if(all(w<3)){
    stop("Error: w must be greater or equal to 3 in at least one dimension")
  }
  
  w_mat <- matrix(1, nrow=w[1], ncol=w[2])
  center_idx<- ceiling(0.5 * length(w_mat))
  w_mat[center_idx] <- NA_real_
  tpi<- r - terra::focal(x = r, w = w_mat, fun = mean, na.rm = na.rm)
  names(tpi)<- "TPI"
  if(include_scale){names(tpi)<- paste0(names(tpi), "_", w[1],"x", w[2])} #Add scale to layer names
  if(og_class=="RasterLayer"){tpi<- raster::raster(tpi)}
  return(tpi)
}