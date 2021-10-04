#' Calculates Relative Difference from Mean Value (RDMV)
#'
#' Calculates Relative Difference from Mean Value (RDMV). RDMV = (focal_value - local_mean)/local_range
#' @param r DEM as a raster layer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number. Default is 3x3.
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @param pad logical value specifying whether rows/columns of NA's should be padded to the edge of the raster to remove edge effects (FALSE by default). If pad is TRUE, na.rm must be TRUE.
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @return a RasterLayer
#' @import raster
#' @export

RDMV<- function(r, w=c(3,3), na.rm=FALSE, pad=FALSE, include_scale=FALSE){
  if(length(w)==1){w<- rep(w, 2)}
  if(any(w<3) | any(0 == (w %% 2))){
    stop("Error: w must be odd and greater than or equal to 3")}
  
  w_mat <- matrix(1, nrow=w[1], ncol=w[2])
  
  localmean<- raster::focal(x = r, w= w_mat, fun=mean, na.rm = na.rm, pad=pad)
  localmax<- raster::focal(x = r, w= w_mat, fun=max, na.rm = na.rm, pad=pad)
  localmin<- raster::focal(x = r, w= w_mat, fun=min, na.rm = na.rm, pad=pad)
  rdmv<- (r - localmean)/(localmax-localmin)
  names(rdmv)<- "rdmv"
  
  if(include_scale){names(rdmv)<- paste0(names(rdmv), "_", w[1],"x", w[2])} #Add scale to layer names
  return(rdmv)
}