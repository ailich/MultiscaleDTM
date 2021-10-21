#' Calculates Topographic Position Index
#'
#' Calculates Topographic Position Index (TPI). This is the value of the focal pixel minus the mean of the surrounding pixels (i.e. local mean but excluding the value of the focal pixel).
#' @param r DEM as a raster layer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number. Default is 3x3.
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @param pad logical value specifying whether rows/columns of NA's should be padded to the edge of the raster to remove edge effects (FALSE by default). If pad is TRUE, na.rm must be TRUE.
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @return a RasterLayer
#' @import raster
#' @references 
#' Weiss, A., 2001. Topographic Position and Landforms Analysis. Presented at the ESRI user conference, San Diego, CA.
#' @export
#' 

TPI<- function(r, w=c(3,3), na.rm=FALSE, pad=FALSE, include_scale=FALSE){
  if(length(w)==1){w<- rep(w, 2)}
  
  if(any(w<3) | any(0 == (w %% 2))){
    stop("Error: w must be odd and greater than or equal to 3")}
  
  w_mat <- matrix(1, nrow=w[1], ncol=w[2])
  center_idx<- ceiling(0.5 * length(w_mat))
  w_mat[center_idx] <- 0 
  if(!na.rm){
    w_mat[w_mat>0]<- 1/sum(w_mat>0)
    tpi<- r - raster::focal(x=r, w=w_mat, fun=sum, na.rm=FALSE, pad=FALSE)
    } else{
      custom_mean_fun<- function(x, na.rm) {
        return(mean(x[-center_idx], na.rm=na.rm))
        } #Using terra::focal with NA weights matrix will be better
      tpi<- r - raster::focal(x = r, w = w_mat, fun = custom_mean_fun, na.rm = TRUE, pad = pad)
    }
  names(tpi)<- "TPI"
  if(include_scale){names(tpi)<- paste0(names(tpi), "_", w[1],"x", w[2])} #Add scale to layer names
return(tpi)
}