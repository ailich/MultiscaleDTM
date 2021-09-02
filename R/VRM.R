#' Implementation of the Sappington et al., (2007) vector ruggedness measure
#'
#' Implementation of the Sappington et al., (2007) vector ruggedness measure. This is a simple wrapper around the vrm function from spatialEco by Jeffrey S Evans, to make the syntax consistent with the rest of the MultiscaleDEM package, and is included in this package since VRM is also a rugosity measure that corrects for slope. If you use this function, please cite the spatialEco package.
#' @param r DEM as a raster layer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @return a RasterLayer
#' @import raster
#' @importFrom spatialEco vrm
#' @export

VRM<- function(r, w, include_scale){
  if(length(w==1)){
    w<- rep(w,2)}
  out<- spatialEco::vrm(x=r, s=w)
  names(out)<- "vrm"
  if(include_scale){names(out)<- paste0(names(out), "_", w[1],"x", w[2])} #Add scale to layer names
  return(out)
}




