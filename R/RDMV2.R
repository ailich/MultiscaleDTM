#' Calculates Relative Difference from Mean Value (RDMV)
#'
#' Calculates Relative Difference from Mean Value (RDMV). RDMV = (focal_value - local_mean)/local_range
#' @param r DEM as a SpatRaster, RasterLayer, RasterStack, or RasterBrick
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number. Default is 3x3.
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @return a SpatRaster or RasterLayer
#' @import terra
#' @importFrom raster raster
#' @references 
#' Lecours, V., Devillers, R., Simms, A.E., Lucieer, V.L., Brown, C.J., 2017. Towards a Framework for Terrain Attribute Selection in Environmental Studies. Environmental Modelling & Software 89, 19-30. https://doi.org/10.1016/j.envsoft.2016.11.027
#' @export

RDMV2<- function(r, w=c(3,3), na.rm=FALSE, pad=FALSE, include_scale=FALSE){
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
  
  localmean<- terra::focal(x = r, w= w, fun=mean, na.rm = na.rm)
  localmax<- terra::focal(x = r, w= w, fun=max, na.rm = na.rm)
  localmin<- terra::focal(x = r, w= w, fun=min, na.rm = na.rm)
  rdmv<- (r - localmean)/(localmax-localmin)
  names(rdmv)<- "rdmv"
  
  if(include_scale){names(rdmv)<- paste0(names(rdmv), "_", w[1],"x", w[2])} #Add scale to layer names
  if(og_class=="RasterLayer"){rdmv<- raster::raster(rdmv)}
  return(rdmv)
}