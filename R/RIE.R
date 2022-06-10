#' Calculates Roughness Index-Elevation
#'
#' Calculates Roughness Index-Elevation. This is the standard deviation of residual topography in a focal window where residual topography is calculated as the focal pixel minus the focal mean.
#' @param r DTM as a SpatRaster or RasterLayer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number. Default is 3x3.
#' @param na.rm A logical indicating whether or not to remove NA values before calculation of SD
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @param filename character Output filename.
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @return a SpatRaster or RasterLayer
#' @examples 
#' r<- rast(volcano, extent= ext(2667400, 2667400 + ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), crs = "EPSG:27200")
#' rie<- RIE(r, w=c(5,5), na.rm = TRUE)
#' plot(rie)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @details Note the original paper by Cavalli et al (2008) uses a fixed 5x5 window and uses 25 as the denominator indicating use of the population standard deviation. This implementation provides a flexible window size and istead calculates the sample standard deviation which uses a denominator of n-1.
#' @references 
#' Cavalli, M., Tarolli, P., Marchi, L., Dalla Fontana, G., 2008. The effectiveness of airborne LiDAR data in the recognition of channel-bed morphology. CATENA 73, 249–260. https://doi.org/10.1016/j.catena.2007.11.001
#' @export
#' 

RIE<- function(r, w=c(3,3), na.rm=FALSE, include_scale=FALSE, filename=NULL, overwrite=FALSE, wopt=list()){
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
  
  resid_topo<- r - terra::focal(x = r, w = w, fun = mean, na.rm = FALSE, wopt=wopt) #Mean is only equal to fitted value if all values present
  rie<- terra::focal(x = resid_topo, w = w, fun = sd, na.rm = na.rm, wopt=wopt)
  names(rie)<- "rie"
  if(include_scale){names(rie)<- paste0(names(rie), "_", w[1],"x", w[2])} #Add scale to layer names
  
  #Return
  if(og_class =="RasterLayer"){
    rie<- raster::raster(rie)
    if(!is.null(filename)){
      return(raster::writeRaster(rie, filename=filename, overwrite=overwrite))
    }
  }
  if(!is.null(filename)){
    return(terra::writeRaster(rie, filename=filename, overwrite=overwrite, wopt=wopt))
  }
  return(rie)
}