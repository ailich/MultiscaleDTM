#' Calculates multiscale slope and aspect using a local planar fit.
#'
#' Calculates multiscale slope and aspect of a DTM over a sliding rectangular window using a using a local planar fit to the surface (Sharpnack and Akin 1969).
#' @param r DTM as a SpatRaster (terra) or RasterLayer (raster) in a projected coordinate system where map units match elevation/depth units (up is assumed to be north for calculations of aspect, northness, and eastness).
#' @param w Vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number. Default is 3x3.
#' @param unit "degrees" or "radians".
#' @param metrics Character vector specifying which terrain attributes to return. The default is to return c("pslope", "paspect", "peastness", and "pnorthness"). These are preceded with a 'p' to differentiate them from the measures calculated by SlpAsp() and 'Qfit' where the 'p' indicates that a planar surface was used for the calculation. Additional measures available include "dz.dx" and "dz.dy" which are the x and y components of slope respectively.
#' @param na.rm Logical indicating whether or not to remove NA values before calculations.
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @param mask_aspect Logical. If TRUE (default), aspect will be set to NA and northness and eastness will be set to 0 when slope = 0. If FALSE, aspect is set to 270 degrees or 3*pi/2 radians ((-pi/2)- atan2(0,0)+2*pi) and northness and eastness will be calculated from this.
#' @param filename character Output filename. Can be a single filename, or as many filenames as there are layers to write a file for each layer
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @return a SpatRaster (terra) or RasterStack/RasterLayer (raster)
#' @examples
#' r<- erupt()
#' pmetrics<- Pfit(r, w = c(5,5), unit = "degrees", na.rm = TRUE)
#' plot(pmetrics)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster stack
#' @importFrom raster writeRaster
#' @details
#' If only first order derivatives are needed, Pfit is faster than Qfit and should provide equivalent results to Qfit for first order derivatives (Jones, 1998) when na.rm=FALSE and approximately the same results otherwise. 
#' 
#' @references
#' Evans, I.S., 1980. An integrated system of terrain analysis and slope mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295.
#' Jones, K. H. (1998). A comparison of algorithms used to compute hill slope as a property of the DEM. Computers & Geosciences, 24(4), 315–323. https://doi.org/10.1016/S0098-3004(98)00032-6
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' @export

Pfit<- function(r, w=c(3,3), unit= "degrees", metrics = c("pslope", "paspect", "peastness", "pnorthness"), na.rm=FALSE, include_scale=FALSE, mask_aspect=TRUE, filename=NULL, overwrite=FALSE, wopt=list()){
  all_metrics<- c("pslope", "paspect", "peastness", "pnorthness", "dz.dx", "dz.dy")
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
  
  # Increase tolerance of inputs to metrics
  metrics<- tolower(metrics) #Make all lowercase
  metrics[metrics=="slope"]<- "pslope" #replace slope with pslope
  metrics[metrics=="aspect"]<- "paspect" #replace aspect with paspect
  metrics[metrics=="eastness"]<- "peastness" #replace eastness with peastness
  metrics[metrics=="northness"]<- "pnorthness" #replace northness with pnorthness
 
  if (any(!(metrics %in% all_metrics))){
    stop("Error: Invlaid metric. Valid metrics include 'pslope', 'paspect', 'peastness', 'pnorthness', 'dz.dx', and 'dz.dy'.")
  }

  needed_metrics<- metrics
  
  if(any(c("peastness", "pnorthness") %in% needed_metrics) & (!("paspect" %in% needed_metrics))){
    needed_metrics<- c(needed_metrics, "paspect")
  }
  
  if(mask_aspect & ("paspect" %in% needed_metrics) & (!("pslope" %in% needed_metrics))){
    needed_metrics<- c(needed_metrics, "pslope")
  }
  
  if(any(c("pslope", "paspect", "peastness", "pnorthness") %in% needed_metrics) & !("dz.dx" %in% needed_metrics)){
    needed_metrics<- c(needed_metrics, "dz.dx")
  }
  
  if(any(c("pslope", "paspect", "peastness", "pnorthness") %in% needed_metrics) & !("dz.dy" %in% needed_metrics)){
    needed_metrics<- c(needed_metrics, "dz.dy")
  }
  
  #Define local coordinate system of window
 
  x_mat<- matrix(data = seq(from = (-xres(r) * floor(w[2]/2)), to = (xres(r) * floor(w[2]/2)), length.out = w[2]), nrow = w[1], ncol=w[2], byrow=TRUE)
  x<- as.vector(t(x_mat)) #Transpose for focal
  
  y_mat<- matrix(data = seq(from = (yres(r) * floor(w[1]/2)), to = (-yres(r) * floor(w[1]/2)), length.out = w[1]), nrow = w[1], ncol=w[2])
  y<- as.vector(t(y_mat)) #Transpose for focal
  
  # forcing through center doesn't affect slope (Jones, 1998) so don't bother estimating the extra coefficient unless na.rm=TRUE
  # If na.rm = TRUE, you want to be able to get the slope even if the central cell is NA
  X<- cbind(x, y) # Z = dx+ey (no intercept)
  
  if(na.rm){
    X<- cbind(X,1)
  } # Z = dx+ey+f (with intercept)
  
  if(!na.rm){
    Xt<- t(X)
    XtX_inv<- solve(Xt %*% X)
  }
  
  # Calculate Regression Parameters
  if(na.rm){
    params<- terra::focalCpp(r, w=w, fun = C_Qfit1_narmT, X_full= X, return_intercept = FALSE, fillvalue=NA, wopt=wopt)
  } else{
    params<- terra::focalCpp(r, w=w, fun = C_Qfit2_narmF, X= X, Xt= Xt, XtX_inv= XtX_inv, fillvalue=NA, wopt=wopt)
  }
  names(params)<- c("dz.dx", "dz.dy")
  
  out<- terra::rast() #Initialize output
  
  if("dz.dx" %in% needed_metrics){
    out<- c(out, params$dz.dx, warn=FALSE)
  }
  
  if("dz.dy" %in% needed_metrics){
    out<- c(out, params$dz.dy, warn=FALSE)
  }
  
  #Use regression parameters to calculate slope and aspect
  if("pslope" %in% needed_metrics){
    slp<- terra::math(terra::math(params$dz.dx^2 + params$dz.dy^2, fun="sqrt", wopt=wopt), fun="atan", wopt=wopt)
    if(unit=="degrees"){
      slp<- slp*(180/pi)
    }
    names(slp)<- "pslope"
    out<- c(out, slp, warn=FALSE)
  }
  
  if("paspect" %in% needed_metrics){
    asp<- (-pi/2) - terra::atan_2(params$dz.dy,params$dz.dx, wopt=wopt) # aspect relative to North
    asp<- ifel(asp < 0, yes = asp+(2*pi), no= asp, wopt=wopt) # Constrain range so between 0 and 2pi
    
    if (mask_aspect & ("paspect" %in% metrics)){
      asp<- terra::mask(asp, mask= slp, maskvalues = 0, updatevalue = NA, wopt=wopt) #Set aspect to undefined where slope is zero (doesn't need to be needed_metrics since also mask northness and eastness as 0).
    }
    
    if("peastness" %in% needed_metrics){
      eastness<- sin(asp)
      if (mask_aspect){
        eastness<- terra::mask(eastness, mask= slp, maskvalues = 0, updatevalue = 0, wopt=wopt) #Set eastness to 0 where slope is zero
      }
      names(eastness)<- "peastness"
      out<- c(out, eastness, warn=FALSE)
    }
    
    if("pnorthness" %in% needed_metrics){
      northness<- cos(asp)
      if (mask_aspect){
        northness<- terra::mask(northness, mask= slp, maskvalues = 0, updatevalue = 0, wopt=wopt) #Set northness to 0 where slope is zero
      }
      names(northness)<- "pnorthness"
      out<- c(out, northness, warn=FALSE)
    }
    if(unit=="degrees"){
      asp<- asp*180/pi
    }
    names(asp)<- "paspect"
    out<- c(out, asp, warn=FALSE)
  }
  
  out<- terra::subset(out, metrics, wopt=wopt) #Subset needed metrics to requested metrics in proper order
  
  if(include_scale){names(out)<- paste0(names(out), "_", w[1],"x", w[2])} #Add scale to layer names
  
  #Return
  if(og_class=="RasterLayer"){
    if(terra::nlyr(out) > 1){
      out<- raster::stack(out) #Convert to RasterStack
      if(!is.null(filename)){
        if(length(filename)==1){
          return(raster::writeRaster(out, filename=filename, overwrite=overwrite, bylayer=FALSE))
        } else{
          return(raster::writeRaster(out, filename=filename, overwrite=overwrite, bylayer=TRUE))
        }
      }
    } else{
      out<- raster::raster(out)
      if(!is.null(filename)){
        return(raster::writeRaster(out, filename=filename, overwrite=overwrite))
      }
    }
  }
  if(!is.null(filename)){
    return(terra::writeRaster(out, filename=filename, overwrite=overwrite, wopt=wopt))
  }
  return(out)
}
