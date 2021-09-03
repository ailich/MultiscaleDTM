#' Calculates multiscale slope, aspect, curvature, and morphometric features
#'
#' Calculates multiscale slope, aspect, curvature, and morphometric features of a DEM over a sliding rectangular window using a quadratic fit to the surface usin the Wood/Evans method according to the equations of Wood 1996. This is an R/C++ implementation of r.param.scale GRASS GIS function. Note however, for aspect, 0 degrees represents north and increases clockwise which differs from the way r.param.scale reports aspect.  Additionally, compared to r.param.scale, curvature has been multiplied by 100 to express curvature as percent gradient per unit length (Albani et al 2004).
#' @param r DEM as a raster layer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number. Default is 3x3.
#' @param unit "degrees" or "radians"
#' @param return_aspect logical indicating whether or not to return aspect in addition to Northness and Eastness (default is FALSE)
#' @param slope_tolerance Slope tolerance that defines a 'flat' surface (degrees; default is 1.0). Relevant for the features layer.
#' @param curvature_tolerance Curvature tolerance that defines 'planar' surface (default is 0.01). Relevant for the features layer.
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @param pad logical value specifying whether rows/columns of NA's should be padded to the edge of the raster to remove edge effects (FALSE by default). If pad is TRUE, na.rm must be TRUE.
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @param mask_aspect A logical. If slope evaluates to 0, aspect (and therefore northness/eastness) will be set to NA when mask_aspect is TRUE (the default). Theoretically, if FALSE this aspect would evaluate to -90 degrees or -pi/2 radians (atan2(0,0)-pi/2) however it can be sporadic due to numerical precision issues related to the matrix math used to derive the regression parameters. Additionally, slope (in radians) will be rounded to 16 digits before checking if it equals zero due to precision issues.
#' @param return_params logical indicating whether to return Wood/Evans regression parameters (default FALSE)
#' @return a RasterStack
#' @import raster
#' @export

WoodEvans<- function(r, w=c(3,3), unit= "degrees", return_aspect= FALSE, slope_tolerance=1, curvature_tolerance=0.01, na.rm=FALSE, pad=FALSE, include_scale=FALSE, mask_aspect=TRUE, return_params= FALSE){
  #Input checks
  if(length(w==1)){
    w<- rep(w,2)}
  if(any(w<3) | any(0 == (w %% 2))){
    stop("Error: w must be odd and greater than or equal to 3")}
  if (!any(unit==c("degrees", "radians"))){
    stop("unit must be 'degrees' or 'radians'")
  }
  if(pad==TRUE){
    if(na.rm==FALSE){
      na.rm<- TRUE
      warning("if pad=TRUE, na.rm must also be TRUE. Changing na.rm to TRUE")
      }
    og_extent<- raster::extent(r)
    r<- raster::extend(r, c(w[1], w[2]), value=NA)
  }
  
  #Process large rasters as smaller chunks
  run_in_blocks<- !raster::canProcessInMemory(r, n = 7)
  if(run_in_blocks==FALSE){
    params<- WoodEvansHelper(r=r, w=w, type= 2, na.rm=na.rm)
  } else{
    block_idx<- raster::blockSize(r, n = 6, minblocks = 2, minrows = w[1])
    out_blocks<- vector(mode = "list", length = block_idx$n)
    block_overlap<- w[1]-1
    for (i in 1:block_idx$n) {
      min_row<- block_idx$row[[i]]
      max_row<- min(min_row + block_idx$nrows[[i]] - 1 + block_overlap, nrow(r))
      block_extent<- raster::extent(r, min_row, max_row, 1, ncol(r))
      curr_block<- raster::crop(r, block_extent)
      out_blocks[[i]]<- WoodEvansHelper(r=curr_block, w=w, type=2, na.rm = na.rm)
    }
    params<- do.call(raster::merge, out_blocks)
    names(params)<- c("a", "b", "c", "d", "e", "f", "SD_resid")
  }
  if(pad==TRUE){
    params<- raster::crop(params, og_extent)
  }
  
  params<- raster::dropLayer(params, 7) #Drop SD_resid layer
  #Use regression parameters to calculate slope and aspect
  slp<- atan(sqrt(params$d^2 + params$e^2))
  names(slp)<- "QuadSlope"
  asp<- atan2(params$e,params$d) - pi/2 #Shift aspect so north is zero
  names(asp)<- "QuadAspect"
  
  asp[asp < 0]<- asp[asp < 0] + 2*pi
  asp[asp >= 2*pi]<- asp[asp >= 2*pi] - 2*pi # Constrain aspect from 0 to 2pi
  if (mask_aspect){
    asp[round(slp,16)==0]<- NA #Mask out outspect values where slope is 0
  }
  
  eastness<- sin(asp)
  names(eastness)<- "QuadEastness"
  northness<- cos(asp)
  names(northness)<- "QuadNorthness"
  
  if(unit=="degrees"){
    slp<- slp*180/pi
    asp<- asp*180/pi
  } else{
    slope_tolerance<- slope_tolerance * (pi/180)
    } #Convert slope tolerance to radians if unit is not degrees
  
  #Curvature
  #Note curvature has been multiplied by 100 to express curvature as percent gradient per unit length (Albani et al 2004)
  
  #  Wood 1996 page 86 & Wilson 2007 page 9-10
  profc<- raster::overlay(params, fun = function(a,b,c,d,e,f)(-200 * (a*d^2 + b*e^2 + c*d*e)) / ((e^2 + d^2)*(1 + e^2 + d^2)^1.5))
  names(profc)<- "ProfCurv"
  
  planc<- raster::overlay(params, fun = function(a,b,c,d,e,f)(200 * (b*d^2 +a*e^2 - c*d*e)) / ((e^2+d^2)^1.5))
  names(planc)<- "PlanCurv"
  
  profc_max<- raster::overlay(params, fun = function(a,b,c,d,e,f) 100*(-a - b + sqrt((a-b)^2+c^2)))
  names(profc_max)<- "ProfCurvMax"
  
  profc_min<- raster::overlay(params, fun = function(a,b,c,d,e,f) 100*(-a - b - sqrt((a-b)^2+c^2)))
  names(profc_min)<- "ProfCurvMin"
  
  mean_curv<- raster::overlay(params, fun = function(a,b,c,d,e,f) 100*(-a - b))
  names(mean_curv)<- "MeanCurv"
  
  #Wood 1996 page 88
  longc<- raster::overlay(params, fun = function(a,b,c,d,e,f) -200 * ((a*d^2 + b*e^2 + c*d*e) / (d^2 + e^2)))
  names(longc)<- "LongCurv"
  
  crosc<- raster::overlay(params, fun = function(a,b,c,d,e,f) -200 * ((b*d^2 + a*e^2 - c*d*e) / (d^2 + e^2)))
  names(crosc)<- "CrossCurv"
  
  #Wood 1996 page 115
  max_curv<- raster::overlay(params, fun = function(a,b,c,d,e,f) 100 * (-a - b + sqrt((a-b)^2+c^2)))
  names(max_curv)<- "MaxCurv"
  
  min_curv<- raster::overlay(params, fun = function(a,b,c,d,e,f) 100 * (-a - b - sqrt((a-b)^2+c^2)))
  names(min_curv)<- "MinCurv"
 
  #Morphometric Features (Wood 1996, Page 120)
  
  classify_features<- classify_features_ff(slope_tolerance, curvature_tolerance) #Define classification function based on slope and curvature tolerance
  features<- raster::overlay(slp, crosc, max_curv, min_curv, fun = classify_features)
  features<- as.factor(features)
  suppressWarnings(levels(features)[[1]]<- data.frame(ID=1:6)) #Make sure all factor levels are present even if it wasn't in the original raster
  levels(features)[[1]]$Feature<- c("Planar", "Pit", "Channel", "Pass", "Ridge", "Peak")
  names(features)<- "Features"
  
  out<- stack(slp, eastness, northness, profc, planc, profc_max, profc_min, mean_curv, max_curv, min_curv, longc, crosc, features)
  if(return_aspect){out<- stack(out, asp)}
  if(return_params){out<- stack(out, params)}
  if(include_scale){names(out)<- paste0(names(out), "_", w[1],"x", w[2])} #Add scale to layer names
  return(out)
}
