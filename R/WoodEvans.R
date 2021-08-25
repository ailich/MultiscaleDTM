#' Calculates multiscale slope, aspect, curvature, and morphometric features
#'
#' Calculates Calculates multiscale slope, aspect, curvature, and morphometric features of a DEM over a sliding rectangular window using a quadratic fit to the surface usin the Wood/Evans method according to the equations of Wood 1996. This is an R/C++ implementation of r.param.scale GRASS GIS function. Note however, for aspect, 0 degrees represents north and increases clockwise which differs from the way r.param.scale reports aspect.  Additionally, compared to r.param.scale, curvature has been multiplied by 100 to express curvature as percent gradient per unit length (Albani et al 2004).
#' @param r DEM as a raster layer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param unit "degrees" or "radians"
#' @param return_aspect logical indicating whether or not to return aspect in addition to Northness and Eastness (default is FALSE)
#' @param slope_tolerance Slope tolerance that defines a 'flat' surface (degrees; default is 1.0). Relevant for the features layer.
#' @param curvature_tolerance Curvature tolerance that defines 'planar' surface (default is 0.01). Relevant for the features layer.
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @param pad logical value specifying whether rows/columns of NA's should be padded to the edge of the raster to remove edge effects (FALSE by default). If pad is TRUE, na.rm must be TRUE.
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @param return_params logical indicating whether to return Wood/Evans regression parameters (default FALSE)
#' @return a RasterStack
#' @import raster
#' @export

WoodEvans<- function(r, w, unit= "degrees", return_aspect= FALSE, slope_tolerance=1, curvature_tolerance=0.01, na.rm=FALSE, pad=FALSE, include_scale=FALSE, return_params= FALSE){
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
  profc<- (-200*(params$a*params$d^2 + params$b*params$e^2 + params$c*params$d*params$e))/((params$e^2+params$d^2)*(1+params$e^2+params$d^2)^1.5)
  names(profc)<- "ProfCurv"
  
  planc<- (200 * (params$b*params$d^2 + params$a*params$e^2 - params$c*params$d*params$e))/((params$e^2 + params$d^2)^1.5)
  names(planc)<- "PlanCurv"
  
  profc_max<- 100 * (-params$a - params$b + sqrt((params$a-params$b)^2 + params$c^2))
  names(profc_max)<- "ProfCurvMax"
  
  profc_min<- 100* (-params$a - params$b - sqrt((params$a-params$b)^2 + params$c^2))
  names(profc_min)<- "ProfCurvMin"
  
  mean_curv<- 100* (-params$a - params$b)
  names(mean_curv)<- "MeanCurv"
  
  max_curv<- 100*(-params$a - params$b + sqrt((params$a-params$b)^2 + params$c^2))
  names(max_curv)<- "MaxCurv"
  
  min_curv<- 100*(-params$a - params$b - sqrt((params$a-params$b)^2 + params$c^2))
  names(min_curv)<- "MinCurv"
  
  #Wood 1996 page 88
  longc<- -200*((params$a*params$d^2 + params$b*params$e^2 + params$c*params$d*params$e)/(params$d^2+params$e^2)) 
  names(longc)<- "LongCurv"
  crosc<- -200*((params$b*params$d^2 + params$a*params$e^2 - params$c*params$d*params$e)/(params$d^2+params$e^2))
  names(crosc)<- "CrossCurv"
  
  #Morphometric Features (Wood 1996, Page 120)
  features<- r
  values(features)<- 1 #Planar
  names(features)<- "Features"
  
  features[min_curv < -curvature_tolerance]<- 3 #Channel
  features[(min_curv < -curvature_tolerance) & (max_curv < -curvature_tolerance)]<- 2 #Pit
  
  features[max_curv > curvature_tolerance]<- 5 #Ridge
  features[(max_curv > curvature_tolerance) & (min_curv < -curvature_tolerance)]<- 4 #Pass
  features[(max_curv > curvature_tolerance) & (min_curv > curvature_tolerance)]<- 6 #Peak
  
  features[slp > slope_tolerance]<- 1 #Planar
  features[(slp > slope_tolerance) & (crosc < -curvature_tolerance)]<- 3 #Channel
  features[(slp > slope_tolerance) & (crosc > curvature_tolerance)]<- 5 #Ridge
  
  features[is.na(slp)]<- NA
  
  features<- as.factor(features)
  levels(features)[[1]]$Feature<- c("Planar", "Pit", "Channel", "Pass", "Ridge", "Peak")
  
  
  out<- stack(slp, eastness, northness, profc, planc, profc_max, profc_min, mean_curv, max_curv, min_curv, longc, crosc, features)
  if(return_aspect){out<- stack(out, asp)}
  if(return_params){out<- stack(out, params)}
  if(include_scale){names(out)<- paste0(names(out), "_", w[1],"x", w[2])} #Add scale to layer names
  return(out)
}
