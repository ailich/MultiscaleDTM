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
#' @param mask_aspect A logical. If TRUE (default), when slope = 0 values, aspect (and therefore northness/eastness) will be set to NA. If FALSE, aspect will to 270 degrees or 3*pi/2 radians (atan2(0,0)-pi/2+2*pi).
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
  
  #Define local coordinate system of window
  x_mat<- matrix(res(r)[2], nrow = w[1], ncol=w[2])
  for (C in 1:w[2]) {
    x_mat[,C]<- x_mat[,C]*C
  }
  x_mat<- x_mat - mean(x_mat)
  x<- as.vector(x_mat)
  
  y_mat<- matrix(res(r)[1], nrow = w[1], ncol=w[2])
  for (R in 1:w[1]) {
    y_mat[R,]<- y_mat[R,]*R
  }
  y_mat<- y_mat - mean(y_mat)
  y<- as.vector(y_mat)
  
  #Explanatory Variable matrix X for quadratic fit
  X<- cbind(1, x^2,y^2, x*y, x, y)
  
  if(pad==TRUE){
    if(na.rm==FALSE){
      na.rm<- TRUE
      warning("if pad=TRUE, na.rm must also be TRUE. Changing na.rm to TRUE")
      }
    og_extent<- raster::extent(r)
    r<- raster::extend(r, c((w[1]-1)/2, (w[2]-1)/2), value=NA)
  }
  
  #Process large rasters as smaller chunks
  run_in_blocks<- !raster::canProcessInMemory(r, n = 8)
  if(run_in_blocks==FALSE){
    params<- raster::brick(r, nl=7, values=FALSE)
    values(params)<- C_multiscale2(r = as.matrix(r), w = w, X=X , na_rm=na.rm) 
    } else{
    
    f_out <- raster::rasterTmpFile()
    params<- raster::brick(r, nl=7, values=FALSE)
    params <- raster::writeStart(params, filename = f_out)
    names(params)<- c("a", "b", "c", "d", "e", "f", "mask")
    
    block_idx<- raster::blockSize(r, n = 8, minblocks = 2, minrows = w[1])
    block_overlap<- (w[1]-1)/2
    nr<- nrow(r)
    nc<- ncol(r)
    for (i in 1:block_idx$n) {
      min_row<- max(c(block_idx$row[[i]] - block_overlap), 1)
      max_row<- min(c(block_idx$row[[i]] + block_idx$nrows[[i]] - 1 + block_overlap, nr))
      curr_block <- raster::getValues(r, row = min_row, nrows = max_row-min_row+1, format="matrix")
      
      out_block<-  C_multiscale2(r = curr_block, w = w, X=X , na_rm=na.rm)  
      #out_block is a formatted as matrix where each column corresponds to a raster layer (this is how writeRaster needs it to be formatted)
      #As you go down rows in this matrix you move across rows in the raster object
      if(i==1){
        out_block<- out_block[1:(nrow(out_block)-(block_overlap*nc)),] #Trim bottom edge of raster
      } else if (i != block_idx$n){
        out_block[(1+block_overlap*nc):(nrow(out_block)-(block_overlap*nc)),] #Trim top and bottom edge of raster
      } else {
        out_block<- out_block[(1+block_overlap*nc):nrow(out_block),] #Trim top edge of raster
      }
      raster::writeValues(params, v= out_block, start= block_idx$row[i])
    }
    params<- raster::writeStop(params)
    }
  names(params)<- c("a", "b", "c", "d", "e", "f", "mask")
  if(pad==TRUE){
    params<- raster::crop(params, og_extent)
  }
  
  mask_raster<- params$mask #mask indicating when all values are the same
  params<- dropLayer(params, 7) #drop mask
  
  #Use regression parameters to calculate slope and aspect
  slp<- atan(sqrt(params$d^2 + params$e^2))
  names(slp)<- "QuadSlope"
  asp<- atan2(params$e,params$d) - pi/2 #Shift aspect so north is zero
  names(asp)<- "QuadAspect"
  
  asp[asp < 0]<- asp[asp < 0] + 2*pi
  asp[asp >= 2*pi]<- asp[asp >= 2*pi] - 2*pi # Constrain aspect from 0 to 2pi
  if (mask_aspect){
    asp[slp==0]<- NA_real_
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
  profc[mask_raster]<- 0
  names(profc)<- "ProfCurv"

  planc<- raster::overlay(params, fun = function(a,b,c,d,e,f)(200 * (b*d^2 +a*e^2 - c*d*e)) / ((e^2+d^2)^1.5))
  planc[mask_raster]<- 0
  names(planc)<- "PlanCurv"

  profc_max<- raster::overlay(params, fun = function(a,b,c,d,e,f) 100*(-a - b + sqrt((a-b)^2+c^2)))
  names(profc_max)<- "ProfCurvMax"
  
  profc_min<- raster::overlay(params, fun = function(a,b,c,d,e,f) 100*(-a - b - sqrt((a-b)^2+c^2)))
  names(profc_min)<- "ProfCurvMin"
  
  mean_curv<- raster::overlay(params, fun = function(a,b,c,d,e,f) 100*(-a - b))
  names(mean_curv)<- "MeanCurv"
  
  #Wood 1996 page 88
  longc<- raster::overlay(params, fun = function(a,b,c,d,e,f) -200 * ((a*d^2 + b*e^2 + c*d*e) / (d^2 + e^2)))
  longc[mask_raster]<- 0
  names(longc)<- "LongCurv"
  
  crosc<- raster::overlay(params, fun = function(a,b,c,d,e,f) -200 * ((b*d^2 + a*e^2 - c*d*e) / (d^2 + e^2)))
  crosc[mask_raster]<- 0
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
