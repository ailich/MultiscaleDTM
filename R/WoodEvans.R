#' Calculates multiscale slope, aspect, curvature, and morphometric features
#'
#' Calculates multiscale slope, aspect, curvature/convexity, and morphometric features of a DEM over a sliding rectangular window using a quadratic fit to the surface (Evans, 1980; Wood, 1996).
#' @param r DEM as a raster layer in a projected coordinate system where map units match elevation/depth units.
#' @param w Vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number. Default is 3x3.
#' @param unit "degrees" or "radians".
#' @param metrics Character vector specifying which terrain attributes to return. The default is to return all available metrics, c("qslope", "qaspect", "qeastness", "qnorthness", "profc", "planc", "meanc", "maxc", "minc", "longc", "crosc", "features"). Slope, aspect, eastness, and northness are preceded with a 'q' to differentiate them from the measures calculated by SlpAsp() where the 'q' indicates that a quadratic surface was used for the calculation. 'profc' is the profile convexity, 'planc' is the plan convexity, 'meanc' is the mean curvature, 'minc' is minimum curvature, 'longc' is longitudinal curvature, crosc is cross-sectional curvature, and 'features' are morphometric features. See details.
#' @param slope_tolerance Slope tolerance that defines a 'flat' surface (degrees; default = 1.0). Relevant for the features layer.
#' @param curvature_tolerance Curvature tolerance that defines 'planar' surface (default = 0.0001). Relevant for the features layer.
#' @param na.rm Logical vector indicating whether or not to remove NA values before calculations.
#' @param pad Logical value specifying whether rows/columns of NA's should be padded to the edge of the raster to remove edge effects (default = FALSE). If pad is TRUE, na.rm must be TRUE.
#' @param include_scale Logical indicating whether to append window size to the layer names (default = FALSE).
#' @param mask_aspect Logical. If TRUE (default), aspect will be set to NA and northness and eastness will be set to 0 when slope = 0. If FALSE, aspect is set to 270 degrees or 3*pi/2 radians (atan2(0,0)-pi/2+2*pi) and northness and eastness will be calculated from this.
#' @param return_params Logical indicating whether to return Wood/Evans regression parameters (default = FALSE).
#' @return a RasterStack
#' @details This function calculates slope, aspect, eastness, northness, profile curvature, planform curvature, mean curvature, maximum curvature, minimum curvature, longitudinal curvature, cross-sectional curvature, and morphometric features using a quadratic surface fit from Z = aX^2+bY^2+cXY+dX+eY+f, where Z is the elevation or depth values, X and Y are the xy coordinates relative to the central cell in the focal window, and a-f are parameters to be estimated (Evans, 1980; Wood, 1996). This is an R/C++ implementation of r.param.scale GRASS GIS function. Note, for aspect, 0 degrees represents north and increases clockwise which differs from the way r.param.scale reports aspect. Additionally, mean curvature is included, which is not available in r.param.scale. All formulas with the exception of mean curvature are from Wood 1996. Mean curvature is calculated according to Wilson et al 2007. All multiplicative constants were removed from curvature formulas so that they are all reported in units of 1/length (Minár et al, 2020). Naming convention for curvatures is not consistent across the literature, however Minár et al (2020) has suggested a framework in which the reported measures of curvature translate to profile convexity = (kn)s), plan convexity = (kp)c), mean curvature = z''mean, maximum curvature = z''min, minimum curvature = z''max, longitudinal curvature = zss, and cross-sectional curvature = zcc. For curvatures, we have adopted a geographic sign convention where convex is positive and concave is negative (i.e., hills are considered convex with positive curvature values; Minár et al, 2020; Wood, 1996).
#' @import raster
#' @references
#' Evans, I.S., 1980. An integrated system of terrain analysis and slope mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295.
#' 
#' Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
#' 
#' Wilson, M.F., O’Connell, B., Brown, C., Guinan, J.C., Grehan, A.J., 2007. Multiscale Terrain Analysis of Multibeam Bathymetry Data for Habitat Mapping on the Continental Slope. Marine Geodesy 30, 3-35. https://doi.org/10.1080/01490410701295962
#' 
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' @export


WoodEvans<- function(r, w=c(3,3), unit= "degrees", metrics= c("qslope", "qaspect", "qeastness", "qnorthness", "profc", "planc", "meanc", "maxc", "minc", "longc", "crosc", "features"), slope_tolerance=1, curvature_tolerance=0.0001, na.rm=FALSE, pad=FALSE, include_scale=FALSE, mask_aspect=TRUE, return_params= FALSE){
  
  all_metrics<- c("qslope", "qaspect", "qeastness", "qnorthness", "profc", "planc", "meanc", "maxc", "minc", "longc", "crosc", "features")
  #Input checks
  if(raster::isLonLat(r)){
    stop("Error: Coordinate system is Lat/Lon. Coordinate system must be projected with elevation/depth units matching map units.")
  }
  if(suppressWarnings(raster::couldBeLonLat(r))){
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
  if(prod(w) < 6){
    stop("Error: Window size must have at least 6 cells to fit surface")
    }
  if (!any(unit==c("degrees", "radians"))){
    stop("unit must be 'degrees' or 'radians'")
  }
  if (any(!(metrics %in% all_metrics))){
    stop("Error: Invlaid metric. Valid metrics include 'qslope', 'qaspect', 'qeastness', 'qnorthness', 'profc', 'planc', 'meanc', 'maxc', 'minc', 'longc', 'crosc', and `features`.")
    }
  
  needed_metrics<- metrics
  
  if(any(c("qeastness", "qnorthness") %in% needed_metrics) & (!("qaspect" %in% needed_metrics))){
    needed_metrics<- c(needed_metrics, "qaspect")
  }
  
  if("features" %in% needed_metrics){
    if(!("qslope" %in% needed_metrics)){needed_metrics<- c(needed_metrics, "qslope")} 
    if(!("crosc" %in% needed_metrics)){needed_metrics<- c(needed_metrics, "crosc")} 
    if(!("maxc" %in% needed_metrics)){needed_metrics<- c(needed_metrics, "maxc")}
    if(!("minc" %in% needed_metrics)){needed_metrics<- c(needed_metrics, "minc")}
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
  run_in_blocks<- !raster::canProcessInMemory(r, n = 7)
  if(run_in_blocks==FALSE){
    params<- raster::brick(r, nl=6, values=FALSE)
    values(params)<- C_multiscale2(r = as.matrix(r), w = w, X=X , na_rm=na.rm) 
    } else{
    
    f_out <- raster::rasterTmpFile()
    params<- raster::brick(r, nl=6, values=FALSE)
    params <- raster::writeStart(params, filename = f_out)
    names(params)<- c("a", "b", "c", "d", "e", "f")
    
    block_idx<- raster::blockSize(r, n = 7, minblocks = 2, minrows = w[1])
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
        out_block<- out_block[(1+block_overlap*nc):(nrow(out_block)-(block_overlap*nc)),] #Trim top and bottom edge of raster
      } else {
        out_block<- out_block[(1+block_overlap*nc):nrow(out_block),] #Trim top edge of raster
      }
      raster::writeValues(params, v= out_block, start= block_idx$row[i])
    }
    params<- raster::writeStop(params)
    }
  names(params)<- c("a", "b", "c", "d", "e", "f")
  if(pad==TRUE){
    params<- raster::crop(params, og_extent)
  }
  
  mask_raster<- (params$d== 0 &  params$e== 0) #mask indicating when d ane e are 0 (slope is 0)
  
  out<- stack() #Initialize output

  #Use regression parameters to calculate slope and aspect
  if("qslope" %in% needed_metrics){
    slp<- atan(sqrt(params$d^2 + params$e^2))
    if(unit=="degrees"){
      slp<- slp*180/pi
    } else{
      slope_tolerance<- slope_tolerance * (pi/180) #Convert slope tolerance to radians if unit is not degrees
      }
    names(slp)<- "qslope"
    out<- stack(out, slp)
  }
  
  if("qaspect" %in% needed_metrics){
    asp<- atan2(params$e,params$d) - pi/2 #Shift aspect so north is zero
    asp[asp < 0]<- asp[asp < 0] + 2*pi
    asp[asp >= 2*pi]<- asp[asp >= 2*pi] - 2*pi # Constrain aspect from 0 to 2pi
    
    if (mask_aspect){
      asp[mask_raster]<- NA_real_
    }
    
    if("qeastness" %in% needed_metrics){
      eastness<- sin(asp)
      if (mask_aspect){
        eastness[mask_raster]<- 0
        }
      names(eastness)<- "qeastness"
      out<- stack(out, eastness)
    }
    
    if("qnorthness" %in% needed_metrics){
      northness<- cos(asp)
      if (mask_aspect){
        northness[mask_raster]<- 0
      }
      names(northness)<- "qnorthness"
      out<- stack(out, northness)
    }
    if(unit=="degrees"){
      asp<- asp*180/pi
      }
    names(asp)<- "qaspect"
    out<- stack(out, asp)
  }
  
  #Curvature
  #Note curvature has been multiplied by 100 to express curvature as percent gradient per unit length (Albani et al 2004)
  
  #  Wood 1996 page 86
  if("profc" %in% needed_metrics){
    profc<- raster::overlay(params, fun = function(a,b,c,d,e,f)(-2 * (a*d^2 + b*e^2 + c*d*e)) / ((e^2 + d^2)*(1 + e^2 + d^2)^1.5))
    profc[mask_raster]<- 0
    names(profc)<- "profc"
    out<- stack(out, profc)
  }
  
  if("planc" %in% needed_metrics){
    planc<- raster::overlay(params, fun = function(a,b,c,d,e,f)(2 * (b*d^2 +a*e^2 - c*d*e)) / ((e^2+d^2)^1.5))
    planc[mask_raster]<- 0
    names(planc)<- "planc"
    out<- stack(out, planc)
  }
  
  #Wood 1996 page 88
  if("longc" %in% needed_metrics){
    longc<- raster::overlay(params, fun = function(a,b,c,d,e,f) -2 * ((a*d^2 + b*e^2 + c*d*e) / (d^2 + e^2)))
    longc[mask_raster]<- 0
    names(longc)<- "longc"
    out<- stack(out, longc)
    }
  
  if("crosc" %in% needed_metrics){
    crosc<- raster::overlay(params, fun = function(a,b,c,d,e,f) -2 * ((b*d^2 + a*e^2 - c*d*e) / (d^2 + e^2)))
    crosc[mask_raster]<- 0
    names(crosc)<- "crosc"
    out<- stack(out, crosc)
    }
  
  #Wood 1996 page 115
  if("maxc" %in% needed_metrics){
    max_curv<- raster::overlay(params, fun = function(a,b,c,d,e,f) (-a - b + sqrt((a-b)^2+c^2)))
    names(max_curv)<- "maxc"
    out<- stack(out, max_curv)
    }
  if("minc" %in% needed_metrics){
    min_curv<- raster::overlay(params, fun = function(a,b,c,d,e,f) (-a - b - sqrt((a-b)^2+c^2)))
    names(min_curv)<- "minc"
    out<- stack(out, min_curv)
    }
  
  #Wilson 2007 page 10
  if("meanc" %in% needed_metrics){
    mean_curv<- raster::overlay(params, fun = function(a,b,c,d,e,f) (-a - b))
    names(mean_curv)<- "meanc"
    out<- stack(out, mean_curv)
  }
  
  #Morphometric Features (Wood 1996, Page 120)
  if("features" %in% needed_metrics){
    classify_features<- classify_features_ff(slope_tolerance, curvature_tolerance) #Define classification function based on slope and curvature tolerance
    features<- raster::overlay(slp, crosc, max_curv, min_curv, fun = classify_features)
    levels(features)<- data.frame(ID=1:6, Feature = c("Planar", "Pit", "Channel", "Pass", "Ridge", "Peak"))
    names(features)<- "features"
    out<- stack(out, features)
  }
  
  out<- raster::subset(out, metrics, drop=TRUE) #Subset needed metrics to requested metrics in proper order
  if(return_params){out<- stack(out, params)}
  if(include_scale){names(out)<- paste0(names(out), "_", w[1],"x", w[2])} #Add scale to layer names
                             
  #identify extreme outliers that are less than Q1% - 100*IQR or greater than Q99% + 100*IQR, where IQR is the range of 1-99% quantiles
  quant <- raster::quantile(
    out[[c(
      which(!raster::is.factor(out))
    )]], 
    probs=c(0, 0.01, 0.99, 1)
  )
                                
  iqr <- quant[ ,3] - quant[ ,2]
  
  outliers <- names(
    which(
      quant[ ,1] < (quant[ ,2] - 100*iqr)  | quant[ ,4] > (quant[ ,3] + 100*iqr)
    )
  )
  
  if(length(outliers != 0)){
    warning(
      "Extreme outliers detected in: ", paste(outliers, collapse=", ")
    )
  }
                             
  return(out)
}
