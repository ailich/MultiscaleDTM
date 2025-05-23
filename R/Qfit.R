#' Helper function factory to classify morphometric features
#'
#' Helper function factory to classify morphometric features according to a modified version of Wood 1996 page 120
#' @param slope_tolerance Slope tolerance that defines a 'flat' surface (degrees; default is 1.0). Relevant for the features layer.
#' @param curvature_tolerance Curvature tolerance that defines 'planar' surface (default is 0.0001). Relevant for the features layer.
#' @importFrom dplyr case_when
#' @return A function that can be passed to raster::overlay to classify morphometric features
#' @references
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.

classify_features_ff<- function(slope_tolerance=1, curvature_tolerance=0.0001){
  #1=PLANAR_FLAT
  #2=PLANAR_SLOPE
  #3=PIT
  #4=CHANNEL
  #5=PASS
  #6=RIDGE
  #7=PEAK
  out_fun<- function(slope, planc, maxc, minc){
    dplyr::case_when(is.na(slope) ~ NA_real_,
                     (slope > slope_tolerance) & (planc > curvature_tolerance) ~ 6,
                     (slope > slope_tolerance) & (planc < -curvature_tolerance) ~ 4,
                     slope > slope_tolerance ~ 2, 
                     (maxc > curvature_tolerance) & (minc > curvature_tolerance) ~ 7,
                     (maxc > curvature_tolerance) & (minc < -curvature_tolerance) ~ 5,
                     maxc > curvature_tolerance ~ 6,
                     (minc < -curvature_tolerance) & (maxc < -curvature_tolerance) ~ 3,
                     minc < -curvature_tolerance ~ 4,
                     TRUE ~ 1)}
  return(out_fun)
}

#' Helper function to filter outliers from regression parameters using interquartile range
#'
#' @param params regression parameters for fitted surface
#' @param outlier_quantile A numeric vector of length two or three. If two numbers are used it specifies the lower (Q1) and upper (Q2) quantiles used for determining the interquantile range (IQR). These values should be between 0 and 1 with Q2 > Q1. An optional third number can be used to specify a the size of a regular sample to be taken which can be useful if the full dataset is too large to fit in memory. Values are considered outliers if they are less than Q1-(100*IQR) or greater than Q2+(100*IQR), where IQR=Q2-Q1.
#' @param wopt list with named options for writing files as in writeRaster

outlier_filter<- function(params, outlier_quantile, wopt=list()){ 
  if(length(outlier_quantile)==2){outlier_quantile[3]<- Inf}
  quant <- terra::global(params, fun = quantile, probs = c(0, outlier_quantile[1], outlier_quantile[2], 1), na.rm = TRUE, maxcell=outlier_quantile[3])
  iqr <- quant[, 3] - quant[, 2]
  outliers <- row.names(quant)[which(quant[, 1] < (quant[, 2] - 100 * iqr) | quant[, 4] > (quant[, 3] + 100 * iqr))]
  
  if (length(outliers) != 0) {
    iq_lims <- matrix(c(quant[, 2] - 100 * iqr, quant[, 3] + 100 * iqr), ncol = 2)
    
    outlier_mask<- rast(params)
    for (i in 1:nlyr(params)) {
      outlier_mask[[i]]<- ((params[[i]] >= iq_lims[i,1]) & (params[[i]] <= iq_lims[i,2])) #0 indicates an outlier
    }
    outlier_mask<- prod(outlier_mask, wopt=wopt) #product of zero indicates an outlier location
    params<- terra::mask(params, mask = outlier_mask, maskvalues=0, updatevalue=NA, wopt=wopt)
    warning("Outliers filtered")
  }
  return(params)
}

#' Calculates multiscale slope, aspect, curvature, and morphometric features using a local quadratic fit
#'
#' Calculates multiscale slope, aspect, curvature, and morphometric features of a DTM over a sliding rectangular window using a local quadratic fit to the surface (Evans, 1980; Wood, 1996).
#' @param r DTM as a SpatRaster (terra) or RasterLayer (raster) in a projected coordinate system where map units match elevation/depth units (up is assumed to be north for calculations of aspect, northness, and eastness).
#' @param w Vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number. Default is 3x3.
#' @param unit "degrees" or "radians".
#' @param metrics Character vector specifying which terrain attributes to return. The default is to return all available metrics, c("elev", "qslope", "qaspect", "qeastness", "qnorthness", "profc", "planc", "twistc", "meanc", "maxc", "minc", "features"). Slope, aspect, eastness, and northness are preceded with a 'q' to differentiate them from the measures calculated by SlpAsp() where the 'q' indicates that a quadratic surface was used for the calculation. 'elev' is the predicted elevation at the central cell (i.e. the intercept term of the regression) and is only relevant when force_center=FALSE. 'profc' is the profile curvature, 'planc' is the plan curvature, 'meanc' is the mean curvature, 'minc' is minimum curvature, and 'features' are morphometric features. See details.
#' @param slope_tolerance Slope tolerance that defines a 'flat' surface (degrees; default = 1.0). Relevant for the features layer.
#' @param curvature_tolerance Curvature tolerance that defines 'planar' surface (default = 0.0001). Relevant for the features layer.
#' @param outlier_quantile A numeric vector of length two or three. If two numbers are used it specifies the lower (Q1) and upper (Q2) quantiles used for determining the interquantile range (IQR). These values should be between 0 and 1 with Q2 > Q1. An optional third number can be used to specify a the size of a regular sample to be taken which can be useful if the full dataset is too large to fit in memory. Values are considered outliers and replaced with NA if they are less than Q1-(100*IQR) or greater than Q2+(100*IQR), where IQR=Q2-Q1. The outlier filter is performed on the results of the regression parameters ('a'-'e' and 'elev') prior to calculation of subsequent terrain attributes. Note that c(0,1) will skip the outlier filtering step and can speed up computations. The default is c(0.01,0.99).
#' @param na.rm Logical indicating whether or not to remove NA values before calculations.
#' @param force_center Logical specifying whether the constrain the model through the central cell of the focal window
#' @param include_scale Logical indicating whether to append window size to the layer names (default = FALSE).
#' @param mask_aspect Logical. If TRUE (default), aspect will be set to NA and northness and eastness will be set to 0 when slope = 0. If FALSE, aspect is set to 270 degrees or 3*pi/2 radians ((-pi/2)- atan2(0,0)+2*pi) and northness and eastness will be calculated from this.
#' @param return_params Logical indicating whether to return Wood/Evans regression parameters (default = FALSE).
#' @param as_derivs Logical indicating whether parameters should be formatted as partial derivatives instead of regression coefficients (default = FALSE) (Minár et al., 2020).
#' @param filename character Output filename. Can be a single filename, or as many filenames as there are layers to write a file for each layer
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @return a SpatRaster (terra) or RasterStack/RasterLayer (raster)
#' @examples
#' r<- erupt()
#' qmetrics<- Qfit(r, w = c(5,5), unit = "degrees", na.rm = TRUE)
#' plot(qmetrics)
#' 
#' # To get only the regression coefficients, set "metrics=c()" and "return_params=TRUE"
#' reg_coefs<- Qfit(r, w = c(5,5), metrics=c(), unit = "degrees", na.rm = TRUE, return_params=TRUE)
#' plot(reg_coefs)
#' @details This function calculates slope, aspect, eastness, northness, profile curvature, plan curvature, mean curvature, twisting curvature, maximum curvature, minimum curvature, morphometric features, and a smoothed version of the elevation surface using a quadratic surface fit from Z = aX^2+bY^2+cXY+dX+eY+f, where Z is the elevation or depth values, X and Y are the xy coordinates relative to the central cell in the focal window, and a-f are parameters to be estimated (Evans, 1980; Minár et al. 2020; Wood, 1996). For aspect, 0 degrees represents north (or if rotated, the direction that increases as you go up rows in your data) and increases clockwise. For calculations of northness (cos(asp)) and eastness (sin(asp)), up in the y direction is assumed to be north, and if this is not true for your data (e.g. you are using a rotated coordinate system), you must adjust accordingly. All curvature formulas are adapted from Minár et al 2020. Therefore  all curvatures are measured in units of 1/length (e.g. m^-1) except twisting curvature which is measured in radians/length (i.e. change in angle per unit distance), and we adopt a geographic sign convention where convex is positive and concave is negative (i.e., hills are considered convex with positive. Naming convention for curvatures is not consistent across the literature, however Minár et al (2020) has suggested a framework in which the reported measures of curvature translate to profile curvature = (kn)s, plan curvature = (kn)c, twisting curvature (Tg)c, mean curvature = kmean, maximum curvature = kmax, minimum curvature = kmin. For morphometric features cross-sectional curvature (zcc) was replaced by planc (kn)c, z''min was replaced by kmax, and z''max was replaced by kmin as these are more robust ways to measures the same types of curvature (Minár et al., 2020). Additionally, the planar feature from Wood (1996) was split into planar flat and slope depending on whether the slope threshold is exceeded or not.
#' @import terra
#' @importFrom raster raster
#' @importFrom raster stack
#' @importFrom raster writeRaster
#' @references
#' Evans, I.S., 1980. An integrated system of terrain analysis and slope mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295.
#' 
#' Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
#' 
#' Wilson, M.F., O’Connell, B., Brown, C., Guinan, J.C., Grehan, A.J., 2007. Multiscale Terrain Analysis of Multibeam Bathymetry Data for Habitat Mapping on the Continental Slope. Marine Geodesy 30, 3-35. https://doi.org/10.1080/01490410701295962
#' 
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' @export

Qfit<- function(r, w=c(3,3), unit= "degrees", metrics= c("elev", "qslope", "qaspect", "qeastness", "qnorthness", "profc", "planc", "twistc", "meanc", "maxc", "minc", "features"), slope_tolerance=1, curvature_tolerance=0.0001, outlier_quantile=c(0.01, 0.99), na.rm=FALSE, force_center=FALSE, include_scale=FALSE, mask_aspect=TRUE, return_params= FALSE, as_derivs= FALSE, filename=NULL, overwrite=FALSE, wopt=list()){
  
  all_metrics<- c("elev", "qslope", "qaspect", "qeastness", "qnorthness", "profc", "planc", "twistc", "meanc", "maxc", "minc", "features")
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
  if((prod(w) < 6) & (!force_center)){
    stop("Error: Window size must have at least 6 cells to fit surface")
  }
  if((prod(w) < 5) & (force_center)){
    stop("Error: Window size must have at least 5 cells to fit surface")
  }
  unit<- tolower(unit) #make lowercase
  if (!any(unit==c("degrees", "radians"))){
    stop("unit must be 'degrees' or 'radians'")
  }
  
  # Increase tolerance of inputs to metrics
  metrics<- tolower(metrics) #Make all lowercase
  metrics[metrics=="slope"]<- "qslope" #replace slope with qslope
  metrics[metrics=="aspect"]<- "qaspect" #replace aspect with qaspect
  metrics[metrics=="eastness"]<- "qeastness" #replace eastness with qeastness
  metrics[metrics=="northness"]<- "qnorthness" #replace northness with qnorthness
  
  if (any(!(metrics %in% all_metrics))){
    stop("Error: Invlaid metric. Valid metrics include 'elev', 'qslope', 'qaspect', 'qeastness', 'qnorthness', 'profc', 'planc', 'twistc', 'meanc', 'maxc', 'minc', and 'features'.")
    }
  
  if(force_center & ("elev" %in% metrics)){
    metrics<- metrics[metrics!="elev"]
    warning("Warning: dropping 'elev' from metrics since force_center is TRUE")
    }
  needed_metrics<- metrics
  
  if(any(c("qeastness", "qnorthness") %in% needed_metrics) & (!("qaspect" %in% needed_metrics))){
    needed_metrics<- c(needed_metrics, "qaspect")
  }
  
  if(mask_aspect & ("qaspect" %in% needed_metrics) & (!("qslope" %in% needed_metrics))){
    needed_metrics<- c(needed_metrics, "qslope")
  }
  
  if("features" %in% needed_metrics){
    if(!("qslope" %in% needed_metrics)){needed_metrics<- c(needed_metrics, "qslope")} 
    if(!("planc" %in% needed_metrics)){needed_metrics<- c(needed_metrics, "planc")} 
    if(!("maxc" %in% needed_metrics)){needed_metrics<- c(needed_metrics, "maxc")}
    if(!("minc" %in% needed_metrics)){needed_metrics<- c(needed_metrics, "minc")}
  }
  
  #Define local coordinate system of window
  x_mat<- matrix(data = seq(from = (-xres(r) * floor(w[2]/2)), to = (xres(r) * floor(w[2]/2)), length.out = w[2]), nrow = w[1], ncol=w[2], byrow=TRUE)
  x<- as.vector(t(x_mat)) #Transpose for focal
  y_mat<- matrix(data = seq(from = (yres(r) * floor(w[1]/2)), to = (-yres(r) * floor(w[1]/2)), length.out = w[1]), nrow = w[1], ncol=w[2])
  y<- as.vector(t(y_mat)) #Transpose for focal
  
  #Explanatory Variable matrix X for quadratic fit
  X<- cbind(x^2,y^2, x*y, x, y) # Z = ax^2+by^2+cxy+dx+ey (no intercept)
  if(!force_center){
    X<- cbind(X,1)
    } # Z = ax^2+by^2+cxy+dx+ey+f (with intercept)
  
  if(!na.rm){
    Xt<- t(X)
    XtX_inv<- solve(Xt %*% X)
    }
  
  # Calculate Regression Parameters
  return_intercept<- "elev" %in% needed_metrics
  
  if((!na.rm) & (!force_center)){
    params<- terra::focalCpp(r, w=w, fun = C_Qfit1_narmF, X= X, Xt= Xt, XtX_inv= XtX_inv, return_intercept = return_intercept, fillvalue=NA, wopt=wopt)
    } else if(na.rm & (!force_center)){
      params<- terra::focalCpp(r, w=w, fun = C_Qfit1_narmT, X_full= X, return_intercept = return_intercept, fillvalue=NA, wopt=wopt)
    } else if((!na.rm) & force_center){
      params<- terra::focalCpp(r, w=w, fun = C_Qfit2_narmF, X= X, Xt= Xt, XtX_inv= XtX_inv, fillvalue=NA, wopt=wopt)} else{
        params<- terra::focalCpp(r, w=w, fun = C_Qfit2_narmT, X_full= X, fillvalue=NA, wopt=wopt)
      }
  
  if(return_intercept){
    names(params)<- c("a","b","c","d","e","elev")
    elev<- params$elev
    params<- params[[-6]] #drop intercept
  } else{
    names(params)<- c("a","b","c","d","e")
  }
  
  # Filter outliers
  if(!all(outlier_quantile[c(1,2)]==c(0,1))){
    params <- outlier_filter(params, outlier_quantile, wopt=wopt)
    }
  
  mask_raster<- terra::app(terra::math(params, fun= "abs", wopt=wopt), fun= "sum", na.rm=FALSE, wopt=wopt) == 0 # Mask of when all params are 0

     
  out<- terra::rast() #Initialize output
  if("elev" %in% needed_metrics){
    out<- c(out, elev, warn=FALSE)
    }

  #Use regression parameters to calculate slope and aspect
  if("qslope" %in% needed_metrics){
    slp<- terra::math(terra::math(params$d^2 + params$e^2, fun="sqrt", wopt=wopt), fun="atan", wopt=wopt)
    if(unit=="degrees"){
      slp<- slp*(180/pi)
    } else{
      slope_tolerance<- slope_tolerance * (pi/180) #Convert slope tolerance to radians if unit is not degrees
      }
    names(slp)<- "qslope"
    out<- c(out, slp, warn=FALSE)
  }
  
  if("qaspect" %in% needed_metrics){
    asp<- (-pi/2) - terra::atan_2(params$e,params$d, wopt=wopt) # aspect relative to North
    asp<- ifel(asp < 0, yes = asp+(2*pi), no= asp, wopt=wopt) # Constrain range so between 0 and 2pi
    
    if (mask_aspect & ("qaspect" %in% metrics)){
      asp<- terra::mask(asp, mask= slp, maskvalues = 0, updatevalue = NA, wopt=wopt) #Set aspect to undefined where slope is zero (doesn't need to be needed_metrics since also mask northness and eastness as 0).
    }
    
    if("qeastness" %in% needed_metrics){
      eastness<- sin(asp)
      if (mask_aspect){
        eastness<- terra::mask(eastness, mask= slp, maskvalues = 0, updatevalue = 0, wopt=wopt) #Set eastness to 0 where slope is zero
        }
      names(eastness)<- "qeastness"
      out<- c(out, eastness, warn=FALSE)
    }
    
    if("qnorthness" %in% needed_metrics){
      northness<- cos(asp)
      if (mask_aspect){
        northness<- terra::mask(northness, mask= slp, maskvalues = 0, updatevalue = 0, wopt=wopt) #Set northness to 0 where slope is zero
        }
      names(northness)<- "qnorthness"
      out<- c(out, northness, warn=FALSE)
    }
    if(unit=="degrees"){
      asp<- asp*180/pi
      }
    names(asp)<- "qaspect"
    out<- c(out, asp, warn=FALSE)
  }
  
  #Curvature
  # Formulas adapted from Minar et al 2020.
  #zxx= 2*a
  #zyy = 2*b
  #zxy=c
  #zx=d
  #zy=e
  
  if("profc" %in% needed_metrics){
    profc<- terra::lapp(params, fun = kns, wopt=wopt)
    profc<- terra::mask(profc, mask= mask_raster, maskvalues = 1, updatevalue = 0, wopt=wopt) #Set curvature to 0 where all parameters are 0
    names(profc)<- "profc"
    out<- c(out, profc, warn=FALSE)
  }
  
  if("planc" %in% needed_metrics){
    planc<- terra::lapp(params, fun = knc, wopt=wopt)
    planc<- terra::mask(planc, mask= mask_raster, maskvalues = 1, updatevalue = 0, wopt=wopt) #Set curvature to 0 where all parameters are 0
    names(planc)<- "planc"
    out<- c(out, planc, warn=FALSE)
  }
  
  if("twistc" %in% needed_metrics){
    twistc<- terra::lapp(params, fun = tgc, wopt=wopt)
    twistc<- terra::mask(twistc, mask= mask_raster, maskvalues = 1, updatevalue = 0, wopt=wopt) #Set curvature to 0 where all parameters are 0
    names(twistc)<- "twistc"
    out<- c(out, twistc, warn=FALSE)
    }
  
  if("maxc" %in% needed_metrics){
    maxc<- terra::lapp(params, fun = kmax, wopt=wopt)
    names(maxc)<- "maxc"
    out<- c(out, maxc, warn=FALSE)
  }
  
  if("minc" %in% needed_metrics){
    minc<- terra::lapp(params, fun = kmin, wopt=wopt)
    names(minc)<- "minc"
    out<- c(out, minc, warn=FALSE)
    }
  
  if("meanc" %in% needed_metrics){
    mean_curv<- terra::lapp(params, fun = kmean, wopt=wopt)
    names(mean_curv)<- "meanc"
    out<- c(out, mean_curv, warn=FALSE)
  }
  
  #Modified version of Morphometric Features (Wood 1996, Page 120)
  if("features" %in% needed_metrics){
    classify_features<- classify_features_ff(slope_tolerance, curvature_tolerance) #Define classification function based on slope and curvature tolerance
    features<- terra::lapp(c(slp, planc, maxc, minc), fun = classify_features, wopt=wopt)
    levels(features)<- data.frame(ID=1:7, features = c("Planar_Flat","Planar_Slope", "Pit", "Channel", "Pass", "Ridge", "Peak"))
    names(features)<- "features"
    out<- c(out, features, warn=FALSE)
  }
  
  if(length(metrics) > 0){out<- terra::subset(out, metrics, wopt=wopt)} #Subset needed metrics to requested metrics in proper order
  if(as_derivs){
    params$a<- 2*params$a
    params$b<- 2*params$b
    names(params)<- c("zxx", "zyy", "zxy", "zx", "zy") #partial derivatives with respect to x^2, y^2, xy, x, and y (Minar et al, 2020)
    }
  
  if(return_params){out<- c(out, params, warn=FALSE)}
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
