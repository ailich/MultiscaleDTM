#' Creates circular focal window
#'
#' Creates circular focal window around central pixel.
#' @param radius radius of circular window
#' @param resolution resolution of intended raster layer (one number or a vector of length 2). Only necessary if unit= "map"
#' @param unit unit for radius. Either "cell" (number of cells, the default) or "map" for map units (e.g. meters).
#' @param return_dismat logical, if TRUE return a matrix of distances from focal cell instead of a matrix to pass to terra::focal.
#' @return a matrix of 1's and NA's showing which cells to include and exclude respectively in focal calculations, or if return_dismat=TRUE, a matrix indicating the distance from the focal cell.
#' @export
circle_window<- function(radius, unit= "cell", resolution, return_dismat = FALSE){
  unit<- tolower(unit)
  if (unit != "map" & unit != "cell"){
    stop("Error: unit must equal 'map' or 'cell'")
  }
  if(length(radius)!=1){
    stop("Error: radius must be a single integer")
  }
  if(unit=="cell"){
    resolution<- c(1,1)
  }
  if (length(resolution)==1){
    resolution<- rep(resolution,2)
  }
  nrows<- floor((radius/resolution[2])*2+1)
  ncols<- floor((radius/resolution[1])*2+1)
  if(nrows %% 2 ==0){
    nrows<- nrows+1
  }
  if(ncols %% 2 ==0){
    ncols<- ncols+1
  } #nrow and ncol must be odd to have central pixel
  x<- matrix(seq(1:ncols), nrow = nrows, ncol =ncols, byrow = TRUE)  - ((ncols+1)/2)
  y<- matrix(seq(1:nrows), nrow=nrows, ncol=ncols, byrow = FALSE) - ((nrows+1)/2)
  x<- x * resolution[1]
  y<- y * resolution[2]
  dis_mat<- sqrt((y^2)+(x^2)) #Distance from center of window
  if(return_dismat){
    return(dis_mat)
  }
  w<- matrix(NA, nrow = nrow(dis_mat), ncol= ncol(dis_mat))
  w[dis_mat <= radius]<- 1
  if(unit=="map"){
    attributes(w)<- append(attributes(w), list(scale=paste0(radius,"MU")))
  } else{
    attributes(w)<- append(attributes(w), list(scale=as.character(radius)))
  }
  attributes(w)<- append(attributes(w), list(shape="circle"))
  return(w)
}

#' Creates annulus focal window
#'
#' Creates annulus focal window around central pixel.
#' @param radius radius of inner annulus c(inner,outer). Inner radius must be less than or equal to outer radius.
#' @param unit unit for radius. Either "cell" (number of cells, the default) or "map" for map units (e.g. meters).
#' @param resolution resolution of intended raster layer (one number or a vector of length 2). Only necessary if unit= "map"
#' @param return_dismat logical, if TRUE return a matrix of distances from focal cell instead of a matrix to pass to terra::focal (default FALSE)
#' @return a matrix of 1's and NA's showing which cells to include and exclude respectively in focal calculations, or if return_dismat=TRUE, a matrix indicating the distance from the focal cell.
#' @export
annulus_window<- function(radius, unit= "cell", resolution, return_dismat=FALSE){
  unit<- tolower(unit)
  if(length(radius) != 2){
    stop("radius must be a vector of length 2")
  }
  if(radius[1] > radius[2]){
    stop("Error: inner radius must be less than or equal to outer radius")
  }
  if((radius[1] < 1) & (unit == "cell")){
    stop("Error: inner radius must be at least 1")
  }
  if(unit=="cell"){
    resolution<- c(1,1)
  }
  if((radius[1] < max(resolution)) & (unit == "map")){
    stop("Error: inner radius must be >= resolution")
  }
  if(length(resolution) > 2){
    stop("Specified inner radius exceeds 2 dimensions")
  }
  
  dis_mat<- circle_window(radius = radius[2], unit = unit, resolution = resolution, return_dismat = TRUE)
  if(return_dismat){
    return(dis_mat)
  }
  w<- matrix(NA, nrow=nrow(dis_mat), ncol=ncol(dis_mat))
  w[(dis_mat >= radius[1]) & (dis_mat <= radius[2])]<- 1
  if(unit=="map"){
    attributes(w)<- append(attributes(w), list(scale=paste0(radius[1],"MUx", radius[2], "MU")))
  } else{
    attributes(w)<- append(attributes(w), list(scale=paste0(radius[1],"x", radius[2])))
  }
  attributes(w)<- append(attributes(w), list(shape="annulus"))
  return(w)
}

#' Calculates Relative Position of a focal cell
#'
#' Calculates the relative position of a focal cell, which represents whether an area is a local high or low. Relative position is the value of the focal cell minus the value of the mean of included values in the focal window. Positive values indicate local highs (i.e. peaks) and negative values indicate local lows (i.e. depressions). Relative Position can be expressed in units of the input DTM raster or can standardized relative to the local topography by dividing by the standard deviation or range of included elevation values in the focal window.
#' @param r DTM as a SpatRaster or RasterLayer.
#' @param w For a "rectangle" focal window, a vector of length 2 specifying dimensions where the first number is the number of rows and the second is the number of columns (or a single number if the number of rows and columns is equal). Window size must be an odd number, and the default is 3x3. For circle and annulus shaped windows, w can be set to NA or NULL and radius can be used instead, or w can be specified using focal weights matrix created by MultiscaleDTM::circle_window or MultiscaleDTM::annulus_window. If "custom" focal window shape is used, w should be a focal weights matrix with 1's for included values and NAs for excluded values.
#' @param shape Character representing the shape of the focal window. Either "rectangle" (default), "circle", or "annulus", or "custom". If a "custom" shape is used, w must be a focal weights matrix.
#' @param radius For "circle" shaped focal windows, a single integer representing the radius. For "annulus" shaped focal windows, a vector of length 2 specifying c(inner, outer) radii of the annulus in "cell" or "map" units. Inner radius must be less than or equal to outer radius. For a circle, the default radius is 1 cell if unit= "cell" or the maximum of the x and y cell resolution if unit="map". There is no default for an annulus window.
#' @param stand Standardization method. Either "none" (the default), "range" or "sd" indicating whether the relative position should be standardized by dividing by the standard deviation or range of included values in the focal window. If stand is 'none' the layer name will be "rpos", otherwise it will be "srpos" to indicate that the layer has been standardized.
#' @param exclude_center Logical indicating whether to exclude the central value from focal calculations (Default=FALSE). Use FALSE for DMV and TRUE for TPI. Note, if a focal weights matrix is supplied to w, setting exclude_center=TRUE will overwrite the center value of w to NA, but setting exclude_center=FALSE will not overwrite the central value to be 1. 
#' @param unit Unit for radius. If a circle or annulus shaped window is specified with the radius parameter either "cell" (number of cells) or "map" for map units (e.g. meters). Otherwise it should be NA or NULL. If radius is specified the default is "cell", otherwise if w is specified, the default is NA.
#' @param na.rm Logical indicating whether or not to remove NA values before calculations.
#' @param include_scale Logical indicating whether to append window size to the layer names (default = FALSE) or a character vector specifying the name you would like to append or a number specifying the number of significant digits. If include_scale = TRUE the number of rows and number of columns will be appended for rectangular or custom windows. For circular windows it will be a single number representing the radius. For annulus windows it will be the inner and outer radius. If unit="map" then window size will have "MU" after the number indicating that the number represents the scale in map units (note units can be extracted from w created with MultiscaleDTM::circle_window and MultiscaleDTM::annulus_window).
#' @param filename Character output filename.
#' @param overwrite Logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt List with named options for writing files as in writeRaster.
#' @return A SpatRaster or RasterLayer.
#' @examples
#' r<- rast(volcano, extent= ext(2667400, 2667400 + 
#' ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
#' crs = "EPSG:27200")
#' rpos<- RelPos(r, w = c(5,5), shape= "rectangle", exclude_center = TRUE, na.rm = TRUE)
#' plot(rpos)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @importFrom dplyr case_when
#' @references 
#' Lecours, V., Devillers, R., Simms, A.E., Lucieer, V.L., Brown, C.J., 2017. Towards a Framework for Terrain Attribute Selection in Environmental Studies. Environmental Modelling & Software 89, 19-30. https://doi.org/10.1016/j.envsoft.2016.11.027
#'
#' Lundblad, E.R., Wright, D.J., Miller, J., Larkin, E.M., Rinehart, R., Naar, D.F., Donahue, B.T., Anderson, S.M., Battista, T., 2006. A benthic terrain classification scheme for American Samoa. Marine Geodesy 29, 89–111. https://doi.org/10.1080/01490410600738021
#'
#' Weiss, A., 2001. Topographic Position and Landforms Analysis. Presented at the ESRI user conference, San Diego, CA.
#'
#' Wilson, J.P., Gallant, J.C. (Eds.), 2000. Terrain Analysis: Principles and Applications. John Wiley & Sons, Inc.
#' @export
#' 
RelPos<- function(r, w=ifelse(tolower(shape)=="rectangle", c(3,3), NA_real_), shape= "rectangle", 
                  radius=dplyr::case_when(is.matrix(w) ~ NA_real_,
                                          tolower(shape)=="circle" & tolower(unit)=="cell" ~ 1,
                                          tolower(shape)=="circle" & tolower(unit)=="map" ~ max(terra::res(r)),
                                          TRUE ~ NA_real_),
                  stand="none", exclude_center= FALSE, 
                  unit=ifelse((!is.matrix(w)) & all(is.null(w) | is.na(w)), "cell", NA_character_), 
                  na.rm=FALSE, include_scale =FALSE, filename=NULL, overwrite=FALSE, wopt=list()){
#Input checks
og_class<- class(r)[1]
if(og_class=="RasterLayer"){
  r<- terra::rast(r) #Convert to SpatRaster
}
if(!(og_class %in% c("RasterLayer", "SpatRaster"))){
  stop("Error: Input must be a 'SpatRaster' or 'RasterLayer'")
}
if(terra::nlyr(r)!=1){
  stop("Error: Input raster must be one layer.")
}
stand<- tolower(stand)
shape<- tolower(shape)
resolution<- terra::res(r)

## Focal Window Checks
w[is.null(w)]<- NA_real_
radius[is.null(radius)]<- NA_real_ 
unit[is.null(unit)]<- NA_character_ #Change NULL values to NA for consistency

if (!(shape %in% c("rectangle", "circle", "annulus", "custom"))){
  stop("Error: shape must be 'rectangle', 'circle', 'annulus', or 'custom'")
}

if((!(shape %in% c("circle", "annulus"))) & any(!is.na(radius))){
  stop("Error: radius should only be specified for cicle and annulus shaped windows")
}

if(is.matrix(w) & any(!is.na(radius))){
  stop("Error: Radius should not be specified if w is a matrix")
}

if(is.matrix(w) & shape=="rectangle"){
  stop("Error: w should be a vector or single number for rectangular focal windows")
}

if(is.matrix(w) & shape %in% c("circle", "annulus")){
  if(!isTRUE(attributes(w)$shape == shape)){
    stop(paste0("Error: shape = ", shape, " but w was not created using ", "MultiscaleDTM::", shape, "_window"))
  }
}


if(length(w)==1){w<- rep(w, 2)}
if(is.vector(w) & any(!is.na(w))){
  if(length(w) > 2){  
    stop("Specified window exceeds 2 dimensions")}
  if(any(0 == (w %% 2))){
    stop("Error: w must be odd")}
  if(all(w<3)){
    stop("Error: w must be greater or equal to 3 in at least one dimension")}
} else if(is.matrix(w)){
  if(any(0 == (dim(w) %% 2))){
    stop("Error: dimensions of w must be odd")}
} else if(any(!is.na(w))){
  stop("Error: w must be a whole number, vector, matrix, or NA/NULL")
}

## unit check
if(unit=="map" & (!shape %in% c("circle", "annulus"))){
  stop("Error: unit can only be 'map' if shape is circle or annulus")
}

if(is.na(radius) & !is.na(unit)){
  stop("Error: unit be NULL or NA if no radius is specified")
} else if(!is.na(radius) & is.na(unit)){
  stop("Error: unit be 'cell' or 'map' if radius is specified")
}

## stand checks
if (length(stand)!= 1){
  stop("Error: stand must be of a character vector of length 1")
}
if (!(stand %in% c("none", "range", "sd"))){
  stop("Error: stand must be 'none', 'range', or 'sd'")
}


# Set up focal window
## if focal weights matrix if w is a vector
if(is.vector(w)){
  if(shape=="rectangle"){
    w_mat <- matrix(1, nrow=w[1], ncol=w[2])
    } else if(shape=="circle"){
      w_mat<- MultiscaleDTM::circle_window(radius = radius, unit = unit, resolution = resolution, return_dismat = FALSE)
    } else{
      w_mat<- MultiscaleDTM::annulus_window(radius = radius, unit = unit, resolution = resolution, return_dismat = FALSE)
    }
  } else{
    w_mat<- w
    if(!all((is.na(w_mat) | w==1))){
      stop("Error: Focal weights matrix, w, must only include NA's and 1's")
    }
  }

if((shape %in% c("circle", "annulus")) & is.na(unit)){
  if(grepl(pattern = "MU", attributes(w_mat)$scale)){
    unit<- "map" 
    } else{
      unit<- "cell"
    }
}

## Exclude center cell
if(exclude_center){
  center_idx<- ceiling(0.5 * length(w_mat))
  w_mat[center_idx] <- NA_real_ #Exclude central cell
}

# Calculate Relative Position
rpos<- r - terra::focal(x = r, w = w_mat, fun = mean, na.rm = na.rm, wopt=wopt)
names(rpos)<- "rpos"

# standardize
if(stand!="none"){
  if(stand=="range"){
    localmax<- terra::focal(x = r, w= w_mat, fun=max, na.rm = na.rm, wopt=wopt)
    localmin<- terra::focal(x = r, w= w_mat, fun=min, na.rm = na.rm, wopt=wopt)
    rpos<- (rpos)/(localmax-localmin)
    names(rpos)<- "srpos"
    } else if(stand=="sd"){
      localsd<- terra::focal(x = r, w= w_mat, fun="sd", na.rm = na.rm, wopt=wopt)
      rpos<- rpos/(localsd)
      names(rpos)<- "srpos"
    }
}

# Extract scale
if(is.character(include_scale)){
  spatial_scale<- include_scale
  } else if(shape=="rectangle"){
    spatial_scale<- paste0(w[1], "x", w[2])
    } else{
  spatial_scale<- attributes(w_mat)$scale
  }
if(is.null(spatial_scale) & shape=="custom"){
  spatial_scale<- paste0(nrow(w_mat), "x", ncol(w_mat))
}
if(is.numeric(include_scale)){ 
  matches<- unlist(regmatches(spatial_scale, gregexpr(pattern = "\\d+\\.?\\d*e?-?\\d*", text = spatial_scale)))
  replacement<- as.character(signif(as.numeric(matches), include_scale))
  if(unit=="map"){
    replacement<- paste0(replacement,"MU")
  }
  if(length(replacement)==1){
    spatial_scale<- replacement
    } else if(length(replacement)==2){
      spatial_scale<- paste0(replacement[1], "x", replacement[2])
      } else{
        stop("Error: Something went wrong in extracting spatial scale")
      }
}
if(isTRUE(include_scale) | is.character(include_scale) | is.numeric(include_scale)){names(rpos)<- paste0(names(rpos), "_", spatial_scale)} #Add scale to layer names

#Return
if(og_class =="RasterLayer"){
  rpos<- raster::raster(rpos)
  if(!is.null(filename)){
    return(raster::writeRaster(rpos, filename=filename, overwrite=overwrite))
  }
}
if(!is.null(filename)){
  return(terra::writeRaster(rpos, filename=filename, overwrite=overwrite, wopt=wopt))
}
return(rpos)
}
