#' Calculates surface area to planar area rugosity
#'
#' Calculates surface area (Jenness, 2004) to planar area rugosity and by default corrects planar area for slope using the arc-chord ratio (Du Preez, 2015). Additionally, the method has been modified to allow for calculations at multiple different window sizes (see details).
#' @param r DEM as a raster layer
#' @param w A single number or a vector of length 2 (row, column) specifying the dimensions of the rectangular window over which surface area will be summed. Window size must be an odd number. 1 refers to "native" scale and surface area and planar area will be calculated per cell (the traditional implementation).
#' @param slope_layer Optionally specify an appropriate slope layer IN RADIANS to use. If not supplied, it will be calculated using the SlpAsp function based on Misiuk et al (2021). The slope layer should have a window size that is 2 larger than the w specified for SAPA.
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @details Planar area is calculated as the x_dis * y_dis if uncorrected for slope and (x_dis * y_dis)/cos(slope) if corrected for slope. When w=1, this is called "native" scale and is equivalent to what is presented in Du Preez (2015) and available in the ArcGIS Benthic Terrain Modeller add-on. In this case operations are performed on a per cell basis where x_dis is the resolution of the raster in the x direction (left/right) and y_dis is the resolution of the raster in the y direction (up/down) and slope is calculated using the Horn (1981) method. To expand this to multiple scales of analysis, at w > 1 slope is calculated based on Misiuk et al (2021) which provides a modification of the Horn method to extend the matric to multiple spatial scales. Planar area is calculated the same way as for w=1 except that now x_dis is the x resolution of the raster * the number of columns in the focal window, and y_dis is y resolution of the raster * the number of rows. For w > 1, surface area is calculated as the sum of surface areas within the focal window. Although the (modified) Horn slope is used by default to be consistent with Du Preez (2015), slope calculated using a different algorithm (e.g. Wood 1996) could be supplied using the slope_layer argument. Additionally, a slope raster can be supplied if you have already calculated it and do not wish to recalculate it. However, be careful to supply a slope layer measured in radians and calculated at the relevant scale (2 larger than the w of SAPA). 
#' @return a RasterLayer
#' @import raster
#' @export

SAPA<- function(r, w = 1, slope_layer= NULL, slope_correction=TRUE, include_scale=FALSE){
  
  if(length(w)==1){w<- rep(w,2)}
  if(any(w %% 2 != 1)){
    stop("w must be odd")
  }
  # if(any(w == 1) & any(w != 1)){
  #     stop("w cannot be 1 in only one direction")
  #   }
  if(all(w==c(1,1))){
    is_native<- TRUE} else{
      is_native<- FALSE} #Indicate whether SAPA is calculated at native scale
  
  x_res<- res(r)[1]
  y_res<- res(r)[2]

  SA<- SurfaceArea(r)
  
  if(!is_native) {
    SA<- focal(SA, w= matrix(data=1, nrow=w[1], ncol=w[2]), fun=sum, na.rm=FALSE)
  } # Sum up surface area in focal window
  
  if(is.null(slope_layer)){
    if (all(w==1)){
      slope_layer<- terrain(r, opt="slope", unit="radians", neighbors=8)
      } else{
        slope_layer<- SlpAsp(r, w=w+2, unit="radians", method="queen", metrics= "slope", include_scale=FALSE)
      }
    }
  PA<- (x_res*w[1]) * (y_res*w[2])
  if(slope_correction){PA<- PA/cos(slope_layer)}#Planar area corrected for slope
  
  sapa<- SA/PA
  names(sapa)<- "sapa"
  
  if(include_scale){
    if(is_native){
      names(sapa)<- paste0(names(sapa), "_native")
    } else{
      names(sapa)<- paste0(names(sapa), "_", w[1], "x", w[2])
      }}
  return(sapa)
}
