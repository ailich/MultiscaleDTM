#' Calculates surface area to planar area rugosity
#'
#' Calculates surface area (Jenness, 2004) to planar area rugosity and by default corrects planar area for slope using the arc-chord ratio (Du Preez, 2015). Additionally, the method has been modified to allow for calculations at multiple different window sizes (see details and Ilich et al. (2023)).
#' @param r DTM as a SpatRaster or RasterLayer in a projected coordinate system where map units match elevation/depth units
#' @param w A single number or a vector of length 2 (row, column) specifying the dimensions of the rectangular window over which surface area will be summed. Window size must be an odd number. 1 refers to "native" scale and surface area and planar area will be calculated per cell (the traditional implementation).
#' @param slope_correction Whether to use the arc-chord ratio to correct planar area for slope (default is TRUE)
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @param slope_layer Optionally specify an appropriate slope layer IN RADIANS to use. If not supplied, it will be calculated using the SlpAsp function based on Misiuk et al (2021). The slope layer should have a window size that is 2 larger than the w specified for SAPA.
#' @param filename character Output filename.
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @examples
#' r<- rast(volcano, extent= ext(2667400, 2667400 + 
#' ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
#' crs = "EPSG:27200")
#' sapa<- SAPA(r, w=c(5,5), slope_correction = TRUE)
#' plot(sapa)
#' @details Planar area is calculated as the x_dis * y_dis if uncorrected for slope and (x_dis * y_dis)/cos(slope) if corrected for slope. When w=1, this is called "native" scale and is equivalent to what is presented in Du Preez (2015) and available in the ArcGIS Benthic Terrain Modeller add-on. In this case operations are performed on a per cell basis where x_dis is the resolution of the raster in the x direction (left/right) and y_dis is the resolution of the raster in the y direction (up/down) and slope is calculated using the Horn (1981) method. To expand this to multiple scales of analysis, at w > 1 slope is calculated based on Misiuk et al (2021) which provides a modification of the Horn method to extend the matric to multiple spatial scales. Planar area is calculated the same way as for w=1 except that now x_dis is the x resolution of the raster * the number of columns in the focal window, and y_dis is y resolution of the raster * the number of rows. For w > 1, surface area is calculated as the sum of surface areas within the focal window. Although the (modified) Horn slope is used by default to be consistent with Du Preez (2015), slope calculated using a different algorithm (e.g. Wood 1996) could be supplied using the slope_layer argument. Additionally, a slope raster can be supplied if you have already calculated it and do not wish to recalculate it. However, be careful to supply a slope layer measured in radians and calculated at the relevant scale (2 larger than the w of SAPA). 
#' @return a SpatRaster or RasterLayer
#' @import terra
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @references
#' Du Preez, C., 2015. A new arc–chord ratio (ACR) rugosity index for quantifying three-dimensional landscape structural complexity. Landscape Ecol 30, 181–192. https://doi.org/10.1007/s10980-014-0118-8
#' 
#' Horn, B.K., 1981. Hill Shading and the Reflectance Map. Proceedings of the IEEE 69, 14-47.
#' 
#' Ilich, A. R., Misiuk, B., Lecours, V., & Murawski, S. A. (2023). MultiscaleDTM: An open-source R package for multiscale geomorphometric analysis. Transactions in GIS, 27(4). https://doi.org/10.1111/tgis.13067
#' 
#' Jenness, J.S., 2004. Calculating landscape surface area from digital elevation models. Wildlife Society Bulletin 32, 829-839.
#' 
#' Misiuk, B., Lecours, V., Dolan, M.F.J., Robert, K., 2021. Evaluating the Suitability of Multi-Scale Terrain Attribute Calculation Approaches for Seabed Mapping Applications. Marine Geodesy 44, 327-385. https://doi.org/10.1080/01490419.2021.1925789
#' @export

SAPA<- function(r, w = 1, slope_correction=TRUE, include_scale=FALSE, slope_layer= NULL, filename=NULL, overwrite=FALSE, wopt=list()){
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
  if(isTRUE(terra::is.lonlat(r, perhaps=FALSE))){
    stop("Error: Coordinate system is Lat/Lon. Coordinate system must be projected with elevation/depth units matching map units.")
  }
  if(terra::is.lonlat(r, perhaps=TRUE, warn=FALSE)){
    warning("Coordinate system may be Lat/Lon. Please ensure that the coordinate system is projected with elevation/depth units matching map units.")
  }
  if(length(w)==1){w<- rep(w,2)}
  if(length(w) > 2){
    stop("Specified window exceeds 2 dimensions")}
  if(any(0 == (w %% 2))){
    stop("w must be odd")
  }
  if(all(w==c(1,1))){
    is_native<- TRUE} else{
      is_native<- FALSE} #Indicate whether SAPA is calculated at native scale
  
  x_res<- terra::xres(r)
  y_res<- terra::yres(r)

  SA<- SurfaceArea(r, wopt=wopt)
  
  if(!is_native) {
    SA<- terra::focal(SA, w= w, fun=sum, na.rm=FALSE, wopt=wopt)
  } # Sum up surface area in focal window
  
  if(is.null(slope_layer) & slope_correction){
    if (all(w==1)){
      slope_layer<- terra::terrain(r, v="slope", unit="radians", neighbors=8, wopt=wopt)
      } else{
        slope_layer<- SlpAsp(r, w=w+2, unit="radians", method="queen", metrics= "slope", na.rm=FALSE, include_scale=FALSE, wopt=wopt)
      }
    }
  PA<- (x_res*w[2]) * (y_res*w[1])
  if(slope_correction){PA<- PA/terra::math(slope_layer, fun="cos", wopt=wopt)}#Planar area corrected for slope
  
  sapa<- SA/PA
  names(sapa)<- "sapa"
  
  if(include_scale){
    if(is_native){
      names(sapa)<- paste0(names(sapa), "_native")
    } else{
      names(sapa)<- paste0(names(sapa), "_", w[1], "x", w[2])
    }}
  
  #Return
  if(og_class =="RasterLayer"){
    sapa<- raster::raster(sapa)
    if(!is.null(filename)){
      return(raster::writeRaster(sapa, filename=filename, overwrite=overwrite))
    }
  }
  if(!is.null(filename)){
    return(terra::writeRaster(sapa, filename=filename, overwrite=overwrite, wopt=wopt))
  }
  return(sapa)
}
