#' Multiscale Slope and Aspect
#'
#' Calculates multiscale slope and aspect based on the slope.k/aspect.k algorithm from Misiuk et al (2021) which extends classical formulations of slope restricted to a 3x3 window to multiple scales. The code from Misiuk et al (2021) was modified to allow for rectangular rather than only square windows.
#' @param r DEM as a raster layer in a projected coordinate system where map units match elevation/depth units
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param unit "degrees" or "radians"
#' @param method "queen" or "rook", indicating how many neighboring cells to use to compute slope for any cell. queen uses 8 neighbors (up, down, left, right, and diagonals) and rook uses 4 (up, down, left, right).
#' @param metrics a character string or vector of character strings of which terrain atrributes to return ("slope" and/or "aspect"). Default is c("slope", "aspect", "eastness", "northness").
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @param mask_aspect A logical. When mask_aspect is TRUE (the default), if slope evaluates to 0, aspect will be set to NA and both eastness and northness will be set to 0. When mask_aspect is FALSE, when slope is 0 aspect will be pi/2 radians or 90 degrees which is the behavior of raster::terrain, and northness and eastness will be calculated from that.
#' @param mask_aspect A logical. If slope evaluates to 0, aspect will be set to NA when mask_aspect is TRUE (the default). If FALSE, when slope is 0 aspect will be pi/2 radians or 90 degrees which is the behavior of raster::terrain.
#' @return a RasterStack or RasterLayer of slope and/or aspect
#' @details When method="rook", slope and aspect are computed according to Fleming and Hoffer (1979) and Ritter (1987). When method="queen", slope and aspect are computed according to Horn (1981). These are the standard slope algorithms found in many GIS packages but are traditionally restricted to a 3 x 3 window size. Misiuk et al (2021) extended these classical formulations  to multiple window sizes. This function modifies the code from Misiuk et al (2021) to allow for rectangular rather than only square windows and also added aspect.
#' @references
#' Fleming, M.D., Hoffer, R.M., 1979. Machine processing of landsat MSS data and DMA topographic data for forest cover type mapping (No. LARS Technical Report 062879). Laboratory for Applications of Remote Sensing, Purdue University, West Lafayette, Indiana.
#' 
#' Horn, B.K., 1981. Hill Shading and the Reflectance Map. Proceedings of the IEEE 69, 14-47.
#' 
#' Misiuk, B., Lecours, V., Dolan, M.F.J., Robert, K., 2021. Evaluating the Suitability of Multi-Scale Terrain Attribute Calculation Approaches for Seabed Mapping Applications. Marine Geodesy 44, 327-385. https://doi.org/10.1080/01490419.2021.1925789
#' 
#' Ritter, P., 1987. A vector-based slope and aspect generation algorithm. Photogrammetric Engineering and Remote Sensing 53, 1109-1111.
#' @import raster
#' @export

SlpAsp <- function(r, w=c(3,3), unit="degrees", method="queen", metrics= c("slope", "aspect", "eastness", "northness"), include_scale=FALSE, mask_aspect=TRUE){ 
  if(raster::isLonLat(r)){
    stop("Error: Coordinate system is Lat/Lon. Coordinate system must be projected with elevation/depth units matching map units.")
  }
  if(suppressWarnings(raster::couldBeLonLat(r))){
    warning("Coordinate system may be Lat/Lon. Please ensure that the coordinate system is projected with elevation/depth units matching map units.")
  }
  if(length(w)==1){w<- rep(w,2)}
  if(length(w) > 2){
    stop("Specified window exceeds 2 dimensions")}
  if(any(0 == (w %% 2))){
    stop("w must be odd")
  }
  if(any(w<3)){
    stop("Error: w must be greater or equal to 3")
  }
  if(!(unit %in% c("degrees", "radians"))){
    stop("unit must be `degrees` or `radians`")
  }
  if(!(method %in% c("queen", "rook"))){
    stop("method must be `queen` or `rook`")
  }
  
  if(any(!(metrics %in% c("slope", "aspect","northness", "eastness")))){
    stop("metrics must be 'slope', 'aspect', 'northness', and/or 'eastness'")
  }
  needed_metrics<- metrics
  if(any(c("eastness", "northness") %in% needed_metrics) & !("aspect" %in% needed_metrics)){
    needed_metrics<- c(needed_metrics, "aspect")
  }
  if(mask_aspect & ("aspect" %in% needed_metrics)){
    needed_metrics<- c(needed_metrics, "slope")
  }
  
  #k is size of window in a given direction
  kx<- w[2]
  ky<- w[1]
  
  #j is the number of cells on either side of the focal cell; l is used to generate the focal matrix
  jx <- (kx/2)-0.5
  jy <- (ky/2)-0.5
  
  lx <- jx-1
  ly <- jy-1
  
  if(method=="queen"){
    
    #create matrix weights for x-component
    xl.end <- matrix(c(1, rep(0, times=kx-1)), ncol=kx, nrow=1)
    xr.end <- matrix(c(rep(0, times=kx-1), 1), ncol=kx, nrow=1)
    
    x.mids <- matrix(0, ncol=kx, nrow=ly)
    
    xl.mid <- matrix(c(2, rep(0, times=kx-1)), ncol=kx, nrow=1)
    xr.mid <- matrix(c(rep(0, times=kx-1), 2), ncol=kx, nrow=1)
    
    xl.mat <- rbind(xl.end, x.mids, xl.mid, x.mids, xl.end)
    xr.mat <- rbind(xr.end, x.mids, xr.mid, x.mids, xr.end) #When upgrade to terra have NA weights instead of zero in *.mat so that you don't get NA if a cell in window that isn't used in calculation is NA
    
    #create matrix weights for y-component
    yt.end <- matrix(c(1, rep(0, times=ky-1)), ncol=1, nrow=ky)
    yb.end <- matrix(c(rep(0, times=ky-1), 1), ncol=1, nrow=ky)
    
    y.mids <- matrix(0, ncol=lx, nrow=ky)
    
    yt.mid <- matrix(c(2, rep(0, times=ky-1)), ncol=1, nrow=ky)
    yb.mid <- matrix(c(rep(0, times=ky-1), 2), ncol=1, nrow=ky)
    
    yt.mat <- cbind(yt.end, y.mids, yt.mid, y.mids, yt.end)
    yb.mat <- cbind(yb.end, y.mids, yb.mid, y.mids, yb.end)
    
    #use focal statistics for e, w, n, s components of the k-neighbourhood
    dz.dx.l <- focal(r, xl.mat, fun=sum, na.rm=FALSE)
    dz.dx.r <- focal(r, xr.mat, fun=sum, na.rm=FALSE)
    
    dz.dy.t <- focal(r, yt.mat, fun=sum, na.rm=FALSE)
    dz.dy.b <- focal(r, yb.mat, fun=sum, na.rm=FALSE)
    
    #calculate dz/dx and dz/dy using the components. 8*j is the weighted run, or distance between ends: 4*j*2, or (4 values in each row)*(length of the side)*(2 sides)
    dz.dx <- (dz.dx.r-dz.dx.l)/(8*jx*res(r)[1])
    dz.dy <- (dz.dy.b-dz.dy.t)/(8*jy*res(r)[2])
  }
  
  if(method=="rook"){
    
    #create matrix weights for x-component
    x.ends <- matrix(0, ncol=kx, nrow=jy)
    
    xl.mid <- matrix(c(1, rep(0, times=kx-1)), ncol=kx, nrow=1)
    xr.mid <- matrix(c(rep(0, times=kx-1), 1), ncol=kx, nrow=1)
    
    xl.mat <- rbind(x.ends, xl.mid, x.ends)
    xr.mat <- rbind(x.ends, xr.mid, x.ends)
    
    #create matrix weights for y-component
    y.ends <- matrix(0, ncol=jx, nrow=ky)
    
    yt.mid <- matrix(c(1, rep(0, times=ky-1)), ncol=1, nrow=ky)
    yb.mid <- matrix(c(rep(0, times=ky-1), 1), ncol=1, nrow=ky)
    
    yt.mat <- cbind(y.ends, yt.mid, y.ends)
    yb.mat <- cbind(y.ends, yb.mid, y.ends)
    
    #use focal statistics for e, w, n, s components of the k-neighbourhood
    dz.dx.l <- focal(r, xl.mat, fun=sum, na.rm=FALSE)
    dz.dx.r <- focal(r, xr.mat, fun=sum, na.rm=FALSE)
    
    dz.dy.t <- focal(r, yt.mat, fun=sum, na.rm=FALSE)
    dz.dy.b <- focal(r, yb.mat, fun=sum, na.rm=FALSE)
    
    #calculate dz/dx and dz/dy using the components. 2*j is the run: (2 sides)*(length of each side)
    dz.dx <- (dz.dx.r-dz.dx.l)/(2*jx*res(r)[1])
    dz.dy <- (dz.dy.b-dz.dy.t)/(2*jy*res(r)[2])
  }
  
  out<- stack() #initialize output
  
  if("slope" %in% needed_metrics){
    slope.k<- (atan(sqrt((dz.dx^2)+(dz.dy^2))))
    if(mask_aspect){
      slp0_idx<- slope.k==0
      }
    if(unit=="degrees" & ("slope" %in% metrics)){
      slope.k<- slope.k * (180/pi)
      }
    
    names(slope.k)<- "slope"
    out<- stack(out, slope.k)
  }
  
  if("aspect" %in% needed_metrics){
    aspect.k<- atan2(dz.dy, -dz.dx)
    aspect.k<- raster::calc(aspect.k, fun= convert_aspect) #convert aspect to clockwise distance from North
    # aspect.k<- (pi/2)- aspect.k #Seems like this could replace convert_aspect
    # aspect.k[aspect.k >= (2*pi)]<- aspect.k[aspect.k >= (2*pi)] - (2*pi)
    # aspect.k[aspect.k < 0]<- aspect.k[aspect.k < 0] + (2*pi)
    
    if("eastness" %in% needed_metrics){
      eastness.k<- sin(aspect.k)
      if(mask_aspect){
        eastness.k[slp0_idx]<- 0 #Set eastness to 0 where slope is zero
      }
      names(eastness.k)<- "eastness"
      out<- stack(out, eastness.k)
    }
    
    if("northness" %in% needed_metrics){
      northness.k<- cos(aspect.k)
      if(mask_aspect){
        northness.k[slp0_idx]<- 0 #Set northenss to 0 where slope is zero
      }
      names(northness.k)<- "northness"
      out<- stack(out, northness.k)
    }
    
    if(mask_aspect){
      aspect.k[slp0_idx]<- NA_real_ #Set aspect to undefined where slope is zero
    }
    if(unit=="degrees"){
      aspect.k<- aspect.k * (180/pi)
      }
    names(aspect.k)<- "aspect"
    out<- stack(out, aspect.k)
    }
  
  out<- raster::subset(out, metrics, drop=TRUE) #Subset needed metrics to requested metrics in proper order
  if(include_scale){names(out)<- paste0(names(out), "_", w[1], "x", w[2])}
  return(out)
}