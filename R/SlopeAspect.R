#' Multiscale Slope and Aspect
#'
#' Calculates multiscale slope and aspect based on the slope.k/aspect.k algorithm from Misiuk et al (2021) which extends classical formulations of slope restricted to a 3x3 window to multiple scales. The code from Misiuk et al (2021) was modified to allow for rectangular rather than only square windows.
#' @param r DEM as a raster layer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param unit "degrees" or "radians"
#' @param method "queen" or "rook", indicating how many neighboring cells to use to compute slope for any cell. queen uses 8 neighbors (up, down, left, right, and diagonals) and rook uses 4 (up, down, left, right).
#' @param metrics a character string or vector of character strings of which terrain atrributes to return ("slope" and/or "aspect"). Default is c("slope", "aspect").
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @param mask_aspect A logical. If slope evaluates to 0, aspect will be set to NA when mask_aspect is TRUE (the default). If FALSE, when slope is 0 aspect will be pi/2 radians or 90 degrees which is the behavior of raster::terrain.
#' @return a RasterStack or RasterLayer of slope and/or aspect
#' @details When method="rook", slope and aspect are computed according to Fleming and Hoffer (1979) and Ritter (1987). When method="queen", slope and aspect are computed according to Horn (1981). These are the standard slope algorithms found in many GIS packages but are traditionally restricted to a 3 x 3 window size. Misiuk et al (2021) extended these classical formulations  to multiple window sizes. This function modifies the code from Misiuk et al (2021) to allow for rectangular rather than only square windows and also added aspect.
#' @import raster
#' @export

SlopeAspect <- function(r, w=c(3,3), unit="degrees", method="queen", metrics= c("slope", "aspect", "northness", "eastness"), include_scale=FALSE, mask_aspect=TRUE){ 
  if(any(w < 3)){
    stop("w must be >=3 and odd")
  }
  
  if(length(w)==1){w<- rep(w,2)}
  
  if(any(w %% 2 != 1)){
    stop("w must be odd")
  }
  
  if(!(unit %in% c("degrees", "radians"))){
    stop("unit must be 'degrees' or 'radians'")
  }
  
  if(!(method %in% c("queen", "rook"))){
    stop("method must be 'queen' or 'rook'")
  }
  
  if(any(!(metrics %in% c("slope", "aspect","northness", "eastness")))){
    stop("metrics must be 'slope', 'aspect', 'northness', and/or 'eastness'")
  }
  
  #k is size of window in a given direction
  kx<- w[1]
  ky<- w[2]
  
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
    xr.mat <- rbind(xr.end, x.mids, xr.mid, x.mids, xr.end)
    
    #create matrix weights for y-component
    yt.end <- matrix(c(1, rep(0, times=ky-1)), ncol=1, nrow=ky)
    yb.end <- matrix(c(rep(0, times=ky-1), 1), ncol=1, nrow=ky)
    
    y.mids <- matrix(0, ncol=lx, nrow=ky)
    
    yt.mid <- matrix(c(2, rep(0, times=ky-1)), ncol=1, nrow=ky)
    yb.mid <- matrix(c(rep(0, times=ky-1), 2), ncol=1, nrow=ky)
    
    yt.mat <- cbind(yt.end, y.mids, yt.mid, y.mids, yt.end)
    yb.mat <- cbind(yb.end, y.mids, yb.mid, y.mids, yb.end)
    
    #use focal statistics for e, w, n, s components of the k-neighbourhood
    dz.dx.l <- focal(r, xl.mat, fun=sum)
    dz.dx.r <- focal(r, xr.mat, fun=sum)
    
    dz.dy.t <- focal(r, yt.mat, fun=sum)
    dz.dy.b <- focal(r, yb.mat, fun=sum)
    
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
    dz.dx.l <- focal(r, xl.mat, fun=sum)
    dz.dx.r <- focal(r, xr.mat, fun=sum)
    
    dz.dy.t <- focal(r, yt.mat, fun=sum)
    dz.dy.b <- focal(r, yb.mat, fun=sum)
    
    #calculate dz/dx and dz/dy using the components. 2*j is the run: (2 sides)*(length of each side)
    dz.dx <- (dz.dx.r-dz.dx.l)/(2*jx*res(r)[1])
    dz.dy <- (dz.dy.b-dz.dy.t)/(2*jy*res(r)[2])
  }
  
  slope.k<- (atan(sqrt((dz.dx^2)+(dz.dy^2))))
  names(slope.k)<- "slope"
  
  aspect.k<- atan2(dz.dy, -dz.dx)
  aspect.k<- raster::calc(aspect.k, fun= convert_aspect)#convert aspect to clockwise distance from North
  if(mask_aspect){
    aspect.k[slope.k==0]<- NA_real_ #Set aspect to undefined where slope is zero
  }
  names(aspect.k)<- "aspect"
  
  northness.k<- cos(aspect.k)
  names(northness.k)<- "northness"
  
  eastness.k<- sin(aspect.k)
  names(eastness.k)<- "eastness"
  
  out<- stack(slope.k, aspect.k, northness.k, eastness.k)
  if(unit=="degrees"){
    out$slope<- out$slope * (180/pi)
    out$aspect<- out$aspect * (180/pi)
  }
  
  out<- raster::subset(out, metrics, drop=TRUE)
  if(include_scale){names(out)<- paste0(names(out), "_", w[1], "x", w[2])}
  return(out)
}
