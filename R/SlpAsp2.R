#' Multiscale Slope and Aspect
#'
#' Calculates multiscale slope and aspect based on the slope.k/aspect.k algorithm from Misiuk et al (2021) which extends classical formulations of slope restricted to a 3x3 window to multiple scales. The code from Misiuk et al (2021) was modified to allow for rectangular rather than only square windows.
#' @param r DTM as a SpatRaster or RasterLayer in a projected coordinate system where map units match elevation/depth units
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param unit "degrees" or "radians"
#' @param method "queen", "rook", or "boundary", indicating how many neighboring cells to use to compute slope for any cell. queen uses 8 neighbors (up, down, left, right, and diagonals), boundary uses all edge pixels of the focal, and rook uses 4 (up, down, left, right). Alternatively, instead of "queen" or "rook", method can be specified as 8 and 4 respectively. 
#' @param metrics a character string or vector of character strings of which terrain attributes to return. Default is to return all available attributes which are c("slope", "aspect", "eastness", "northness").
#' @param na.rm Logical indicating whether or not to remove NA values before calculations. Only applicable if method is "queen" or "8".
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @param mask_aspect A logical. When mask_aspect is TRUE (the default), if slope evaluates to 0, aspect will be set to NA and both eastness and northness will be set to 0. When mask_aspect is FALSE, when slope is 0 aspect will be pi/2 radians or 90 degrees which is the behavior of raster::terrain, and northness and eastness will be calculated from that.
#' @param filename character Output filename. Can be a single filename, or as many filenames as there are layers to write a file for each layer
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @import terra
#' @importFrom raster raster
#' @importFrom raster stack
#' @importFrom raster writeRaster
#' @return a SpatRaster or RasterStack of slope and/or aspect (and components of aspect)
#' @examples
#' r<- rast(volcano, extent= ext(2667400, 2667400 + 
#' ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
#' crs = "EPSG:27200")
#' slp_asp<- SlpAsp(r = r, w = c(5,5), unit = "degrees", 
#' method = "queen", metrics = c("slope", "aspect", 
#' "eastness", "northness"))
#' plot(slp_asp)
#' @details When method="rook", slope and aspect are computed according to Fleming and Hoffer (1979) and Ritter (1987). When method="queen", slope and aspect are computed according to Horn (1981). These are the standard slope algorithms found in many GIS packages but are traditionally restricted to a 3 x 3 window size. Misiuk et al (2021) extended these classical formulations  to multiple window sizes. This function modifies the code from Misiuk et al (2021) to allow for rectangular rather than only square windows and also added aspect.
#' @references
#' Fleming, M.D., Hoffer, R.M., 1979. Machine processing of landsat MSS data and DMA topographic data for forest cover type mapping (No. LARS Technical Report 062879). Laboratory for Applications of Remote Sensing, Purdue University, West Lafayette, Indiana.
#' 
#' Horn, B.K., 1981. Hill Shading and the Reflectance Map. Proceedings of the IEEE 69, 14-47.
#' 
#' Misiuk, B., Lecours, V., Dolan, M.F.J., Robert, K., 2021. Evaluating the Suitability of Multi-Scale Terrain Attribute Calculation Approaches for Seabed Mapping Applications. Marine Geodesy 44, 327-385. https://doi.org/10.1080/01490419.2021.1925789
#' 
#' Ritter, P., 1987. A vector-based slope and aspect generation algorithm. Photogrammetric Engineering and Remote Sensing 53, 1109-1111.
#' @import terra
#' @importFrom raster raster
#' @importFrom raster stack
#' @export

SlpAsp2 <- function(r, w=c(3,3), unit="degrees", method="queen", metrics= c("slope", "aspect", "eastness", "northness"), na.rm=FALSE, include_scale=FALSE, mask_aspect=TRUE, filename=NULL, overwrite=FALSE, wopt=list()){ 
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
  if(any(w<3)){
    stop("Error: w must be greater or equal to 3")
  }
  unit<- tolower(unit) #make lowercase
  if(!(unit %in% c("degrees", "radians"))){
    stop("unit must be 'degrees' or 'radians'")
  }
  
  if(method==4){method<- "rook"}
  if(method==8){method<- "queen"}
  
  # if(!(method %in% c("queen", "rook", "boundary"))){
  #   stop("method must be 'queen', 'rook', 'boundary', '8', or '4'")
  # }
  metrics<- tolower(metrics) #Make all lowercase
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
  if(method=="rook" & na.rm){
    warning("na.rm=TRUE is only relevant if `method` is 'queen' or '8'")
    na.rm<- FALSE
  }
  
  if(na.rm){
    Non_NA_rast<- !is.na(r)
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
    xl.end <- matrix(c(1, rep(NA_real_, times=kx-1)), ncol=kx, nrow=1)
    xr.end <- matrix(c(rep(NA_real_, times=kx-1), 1), ncol=kx, nrow=1)
    
    x.mids <- matrix(NA_real_, ncol=kx, nrow=ly)
    
    xl.mid <- matrix(c(2, rep(NA_real_, times=kx-1)), ncol=kx, nrow=1)
    xr.mid <- matrix(c(rep(NA_real_, times=kx-1), 2), ncol=kx, nrow=1)
    
    xl.mat <- rbind(xl.end, x.mids, xl.mid, x.mids, xl.end)
    xr.mat <- rbind(xr.end, x.mids, xr.mid, x.mids, xr.end)
    
    #create matrix weights for y-component
    yt.end <- matrix(c(1, rep(NA_real_, times=ky-1)), ncol=1, nrow=ky)
    yb.end <- matrix(c(rep(NA_real_, times=ky-1), 1), ncol=1, nrow=ky)
    
    y.mids <- matrix(NA_real_, ncol=lx, nrow=ky)
    
    yt.mid <- matrix(c(2, rep(NA_real_, times=ky-1)), ncol=1, nrow=ky)
    yb.mid <- matrix(c(rep(NA_real_, times=ky-1), 2), ncol=1, nrow=ky)
    
    yt.mat <- cbind(yt.end, y.mids, yt.mid, y.mids, yt.end)
    yb.mat <- cbind(yb.end, y.mids, yb.mid, y.mids, yb.end)
    
    #use focal statistics for e, w, n, s components of the k-neighbourhood
    dz.dx.l <- terra::focal(r, xl.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    dz.dx.r <- terra::focal(r, xr.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    dz.dy.t <- terra::focal(r, yt.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    dz.dy.b <- terra::focal(r, yb.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    
    if(!na.rm){
      #calculate dz/dx and dz/dy using the components. 8*j is the weighted run, or distance between ends: 4*j*2, or (4 values in each row)*(length of the side)*(2 sides)
      dz.dx <- (dz.dx.r-dz.dx.l)/(8*jx*terra::res(r)[1])
      dz.dy <- (dz.dy.t-dz.dy.b)/(8*jy*terra::res(r)[2])
    } else{
      weights.l<-terra::focal(Non_NA_rast, w=xl.mat, fun=sum, na.rm=TRUE)
      weights.r<-terra::focal(Non_NA_rast, w=xr.mat, fun=sum, na.rm=TRUE)
      weights.t<-terra::focal(Non_NA_rast, w=yt.mat, fun=sum, na.rm=TRUE)
      weights.b<-terra::focal(Non_NA_rast, w=yb.mat, fun=sum, na.rm=TRUE)
      dz.dx <- ((dz.dx.r/weights.r) - (dz.dx.l/weights.l))/(2*jx*terra::xres(r))
      dz.dy <- ((dz.dy.t/weights.t) - (dz.dy.b/weights.b))/(2*jy*terra::yres(r))
    }
  }
  
  
  if(method=="rook"){
    
    #create matrix weights for x-component
    x.ends <- matrix(NA_real_, ncol=kx, nrow=jy)
    
    xl.mid <- matrix(c(1, rep(NA_real_, times=kx-1)), ncol=kx, nrow=1)
    xr.mid <- matrix(c(rep(NA_real_, times=kx-1), 1), ncol=kx, nrow=1)
    
    xl.mat <- rbind(x.ends, xl.mid, x.ends)
    xr.mat <- rbind(x.ends, xr.mid, x.ends)
    
    #create matrix weights for y-component
    y.ends <- matrix(NA_real_, ncol=jx, nrow=ky)
    
    yt.mid <- matrix(c(1, rep(NA_real_, times=ky-1)), ncol=1, nrow=ky)
    yb.mid <- matrix(c(rep(NA_real_, times=ky-1), 1), ncol=1, nrow=ky)
    
    yt.mat <- cbind(y.ends, yt.mid, y.ends)
    yb.mat <- cbind(y.ends, yb.mid, y.ends)
    
    #use focal statistics for e, w, n, s components of the k-neighbourhood
    dz.dx.l <- terra::focal(r, xl.mat, fun=sum, na.rm=FALSE, wopt=wopt)
    dz.dx.r <- terra::focal(r, xr.mat, fun=sum, na.rm=FALSE, wopt=wopt)
    
    dz.dy.t <- terra::focal(r, yt.mat, fun=sum, na.rm=FALSE, wopt=wopt)
    dz.dy.b <- terra::focal(r, yb.mat, fun=sum, na.rm=FALSE, wopt=wopt)
    
    #calculate dz/dx and dz/dy using the components. 2*j is the run: (2 sides)*(length of each side)
    dz.dx <- (dz.dx.r-dz.dx.l)/(2*jx*terra::xres(r))
    dz.dy <- (dz.dy.t-dz.dy.b)/(2*jy*terra::yres(r))
  }
  
  if(method=="boundary"){
    #k is size of window in a given direction
    kx<- w[2]
    ky<- w[1]
    
    #j is the number of cells on either side of the focal cell; l is used to generate the focal matrix
    jx <- (kx/2)-0.5
    jy <- (ky/2)-0.5
    
    lx <- jx-1
    ly <- jy-1
    
    #create matrix weights for x-component
    
    xl.mat <- matrix(data=NA, nrow=ky, ncol=kx)
    xl.mat[,1]<- 2
    xl.mat[1,1]<- 1
    xl.mat[ky,1]<- 1
    print("xl.mat")
    print(xl.mat)
    
    xr.mat <- matrix(data=NA, nrow=ky, ncol=kx)
    xr.mat[,kx]<- 2
    xr.mat[1,kx]<- 1
    xr.mat[ky,kx]<- 1
    print("xr.mat")
    print(xr.mat)
    
    #create matrix weights for y-component
    yt.mat <- matrix(data=NA, nrow=ky, ncol=kx)
    yt.mat[1,]<- 2
    yt.mat[1,1]<- 1
    yt.mat[1,kx]<- 1
    print("yt.mat")
    print(yt.mat)
    
    yb.mat <- matrix(data=NA, nrow=ky, ncol=kx)
    yb.mat[ky,]<- 2
    yb.mat[ky,1]<- 1
    yb.mat[ky,kx]<- 1
    print("yb.mat")
    print(yb.mat)
    
    #use focal statistics for e, w, n, s components of the k-neighbourhood
    dz.dx.l <- terra::focal(r, xl.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    dz.dx.r <- terra::focal(r, xr.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    dz.dy.t <- terra::focal(r, yt.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    dz.dy.b <- terra::focal(r, yb.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    
    if(!na.rm){
      xw<- sum(xl.mat,na.rm=TRUE)+sum(xr.mat,na.rm=TRUE) #Sum of x weights
      yw<- sum(yt.mat,na.rm=TRUE)+sum(yb.mat,na.rm=TRUE) #Sum of y weights
      
      dz.dx <- (dz.dx.r-dz.dx.l)/(xw*jx*terra::res(r)[1])
      dz.dy <- (dz.dy.t-dz.dy.b)/(yw*jy*terra::res(r)[2])
    } else{
      weights.l<-terra::focal(Non_NA_rast, w=xl.mat, fun=sum, na.rm=TRUE, wopt=wopt)
      weights.r<-terra::focal(Non_NA_rast, w=xr.mat, fun=sum, na.rm=TRUE, wopt=wopt)
      weights.t<-terra::focal(Non_NA_rast, w=yt.mat, fun=sum, na.rm=TRUE, wopt=wopt)
      weights.b<-terra::focal(Non_NA_rast, w=yb.mat, fun=sum, na.rm=TRUE, wopt=wopt)
      dz.dx <- ((dz.dx.r/weights.r) - (dz.dx.l/weights.l))/(2*jx*terra::xres(r))
      dz.dy <- ((dz.dy.t/weights.t) - (dz.dy.b/weights.b))/(2*jy*terra::yres(r))
    }
  }
  
  if(method=="boundary2"){
    #create matrix weights
    
    mat_full<- matrix(data=NA, nrow=ky, ncol=kx)
    mat_full[1,]<-1
    mat_full[,1]<-1
    mat_full[ky,]<-1
    mat_full[,kx]<-1
    mat_full[jy+1,1]<-2
    mat_full[jy+1,kx]<-2
    
    mat_full[1, jx+1]<-2
    mat_full[ky, jx+1]<-2
    
    #create matrix weights for x-component
    
    xl.mat <- matrix(data=NA, nrow=ky, ncol=kx)
    xl.mat[, 1:jx]<- mat_full[, 1:jx]
    print("xl.mat")
    print(xl.mat)
    
    xr.mat <- matrix(data=NA, nrow=ky, ncol=kx)
    xr.mat[, (jx+2):kx]<- mat_full[, (kx-jx+1):kx]
    print("xr.mat")
    print(xr.mat)
    
    #create matrix weights for y-component
    yt.mat <- matrix(data=NA, nrow=ky, ncol=kx)
    yt.mat[1:jy,]<- mat_full[1:jy,]
    print("yt.mat")
    print(yt.mat)
    
    yb.mat <- matrix(data=NA, nrow=ky, ncol=kx)
    yb.mat[(jy+2):ky,]<- mat_full[(jy+2):ky,]
    print("yb.mat")
    print(yb.mat)
    
    
    #use focal statistics for e, w, n, s components of the k-neighbourhood
    dz.dx.l <- terra::focal(r, xl.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    dz.dx.r <- terra::focal(r, xr.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    dz.dy.t <- terra::focal(r, yt.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    dz.dy.b <- terra::focal(r, yb.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    
    if(!na.rm){
      nc<- (w[1]*2)+(w[2]*2)-4 # Number of border cells (also is equivalent sum of x or y weights)
      dz.dx <- (dz.dx.r-dz.dx.l)/(nc*jx*terra::res(r)[1])
      dz.dy <- (dz.dy.t-dz.dy.b)/(nc*jy*terra::res(r)[2])
    } else{
      weights.l<-terra::focal(Non_NA_rast, w=xl.mat, fun=sum, na.rm=TRUE)
      weights.r<-terra::focal(Non_NA_rast, w=xr.mat, fun=sum, na.rm=TRUE)
      weights.t<-terra::focal(Non_NA_rast, w=yt.mat, fun=sum, na.rm=TRUE)
      weights.b<-terra::focal(Non_NA_rast, w=yb.mat, fun=sum, na.rm=TRUE)
      dz.dx <- ((dz.dx.r/weights.r) - (dz.dx.l/weights.l))/(2*jx*terra::xres(r))
      dz.dy <- ((dz.dy.t/weights.t) - (dz.dy.b/weights.b))/(2*jy*terra::yres(r))
    }
  }
  
  if(method=="boundary3"){
    #k is size of window in a given direction
    kx<- w[2]
    ky<- w[1]
    
    #j is the number of cells on either side of the focal cell; l is used to generate the focal matrix
    jx <- (kx/2)-0.5
    jy <- (ky/2)-0.5
    
    lx <- jx-1
    ly <- jy-1
    
    #create matrix weights for x-component
    
    xl.mat <- matrix(data=NA, nrow=ky, ncol=kx)
    xl.mat[,1]<- c(seq(1,2, length.out=jy+1), seq(2, 1, length.out=jy+1)[-1])
    print("xl.mat")
    print(xl.mat)
    
    xr.mat <- matrix(data=NA, nrow=ky, ncol=kx)
    xr.mat[,kx]<- c(seq(1,2, length.out=jy+1), seq(2, 1, length.out=jy+1)[-1])
    print("xr.mat")
    print(xr.mat)
    
    #create matrix weights for y-component
    yt.mat <- matrix(data=NA, nrow=ky, ncol=kx)
    yt.mat[1,]<- c(seq(1,2, length.out=jx+1), seq(2, 1, length.out=jx+1)[-1])
    print("yt.mat")
    print(yt.mat)
    
    yb.mat <- matrix(data=NA, nrow=ky, ncol=kx)
    yb.mat[ky,]<- c(seq(1,2, length.out=jx+1), seq(2, 1, length.out=jx+1)[-1])
    print("yb.mat")
    print(yb.mat)
    
    #use focal statistics for e, w, n, s components of the k-neighbourhood
    dz.dx.l <- terra::focal(r, xl.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    dz.dx.r <- terra::focal(r, xr.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    dz.dy.t <- terra::focal(r, yt.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    dz.dy.b <- terra::focal(r, yb.mat, fun=sum, na.rm=na.rm, wopt=wopt)
    
    if(!na.rm){
      xw<- sum(xl.mat,na.rm=TRUE)+sum(xr.mat,na.rm=TRUE) #Sum of x weights
      yw<- sum(yt.mat,na.rm=TRUE)+sum(yb.mat,na.rm=TRUE) #Sum of y weights
      
      dz.dx <- (dz.dx.r-dz.dx.l)/(xw*jx*terra::res(r)[1])
      dz.dy <- (dz.dy.t-dz.dy.b)/(yw*jy*terra::res(r)[2])
    } else{
      weights.l<-terra::focal(Non_NA_rast, w=xl.mat, fun=sum, na.rm=TRUE, wopt=wopt)
      weights.r<-terra::focal(Non_NA_rast, w=xr.mat, fun=sum, na.rm=TRUE, wopt=wopt)
      weights.t<-terra::focal(Non_NA_rast, w=yt.mat, fun=sum, na.rm=TRUE, wopt=wopt)
      weights.b<-terra::focal(Non_NA_rast, w=yb.mat, fun=sum, na.rm=TRUE, wopt=wopt)
      dz.dx <- ((dz.dx.r/weights.r) - (dz.dx.l/weights.l))/(2*jx*terra::xres(r))
      dz.dy <- ((dz.dy.t/weights.t) - (dz.dy.b/weights.b))/(2*jy*terra::yres(r))
    }
  }
  
  out<- terra::rast() #initialize output
  
  if("slope" %in% needed_metrics){
    slope.k<- terra::math(terra::math(dz.dx^2 + dz.dy^2, fun= "sqrt", wopt=wopt), fun= "atan", wopt=wopt)
    if(unit=="degrees" & ("slope" %in% metrics)){
      slope.k<- slope.k * (180/pi)
      }
    names(slope.k)<- "slope"
    out<- c(out, slope.k, warn=FALSE)
  }
  
  if("aspect" %in% needed_metrics){
    aspect.k<- (-pi/2) - terra::atan_2(dz.dy, dz.dx, wopt=wopt) # aspect relative to North
    aspect.k<- ifel(aspect.k < 0, yes = aspect.k+(2*pi), no= aspect.k, wopt=wopt) # Constrain range so between 0 and 2pi
    
    if("eastness" %in% needed_metrics){
      eastness.k<- terra::math(aspect.k, fun="sin", wopt=wopt)
      if(mask_aspect){
        eastness.k<- terra::mask(eastness.k, mask= slope.k, maskvalues = 0, updatevalue = 0, wopt=wopt) #Set eastness to 0 where slope is zero
        }
      names(eastness.k)<- "eastness"
      out<- c(out, eastness.k, warn=FALSE)
    }
    
    if("northness" %in% needed_metrics){
      northness.k<- terra::math(aspect.k, fun="cos", wopt=wopt)
      if(mask_aspect){
        northness.k<- terra::mask(northness.k, mask= slope.k, maskvalues = 0, updatevalue = 0, wopt=wopt) #Set northenss to 0 where slope is zero
        }
      names(northness.k)<- "northness"
      out<- c(out, northness.k, warn=FALSE)
    }
    
    if(mask_aspect){
      aspect.k<- terra::mask(aspect.k, mask= slope.k, maskvalues = 0, updatevalue = NA, wopt=wopt) #Set aspect to undefined where slope is zero
      }
    if(unit=="degrees"){
      aspect.k<- aspect.k * (180/pi)
      }
    names(aspect.k)<- "aspect"
    out<- c(out, aspect.k, warn=FALSE)
    }
  
  if(!is.null(metrics)){out<- terra::subset(out, metrics, wopt=wopt)} #Subset needed metrics to requested metrics in proper order
  if(include_scale){names(out)<- paste0(names(out), "_", w[1], "x", w[2])}
  
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
