#' Calculates multiscale slope and aspect
#'
#' Calculates Calculates multiscale slope and aspect of a DEM over a sliding rectangular window using a quadratic fit to the surface according to the Wood/Evans procedure. Note that aspect has been shifted by pi/2 (90 degrees) so that 0 degrees represents North rather than East as this is how many measures of aspect are reported (e.g. Horn's method), but differs from the way that the Wood/Evans method is traditionally reported.
#' @param r DEM as a raster layer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param unit "degrees" or "radians"
#' @param return_aspect logical indicating whether or not to return aspect in addition to Northness and Eastness (default is FALSE)
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @param pad logical value specifying whether rows/columns of NA's should be padded to the edge of the raster to remove edge effects (FALSE by default). If pad is TRUE, na.rm must be TRUE.
#' @param return_params logical indicating whether to return Wood/Evans regression parameters (default FALSE)
#' @param type Integer specifying whether to fit the surface with a linear (1) or quadratic (2) fit. Default is type 2.
#' @return a RasterStack
#' @import raster
#' @export

WoodEvans_SlopeAspect<- function(r, w, unit= "degrees", return_aspect= FALSE, na.rm=FALSE, pad=FALSE, return_params= FALSE, type=2){
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
  run_in_blocks<- !raster::canProcessInMemory(r, n = 6)
  if(run_in_blocks==FALSE){
    params<- WoodEvansHelper(r=r, w=w, type= type, na.rm=na.rm)
  } else{
    block_idx<- raster::blockSize(r, n = 6, minblocks = 2, minrows = w[1])
    out_blocks<- vector(mode = "list", length = block_idx$n)
    block_overlap<- w[1]-1
    for (i in 1:block_idx$n) {
      min_row<- block_idx$row[[i]]
      max_row<- min(min_row + block_idx$nrows[[i]] - 1 + block_overlap, nrow(r))
      block_extent<- raster::extent(r, min_row, max_row, 1, ncol(r))
      curr_block<- raster::crop(r, block_extent)
      out_blocks[[i]]<- WoodEvansHelper(r=curr_block, w=w, type=type, na.rm = na.rm)
    }
    params<- do.call(raster::merge, out_blocks)
    names(params)<- c("a", "b", "c", "d", "e", "f")
  }
  if(pad==TRUE){
    params<- raster::crop(params, og_extent)
  }
  
  #Use regression parameters to calculate slope and aspect
  slp<- atan(sqrt(params$d^2 + params$e^2))
  names(slp)<- paste0("QuadSlope", "_", w[1], "x", w[2])
  asp<- atan2(params$e,params$d) + pi/2
  names(asp)<- paste0("QuadAspect", "_", w[1], "x", w[2])
  eastness<- sin(asp)
  names(eastness)<- paste0("QuadEastness", "_", w[1], "x", w[2])
  northness<- cos(asp)
  names(northness)<- paste0("QuadNorthness", "_", w[1], "x", w[2])
  if(unit=="degrees"){
    slp<- slp*180/pi
    asp<- asp*180/pi
  }
  
  out<- stack(slp, eastness, northness)
  if(return_aspect){out<- stack(out, asp)}
  if(return_params){
    names(params)<- paste0(names(params), "_", w[1], "x", w[2])
    out<- stack(out, params)}
  return(out)
}



