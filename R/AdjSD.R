#' Calculates standard deviation of bathymetry (a measure of rugosity) adjusted for slope
#'
#' Calculates standard deviation of bathymetry (a measure of rugosity). Using a sliding rectangular window a plane is fit to the data and the standard deviation of the residuals is calculated.
#' @param r DEM as a raster layer
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @param pad logical value specifying whether rows/columns of NA's should be padded to the edge of the raster to remove edge effects (FALSE by default). If pad is TRUE, na.rm must be TRUE.
#' @param include_scale logical indicating whether to append window size to the layer names (default = FALSE)
#' @return a RasterLayer of adjusted rugosity
#' @import raster
#' @export

AdjSD<- function(r, w=c(3,3), na.rm=FALSE, pad=FALSE, include_scale=FALSE){
  #Input checks
  if(length(w==1)){
    w<- rep(w,2)}
  if(any(w<3) | any(0 == (w %% 2))){
    stop("Error: w must be odd and greater than or equal to 3")}
  if(pad==TRUE){
    if(na.rm==FALSE){
      na.rm<- TRUE
      warning("if pad=TRUE, na.rm must also be TRUE. Changing na.rm to TRUE")
    }
    og_extent<- raster::extent(r)
    r<- raster::extend(r, c((w[1]-1)/2, (w[2]-1)/2), value=NA)
  }
  
  #Process large rasters as smaller chunks
  run_in_blocks<- !raster::canProcessInMemory(r, n = 6)
  if(run_in_blocks==FALSE){
    params<- WoodEvansHelper(r=r, w=w, type= 1, na.rm=na.rm)
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
  }
  if(pad==TRUE){
    params<- raster::crop(params, og_extent)
  }
  names(params)<- c("a", "b", "c", "d", "e", "f", "adjSD")
  out<- params$adjSD
  if(include_scale){names(out)<- paste0(names(out), "_", w[1],"x", w[2])} #Add scale to layer names
  return(out)
  }
  



