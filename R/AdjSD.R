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
  if(length(w)==1){
    w<- rep(w,2)}
  if(length(w) > 2){
    stop("Specified window exceeds 2 dimensions")}
  if(any(0 == (w %% 2))){
    stop("Error: w must be odd")}
  if(all(w<3)){
    stop("Error: w must be greater or equal to 3 in at least one dimension")
  }
  if(prod(w) < 4){
    stop("Error: Window size must have at least 4 cells to fit surface and have residuals")
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
  
  #Explanatory Variable matrix X for linear fit
  X<- cbind(1, x, y)
  
  if(pad==TRUE){
    if(na.rm==FALSE){
      na.rm<- TRUE
      warning("if pad=TRUE, na.rm must also be TRUE. Changing na.rm to TRUE")
    }
    og_extent<- raster::extent(r)
    r<- raster::extend(r, c((w[1]-1)/2, (w[2]-1)/2), value=NA)
  }
  
  #Process large rasters as smaller chunks
  run_in_blocks<- !raster::canProcessInMemory(r, n = 2)
  if(run_in_blocks==FALSE){
    out<- raster(r)
    values(out)<- C_multiscale1(r = as.matrix(r), w= w, X=X, na_rm=na.rm)
  } else{
    f_out <- raster::rasterTmpFile()
    out<- raster::raster(r)
    out <- raster::writeStart(out, filename = f_out)
    
    block_idx<- raster::blockSize(r, n = 2, minblocks = 2, minrows = w[1])
    block_overlap<- (w[1]-1)/2
    nr<- nrow(r)
    nc<- ncol(r)
    for (i in 1:block_idx$n) {
      min_row<- max(c(block_idx$row[[i]] - block_overlap), 1)
      max_row<- min(c(block_idx$row[[i]] + block_idx$nrows[[i]] - 1 + block_overlap, nr))
      curr_block <- raster::getValues(r, row = min_row, nrows = max_row-min_row+1, format="matrix")
      
      out_block<- C_multiscale1(r = curr_block, w= w, X=X, na_rm=na.rm)
      #out_block is a formatted as vector where you move across rows in the raster object with each element (required for writeValues)
      if(i==1){
        out_block<- out_block[1:(length(out_block)-(block_overlap*nc))] #Trim bottom edge of raster
      } else if (i != block_idx$n){
        out_block<- out_block[(1+block_overlap*nc):(length(out_block)-(block_overlap*nc))] #Trim top and bottom edge of raster
      } else {
        out_block<- out_block[(1+block_overlap*nc):length(out_block)] #Trim top edge of raster
      }
      raster::writeValues(out, v= out_block, start= block_idx$row[i])
    }
    out<- raster::writeStop(out)
  }
  if(pad==TRUE){
    out<- raster::crop(out, og_extent)
  }
  names(out)<- "adjSD"
  if(include_scale){names(out)<- paste0(names(out), "_", w[1],"x", w[2])} #Add scale to layer names
  return(out)
}