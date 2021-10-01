#' Calculates surface area of a DEM
#'
#' Calculates surface area on a per cell basis of a DEM based on Jenness, 2004. This wrapper for sp::surfaceArea that natively works on rasters.
#' @param r DEM as a raster layer
#' @return a RasterLayer
#' @import raster
#' @importFrom sp surfaceArea
#' @export

SurfaceArea2<- function(r, run_in_blocks){
  #run_in_blocks<- !raster::canProcessInMemory(r, n = 2)
  x_res<- res(r)[1]
  y_res<- res(r)[2]
  
  if(run_in_blocks==FALSE){
    SA<- sp::surfaceArea(as.matrix(r), cellx =x_res, celly=y_res, byCell=TRUE)
    SA<- raster::raster(SA, template=r)
    } else{
      f_out <- raster::rasterTmpFile()
      SA<- raster::raster(r)
      SA <- raster::writeStart(SA, filename = f_out)
      
      block_idx<- raster::blockSize(r, n = 2, minblocks = 2, minrows = 3)
      block_overlap<- 1
      nr<- nrow(r)
    for (i in 1:block_idx$n) {
      min_row<- max(c(block_idx$row[[i]] - block_overlap), 1)
      max_row<- min(c(block_idx$row[[i]] + block_idx$nrows[[i]] - 1 + block_overlap, nr))
      curr_block <- raster::getValues(r, row = min_row, nrows = max_row-min_row+1, format="matrix")
      out_block<- sp::surfaceArea(curr_block, cellx =x_res, celly=y_res, byCell=TRUE)
      if(i==1){
        out_block<- out_block[-nrow(out_block),] #Trim last row
      } else if (i != block_idx$n){
        out_block<- out_block[2:(nrow(out_block)-1),] #Trim first and last row
      } else {
        out_block<- out_block[-1,] #Trim first row
      }
      SA<- raster::writeValues(SA, v= as.vector(t(out_block)), start= block_idx$row[i])
    }
      SA<- raster::writeStop(SA)
      }
  names(SA)<- "SA"
  return(SA)
  }