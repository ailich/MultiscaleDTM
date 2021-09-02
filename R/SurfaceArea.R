#' Calculates surface area of a DEM
#'
#' Calculates surface area on a per cell basis of a DEM based on Jenness, 2004. This wrapper for sp::surfaceArea that natively works on rasters.
#' @param r DEM as a raster layer
#' @return a RasterLayer
#' @import raster
#' @importFrom sp surfaceArea
#' @export

SurfaceArea<- function(r){
  run_in_blocks<- !raster::canProcessInMemory(r, n = 1)
  x_res<- res(r)[1]
  y_res<- res(r)[2]
  
  if(run_in_blocks==FALSE){
    SA<- sp::surfaceArea(as.matrix(r), cellx =x_res, celly=y_res, byCell=TRUE)
    SA<- raster(SA, template=r)
    } else{
    block_idx<- raster::blockSize(r, n = 1, minblocks = 2, minrows = 3)
    out_blocks<- vector(mode = "list", length = block_idx$n)
    block_overlap<- 2
    for (i in 1:block_idx$n) {
      min_row<- block_idx$row[[i]]
      max_row<- min(min_row + block_idx$nrows[[i]] - 1 + block_overlap, nrow(r))
      block_extent<- raster::extent(r, min_row, max_row, 1, ncol(r))
      curr_block<- raster::crop(r, block_extent)
      out_blocks[[i]]<- sp::surfaceArea(as.matrix(curr_block), cellx =x_res, celly=y_res, byCell=TRUE)
      if(i<block_idx$n){
        out_blocks[[i]][nrow(out_blocks[[i]]),]<- NA_real_
      } #Set to last row to NA so merging chunks occurs properly
      out_blocks[[i]]<- raster(out_blocks[[i]], template=curr_block)
    }
    SA<- do.call(raster::merge, out_blocks)
    }
  names(SA)<- "SA"
  return(SA)
  }