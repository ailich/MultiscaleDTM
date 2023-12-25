#' FocalCpp for use in parallel code
#'
#' FocalCpp for use in parallel code
#' @param x a packed SpatRaster
#' @param ... Other arguments for terra::focalCpp
#' @return a packed SpatRaster
#' @import terra
#'
focalCpp_parallel<- function(x,...){
  terra::wrap(terra::focalCpp(terra::unwrap(x),...))
}

#' Break raster into smaller chunks
#'
#' Break raster into smaller chunks
#' @param r SpatRaster
#' @param n_chunks Desired number of chunks to break raster into
#' @param buffer the number of rows above/below the cell value that the calculation needs access to. For cell by cell calculations this should be zero and for standard focal operations this would be (w-1)/2 where w is the number of rows in the focal window.
#' @return a list containing a dataframe specifying how rasters are chunked and a list of the chunked SpatRasters
#' @import terra
#' @importFrom dplyr mutate
#' @importFrom dplyr lead
#' @importFrom dplyr %>%

chunk_raster<- function(r, n_chunks, buffer){
  breaks<- data.frame(write_start = ceiling(seq(1, nrow(r)+1, length.out = n_chunks+1)))
  breaks<- breaks[!duplicated(breaks$write_start), ,drop=FALSE] #remove duplicates
  breaks<- breaks %>% dplyr::mutate(write_end=dplyr::lead(write_start, n=1)-1)
  breaks<- breaks %>% dplyr::mutate(chunk_start=write_start - buffer)
  breaks<- breaks %>% dplyr::mutate(chunk_end = write_end + buffer)
  breaks<- breaks[-nrow(breaks),]
  breaks$chunk_end[breaks$chunk_end > nrow(r)]<- nrow(r)
  breaks$chunk_start[breaks$chunk_start < 1]<- 1

  r_list<- vector(mode = "list", length = nrow(breaks))
  for (i in 1:nrow(breaks)) {
    r_list[[i]]<- terra::wrap(r[breaks$chunk_start[i]:breaks$chunk_end[i], ,drop=FALSE])
  }
  return(list(breaks, r_list))
}

#' Merge raster chunks from a list into a single raster
#'
#' Merge raster chunks from a list into a single raster
#' @param x A list of SpatRasters
#' @param breaks_df Dataframe output from chunk_raster
#' @return a SpatRaster
#' @import terra

combine_raster_chunks<- function(x, breaks_df){
  for (i in 1:length(x)) {
    st<- 1 + breaks_df$write_start[i] - breaks_df$chunk_start[i]
    end<- nrow(x[[i]]) - (breaks_df$chunk_end[i] - breaks_df$write_end[i])
    x[[i]]<- x[[i]][st:end, , , drop=FALSE]  # Trim chunk
  }
  out<- do.call(terra::merge, x)
}
