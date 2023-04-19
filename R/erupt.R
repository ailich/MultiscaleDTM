#' Create georeferenced version of R's built in volcano dataset
#'
#' Create georeferenced version of R's built in volcano dataset. Useful dataset for generating quick examples.
#' @return SpatRaster
#' @examples 
#' r<- erupt()
#' @import terra
#' @export
erupt<- function(){
  return(rast(volcano, 
              extent= ext(2667400, 2667400 + ncol(volcano)*10, 
                          6478700, 6478700 + nrow(volcano)*10), 
              crs = "EPSG:27200"))
}