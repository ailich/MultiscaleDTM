#' Helper function to convert aspect to clockwise distance from North
#'
#' #' Helper function to convert aspect to clockwise distance from North
#' @param aspect a number in radians representing aspect as calculated by atan2(dz.dy, -dz.dx)
#' @return aspect as to clockwise distance in radians from North
#' @details Adapted from https://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-aspect-works.htm

convert_aspect<- function(aspect){
  if(is.na(aspect)){
    return(NA_real_)
  } else if(aspect < 0){
    return(pi/2 - aspect)
  } else if (aspect > (pi/2)){
    return((2*pi) - aspect + (pi/2))
  } else{
    return((pi/2)- aspect)}
}