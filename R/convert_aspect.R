#' Helper function to convert aspect to clockwise distance from North
#'
#' #' Helper function to convert aspect to clockwise distance from North
#' @param aspect a number in radians representing aspect as calculated by atan2(dz.dy, -dz.dx)
#' @importFrom dplyr case_when
#' @return aspect as to clockwise distance in radians from North
#' @details Adapted from https://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-aspect-works.htm

convert_aspect<- function(aspect){
  out<- case_when(is.na(aspect) ~ NA_real_,
                  aspect < 0 ~ (pi/2 - aspect),
                  aspect > (pi/2) ~ (2*pi) - aspect + (pi/2),
                  TRUE ~ (pi/2)- aspect)
  return(out)
}