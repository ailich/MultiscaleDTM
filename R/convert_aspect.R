#' Helper function to convert aspect to clockwise distance from North
#'
#' #' Helper function to convert aspect to clockwise distance from North (0 to 2pi) from counterclockwise distance from east (-pi to pi)
#' @param aspect a number in radians representing aspect as calculated by atan2(dz.dy, -dz.dx)
#' @importFrom dplyr case_when
#' @return aspect as to clockwise distance in radians from North
#' @details Adapted from https://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-aspect-works.htm

convert_aspect<- function(aspect){
  out<- dplyr::case_when(is.na(aspect) ~ NA_real_,
                  aspect > (pi/2) ~ (2*pi) - aspect + (pi/2),
                  TRUE ~ (pi/2)- aspect)
  return(out)
}

#' Helper function to convert aspect to clockwise distance from North
#'
#' #' Helper function to 
#' @param aspect a number in radians representing aspect as calculated by atan2(e,d)
#' @importFrom dplyr case_when
#' @return aspect as to clockwise distance in radians from North

convert_aspect2<- function(aspect){
  out<- dplyr::case_when(is.na(aspect) ~ NA_real_,
                         aspect <= (-pi/2) ~ (-pi/2) - aspect,
                         TRUE ~ (-pi/2) - aspect + (2*pi))
  return(out)
}