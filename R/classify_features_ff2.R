#' Helper function factory to classify morphometric features
#'
#' Helper function factory to classify morphometric features according to Wood 1996 page 120
#' @param slope_tolerance Slope tolerance that defines a 'flat' surface (degrees; default is 1.0). Relevant for the features layer.
#' @param curvature_tolerance Curvature tolerance that defines 'planar' surface (default is 0.01). Relevant for the features layer.
#' @importFrom dplyr case_when
#' @return A function that can be passed to raster::overlay to classify morphometric features
#' @references
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.

classify_features_ff2<- function(slope_tolerance, curvature_tolerance){
  #0=PLANAR
  #1=PIT
  #2=CHANNEL
  #3=PASS
  #4=RIDGE
  #5=PEAK
  out_fun<- function(slope, crosc, maxic, minic){
    dplyr::case_when(is.na(slope) ~ NA_real_,
                     (slope > slope_tolerance) & (crosc > curvature_tolerance) ~ 4,
                     (slope > slope_tolerance) & (crosc < -curvature_tolerance) ~ 2,
                     slope > slope_tolerance ~ 0,
                     (maxic > curvature_tolerance) & (minic > curvature_tolerance) ~ 5,
                     (maxic > curvature_tolerance) & (minic < -curvature_tolerance) ~ 3,
                     maxic > curvature_tolerance ~ 4,
                     (minic < -curvature_tolerance) & (maxic < -curvature_tolerance) ~ 1,
                     minic < -curvature_tolerance ~ 2,
                     TRUE ~ 0)}
  return(out_fun)
}
