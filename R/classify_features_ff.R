#' Helper function factory to classify morphometric features
#'
#' Helper function factory to classify morphometric features according to Wood 1996 page 120
#' @param slope_tolerance Slope tolerance that defines a 'flat' surface (degrees; default is 1.0). Relevant for the features layer.
#' @param curvature_tolerance Curvature tolerance that defines 'planar' surface (default is 0.01). Relevant for the features layer.
#' @importFrom dplyr case_when
#' @return A function that can be passed to raster::overlay to classify morphometric features
#' @references
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' @export

classify_features_ff<- function(slope_tolerance, curvature_tolerance){
  #1=PLANAR
  #2=PIT
  #3=CHANNEL
  #4=PASS
  #5=RIDGE
  #6=PEAK
  out_fun<- function(slope, crosc, maxic, minic){
    dplyr::case_when(is.na(slope) ~ NA_real_,
                     (slope > slope_tolerance) & (crosc > curvature_tolerance) ~ 5,
                     (slope > slope_tolerance) & (crosc < -curvature_tolerance) ~ 3,
                     slope > slope_tolerance ~ 1,
                     (maxic > curvature_tolerance) & (minic > curvature_tolerance) ~ 6,
                     (maxic > curvature_tolerance) & (minic < -curvature_tolerance) ~ 4,
                     maxic > curvature_tolerance ~ 5,
                     (minic < -curvature_tolerance) & (maxic < -curvature_tolerance) ~ 2,
                     minic < -curvature_tolerance ~ 3,
                     TRUE ~ 1)}
  return(out_fun)
}
