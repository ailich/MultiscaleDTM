#' IN DEVELOPMENT. Helper function factory to classify morphometric features
#'
#' IN DEVELOPMENT. Helper function factory to classify morphometric features according to Wood 1996. Function is not vectorized so it does not work with overlay as of now. 
#' @param slope_tolerance Slope tolerance that defines a 'flat' surface (degrees; default is 1.0). Relevant for the features layer.
#' @param curvature_tolerance Curvature tolerance that defines 'planar' surface (default is 0.01). Relevant for the features layer.

classify_features_ff<- function(slope_tolerance, curvature_tolerance){
  #1=PLANAR
  #2=PIT
  #3=CHANNEL
  #4=PASS
  #5=RIDGE
  #6=PEAK
  out_fun<- function(slope, crosc, maxic, minic){
    if(is.na(slope)){return(NA_integer_)}
    if (slope > slope_tolerance) {
      if (crosc > curvature_tolerance) {return (5)}
      if (crosc < -curvature_tolerance) {return (3)} else {
        return (1) }}
    if (maxic > curvature_tolerance){
      if (minic > curvature_tolerance) {return (6)}
      if (minic < -curvature_tolerance) {return (4)} else {
        return (5)}}
    if (minic < -curvature_tolerance) {
      if (maxic < -curvature_tolerance) {return (2)} else {
        return (3)}}
    return(1)}
  return(out_fun)
}



