#' Calculate normal contour curvature
#'
#' Calculate normal contour curvature (kn)c, which is the principal representative of the plan curvature group based on regression coefficients from the equation Z =ax^2+by^2+cxy+dx+ey(+f).
#' @param a regression coefficient
#' @param b regression coefficient
#' @param c regression coefficient
#' @param d regression coefficient
#' @param e regression coefficient
#' @references
#' Evans, I.S., 1980. An integrated system of terrain analysis and slope mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295.
#' 
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' 
#' Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
knc<- function(a,b,c,d,e){
  #Planc: Normal Contour Curvature
  out<- -2*(a*(e^2) - c*d*e + b*d^2)/((d^2+e^2) * sqrt(1+d^2+e^2))
  return(out)
}

#' Calculate normal slope line curvature
#'
#' Calculate normal slope line curvature (kn)s, which is the principal representative of the profile curvature group based on regression coefficients from the equation Z =ax^2+by^2+cxy+dx+ey(+f).
#' @param a regression coefficient
#' @param b regression coefficient
#' @param c regression coefficient
#' @param d regression coefficient
#' @param e regression coefficient
#' @references
#' Evans, I.S., 1980. An integrated system of terrain analysis and slope mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295.
#' 
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' 
#' Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414

kns<- function(a,b,c,d,e){
  #Profc: Normal Slope Line Curvature
  out<- (-2 * (a*d^2 + c*d*e + b*e^2)) / ((d^2 + e^2)*(1 + d^2 + e^2)^1.5)
  return(out)
}

#' Calculate contour geodesic torsion
#'
#' Calculate contour geodesic torsion (tg)c, which is the principal representative of the twisting curvature group based on regression coefficients from the equation Z =ax^2+by^2+cxy+dx+ey(+f).
#' @param a regression coefficient
#' @param b regression coefficient
#' @param c regression coefficient
#' @param d regression coefficient
#' @param e regression coefficient
#' @references
#' Evans, I.S., 1980. An integrated system of terrain analysis and slope mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295.
#' 
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' 
#' Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414

tgc<- function(a,b,c,d,e){
  #TwistC: Contour geodesic torsion
  out<- (2*d*e*(a-b) - c*(d^2-e^2))/((d^2+e^2)*(1+d^2+e^2))
  return(out)
}

#' Calculate mean curvature
#'
#' Calculate mean curvature, kmean, from the equation Z =ax^2+by^2+cxy+dx+ey(+f).
#' @param a regression coefficient
#' @param b regression coefficient
#' @param c regression coefficient
#' @param d regression coefficient
#' @param e regression coefficient
#' @references
#' Evans, I.S., 1980. An integrated system of terrain analysis and slope mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295.
#' 
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' 
#' Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414

kmean<- function(a,b,c,d,e){
  #Mean Curvature
  out<- -(a*(1+e^2) - c*d*e + b*(1+d^2)) / (sqrt((1+d^2+e^2)^3))
  return(out)
}

#' Calculate unsphericity curvature
#'
#' Calculate  unsphericity curvature, ku, from the equation Z =ax^2+by^2+cxy+dx+ey(+f).
#' @param a regression coefficient
#' @param b regression coefficient
#' @param c regression coefficient
#' @param d regression coefficient
#' @param e regression coefficient
#' @references
#' Evans, I.S., 1980. An integrated system of terrain analysis and slope mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295.
#' 
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' 
#' Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414

ku<- function(a,b,c,d,e){
  #unsphericity curvature
  out<- sqrt(((a*(1+e^2) - c*d*e +b*(1+d^2)) / (2*sqrt((1+d^2+e^2)^3)))^2 - ((4*a*b-c^2)/(1+d^2+e^2)^2))
  return(out)
}

#' Calculate min curvature
#'
#' Calculate min curvature, kmin, from the equation Z =ax^2+by^2+cxy+dx+ey(+f).
#' @param a regression coefficient
#' @param b regression coefficient
#' @param c regression coefficient
#' @param d regression coefficient
#' @param e regression coefficient
#' @references
#' Evans, I.S., 1980. An integrated system of terrain analysis and slope mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295.
#' 
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' 
#' Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414

kmin<- function(a,b,c,d,e){
  #Min Curvature
  out<- kmean(a,b,c,d,e)-ku(a,b,c,d,e)
  return(out)
}

#' Calculate max curvature
#'
#' Calculate max curvature, kmax, from the equation Z =ax^2+by^2+cxy+dx+ey(+f).
#' @param a regression coefficient
#' @param b regression coefficient
#' @param c regression coefficient
#' @param d regression coefficient
#' @param e regression coefficient
#' @references
#' Evans, I.S., 1980. An integrated system of terrain analysis and slope mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295.
#' 
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' 
#' Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414

kmax<- function(a,b,c,d,e){
  #Max Curvature
  out<- kmean(a,b,c,d,e)+ku(a,b,c,d,e)
  return(out)
}