#' Interactive Shiny app to look at terrain attributes
#'
#' Interactive Shiny app to look at terrain attributes based on a surface fit using a Wood/Evans Quadratic Equation: Z =ax^2+by^2+cxy+dx+ey(+f)
#' @import shiny
#' @import rgl
#' @references
#' Evans, I.S., 1980. An integrated system of terrain analysis and slope mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295.
#' 
#' Wood, J., 1996. The geomorphological characterisation of digital elevation models (Ph.D.). University of Leicester.
#' 
#' Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
#' @export
explore_terrain<- function(){
  appDir<- system.file("Terrain_Attributes_Explorer_App", package = "MultiscaleDEM")
  shiny::runApp(appDir, display.mode = "normal")
  }