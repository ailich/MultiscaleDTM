#' Helper function to calculates Wood/Evans regression parameters 
#'
#' Helper function to calculates Wood/Evans regression parameters 
#' @param r DEM as a matrix
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param type 1 (linear) or 2 (quadratic). Default is type 2
#' @param na.rm A logical vector indicating whether or not to remove NA values before calculations
#' @return a list containing matrices of regression parameters (a-f) and the standard deviation of residuals

WoodEvansHelper2<- function(r, w, type, resolution, na.rm){
  if(type!=1 & type!=2){
    stop("type must be equal to 1 (linear) or 2 (quadratic)")
  }
  
  #Define local coordinate system of window
  x_mat<- matrix(resolution[2], nrow = w[1], ncol=w[2])
  for (C in 1:w[2]) {
    x_mat[,C]<- x_mat[,C]*C
  }
  x_mat<- x_mat - mean(x_mat)
  x<- as.vector(x_mat)
  
  y_mat<- matrix(resolution[1], nrow = w[1], ncol=w[2])
  for (R in 1:w[1]) {
    y_mat[R,]<- y_mat[R,]*R
  }
  y_mat<- y_mat - mean(y_mat)
  y<- as.vector(y_mat)
  
  #Explanatory Variable matrix X
  if(type==2){
    X<- cbind(1, x^2,y^2, x*y, x, y)} else{
      X<- cbind(1, x*y, x, y)
    }
  
  
  #Run C++ code to calculate regression parameters based on a quadratic fit to elevation values and convert to raster stack
  if(type==1){
    params<- C_multiscale1(r = as.matrix(r), w = w, X = X, na_rm = na.rm)
  } else{
    params<- C_multiscale2(r = as.matrix(r), w = w, X = X, na_rm = na.rm)
    }
  return(params)
}