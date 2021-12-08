knc<- function(a,b,c,d,e,f){
  #Planc: Normal Contour Curvature
  out<- -(2*a*(e^2) - 2*c*d*e + 2*b*d^2)/((d^2+e^2) * sqrt(1+d^2+e^2))
  return(out)
}

kns<- function(a,b,c,d,e,f){
  #Profc: Normal Slope Line Curvature
  out<- (-2 * (a*d^2 + b*e^2 + c*d*e)) / ((e^2 + d^2)*(1 + e^2 + d^2)^1.5)
  return(out)
}

tgc<- function(a,b,c,d,e,f){
  #TwistC: Contour geodesic torsion
  out<- (d*e*(2*a-2*b) - c*(d^2-e^2))/((d^2+e^2)*(1+d^2+e^2))
  return(out)
}

kmean<- function(a,b,c,d,e,f){
  #Mean Curvature
  out<- -((1+e^2)*2*a - 2*c*d*e + (1+d^2)*2*b) / (2*sqrt((1+d^2+e^2)^3))
  return(out)
}

ku<- function(a,b,c,d,e,f){
  #unsphericity curvature
  out<- sqrt((((1+e^2)*2*a - 2*c*d*e +(1+d^2)*2*b) / (2*sqrt((1+d^2+e^2)^3)))^2 - ((2*a*2*b-c^2)/(1+d^2+e^2)^2))
  return(out)
}

kmin<- function(a,b,c,d,e,f){
  #Min Curvature
  out<- kmean(a,b,c,d,e)-ku(a,b,c,d,e)
  return(out)
}

kmax<- function(a,b,c,d,e,f){
  #Max Curvature
  out<- kmean(a,b,c,d,e)+ku(a,b,c,d,e)
  return(out)
}