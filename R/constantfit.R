#####constant fitter

constantobjective = function(N,A){
  value = N*(log(N/A))-N
  value[is.nan(value)]<-0
  return(value)
}