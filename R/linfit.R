#######linear fitter
linfit <-function(N,t){
  
  if (length(t)==0){ # if no data is recieved, set cost to 0
    return(list("a"=0,"b"=0,"cost"=-Inf))
  }  
  if (length(t)==1){ #if only one point is supplied, treat as constant intensity
    return(list("a"=0,"b"=N/t[1],"cost"=log(N/t[1])))
  }
  
  t=rep(t,N)
  
  epsilon  <- 1
  n        <- length(t)
  i        <- 1
  M        <- t[n]
  S        <- M^2/2
  
  #start with some initial values for a and b
  coef     <- matrix(c(0,length(t)/M),nrow=2)
  
  new.cost <- -Inf
  
  while (abs(epsilon) >10e-4 && i<100){
    old.cost <- new.cost
    
    #first derivatives
    f=matrix(c(sum(t/(coef[1]*t+coef[2]))-S,sum(1/(coef[1]*t+coef[2]))-M),nrow=2)
    
    #hessian values
    fa <- -sum((t/(coef[1]*t+coef[2]))^2)
    fb <- -sum(t/(coef[1]*t+coef[2])^2)
    ga <- fb
    gb <- -sum(1/(coef[1]*t+coef[2])^2)
    
    #create the inverse of the hessian
    invhess <- 1/(fa*gb-fb*ga)*matrix(c(gb,-fb,-ga,fa),nrow=2,ncol=2,byrow=TRUE) 
    
    #run newtons method
    coef=coef-invhess%*%f 
    
    if(coef[2]<0 || (M*coef[1]+coef[2])<0){#if best line has a negative intensity, dont use linear
      
      return(list("a"=0,"b"=length(t)/t[length(t)],"cost"=constantobjective(length(t),t[length(t)])))
      
    }
    #calculate the new cost
    logsum=coef[1]*t+coef[2]
    logsum=logsum[logsum>0]
    new.cost <- sum(log(logsum))-coef[1]*S-coef[2]*M
    epsilon  <- new.cost-old.cost
    
    i <- i+1
  }
  list("a"=coef[1],"b"=coef[2],"cost"=new.cost)
}
