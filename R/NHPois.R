
#Non-homogeneous Poisson Process Generator!!
NHPois <- function(time_length, equation){
  
  expr<-function(x){
    x<-x
    eval(equation)
  }
  current <- 0
  i<-1
  data <- vector()
  maxrate = optim(par=time_length/2,fn=expr,method='L-BFGS-B',lower=0,upper=time_length,control = list(fnscale = -1))$value
  while(current < time_length){
    current<- current+ rexp(n=1,rate=maxrate)
    u <- runif(1)
    
    if (u < expr(current) / ( maxrate)){
      data[i] <- current
      i<-i+1
    }
    
  }
  return(data[c(-length(data))])
}

