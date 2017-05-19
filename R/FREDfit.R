###########################FRED fitter

FREDfunc= function(x,t){
  A=exp(1)^x[1]
  lambda=exp(1)^x[2]
  tau1=exp(1)^x[3]
  tau2=exp(1)^x[4]
  b=exp(1)^x[5]
  sum(log(A*lmabda*exp(1)^(-tau1/t-t/tau2)))-(2*a*lambda*(4*tau1*tau2)^0.5)*besselI((4*tau1/tau2)^0.5,1)-b*t[length(t)]
}

FREDfit <- function(N,x){
  if(length(x)<5){
    return(list("A"=0,"lambda"=0,"tau1"=0,"tau2"=0,"b"=0,"cost"=-Inf))
  }
  x=rep(x,N)
  
  #fit the data normally
  fit <- optim(par=c(1,1,1,1,log(length(x)/x[length(x)])),fn=powfunc,control = list(fnscale = -1),t=x)
  pars=fit$par
  A=exp(1)^pars[[1]]
  lambda=exp(1)^pars[[2]]
  tau1=exp(1)^pars[[3]]
  tau2=exp(1)^pars[[4]]
  b=exp(1)^pars[[5]]
  cost=fit$value
  return(list("A"=A,"lambda"=lambda,"tau1"=tau1,"tau2"=tau2,"b"=b,"cost"=cost))
}