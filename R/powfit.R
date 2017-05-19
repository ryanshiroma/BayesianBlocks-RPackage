##############################power fitter
powfunc= function(x,t){
  a=x[1]
  b=x[2]
  c=x[3]
  tn=t[length(t)]
  eb=exp(1)^b
  sum(log(a*t^eb+c))-((a*tn^(eb+1))/(eb+1)+c*tn)
}

powgunc= function(x,t){
  a=x[1]
  b=x[2]
  c=x[3]
  eb=exp(1)^b
  atbc=a*t^(eb)+c
  teb=t^eb
  
  tn=t[length(t)]
  logt=log(t)
  fa=sum(teb/atbc)-((tn^(eb+1))/(eb+1))
  fb=sum((a*eb*t^(eb+1)*(eb*logt+logt-1))/(eb+1)^2)-(a*eb*tn^(eb+1)*((eb+1)*log(tn)-1))/(eb+1)^2
  fc=sum(1/atbc)-tn
  return(c(fa,fb,fc))
}


powfit <- function(N,x){
  if(length(x)<15){
    return(list("a"=0,"b"=0,"c"=length(x)/x[length(x)],"direction"="forward","cost"=constantobjective(length(N),x[length(x)])))
  }
  x=rep(x,N)
  
  #fit the data normally
  fit <- optim(par=c(1,1,length(x)/x[length(x)]),fn=powfunc,gr=powgunc,control = list(fnscale = -1),t=x)
  pars=fit$par
  a=pars[[1]]
  b=exp(1)^pars[[2]]
  c=pars[[3]]
  cost=fit$value
  forward=list("a"=a,"b"=b,"c"=c,"direction"="forward", "cost"=cost)
  if(sum((a*c(0,x)^b+c)<0)>0){
    forward$cost=-Inf
  }
  
  #fit the data in reverse
  x=x[length(x)]-rev(x[1:length(x)-1])
  fit <- optim(par=c(1,1,length(x)/x[length(x)]),fn=powfunc,gr=powgunc,control = list(fnscale = -1),t=x)
  pars=fit$par
  a=pars[[1]]
  b=exp(1)^pars[[2]]
  c=pars[[3]]
  cost=fit$value
  reverse=list("a"=a,"b"=b,"c"=c,"direction" = "reverse","cost"=cost)
  if(sum((a*c(0,x)^b+c)<0)>0){
    reverse$cost=-Inf
  }
  
  if (forward$cost>reverse$cost){
    return(forward)
  }else{
    return(reverse)
  }
}