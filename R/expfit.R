##############################exponential fitter
exp.f.iwls= function(x,Yk,xk){
  a=x[1]
  b=x[2]
  c=x[3]
  f=a*exp(1)^(b*xk)+c
  total=sum(1/f)
  sum((f-Yk)^2/(f*total))
}

exp.g.iwls= function(x,Yk,xk){
  a=x[1]
  b=x[2]
  c=x[3]
  fa=exp(1)^(b*xk)
  fb=a*xk*exp(1)^(b*xk)
  fc=1
  f=a*exp(1)^(b*xk)+c
  f2=f^2
  f3=f2*f
  y2=Yk^2
  total=sum(1/f)
  part1=total*(f2-y2)
  part2=(f3-2*y*f2+y2*f)
  denom=f2*total^2
  da=fa*part1-sum(fa/f2)*part2
  db=fb*part1-sum(fb/f2)*part2
  dc=fc*part1-sum(fc/f2)*part2
  return(c(da,db,dc)/denom)
}

exp.f.ML= function(x,t){
  a=x[1]
  b=x[2]
  c=x[3]
  tn=t[length(t)]
  eb=exp(1)^b
  ec=exp(1)^c
  sum(log(a*exp(1)^(eb*t)+ec))-(a/(eb)*(exp(1)^(eb*tn)-1)+ec*tn)
}

exp.g.ML= function(x,t){
  a=x[1]
  b=x[2]
  c=x[3]
  eb=exp(1)^b
  ec=exp(1)^c
  teb=exp(1)^(t*b)
  tn=t[length(t)]
  eteb=exp(1)^(t*eb)
  denom=(a*eteb+ec)
  etneb=exp(1)^(tn*eb)
  fa=sum(eteb/denom)-(1/eb)*(exp(1)^(eb*tn)-1)
  fb=sum((a*t*eb*eteb)/denom)-(a/eb)*(tn*eb*etneb-etneb+1)
  fc=sum(1/denom)-tn*ec
  return(c(fa,fb,fc))
}

expfit <- function(N,x){
  x=rep(x,N)
  if(length(x)<15){
    return(list("a"=0,"b"=0,"c"=length(x)/x[length(x)],"direction"="forward","cost"=constantobjective(length(N),x[length(x)])))
  }
  #fit the data normally
  fit <- optim(par=c(0,0,log(length(x)/x[length(x)])),fn=exp.f.ML,gr=exp.g.ML,control = list(fnscale = -1),t=x)
  pars=fit$par
  a=pars[[1]]
  b=exp(1)^pars[[2]]
  c=exp(1)^pars[[3]]
  cost=fit$value
  forward=list("a"=a,"b"=b,"c"=c,"direction"="forward", "cost"=cost)
  
  #fit the data in reverse
  x=x[length(x)]-rev(x[1:length(x)-1])
  fit <- optim(par=c(0,0,log(length(x)/x[length(x)])),fn=exp.f.ML,gr=exp.g.ML,control = list(fnscale = -1),t=x)
  pars=fit$par
  a=pars[[1]]
  b=exp(1)^pars[[2]]
  c=exp(1)^pars[[3]]
  #cost=sum(log(a*exp(1)^(b*x)+c))-((a/b)*(exp(1)^(b*x[length(x)])-1)+c*x[length(x)])
  cost=fit$value
  reverse=list("a"=a,"b"=b,"c"=c,"direction" = "reverse","cost"=cost)
  
  if (forward$cost>reverse$cost){
    return(forward)
  }else{
    return(reverse)
  }
  
}

iwls = function(N,t) {
  if (length(t)<40){ # if no data is recieved, set cost to 0
    return(list("a"=0,"b"=0,"cost"=-Inf))
  }  
  #if (length(t)<){ #if only one point is supplied, treat as constant intensity
  #  return(list("a"=0,"b"=N/t[1],"cost"=log(N/t[1])))
  #}
  x=rep(t,N)
  bins=as.integer(max(length(x)/7,2))
  # number of breaks needs to depend on the time length to give accurate estimates
  histogram = hist(x, breaks =seq(0,x[length(x)],length.out=bins+1), plot = FALSE)
  L = max(histogram$breaks)         # total interval length
  N = length(histogram$counts)      # number of blocks on this interval
  binwidth=L/N
  k = 1:N                           # k'th block index (1<=k<=N)
  xk = (k - 1/2)*binwidth            # midpoints of blocks (from page 367)
  Yk = histogram$counts/binwidth            # event counts in bins
  epsilon=1
  est.weights=rep(1,N)
  M=x[length(x)]
  S=M^2/2
  old.cost=-Inf
  while(epsilon>1e-4){
    m = lm(Yk ~ xk,weights=est.weights)
    est.weights = 1/m$fitted.values
    coef=m$coefficients
    new.cost=(sum(log(coef[[2]]*x+coef[[1]]))-coef[[2]]*S-coef[[1]]*M)
    if (is.nan(new.cost) || sum(est.weights<0)>0){
      return(list("a"=0,"b"=0,"cost"=-Inf))
    }
    epsilon = new.cost-old.cost
    old.cost=new.cost
  }
  
  list("a"=coef[2],"b"=coef[1],"cost"=new.cost)
}

