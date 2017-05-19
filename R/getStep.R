

#this function returns a BB object for the optimal blocks at any given cell location "n" where 1 < n < N
getStep <- function(BBobj,n){
  
  lasttype      <- BBobj$lasttype[1:n]
  lastchange    <- BBobj$lastchange[1:n]
  lastparams    <- BBobj$lastparams[1:n]
  
  
  left          <- vector()
  params        <- list()
  type          <- vector()
  type[1]       <- lasttype[n]
  left[1]       <- lastchange[n]
  params[[1]]   <- lastparams[[n]]
  i=1
  while (left[i] > 1) {
    
    left[i+1]      <- lastchange[left[i]]
    type[i+1]      <- lasttype[left[i]]
    params[[i+1]]  <- lastparams[[ left[i] ]]
    
    i <- i+1
  }
  left      <- rev(left)+1
  type      <- rev(type)
  params    <- rev(params)
  
  if (length(left) == 1){
    right <- n
  }
  else {
    right <- c(left[2:length(left)]-1,n)
  }
  
  BBdata   <- list("data"         = BBobj$data[1:n],
                   "N"            = BBobj$N[1:n],
                   "A"            = BBobj$A[1:n],
                   "left"         = left,
                   "right"        = right,
                   "type"         = type,
                   "params"       = params,
                   "opt"          = BBobj$opt[1:n],
                   "lastchange"   = lastchange,
                   "lastparams"   = lastparams,
                   "lasttype"     = lasttype,
                   "pruned"       = BBobj$pruned)
  
  BBobject <- structure(BBdata, class = "BB")
  
  return(BBobject)
}