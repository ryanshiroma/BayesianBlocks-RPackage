#' Optimum Piecewise Segmentation of event arrival data
#'
#' This function runs the PELT algorithm to partition photon event data into 
#' piecewise models. It accepts time-to-event, time-to-spill, and binned data.
#'
#' @param data A vector of time-to-event cell arrival times \code{nx1}
#' @param N A vector containing the number of photons per cell \code{nx1}
#' @param c Block penalty scalar value \code{c>0} 
#' @param type A vector containing the block types to be conisdered in the model 
#' \itemize{
#'   \item "constant" for constant blocks \code{y=a},
#'   \item "linear" for linear function blocks \code{y=ax+b},
#'   \item "exponential" for exponential function blocks \code{y=ae^(bx)+c},
#'   \item "linear" for power function blocks \code{y=at^b+c}
#' }
#' @param alpha The significance level desired for model comparison
#' @param verbose Logical indicating whether you want algorithm progress output
#' @return A BB object containing all segmentation information
#' @examples
#' #create some random data
#' set.seed(1)
#' equation1 <- expression(40)
#' equation2 <- expression(40+100*x)
#' equation3 <- expression(240)
#' x1 <- NHPois(time_length =3,equation = equation1)
#' x2 <- NHPois(time_length =2,equation = equation2)+x1[length(x1)]
#' x3 <- NHPois(time_length =2,equation = equation3)+x2[length(x2)]
#' x=c(x1,x2,x3)
#' N <- rep(1,length(x))#let all cells be of size 1
#' 
#' #run model
#' model <- optinterval(x, N,3,type=c("constant","linear"),pvalue=0.05,verbose=TRUE)
#' 
#' #plot model
#' plot(model,show=c("blocks","points"),xlim=c(0,x[length(x)]),ylim=c(0,400))
optinterval = function(data,N,c,type=c("constant"),alpha=0.05,verbose=FALSE){
  start.time <- Sys.time()
  
  chi.lin    <- qchisq(1-pvalue,df=1)/2
  chi.pow    <- qchisq(1-pvalue,df=2)/2
  chi.exp    <- qchisq(1-pvalue,df=2)/2
  if("power" %in% type | "exponential" %in% type){
    maxpen   <- chi.pow +c  
  }else if("linear" %in% type){
    maxpen   <- chi.lin+c
  }else{
    maxpen   <- c
  }
  
  n          <- length(N)
  percent    <- floor(n/100)
  lowerbound <- data[1]-(data[2]-data[1])/2
  upperbound <- data[n]+(data[n]-data[n-1])/2
  xx         <- c(lowerbound,(data[1:(length(data)-1)] + data[2:(length(data))])/2,upperbound) # voronoi cell vertices
  A          <- diff(xx) # length of each voronoi cell
  
  opt           <- rep(0,n+1)
  lastchange    <- rep(1,n)
  changeA       <- rep(0,n) 
  changeN       <- rep(0,n)
  optint        <- matrix(-Inf,nrow=n,ncol=length(type))
  last          <- rep(0,length(type))
  lasttype      <- rep("None",n)
  lastparams    <- list()
  unpruned <- NULL
  endobj   <- rep(0,n)
  
  
  
  #begin looping through each point
  for (i in 1:n){
    
    unpruned          <- c(unpruned,i)
    
    ##### constant blocks ##########
    if ("constant" %in% type){
      changeA[unpruned]  <- changeA[unpruned] + A[i]
      changeN[unpruned]  <- changeN[unpruned] + N[i]
      optint[unpruned,which(type=="constant")]  <- opt[unpruned] + constantobjective(changeN[unpruned],changeA[unpruned])-c
      last[which(type=="constant")]             <- which.max(optint[unpruned,which(type=="constant")])
    }
    ################################
    
    ##### linear blocks ############
    if ("linear" %in% type){
      linblocks=list()
      for(j in unpruned){
        x=data[j:i]-max(data[j-1],0)
        linblocks[[j]]   <- linfit(N[j:i],x)
        #linblocks[[j]]   <- iwls(N[j:i],x)
        optint[j,which(type=="linear")]       <- opt[j] + linblocks[[j]][["cost"]]-c-chi.lin
      }
      last[which(type=="linear")]              <- which.max(optint[unpruned,which(type=="linear")])
    }
    
    ##### exp blocks ############
    if ("exponential" %in% type){
      expblocks=list()
      #print(unpruned)
      for(j in unpruned){
        x=data[j:i]-max(data[j-1],0)
        expblocks[[j]]   <- expfit(N[j:i],x)
        optint[j,which(type=="exponential")]       <- opt[j] + expblocks[[j]][["cost"]]-c-chi.exp
      }
      last[which(type=="exponential")]                  <- which.max(optint[unpruned,which(type=="exponential")])
    }
    ################################
    
    ##### power blocks ############
    if ("power" %in% type){
      powblocks=list()
      for(j in unpruned){
        x=data[j:i]-max(data[j-1],0)
        powblocks[[j]]   <- powfit(N[j:i],x)
        optint[j,which(type=="power")]       <- opt[j] + powblocks[[j]][["cost"]]-c-chi.pow
      }
      last[which(type=="power")]               <- which.max(optint[unpruned,which(type=="power")])
    }
    ################################  
    
    ##### FRED blocks ############
    if ("FRED" %in% type){
      FREDblocks=list()
      for(j in unpruned){
        x=data[j:i]-max(data[j-1],0)
        FREDblocks[[j]]   <- FREDfit(N[j:i],x)
        optint[j,which(type=="FRED")]       <- opt[j] + FREDblocks[[j]][["cost"]]-c-chi.exp
      }
      last[which(type=="FRED")]              <- which.max(optint[unpruned,which(type=="FRED")])
    }
    ################################
    
    ################################
    
    bestshape<-which.max(apply(optint,2,max))
    
    lastchange[i]     <- unpruned[bestshape]-1
    if(("FRED" %in% type) && (bestshape==which(type=="FRED"))){ #pow block is best
      lasttype[i]       <- "FRED"
      lastparams[[i]]   <- FREDblocks[[unpruned[which(type=="FRED")]]]
    }
    else if (("linear" %in% type) && (bestshape==which(type=="linear"))){#linear block is best
      lasttype[i]       <- "linear"
      lastparams[[i]]   <- linblocks[[unpruned[which(type=="linear")]]]
    }
    
    else if(("exponential" %in% type) && (bestshape==which(type=="exponential"))){ #exp block is best
      lasttype[i]       <- "exponential"
      lastparams[[i]]   <- expblocks[[unpruned[which(type=="exponential")]]]
    }
    
    else if(("power" %in% type) && (bestshape==which(type=="power"))){ #pow block is best
      lasttype[i]       <- "power"
      lastparams[[i]]   <- powblocks[[unpruned[which(type=="power")]]]
    }
    
    
    if(("constant" %in% type) && (bestshape==which(type=="constant"))){ #constant block is best
      lasttype[i]       <- "constant"
      lastparams[[i]]   <- list("b"    = sum(N[lastchange[i]:i])/sum(A[lastchange[i]:i]),
                                "cost" = constantobjective(sum(N[lastchange[i]:i]),sum(A[lastchange[i]:i])))
    }
    
    opt[i+1]          <- max(optint[unpruned,bestshape])
    unpruned          <- unpruned[((optint[unpruned,bestshape]+maxpen-opt[i+1])>0)]
    
    if((verbose==TRUE) && (i %% percent==0)){#print out the progress of the algorithm
      cat("\n",round(100*i/n,0),"% of points completed",round(100*(i-length(unpruned))/i,0),"% pruned ")
    }
    
  }
  BBdata   <- list("data"         = data, 
                   "N"            = N,
                   "A"            = A,
                   "opt"          = opt,
                   "lastchange"   = lastchange,
                   "lastparams"   = lastparams,
                   "lasttype"     = lasttype,
                   "pruned"       = setdiff(1:i,unpruned))
  
  
  
  print(Sys.time()-start.time)
  model <- getStep(BBdata,n)
  summary(model)
  return(model)
}

