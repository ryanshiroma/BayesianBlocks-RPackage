#S3 generic plot function 

"plot.BB" <- function(x,show=c("hist","blocks"),binwidth=NULL,bins=NULL,ylim=NULL,xlim=NULL,xact=NULL,main="Bayesian Blocks",xlab="Time",ylab="Intensity") {
  
  data       <- x$data
  n          <- length(data)
  intensity  <- round(x$N/x$A,4)
  legend     <- vector()
  legend.col <- vector()
  
  #check if xlim is supplied
  if (!is.null(xlim)){
    lowerbound=xlim[1]
    upperbound=xlim[2]
  }
  else{
    lowerbound <- data[1]-(data[2]-data[1])/2
    upperbound <- data[n]+(data[n]-data[n-1])/2  
  }
  
  #check if ylim is supplied
  if(is.null(ylim)){
    ylim=c(0,max(intensity)*0.75)
  }
  
  #check if the hist bins are supplied
  if (is.null(binwidth) & is.null(bins)){
    binwidth=data[n]/100
    bins=100
  }
  else if (is.null(bins)){
    bins   <- round(data[n]/binwidth,0)
  }
  else{
    binwidth=data[n]/bins
  }
  
  
  #initialize plot
  plot(NA, ylim = ylim, xlim = c(lowerbound,upperbound),xaxt=xact,bty="n",
       main=main,xlab=xlab,ylab=ylab)  
  grid(NA,NULL)#add grid for y-axis
  
  
  
  #plot histogram if requested
  if ("hist" %in% show){
    histdata <- rep(data,x$N)
    end      <- histdata[length(histdata)]
    i        <- 0:(end/binwidth)
    
    height   <- hist(histdata, breaks=c(i,bins+1)*binwidth,plot=FALSE)$counts/binwidth
    rect(xleft=i*binwidth,ybottom=0,xright= (i+1)*binwidth, ytop=height,col="grey70",border="grey40")
    legend=c(legend,"Binned Data")
    legend.col=c(legend.col,"grey70")
  }
  
  
  #plot individual intensities if requested
  if ("points" %in% show){ 
    segments(x0=c(0,data[1:(length(data)-1)]), y0=intensity[1:length(intensity)], 
             x1=data[1:length(data)], y1=intensity[1:length(intensity)],
             col = "cornflowerblue", lwd = 3)
    legend=c(legend,"Data Points")
    legend.col=c(legend.col,"cornflowerblue")
  }
  
  
  #plot blocks if requested
  if ("blocks" %in% show){
    #plot constant blocks
    indices <- which(x$type=="constant")
    if (length(indices)>0){
      b=sapply( x$params[indices], "[[" , "b" )
      segments(x0=data[x$left[indices]], y0=b, 
               x1=data[x$right[indices]], y1=b,
               col = "red", lwd = 3)
    }
    #plot linear blocks
    indices <- which(x$type=="linear")
    if(length(indices)>0){
      a=sapply( x$params[indices], "[[" , "a" )
      b=sapply( x$params[indices], "[[" , "b" )
      lefts=data[x$left[indices]]
      rights=data[x$right[indices]]
      segments(x0=lefts, y0=b, 
               x1=rights, y1=a*(rights-lefts)+b,
               col = "red", lwd = 3)
    }
    
    indices <- which(x$type=="exponential")
    if(length(indices)>0){
      a=sapply( x$params[indices], "[[" , "a" )
      b=sapply( x$params[indices], "[[" , "b" )
      c=sapply( x$params[indices], "[[" , "c" )
      dir=sapply( x$params[indices], "[[" , "direction" )
      lefts=data[x$left[indices]]
      rights=data[x$right[indices]]
      
      for (i in 1:length(a)){
        xv=seq(lefts[i],rights[i],length.out=100)
        if(dir[i]=="forward"){
          lines(xv,a[i]*exp(1)^(b[i]*(xv-lefts[i]))+c[i],col="red",lwd=3)
        }else{
          lines(rev(xv),a[i]*exp(1)^(b[i]*(rights[i]-rev(xv)))+c[i],col="red",lwd=3)
        }
      }
    }
    
    indices <- which(x$type=="power")
    if(length(indices)>0){
      a=sapply( x$params[indices], "[[" , "a" )
      b=sapply( x$params[indices], "[[" , "b" )
      c=sapply( x$params[indices], "[[" , "c" )
      dir=sapply( x$params[indices], "[[" , "direction" )
      lefts=data[x$left[indices]]
      rights=data[x$right[indices]]
      
      for (i in 1:length(a)){
        xv=seq(lefts[i],rights[i],length.out=100)
        if(dir[i]=="forward"){
          lines(xv,a[i]*(xv-lefts[i])^b[i]+c[i],col="red",lwd=3)
        }else{
          lines(rev(xv),a[i]*(rights[i]-rev(xv))^b[i]+c[i],col="red",lwd=3)
        }
      }
      
    }
    #plot FRED blocks
    indices <- which(x$type=="FRED")
    if(length(indices)>0){
      A=sapply( x$params[indices], "[[" , "A" )
      lambda=sapply( x$params[indices], "[[" , "lambda" )
      tau1=sapply( x$params[indices], "[[" , "tau1" )
      tau2=sapply( x$params[indices], "[[" , "tau2" )
      b=sapply( x$params[indices], "[[" , "b" )
      lefts=data[x$left[indices]]
      rights=data[x$right[indices]]
      for (i in 1:length(indices)){
        xv=seq(lefts[i],rights[i],length.out=100)
        tb=xv-lefts[i]
        lines(xv,a[i]*lambda[i]*exp(1)^(-tau1/tb-tb/tau2)+b[i],col="red",lwd=3)
      }
    }
    
    
    
    legend=c(legend,"Blocks")
    legend.col=c(legend.col,"red")
    
  }
  
  #plot legend
  legend(x=x$data[1],y=ylim[2]*0.98,legend= legend, lty=rep(1,length(legend)), 
         lwd=rep(3,length(legend)), col=legend.col) 
}
