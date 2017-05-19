
#S3 Generic summary function

summary.BB <- function(x) {
  cat(length(x$data), 'data points\n\n',
      length(x$left),'total blocks:\n\t',
      sum(x$type=='constant'),'constant blocks\n\t',
      sum(x$type=='linear'),'linear blocks\n\t',  
      sum(x$type=='exponential'),'exponential blocks\n\t',
      sum(x$type=='power'),'power function blocks\n\t',
      sum(x$type=='FRED'),'FRED blocks\n\n',
      floor((length(x$pruned)*100)/length(x$data)),'% of the points were pruned\n\n')
}
