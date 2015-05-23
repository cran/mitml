mitmlComplete <- function(x, print=0, force.list=FALSE){

  if(sum(print<=0)>1) stop("Only one negative or zero value is allowed in 'print'.")

  dat <- x$data
  srt <- attr(x$data,"sort")
  attr(dat,"sort") <- NULL
  attr(dat,"group") <- NULL

  m <- x$iter$m
  ind <- x$index.mat
  rpm <- x$replacement.mat

  if(class(print)%in%c("integer","numeric")){

    if(length(print)==1){
      com <- dat
      if(print>0){ com[ind] <- rpm[,print] }
      out <- com[srt,]
      if(force.list) out <- list(out)
    }else{
      out <- list()
      for(ii in print){
        com <- dat
        if(ii>0){ com[ind] <- rpm[,ii] }
        out <- c(out,list(com[srt,]))
      }
    }

  }else{

    if(!print%in%c("list","all")) stop("Invalid 'print' argument.")
    out <- list()
    for(ii in 1:m){
      com <- dat
      com[ind] <- rpm[,ii]
      out <- c(out,list(com[srt,]))
    }

  }
  
  if(class(out)=="list") class(out) <- c("mitml.list","list")
  out

}

