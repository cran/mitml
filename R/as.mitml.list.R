as.mitml.list <- function(x){
# adds a class attribute "mitml.list" to its argument

  if(!is.list(x)) stop("Argument must be a 'list'.")

  if(any(!sapply(x,is.data.frame))) x <- lapply(x,as.data.frame)
  class(x) <- c(class(x),"mitml.list")
  x
    
}
