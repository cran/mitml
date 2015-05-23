print.mitml <- function(x,...){
# print method for objects of class "mitml"

  cl <- x$call
  vrs <-x$model 
  itr <- x$iter
  ngr <- length(unique(attr(x$data,"group")))
  cat("\nCall:\n", paste(deparse(cl)), sep="\n")

  cat("\nCluster:\t\t\t", vrs$clus, sep=" ", collapse="\n")
  cat("Target variables:\t\t", vrs$yvrs, collapse="\n")
  cat("Fixed effect predictors:\t", vrs$pvrs, collapse="\n")
  cat("Random effect predictors:\t", vrs$qvrs, collapse="\n")

  cat("\nPerformed", itr$burn, "burn-in iterations, and generated", itr$m, 
      "imputed data sets,\neach", itr$iter,"iterations apart.",
      if(ngr>1){c("\nImputations were carried out seperately within", ngr, "groups.\n")},"\n")

  invisible()
}

