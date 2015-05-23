print.mitml.summary <- function(x,...){
# print method for objects of class "summary.mitml"

  cl <- x$call
  vrs <- x$model
  itr <- x$iter
  ngr <- x$ngr
  mdr <- x$missing.rates
  Rhat <- x$Rhat
  cat("\nCall:\n", paste(deparse(cl)), sep="\n")

  cat("\nCluster:\t\t\t", vrs$clus, sep=" ", collapse="\n")
  cat("Target variables:\t\t", vrs$yvrs, collapse="\n")
  cat("Fixed effect predictors:\t", vrs$pvrs, collapse="\n")
  cat("Random effect predictors:\t", vrs$qvrs, collapse="\n")

  cat("\nPerformed", itr$burn, "burn-in iterations, and generated", itr$m, 
      "imputed data sets,\neach", itr$iter,"iterations apart.",
      if(ngr>1){c("\nImputations were carried out seperately within", ngr, "groups.")},"\n")

  Rhatout <- matrix(c( sapply(Rhat, function(z) min(z[,"Rhat"])),
             sapply(Rhat, function(z) quantile(z[,"Rhat"],.25)),
             sapply(Rhat, function(z) mean(z[,"Rhat"])),
             sapply(Rhat, function(z) median(z[,"Rhat"])),
             sapply(Rhat, function(z) quantile(z[,"Rhat"],.75)),
             sapply(Rhat, function(z) max(z[,"Rhat"])) ), ncol=6 )
  rownames(Rhatout) <- names(Rhat)
  colnames(Rhatout) <- c("Min","25%","Mean","Median","75%","Max")
  cat("\nPotential scale reduction after imputation (Rhat):\n")
  print.table(round(Rhatout,5))

  mdrout <- t(as.matrix(mdr))
  rownames(mdrout) <- "MD%"
  cat("\nPercent missing data per variable:\n")
  print.table(mdrout)

  cat("\n")

  invisible()
}

