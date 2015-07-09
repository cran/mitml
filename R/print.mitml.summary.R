print.mitml.summary <- function(x,...){
# print method for objects of class "summary.mitml"

  cl <- x$call
  vrs <- x$model
  itr <- x$iter
  ngr <- x$ngr
  mdr <- x$missing.rates
  Rhat <- x$Rhat
  cat("\nCall:\n", paste(deparse(cl)), sep="\n")
  cat("\n")

  cat(formatC("Cluster variable:",width=-25), vrs$clus, sep=" ", collapse="\n")
  cat(formatC("Target variables:",width=-25), vrs$yvrs, collapse="\n")
  cat(formatC("Fixed effect predictors:",width=-25), vrs$pvrs, collapse="\n")
  cat(formatC("Random effect predictors:",width=-25), vrs$qvrs, collapse="\n")

  cat("\nPerformed", sprintf("%.0f",itr$burn), "burn-in iterations, and generated", sprintf("%.0f",itr$m),
      "imputed data sets,\neach", sprintf("%.0f",itr$iter), "iterations apart.",
      if(ngr>1){c("\nImputations were carried out seperately within", sprintf("%.0f",ngr), "groups.")},"\n")

  Rhatout <- matrix(c( sapply(Rhat, function(z) min(z[,"Rhat"])),
             sapply(Rhat, function(z) quantile(z[,"Rhat"],.25)),
             sapply(Rhat, function(z) mean(z[,"Rhat"])),
             sapply(Rhat, function(z) median(z[,"Rhat"])),
             sapply(Rhat, function(z) quantile(z[,"Rhat"],.75)),
             sapply(Rhat, function(z) max(z[,"Rhat"])) ), ncol=6 )
  rownames(Rhatout) <- c("Beta:","Psi:","Sigma:")
  colnames(Rhatout) <- c("Min","25%","Mean","Median","75%","Max")
  cat("\nPotential scale reduction (Rhat, imputation phase):\n")
  print.table(round(Rhatout,3))

  cat("\nLargest potential scale reduction (Rhat):\n")
  maxb <- Rhat$beta[which.max(Rhat$beta[,4]),]
  maxp <- Rhat$psi[which.max(Rhat$psi[,4]),]
  maxs <- Rhat$sigma[which.max(Rhat$sigma[,4]),]
  cat("Beta: [", paste(maxb[1:2],collapse=",") ,"], ",
      "Psi: [", paste(maxp[1:2],collapse=",") ,"], ",
      "Sigma: [", paste(maxs[1:2],collapse=",") ,"]\n", sep="")

  mdrout <- t(as.matrix(mdr))
  rownames(mdrout) <- "MD%"
  cat("\nMissing data per variable:\n")
  print.table(mdrout)

  cat("\n")

  invisible()
}

