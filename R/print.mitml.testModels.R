print.mitml.testModels <- function(x,...){
# print method for MI estimates

  cl <- x$call
  test <- x$test
  mth <- x$method
  use <- x$use
  m <- x$m
  adj <- x$adj.df
  dfc <- x$df.com

  cat("\nCall:\n", paste(deparse(cl)), sep="\n")

  cat("\nModel comparison calculated from",m,"imputed data sets.")
  cat("\nCombination method:",mth,
    if(mth=="D2"){paste("(",use,")",sep="")},"\n")

  cat("\n")
  out <- sprintf(c("%.3f","%.0f","%.3f","%.3f","%.3f"),test)
  w <- max(sapply(c(out,colnames(test)),nchar))
  cat("  ",format(colnames(test),justify="right",width=w),"\n")
  cat("  ",format(out,justify="right",width=w),"\n")

  if(mth=="D1"){
  cat(if(adj){c("\nHypothesis test adjusted for small samples with",
              paste("df=[",paste(dfc,collapse=","),"]\ncomplete-data degrees of freedom.",sep=""))
      }else{"\nUnadjusted hypothesis test as appropriate in larger samples."},"\n")
  }

  cat("\n")

  invisible()
}

