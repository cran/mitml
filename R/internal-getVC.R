# ***
# Functions to extract variance components from supported classes
# of statistical models
#

# *** lmer method
.getVC.lmer <- function(model){

  if(!requireNamespace("lme4", quietly=TRUE)) stop("The 'lme4' package must be installed in order to use this function.")
  m <- length(model)
  vlist <- addp <- NULL

  # variance components
  vc <- lapply(model,lme4::VarCorr)
  clus <- names(vc[[1]])
  vlist <- list()
  for(vv in clus){
    q <- dim(vc[[1]][[vv]])[1]
    v.cl <- vapply(vc, function(z) z[[vv]], FUN.VALUE=matrix(0,q,q))
    if(is.null(dim(v.cl))) dim(v.cl) <- c(1,1,m)
    dimnames(v.cl)[1:2] <- lapply(dimnames(vc[[1]][[vv]]), function(z) sub("^[(]Intercept[)]$","Intercept",z))
    vlist[[paste("|",vv,sep="")]] <- v.cl
  }
  rv <- sapply(vc, function(z) attr(z,"sc")^2)
  dim(rv) <- c(1,1,m)
  dimnames(rv) <- list("Residual","Residual",NULL)
  vlist <- c(vlist, list(rv))

  # additional parameters
  for(vv in clus){
    if("(Intercept)"%in%colnames(vc[[1]][[vv]])){
      iv <- sapply(vc, function(z) z[[vv]]["(Intercept)","(Intercept)"])
      icc <- iv / (iv + rv[1,1,])
      addp <- c(addp, mean(icc))
      names(addp) <- paste("ICC|",vv,sep="")
    }
  }

  out <- list(vlist=vlist,addp=addp)
  out
}

# *** nlme method
.getVC.nlme <- function(model){

  if(!requireNamespace("nlme", quietly=TRUE)) stop("The 'nlme' package must be installed in order to use this function.")
  m <- length(model)
  vlist <- addp <- NULL

  # variance components (single clustering)
  vc <- lapply(model,nlme::getVarCov)
  clus <- attr(vc[[1]],"group.levels")
  q <- dim(vc[[1]])[1]
  v.cl <- vapply(vc, identity, FUN.VALUE=matrix(0,q,q))
  if(is.null(dim(v.cl))) dim(v.cl) <- c(1,1,m)
  dimnames(v.cl)[1:2] <- lapply(dimnames(vc[[1]]), function(z) sub("^[(]Intercept[)]$","Intercept",z))
  vlist <- list()
  vlist[[paste("|",clus,sep="")]] <- v.cl
  # residual variance
  rv <- sapply(model, function(z) z$sigma^2)
  dim(rv) <- c(1,1,m)
  dimnames(rv) <- list("Residual","Residual",NULL)
  vlist <- c(vlist, list(rv))

  # additional parameters (pre-processed)
  if("(Intercept)"%in%colnames(vc[[1]])){
    iv <- sapply(vc, function(z) z["(Intercept)","(Intercept)"])
    icc <- iv / (iv + rv[1,1,])
    addp <- c(addp, mean(icc))
    names(addp) <- paste("ICC|",clus,sep="")
  }

  out <- list(vlist=vlist,addp=addp)
  out
}

# *** lm method
.getVC.lm <- function(model){

  m <- length(model)
  vlist <- addp <- NULL

  rv <- sapply(model, function(z) summary(z)$sigma^2 )
  dim(rv) <- c(1,1,m)
  dimnames(rv) <- list("Residual","Residual",NULL)
  vlist <- c(vlist, list(rv))

  out <- list(vlist=vlist,addp=addp)
  out
}

