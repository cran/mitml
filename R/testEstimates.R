testEstimates <- function(model, var.comp=FALSE, df.com=NULL){
# combine scalar estimates from the analysis of multiply imputed data

  cls <- class(model[[1]])
  coef.method <- vc.method <- "default"

  # identify procedures for model class
  if(cls=="lm") vc.method <- "lm"
  if(length(grep("merMod",cls)) > 0 & coef.method=="default"){
    if(!requireNamespace("lme4", quietly=TRUE)) stop("The 'lme4' package must be installed in order to handle 'merMod' class objects.")
    coef.method <- vc.method <- "lmer"
  }
  if(length(grep("lme",cls)) > 0 & coef.method=="default"){
    if(!requireNamespace("nlme", quietly=TRUE)) stop("The 'nlme' package must be installed in order to handle 'lme' class objects.")
    coef.method <- vc.method <- "nlme"
  }

  # *** fixed coefficients

  fe <- switch(coef.method,
    lmer=.getCOEF.lmer(model,diagonal=TRUE),
    nlme=.getCOEF.nlme(model,diagonal=TRUE),
    default=.getCOEF.default(model,diagonal=TRUE)
  )

  m <- length(model)
  Qhat <- fe$Qhat
  Uhat <- fe$Uhat
  if(is.null(dim(Qhat))){ 
    dim(Qhat) <- c(1,m)
    dim(Uhat) <- c(1,m)
    dimnames(Qhat) <- dimnames(Uhat) <- list(fe$nms, NULL)
  }

  Qbar <- apply(Qhat,1,mean)
  Ubar <- apply(Uhat,1,mean)
  B <- apply(Qhat,1,var)
  T <- Ubar + (1+m^(-1)) * B

  r <- (1+m^(-1))*B/Ubar 
  v <- vm <- (m-1)*(1+r^(-1))^2
  fmi <- (r+2/(v+3))/(r+1)

  se <- sqrt(T)
  t <- Qbar/se

  if(!is.null(df.com)){
    lam <- r/(r+1)
    vobs <- (1-lam)*((df.com+1)/(df.com+3))*df.com
    v <- (vm^(-1)+vobs^(-1))^(-1)
  }

  p <- 1-pt(abs(t),df=v)
  out <- matrix(c(Qbar,se,t,v,p,r,fmi),ncol=7)
  colnames(out) <- c("Estimate","Std.Error","t.value","df","p.value","RIV","FMI")
  rownames(out) <- names(Qbar)

  # *** variance components
  vout <- NULL
  if(var.comp){

    vc <- switch(vc.method,
      lmer=.getVC.lmer(model),
      nlme=.getVC.nlme(model),
      lm=.getVC.lm(model),
      default=warning("Computation of variance components not supported for objects of class '",cls,"' (see ?with.mitml.list for examples to calculate these manually).")
    )
   
    vlist <- vc$vlist
    addp <- vc$addp

    if(!is.null(vlist)){
      vlist <- lapply(vlist, function(z) apply(z,1:2,mean) )
      ln <- names(vlist)
      nms <- vout <- c()
      for(vv in 1:length(vlist)){
        vc <- vlist[[vv]]
        rn <- rownames(vc)
        cn <- colnames(vc)
        for(rr in 1:nrow(vc)){
        for(cc in 1:ncol(vc)){
          if(cc>=rr){
            vout <- c(vout, vc[rr,cc])
            nms <- c(nms, paste(rn[rr],"~~",cn[cc],ln[vv],sep=""))
          }
        }}
      }
    }
    vout <- matrix(vout,ncol=1)
    colnames(vout) <- "Estimate"
    rownames(vout) <- nms
    if(!is.null(addp)) vout <- rbind(vout, as.matrix(addp))

  }
  
  out <- list(
    call=match.call(),
    estimates=out,
    var.comp=vout,
    m=m,
    adj.df=!is.null(df.com),
    df.com=df.com,
    cls.method=coef.method
  )
  class(out) <- "mitml.testEstimates"
  out

}

