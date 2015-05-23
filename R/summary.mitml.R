summary.mitml <- function(object,n.Rhat=3,...){
# summary method for objects of class "mitml"

  inc <- object$data
  ngr <- length(unique(attr(object$data,"group")))
  prm <- object$par.imputation

  # percent missing
  mdr <- sapply(inc, FUN=function(x){mean(is.na(x))})
  mdr[] <- sprintf(mdr*100,fmt="%.1f")
  mdr <- gsub("^0.0$","0",mdr)

  # potential scale reduction
  iter <- dim(prm[[1]])[3]
  mod <- iter %% n.Rhat
  n <- rep( (iter-mod) / n.Rhat , n.Rhat)
  nmat <- matrix(c(cumsum(n) - n + 1, cumsum(n)), nrow=n.Rhat)
  m <- object$iter$m
  n <- n[1]

  Rhat.list <- list(beta=NULL,psi=NULL,sigma=NULL)
  for(pp in c("beta","psi","sigma")){
    ni <- dim(prm[[pp]])[1]
    nj <- dim(prm[[pp]])[2]
    nl <- dim(prm[[pp]])[4]
    Rhat.mat <- matrix(NA,0,4)

    for(ii in 1:ni){
    for(jj in 1:nj){
    for(ll in 1:nl){

      chains <- apply(nmat, 1, function(x) prm[[pp]][ii,jj,x[1]:x[2],ll])

      # values per chain
      mns <- apply(chains,2,mean)
      vrs <- apply(chains,2,var)
      Bdivn <- var(mns)
      W <- mean(vrs)
      muhat <- mean(chains)
      sighat2 <- (n-1)/n * W + Bdivn
      # sampling distribution
      Vhat <- sighat2 + Bdivn/m
      var.Vhat <- ((n-1)/n)^2*(1/m)*var(vrs) + ((m+1)/(m*n))^2*2/(m-1)*(Bdivn*n)^2 +
                  2*((m+1)*(n-1)/(m*n^2)) * (n/m)*(cov(vrs,mns^2)-2*muhat*cov(vrs,mns))
      df <- 2*Vhat^2 / var.Vhat
      # compute Rhat
      Rhat <- sqrt( (Vhat/W)*df/(df-2) )
      Rhat.mat <- rbind(Rhat.mat, c(ii,jj,ll,Rhat))

    }}}
    colnames(Rhat.mat) <- c("i1","i2","grp","Rhat")
    Rhat.list[[pp]] <- Rhat.mat
  }

  smr <- list(
    call=object$call,
    model=object$model,
    prior=object$prior,
    iter=object$iter,
    ngr=ngr,
    missing.rates=mdr,
    Rhat=Rhat.list
  )
  
  class(smr) <- "mitml.summary"
  smr
}

