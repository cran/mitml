.onAttach <- function(libname, pkgname) {
    packageStartupMessage("*** This is beta software. Please report any bugs!")
}

# moving window average for time series
.movingAverage <- function(x, B, fill=TRUE){

    x1 <- cumsum(x)
    N <- length(x)
    y <- rep(NA,N)
    i <- seq(B+1 , N-B)
    xdiff <- x1[ -seq(1,B) ] - x1[ -seq(N-B+1,N) ] 
    xdiff <- xdiff[ - seq(1,B) ] 

    y[i]  <- ( x1[i] + xdiff - c(0,x1[ -seq(N-2*B,N) ]) ) / (2*B+1)

  # fill NAs at beginning and end of time series
  if(fill){
    j <- seq(0,B-1)
    ybeg <- sapply(j, function(z) sum( x[ seq(1,(2*z+1)) ]) / (2*z+1) )
    yend <- sapply(rev(j), function(z) sum( x[ seq(N-2*z,N) ] ) / (2*z+1) )
    y[j+1] <- ybeg
    y[rev(N-j)] <- yend
  }

  y
}

# Gelman-Rubin (1992) criterion for convergence (Rhat)
.GelmanRubin <- function(x,m){

  iter <- ncol(x)
  mod <- iter %% m
  n <- rep( (iter-mod)/m , m )
  nmat <- matrix(c(cumsum(n)-n+1, cumsum(n)), nrow=m)
  n <- n[1]

  Rhat <- numeric(nrow(x))
  for(ii in 1:nrow(x)){

    # values per chain
    chs <- apply(nmat, 1, function(j) x[ii,j[1]:j[2]])
    mns <- apply(chs,2,mean)
    vrs <- apply(chs,2,var)
    Bdivn <- var(mns)
    W <- mean(vrs)
    muhat <- mean(chs)
    sighat2 <- (n-1)/n * W + Bdivn
    # sampling distribution
    Vhat <- sighat2 + Bdivn/m
    var.Vhat <- ((n-1)/n)^2*(1/m)*var(vrs) + ((m+1)/(m*n))^2*2/(m-1)*(Bdivn*n)^2 +
                2*((m+1)*(n-1)/(m*n^2)) * (n/m)*(cov(vrs,mns^2)-2*muhat*cov(vrs,mns))
    df <- 2*Vhat^2 / var.Vhat
    # compute Rhat
    Rhat[ii] <- sqrt( (Vhat/W)*df/(df-2) )

  }
  Rhat

}

# criterion for goodness of approximation (Hoff, 2009)
.SDprop <- function(x){

  np <- nrow(x)
  sp0 <- numeric(np)
  for(i in 1:np){
    arp <- ar(x[i,], aic=TRUE)
    sp0[i] <- arp$var.pred/(1 - sum(arp$ar))^2   # spectral density at frequency 0
  }
  v <- apply(x, 1, var)   # variance of chain
  mcmc.v <- sp0/ncol(x)   # mcmc-variance (correcting for autocorrelation)
  sqrt(mcmc.v / v)   # proportion of variance due to sampling inefficiency

}

