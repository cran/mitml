plot.mitml <- function(x, print=c("beta","psi","sigma"), pos=NULL, group="all",
                       trace=c("imputation","burnin","all"), thin=1, smooth=3,
                       n.Rhat=3, export=c("none","png","pdf"), dev.args=list(),
                       ...){

# plot method for objects of class "mitml"

  # .movingAverage <- mitml:::.movingAverage
  # .smoothedACF <- mitml:::.smoothedACF

  # retrieve data and variable names
  vrs <- x$model
  clus <- x$model$clus
  pvrs <- seq_along(attr(vrs,"full.names")$pvrs)
  qvrs <- seq_along(attr(vrs,"full.names")$qvrs)
  names(pvrs) <- attr(vrs,"full.names")$pvrs
  names(qvrs) <- attr(vrs,"full.names")$qvrs

  # check for random L1
  rl1 <- x$random.L1=="full"

  print <- match.arg(print,several.ok=TRUE)
  trace <- match.arg(trace)
  export <- match.arg(export)

  # check print and position for selected parameters
  if(!is.null(pos) & length(print)>1){
    pos <- NULL
    warning("The 'pos' argument may only be used when 'print' is cleary defined as one of 'beta', 'psi', or 'sigma' (see '?plot').")
  }

  # grouping
  grp.labels <- unique(attr(x$data,"group"))
  grp.imp <- length(grp.labels)
  if(class(group)=="numeric") grp.labels <- grp.labels[group]
  grp <- length(grp.labels)

  # export, graphical parameters
  if(export!="none"){
    wd <- getwd()
    out <- file.path(wd,"mitmlPlots")
    if(!file.exists(out)) dir.create(out)
  }else{
    do.call(dev.new,dev.args)
    devAskNewPage(ask=FALSE)
  }
  oldpar <- par(no.readonly=TRUE)

  # ***
  # start plotting
  #

  for(gg in 1:grp){

  # grouping
  if(grp.imp>1){
    glab <- paste(",Group:",grp.labels[gg],sep="")
    gfile <- paste("Group-",grp.labels[gg],"_",sep="")
  }else{
    glab <- gfile <- ""
  }

  # expand yvrs for multiple categories
  yvrs <- vrs$yvrs
  cvrs <- attr(x$data,"cvrs")
  # expand categorical target variables
  if(!is.null(cvrs)){
    yvrs <- c(yvrs[!yvrs%in%cvrs], cvrs)
    for(cc in 1:length(cvrs)){
      cv <- cvrs[cc]
      ci <- which(yvrs==cv)
      yi <- 1:length(yvrs)
      nlev <- attr(x$data,"levels")[gg,cc]
      if(nlev>2){
        newy <- paste0(cv, 1:(nlev-1))
      }else{
        newy <- cv
      }
      sel0 <- yi[yi<ci]
      sel1 <- yi[yi>ci]
      yvrs <- c(yvrs[sel0],newy,yvrs[sel1])
    }
  }
  ynam <- yvrs
  yvrs <- seq_along(yvrs)
  names(yvrs) <- ynam

  # number of iterations
  n <- dim(x$par.burnin[["beta"]])[3]+dim(x$par.imputation[["beta"]])[3]
  nb <- dim(x$par.burnin[["beta"]])[3]
  ni <- dim(x$par.imputation[["beta"]])[3]
  niter <- x$iter[["iter"]]
  # thinned-sample indicators
  s <- seq.int(thin,n,by=thin)
  sb <- seq.int(thin,nb,by=thin)
  si <- seq.int(thin,ni,by=thin)
  lag <- ceiling(niter/thin)

  # *** plots for fixed regression coefficients
  if("beta" %in% print){

  # check if pos is badly defined
  if(!is.null(pos)){
    if(pos[1] > max(pvrs) | pos[1] < min(pvrs) | pos[2] > max(yvrs) | pos[2] < min(yvrs)){
      .restoreDevice(oldpar,export,close=TRUE)
      stop("There is no entry [",pos[1],",",pos[2],"] in 'beta'.")
    }
  }

  for(ii in yvrs){
    for(jj in pvrs){

      # skip if individual parameters requested
      if(!is.null(pos)){
        if(!(pos[1]==jj & pos[2]==ii)) next
      }
    
      if(export!="none"){
        filename <- paste("BETA_",gfile,names(yvrs[ii]),"_ON_",names(pvrs[jj]),".",export,sep="")
        filename <- gsub("[(),]","",filename)
        filename <- gsub("[[:space:]]","-",filename)
        out.args <- c(list(file=file.path(out,filename)),dev.args)
        do.call(export, out.args)
      }
      
      layout(matrix(c(1,2,3,4),2,2), c(5,1), c(1.13,1))

      # choose section of trace
      switch(trace,
        imputation={
          trc <- x$par.imputation[["beta"]][jj,ii,,gg][si]
        },
        burnin={
          trc <- x$par.burnin[["beta"]][jj,ii,,gg][sb]
        },
        all={
          trc <- c(x$par.burnin[["beta"]][jj,ii,,gg],
                   x$par.imputation[["beta"]][jj,ii,,gg])[s]
        }
      )

      # trace plot
      par(mar=c(3,3,2,0)+0.5, mgp=c(2,1,0), font.lab=2)
      ymin <- min(trc)
      ymax <- max(trc)
      yr <- ymax-ymin
      plot(trc, type=ifelse(trace=="all","n","l"), ylab="Trace", xlab="Iteration",
           xaxt="n", ylim=c(ymin-yr*.03, ymax+yr*.03),
           ...)
      if(trace=="all"){
        lines(which(s<=nb), trc[s<=nb], col="grey75")
        lines(which(s>=nb), trc[s>=nb], col="black")
      }
      axt <- axTicks(1)
      title(main=paste("Beta [",jj,",",ii,glab,"]: ",names(yvrs[ii])," ON ",names(pvrs[jj]),sep=""), cex.main=1)
      if(trace=="imputation"){
        axl <- sprintf("%d", thin*(axt+nb))
      }else{
        axl <- sprintf("%d", thin*axt)
      }
      axis(side=1, at=axt, labels=axl)

      # trend line for trace (moving window average)
      if(all(is.numeric(smooth),smooth>0)){
        B <- floor(niter/(smooth*thin))
        mwa <- .movingAverage(trc,B,fill=TRUE)
        lines(mwa, col="grey60")
      }

      # blue line
      if(trace=="all") abline(v=ceiling(nb/thin), col="blue")

      # further plots
      if(trace=="burnin"){
        drw <- x$par.burnin[["beta"]][jj,ii,sb,gg]
      }else{
        drw <- x$par.imputation[["beta"]][jj,ii,si,gg]
      }

      # autocorrelation plot
      par(mar=c(3,3,1,0)+0.5)
      ac <- acf(drw, lag.max=lag+2, plot=F)
      plot(ac[1:lag], ylim=c(-.1,1), yaxt="n", main=NULL, ylab="ACF", ci=0,
           ...)
      axis(side=2, at=c(0,.5,1))
      abline(h=c(-.1,.1), col="blue")

      # kernel density plot
      par(mar=c(3,0,2,0)+0.5, mgp=c(2,0,0))
      ddrw <- density(drw)
      plot(x=ddrw$y, y=ddrw$x, type="l", xaxt="n", yaxt="n", xlab="", ylab="",
           ylim=c(ymin-yr*.03, ymax+yr*.03),
           ...)

      # posterior summary
      par(mar=c(1,-0.5,0,-0.5)+0.5)
      plot.new()
      text(0,0.5,paste("EAP:   ", sprintf(fmt="%.3f", mean(drw)), "\n",
                       "MAP:   ", sprintf(fmt="%.3f", ddrw$x[which.max(ddrw$y)]), "\n",
                       "SD:    ", sprintf(fmt="%.3f", sd(drw)), "\n",
                       "2.5%:  ", sprintf(fmt="%.3f", quantile(drw,.025)), "\n",
                       "97.5%: ", sprintf(fmt="%.3f", quantile(drw,.975)), "\n",
                       "Rhat:  ", sprintf(fmt="%.3f", .GelmanRubin(t(drw),n.Rhat)), "\n",
                       "ACF-k: ", sprintf(fmt="%.3f", .smoothedACF(ac,k=lag,sd=.5)), "\n",
                       sep=""), adj=c(0,.5), cex=.8, family="mono", font=2, ...)
    
      if(export!="none"){
        dev.off()
      }else{
        devAskNewPage(ask=TRUE)
      }
  }}}

  # *** plots for random effects' variance components
  if("psi" %in% print){

  # inde matrix
  bvec <- t(expand.grid(qvrs, yvrs))

  # attempt to fix pos if badly defined
  if(!is.null(pos)){
    pos0 <- pos
    if(pos[2]>pos[1]){   # fix if pos is redundant/transposed
      pos[1] <- pos0[2]
      pos[2] <- pos0[1]
    }
    if(any(pos0 > max(yvrs)) | any(pos0 < min(yvrs))){
      .restoreDevice(oldpar,export,close=TRUE)
      stop("There is no entry [",pos0[1],",",pos0[2],"] in 'psi'.")
    }
    if(!identical(pos,pos0))
      warning("Could not use entry [",pos0[1],",",pos0[2],"] in 'psi'. Used [",pos[1],",",pos[2],"] instead.")
  }

  for(ii in 1:(length(yvrs)*length(qvrs))){
    for(jj in ii:(length(yvrs)*length(qvrs))){

      # skip if individual parameters requested
      if(!is.null(pos)){
        if(!(pos[1]==jj & pos[2]==ii)) next
      }
    
      if(export!="none"){
        filename <- paste("PSI_",gfile,names(yvrs[bvec[2,ii]]),"_ON_",names(qvrs[bvec[1,ii]]),
          "_WITH_",names(yvrs[bvec[2,jj]]),"_ON_",names(qvrs[bvec[1,jj]]),".",export,sep="")
        filename <- gsub("[(),]","",filename)
        filename <- gsub("[[:space:]]","-",filename)
        out.args <- c(list(file=file.path(out,filename)),dev.args)
        do.call(export, out.args)
      }

      layout(matrix(c(1,2,3,4),2,2), c(5,1), c(1.13,1))

      switch(trace,
        imputation={
          trc <- x$par.imputation[["psi"]][jj,ii,,gg][si]
        },
        burnin={
          trc <- x$par.burnin[["psi"]][jj,ii,,gg][sb]
        },
        all={
          trc <- c(x$par.burnin[["psi"]][jj,ii,,gg],
                   x$par.imputation[["psi"]][jj,ii,,gg])[s]
        }
      )

      # trace plot
      par(mar=c(3,3,2,0)+0.5, mgp=c(2,1,0), font.lab=2)
      ymin <- min(trc)
      ymax <- max(trc)
      yr <- ymax-ymin
      plot(trc, type=ifelse(trace=="all","n","l"), ylab="Trace", xlab="Iteration",
           xaxt="n", ylim=c(ymin-yr*.03, ymax+yr*.03),
           ...)
      if(trace=="all"){
        lines(which(s<=nb), trc[s<=nb], col="grey75")
        lines(which(s>=nb), trc[s>=nb], col="black")
      }
      title(main=paste("Psi [",jj,",",ii,glab,"]: (",names(yvrs[bvec[2,ii]])," ON ",names(qvrs[bvec[1,ii]]),
            ") WITH (",names(yvrs[bvec[2,jj]])," ON ",names(qvrs[bvec[1,jj]]),")",sep=""), cex.main=1)
      axt <- axTicks(1)
      if(trace=="imputation"){
        axl <- sprintf("%d", thin*(axt+nb))
      }else{
        axl <- sprintf("%d", thin*axt)
      }
      axis(side=1, at=axt, labels=axl)

      # trend line for trace (moving window average)
      if(all(is.numeric(smooth),smooth>0)){
        B <- floor(niter/(smooth*thin))
        mwa <- .movingAverage(trc,B,fill=TRUE)
        lines(mwa, col="grey60")
      }

      # blue line
      if(trace=="all") abline(v=ceiling(nb/thin), col="blue")

      # further plots
      if(trace=="burnin"){
        drw <- x$par.burnin[["psi"]][jj,ii,sb,gg]
      }else{
        drw <- x$par.imputation[["psi"]][jj,ii,si,gg]
      }

      # autocorrelation plot
      par(mar=c(3,3,1,0)+0.5)
      ac <- acf(drw, lag.max=lag+2, plot=F)
      plot(ac[1:lag], ylim=c(-.1,1), yaxt="n", main=NULL, ylab="ACF", ci=0,
           ...)
      axis(side=2, at=c(0,.5,1))
      abline(h=c(-.1,.1), col="blue")

      # kernel density plot
      par(mar=c(3,0,2,0)+0.5, mgp=c(2,0,0))
      ddrw <- density(drw)
      plot(x=ddrw$y, y=ddrw$x, type="l", xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(ymin-yr*.03, ymax+yr*.03), ...)

      # posterior summary
      par(mar=c(1,-0.5,0,-0.5)+0.5)
      plot.new()
      text(0,0.5,paste("EAP:   ", sprintf(fmt="%.3f", mean(drw)), "\n",
                       "MAP:   ", sprintf(fmt="%.3f", ddrw$x[which.max(ddrw$y)]), "\n",
                       "SD:    ", sprintf(fmt="%.3f", sd(drw)), "\n",
                       "2.5%:  ", sprintf(fmt="%.3f", quantile(drw,.025)), "\n",
                       "97.5%: ", sprintf(fmt="%.3f", quantile(drw,.975)), "\n",
                       "Rhat:  ", sprintf(fmt="%.3f", .GelmanRubin(t(drw),n.Rhat)), "\n",
                       "ACF-k: ", sprintf(fmt="%.3f", .smoothedACF(ac,k=lag,sd=.5)), "\n",
                       sep=""), adj=c(0,.5), cex=.8, family="mono", font=2, ...)

    if(export!="none"){
      dev.off()
    }else{
      devAskNewPage(ask=TRUE)
    }
  }}}

  # *** plots for residual variance components
  if("sigma" %in% print){

  # cluster-specific covariance matrices stacked in rows
  clus2 <- unique(x$data[,clus])
  clus3 <- if(rl1) seq_along(clus2) else 1

  # attempt to fix pos if badly defined
  if(!is.null(pos)){
    pos0 <- pos
    dims <- dim(x$par.imputation$sigma)
    if(pos[2] > length(yvrs)){   # fix if pos is transposed
      pos[1] <- pos0[2]
      pos[2] <- pos0[1]
      pos0 <- pos
    }
    if(pos[2] > ((pos[1]-1)%%length(yvrs))+1){   # fix if pos is redundant
      pos[1] <- pos0[1] - pos0[1]%%length(yvrs) + pos0[2]
      pos[2] <- pos0[1]%%length(yvrs)
    }
    if(all(pos0 > max(yvrs)) | any(pos0 < min(yvrs)) | max(pos0) > dims[1]){
      .restoreDevice(oldpar,export,close=TRUE)
      stop("There is no entry [",pos0[1],",",pos0[2],"] in 'sigma'.")
    }
    if(!identical(pos,pos0))
      warning("Could not use entry [",pos0[1],",",pos0[2],"] in 'sigma'. Used [",pos[1],",",pos[2],"] instead.")
  }

  for(icl in clus3){

  for(ii in 1:length(yvrs)){
    for(jj in ii:length(yvrs)){

      # adjust row index for cluster-specific covariance matrices
      jj2 <- jj+(icl-1)*length(yvrs)

      # skip if individual parameters requested
      if(!is.null(pos)){
        if(!(pos[1]==jj2 & pos[2]==ii)) next
      }
    
      if(export!="none"){
        filename <- paste0("SIGMA_",gfile,names(yvrs[ii]),"_WITH_",names(yvrs[jj]), 
                           if(rl1) paste0(clus,clus2[icl]),".",export)
        filename <- gsub("[(),]","",filename)
        filename <- gsub("[[:space:]]","-",filename)
        out.args <- c(list(file=file.path(out,filename)),dev.args)
        do.call(export, out.args)
      }

      layout(matrix(c(1,2,3,4),2,2), c(5,1), c(1.13,1))

      switch(trace,
        imputation={
          trc <- x$par.imputation[["sigma"]][jj2,ii,,gg][si]
        },
        burnin={
          trc <- x$par.burnin[["sigma"]][jj2,ii,,gg][sb]
        },
        all={
          trc <- c(x$par.burnin[["sigma"]][jj2,ii,,gg],
                   x$par.imputation[["sigma"]][jj2,ii,,gg])[s]
        }
      )

      # trace plots
      par(mar=c(3,3,2,0)+0.5, mgp=c(2,1,0), font.lab=2)
      ymin <- min(trc)
      ymax <- max(trc)
      yr <- ymax-ymin
      plot(trc, type=ifelse(trace=="all","n","l"), ylab="Trace", xlab="Iteration",
           xaxt="n", ylim=c(ymin-yr*.03, ymax+yr*.03),
           ...)
      if(trace=="all"){
        lines(which(s<=nb), trc[s<=nb], col="grey75")
        lines(which(s>=nb), trc[s>=nb], col="black")
      }
      title(main=paste0("Sigma [",jj2,",",ii,glab,"]: ",names(yvrs[ii])," WITH ",names(yvrs[jj]),
                        if(rl1) paste0(" [",clus,":",clus2[icl],"]")), cex.main=1)
      axt <- axTicks(1)
      if(trace=="imputation"){
        axl <- sprintf("%d", thin*(axt+nb))
      }else{
        axl <- sprintf("%d", thin*axt)
      }
      axis(side=1, at=axt, labels=axl)

      # trend line for trace (moving window average)
      if(all(is.numeric(smooth),smooth>0)){
        B <- floor(niter/(smooth*thin))
        mwa <- .movingAverage(trc,B,fill=TRUE)
        lines(mwa, col="grey60")
      }

      # blue line
      if(trace=="all") abline(v=ceiling(nb/thin), col="blue")

      # further plots
      if(trace=="burnin"){
        drw <- x$par.burnin[["sigma"]][jj2,ii,sb,gg]
      }else{
        drw <- x$par.imputation[["sigma"]][jj2,ii,si,gg]
      }

      # autocorrelation plot
      par(mar=c(3,3,1,0)+0.5)
      ac <- acf(drw, lag.max=lag+2, plot=F)
      plot(ac[1:lag], ylim=c(-.1,1), yaxt="n", main=NULL, ylab="ACF", ci=0,
           ...)
      axis(side=2, at=c(0,.5,1))
      abline(h=c(-.1,.1), col="blue")

      # kernel density plot
      par(mar=c(3,0,2,0)+0.5, mgp=c(2,0,0))
      ddrw <- density(drw)
      plot(x=ddrw$y, y=ddrw$x, type="l", xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(ymin-yr*.03, ymax+yr*.03), ...)

      # posterior summary
      par(mar=c(1,-0.5,0,-0.5)+0.5)
      plot.new()
      text(0,0.5,paste("EAP:   ", sprintf(fmt="%.3f", mean(drw)), "\n",
                       "MAP:   ", sprintf(fmt="%.3f", ddrw$x[which.max(ddrw$y)]), "\n",
                       "SD:    ", sprintf(fmt="%.3f", sd(drw)), "\n",
                       "2.5%:  ", sprintf(fmt="%.3f", quantile(drw,.025)), "\n",
                       "97.5%: ", sprintf(fmt="%.3f", quantile(drw,.975)), "\n",
                       "Rhat:  ", sprintf(fmt="%.3f", .GelmanRubin(t(drw),n.Rhat)), "\n",
                       "ACF-k: ", sprintf(fmt="%.3f", .smoothedACF(ac,k=lag,sd=.5)), "\n",
                       sep=""), adj=c(0,.5), cex=.8, family="mono", font=2, ...)

    if(export!="none"){
      dev.off()
    }else{
      devAskNewPage(ask=TRUE)
    }
  }}}}

  }

  plot.new()
  par(oldpar)
  if(export=="none") devAskNewPage(ask=FALSE)
  dev.off()
  invisible()
}

# restore and shut down parameters upon error
.restoreDevice <- function(pars, export, close=TRUE){

  par(pars)
  if(export=="none") devAskNewPage(ask=FALSE)
  if(close) dev.off()
  invisible()

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

# lag-k autocorrelation smoothed by values of a normal density
.smoothedACF <- function(x, k, sd=.5){

  x0 <- x$ac[-1,1,1]
  n <- length(x0)
  add <- n-k
  x0 <- x0[(k-add):n]

  # weights based on normal density
  w <- dnorm(-add:add, 0, sd)
  y <- sum( x0 * (w/sum(w)) )

  y

}

