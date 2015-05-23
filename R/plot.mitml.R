plot.mitml <- function(x, print=c("beta","psi","sigma"), group="all", export=c("none","png","pdf"), dev.args=list(), ...){
# plot method for objects of class "mitml"

  vrs <- x$model
  yvrs <- 1:length(vrs$yvrs); names(yvrs) <- vrs$yvrs
  pvrs <- 1:length(vrs$pvrs); names(pvrs) <- vrs$pvrs
  qvrs <- 1:length(vrs$qvrs); names(qvrs) <- vrs$qvrs

  print <- match.arg(print,several.ok=TRUE)
  export <- match.arg(export)

  grp.labels <- unique(attr(x$data,"group"))
  grp.imp <- length(grp.labels)
  if(class(group)=="numeric") grp.labels <- grp.labels[group]
  grp <- length(grp.labels)

  if(export!="none"){
    wd <- getwd()
    out <- file.path(wd,"panPlots")
    if(!file.exists(out)) dir.create(out)
  }else{
    do.call(dev.new,dev.args)
    devAskNewPage(ask=FALSE)
  }
  oldpar <- par(no.readonly=TRUE)

  for(gg in 1:grp){
  if(grp.imp>1){
    glab <- paste(",Group:",grp.labels[gg],sep="")
    gfile <- paste("Group-",grp.labels[gg],"_",sep="")
  }else{
    glab <- gfile <- ""
  }

  # plots for fixed regression coefficients
  if("beta" %in% print){
  for(ii in yvrs){
    for(jj in pvrs){
    
    if(export!="none"){
      filename <- paste("BETA_",gfile,names(yvrs[ii]),"_ON_",names(pvrs[jj]),".",export,sep="")
      filename <- gsub("[()]","",filename)
      out.args <- c(list(file=file.path(out,filename)),dev.args)
      do.call(export, out.args)
    }
      
    layout(matrix(c(1,2,3,3),2,2), c(8,1), c(1.13,1))

    par(mar=c(3,3,2,0)+0.5, mgp=c(2,1,0), font.lab=2)
    trc <- x$par.burnin[["beta"]][jj,ii,,gg]
    plot(trc, type="l", ylab="Trace", xlab="Burn-in iteration", 
         main=paste("Beta [",ii,",",jj,glab,"]: ",names(yvrs[ii])," ON ",names(pvrs[jj]),sep=""),
         ...)

    par(mar=c(3,3,1,0)+0.5)
    drw <- x$par.imputation[["beta"]][jj,ii,,gg]
    acf(drw, lag.max=x$iter[["iter"]], ylim=c(-.07,1), yaxp=c(0,1,2), main=NULL, ylab="ACF", ci=0, ...)
    abline(h=0.1)

    par(mar=c(3,0,2,0)+0.5, mgp=c(2,0,0))
    ddrw <- density(drw)
    plot(x=ddrw$y, y=ddrw$x, type="l", xaxt="n", yaxt="n", xlab="", ylab="", ...)

    if(export!="none"){
      dev.off()
    }else{
      devAskNewPage(ask=TRUE)
    }}
  }}

  # plots for random effects
  if("psi" %in% print){
  bvec <- t(expand.grid(yvrs, qvrs))
  for(ii in 1:(length(yvrs)*length(qvrs))){
    for(jj in ii:(length(yvrs)*length(qvrs))){

      if(export!="none"){
        filename <- paste("PSI_",gfile,names(yvrs[bvec[1,ii]]),"_ON_",names(qvrs[bvec[2,ii]]),
          "_WITH_",names(yvrs[bvec[1,jj]]),"_ON_",names(qvrs[bvec[2,jj]]),".",export,sep="")
        filename <- gsub("[()]","",filename)
        out.args <- c(list(file=file.path(out,filename)),dev.args)
        do.call(export, out.args)
      }

      layout(matrix(c(1,2,3,3),2,2), c(8,1), c(1.13,1))

      par(mar=c(3,3,2,0)+0.5, mgp=c(2,1,0), font.lab=2)
      trc <- x$par.burnin[["psi"]][ii,jj,,gg]
      plot(trc, type="l", ylab="Trace", xlab="Burn-in iteration", 
           main=paste("Psi [",ii,",",jj,glab,"]: (",names(yvrs[bvec[1,ii]])," ON ",names(qvrs[bvec[2,ii]]),
                      ") WITH (",names(yvrs[bvec[1,jj]])," ON ",names(qvrs[bvec[2,jj]]),")",sep=""), 
           ...)

      par(mar=c(3,3,1,0)+0.5)
      drw <- x$par.imputation[["psi"]][ii,jj,,gg]
      acf(drw, lag.max=x$iter[["iter"]], ylim=c(-.07,1), yaxp=c(0,1,2), main=NULL, ylab="ACF", ci=0, ...)
      abline(h=0.1)

      par(mar=c(3,0,2,0)+0.5, mgp=c(2,0,0))
      ddrw <- density(drw)
      plot(x=ddrw$y, y=ddrw$x, type="l", xaxt="n", yaxt="n", xlab="", ylab="", ...)

    if(export!="none"){
      dev.off()
    }else{
      devAskNewPage(ask=TRUE)
    }}
  }}

  # plots for residuals
  if("sigma" %in% print){
  for(ii in 1:length(yvrs)){
    for(jj in ii:length(yvrs)){

      if(export!="none"){
        filename <- paste("SIGMA_",gfile,names(yvrs[ii]),"_WITH_",names(yvrs[jj]),".",export,sep="")
        filename <- gsub("[()]","",filename)
        out.args <- c(list(file=file.path(out,filename)),dev.args)
        do.call(export, out.args)
      }

      layout(matrix(c(1,2,3,3),2,2), c(8,1), c(1.13,1))

      par(mar=c(3,3,2,0)+0.5, mgp=c(2,1,0), font.lab=2)
      trc <- x$par.burnin[["sigma"]][ii,jj,,gg]
      plot(trc, type="l", ylab="Trace", xlab="Burn-in iteration", 
           main=paste("Sigma [",ii,",",jj,glab,"]: ",names(yvrs[ii])," WITH ",names(yvrs[jj]),sep=""), 
           ...)

      par(mar=c(3,3,1,0)+0.5)
      drw <- x$par.imputation[["sigma"]][ii,jj,,gg]
      acf(drw, lag.max=x$iter[["iter"]], ylim=c(-.07,1), yaxp=c(0,1,2), main=NULL, ylab="ACF", ci=0, ...)
      abline(h=0.1)

      par(mar=c(3,0,2,0)+0.5, mgp=c(2,0,0))
      ddrw <- density(drw)
      plot(x=ddrw$y, y=ddrw$x, type="l", xaxt="n", yaxt="n", xlab="", ylab="", ...)

    if(export!="none"){
      dev.off()
    }else{
      devAskNewPage(ask=TRUE)
    }}
  }}

  }

  par(oldpar)
  if(export=="none"){
    plot.new()
    devAskNewPage(ask=FALSE)
  }
  dev.off()
  invisible()
}

