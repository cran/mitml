panImpute <- function(data, type, formula, n.burn=5000, n.iter=100, m=10, group=NULL, prior=NULL, seed=NULL, save.pred=FALSE){
# wrapper function for the Gibbs sampler in the pan package

  # *** checks
  if(!missing(type) & !missing(formula)) stop("Only one of 'type' or 'formula' may be specified.")

  # preserve original order
  if(!is.data.frame(data)) as.data.frame(data)
  data <- cbind(data, original.order=1:nrow(data))
  if(!missing(type)) type <- c(type,0)

  # address additional grouping
  if(!missing(formula)){
    if(is.null(group)){
      group <- rep(1,nrow(data))
    }else{
      group <- data[,group]
      if(length(group)!=nrow(data)) stop("Argument 'group' is not correctly specified.")
    }
  }
  if(!missing(type)){
    if(sum(type==-1)>1) stop("Argument 'group' is not correctly specified,")
    if(sum(type==-1)==0){
      group <- rep(1,nrow(data))
    }else{
      group <- data[,type==-1]
      type[type==-1] <- 0
    }
  }
  group.original <- group
  group <- as.numeric(factor(group,levels=unique(group)))


  # *** prepare input (formula)
  if(!missing(formula)){

    ft <- terms(formula)
    tl <- attr(ft,"term.labels")
    vrs <- attr(ft,"variables")[-1]
  
    # responses
    yvrs <- as.character(vrs)[attr(ft,"response")]
    yvrs <- strsplit(yvrs, split="[[:blank:]]*[[:punct:]][[:blank:]]*")[[1]]
    # reorder: needed for indicator matrix
    yvrs <- yvrs[na.omit(match(colnames(data),yvrs))]
  
    # cluster id
    clt <- tl[grep("\\|",tl)]
    if(length(clt)==0) stop("Cluster indicator not found in formula\n\n",formula,"\n\nPlease specify the cluster indicator and at least one random term using the '|' operator.")
    clt <- strsplit( clt, split="[[:blank:]]*\\|[[:blank:]]*" )[[1]]
    clus <- clname <- clt[2]

    # order data
    data <- data[ order(group,data[,clname]), ]
    group.original <- group.original[ order(group) ]
    group <- group[ order(group) ]

    # predictors: fixed
    pvrs <- c(if(attr(ft,"intercept")){"(Intercept)"}, tl[-grep("\\|",tl)])
    # random
    cl.fml <- as.formula(paste("~",clt[1])) 
    cl.ft <- terms(cl.fml)
    qvrs <- c(if(attr(cl.ft,"intercept")){"(Intercept)"}, attr(cl.ft,"term.labels"))

    # model matrix
    attr(data,"na.action") <- identity
    mm <- model.matrix(formula, data=data)
    mm <- mm[,-grep("\\|",colnames(mm)),drop=FALSE]
    mm <- cbind(mm, data[qvrs[!qvrs%in%c(pvrs,"(Intercept)")]])
    # add random-only intercept
    if("(Intercept)"%in%qvrs & !"(Intercept)"%in%pvrs){
      add.int <- matrix(1,nrow(data),1)
      colnames(add.int) <- "(Intercept)"
      mm <- cbind(add.int,mm)
    }
    all.pred <- colnames(mm)
    pnames <- colnames(data)[colnames(data)%in%all.pred]
    psave <- all.pred[!all.pred%in%c(pnames,"(Intercept)")]
   
    y <- as.matrix(data[yvrs])
    clus <- data[,clus]
    pred <- as.matrix(mm)
    xcol <- which(colnames(mm)%in%pvrs)
    zcol <- which(colnames(mm)%in%qvrs)

  }


  # *** prepare input (type)
  if(!missing(type)){
  
    if(ncol(data)!=length(type)) stop("Length of 'type' must be equal to the number of colums in 'data'.")
    if(sum(type==-2)<1) stop("Cluster indicator not found.")
    if(sum(type==-2)>1) stop("Only one cluster indicator may be specified.")
    if(save.pred) warning("Option 'save.pred' is ignored if 'type' is specified")
    save.pred=FALSE 
 
    data <- data[ order(group,data[,type==-2]), ]
    group.original <- group.original[ order(group) ]
    group <- group[ order(group) ]

    clname <- colnames(data)[type==-2]
    clus <- data[,clname]
    yvrs <- colnames(data)[type==1]
    y <- as.matrix(data[yvrs])
  
    pred <- cbind(1,as.matrix(data[type%in%c(2,3)]))
    pnames <- colnames(data)[type%in%c(2,3)]
    pvrs <- c("(Intercept)",pnames)
    qvrs <- c("(Intercept)",colnames(data)[type==3])
    colnames(pred) <- pvrs
    
    xcol <- 1:length(pvrs)
    zcol <- xcol[pvrs%in%qvrs]

  }

  # * * * * * * * * * * * * * * * * * * * *

  if(sum(is.na(y))==0) stop("Target variables do not contain any missing data.")

  if(is.null(prior)){
    prior <- list( a=ncol(y), Binv=diag(1,ncol(y)),
      c=ncol(y)*length(zcol), Dinv=diag(1,ncol(y)*length(zcol)) )
  }

  if(is.null(seed)){
    set.seed(as.integer(runif(1,0,10^6)))
  }else{
    set.seed(as.integer(seed))
  }
  rns <- sapply(unique(group), function(x,m) as.integer(runif(m+1,0,10^6)), m=m)

  # prepare output
  ind <- which(is.na(data), arr.ind=TRUE, useNames=FALSE)
  ind <- ind[ ind[,2] %in% which(colnames(data)%in%colnames(y)), ]
  rpm <- matrix(NA, nrow(ind), m)

  ng <- length(unique(group))
  np <- length(xcol)
  nq <- length(zcol)
  nr <- ncol(y)
  bpar <- list(beta=array( NA, c(np,nr,n.burn,ng) ),
               psi=array( NA, c(nr*nq,nr*nq,n.burn,ng) ),
               sigma=array( NA, c(nr,nr,n.burn,ng) ))
  ipar <- list(beta=array( NA, c(np,nr,n.iter*m,ng) ),
               psi=array( NA, c(nr*nq,nr*nq,n.iter*m,ng) ),
               sigma=array( NA, c(nr,nr,n.iter*m,ng) ))
  
  # burn-in
  cat("Runnin burn-in phase ...\n")
  flush.console()
  glast <- as.list(unique(group))
  for(gg in unique(group)){

    gi <- group==gg
    gy <- y[gi,]
    gpred <- pred[gi,]
    gclus <- clus[gi]

    cur <- pan::pan(gy, subj=gclus, gpred, xcol, zcol, prior, seed=rns[1,gg], iter=n.burn)
    glast[[gg]] <- cur$last

    bpar[["beta"]][,,,gg] <- cur$beta
    bpar[["psi"]][,,,gg] <- cur$psi
    bpar[["sigma"]][,,,gg] <- cur$sigma

  }
  
  # imputation
  for(ii in 1:m){
    cat("Creating imputed data set (",ii,"/",m,") ...\n")
    flush.console()

    gy.imp <- as.list(unique(group))
    for(gg in unique(group)){

      gi <- group==gg
      gy <- y[gi,]
      gpred <- pred[gi,]
      gclus <- clus[gi]
  
      cur <- pan::pan(gy, subj=gclus, gpred, xcol, zcol, prior, seed=rns[ii+1,gg], iter=n.iter, 
        start=glast[[gg]])
      glast[[gg]] <- cur$last
  
      # populate output
      gy.imp[[gg]] <- cur$y
      ipar[["beta"]][,,(n.iter*(ii-1)+1):(n.iter*ii),gg] <- cur$beta
      ipar[["psi"]][,,(n.iter*(ii-1)+1):(n.iter*ii),gg] <- cur$psi
      ipar[["sigma"]][,,(n.iter*(ii-1)+1):(n.iter*ii),gg] <- cur$sigma

    }
    y.imp <- do.call(rbind,gy.imp)
    rpm[,ii] <- y.imp[is.na(y)]

  }
  cat("Done!\n")

  # clean up
  srt <- data[,ncol(data)]
  data=data[,-ncol(data)]

  # prepare output data
  if(save.pred) data <- cbind(data,mm[psave])
  attr(data,"sort") <- srt
  attr(data,"group") <- group.original

  out <- list(
    data=data,
    replacement.mat=rpm,
    index.mat=ind,
    call=match.call(),
    model=list(clus=clname, yvrs=yvrs, pvrs=pvrs, qvrs=qvrs),
    prior=prior,
    iter=list(burn=n.burn, iter=n.iter, m=m),
    par.burnin=bpar,
    par.imputation=ipar
  )
  class(out) <- "mitml"
  out
  
}

