clusterMeans <- function(x, cluster, adj=FALSE, group=NULL){
# calculate cluster means

  # get objects if names are given
  isname <- c(length(x)==1, length(cluster)==1, length(group)==1) &
    c(is.character(x), is.character(cluster), is.character(group))
  if(any(isname)){
    parent <- parent.frame()
    if(isname[1]) x <- eval(parse(text=x),parent)
    if(isname[2]) cluster <- eval(parse(text=cluster),parent)
    if(isname[3]) group <- eval(parse(text=group),parent)
  }

  if(!is.null(group)) {
    if(!is.factor(group)) group <- factor(group)
    glv <- length(unique(group))+1
    cluster <- as.numeric(cluster) + as.numeric(group)/glv
  }else{
    cluster <- as.numeric(cluster)
  }
  cluster <- as.numeric(as.factor(cluster))
  n.obs <- rowsum(as.integer(!is.na(x)), cluster)
  gm <- rowsum(x, cluster, na.rm = T)/n.obs
  gm[is.nan(gm)] <- NA
  gm <- gm[cluster]
  if(adj){
    n.obs <- n.obs[cluster]
    ((n.obs * gm) - x)/(n.obs - 1)
  }else{
    gm
  }
}
