e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

getArea <- function(X=X,buff=buff){
  N.session=length(X)
  area=rep(NA,N.session)
  for(a in 1:N.session){
    xlim=c(min(X[[a]][,1]),max(X[[a]][,1]))+c(-buff[[a]],buff[[a]])
    ylim=c(min(X[[a]][,2]),max(X[[a]][,2]))+c(-buff[[a]],buff[[a]])
    area[a]=diff(xlim)*diff(ylim)
  }
  return(area)
}

sim.SMR.multisession <-
  function(N.session=NA,lambda=NA,n.marked=NA,lam0=NA,theta.d=NA,sigma=NA,K=NA,X=X,buff=NA,
           theta.marked=NA,theta.unmarked=NA,K1D=NA,tlocs=NA,marktype="natural",obstype="poisson"){
    if(!marktype%in%c("natural","premarked")){
      stop("marktype must be 'natural' or 'premarked'")
    }
    if(length(lambda)!=N.session)stop("lambda must be of length N.session")
    if(length(n.marked)!=N.session)stop("n.marked must be of length N.session")
    if(length(lam0)!=N.session)stop("lam0 must be of length N.session")
    # if(length(theta.d)!=N.session)stop("theta.d must be of length N.session")
    if(length(sigma)!=N.session)stop("sigma must be of length N.session")
    if(length(K)!=N.session)stop("K must be of length N.session")
    if(length(X)!=N.session)stop("X must be of length N.session")
    if(length(buff)!=N.session)stop("buff must be of length N.session")
    if(!is.matrix(theta.marked))stop("theta.marked must be a matrix")
    if(nrow(theta.marked)!=N.session)stop("theta.marked must have N.session rows")
    if(ncol(theta.marked)!=3)stop("theta.marked must have N.session columns")
    if(!all(rowSums(theta.marked)==1))stop("theta.marked rows must all sum to 1.")
    if(length(theta.unmarked)!=N.session)stop("theta.unmarked must be of length N.session")
    if(obstype=="negbin"){
      if(length(theta.d)!=N.session)stop("theta.d must be of length N.session")
    }else{
      theta.d=rep(theta.d,N.session)
    }
    if(length(theta.unmarked)!=N.session)stop("theta.unmarked must be of length N.session")
    
    #realized N
    N=rpois(N.session,lambda)
    
    library(abind)
    xlim=ylim=matrix(NA,N.session,2)
    s=D=vector("list",N.session)
    J=rep(NA,N.session)
    
    for(g in 1:N.session){
      X[[g]]=as.matrix(X[[g]])
      xlim[g,]=c(min(X[[g]][,1]),max(X[[g]][,1]))+c(-buff[g],buff[g])
      ylim[g,]=c(min(X[[g]][,2]),max(X[[g]][,2]))+c(-buff[g],buff[g])
      s[[g]]<- cbind(runif(N[g], xlim[g,1],xlim[g,2]), runif(N[g],ylim[g,1],ylim[g,2]))
      D[[g]]<- e2dist(s[[g]],X[[g]])
      J[g]=nrow(X[[g]])
    }
    
    #trap operation
    if(!any(is.na(K1D))){
      for(g in 1:N.session){
        if(any(K1D[[g]]>K[g])){
          stop("Some entries in K1D[[g]] are greater than K[g].")
        }
        if(is.null(dim(K1D[[g]]))){
          if(length(K1D[[g]])!=J[g]){
            stop("K1D[[g]] vector must be of length J[g].")
          }
        }
      }
    }else{
      K1D=vector("list",N.session)
      for(g in 1:N.session){
        K1D[[g]]=rep(K[g],J[g])
      }
    }
    
    #simulate sessions one at a time
    data=vector("list",N.session)
    for(g in 1:N.session){
      data[[g]]=sim.SMR(N=N[g],n.marked=n.marked[g],marktype=marktype,
                   theta.marked=theta.marked[g,],theta.unmarked=theta.unmarked[g],
                   lam0=lam0[g],sigma=sigma[g],K=K[g],X=X[[g]],buff=buff[g],tlocs=tlocs[g],
                   obstype=obstype,theta.d=theta.d[g])
    }
    
    #combine session data
    n.samples=rep(NA,N.session)
    for(g in 1:N.session){
      n.samples[g]=length(data[[g]]$this.j)
    }
    n.samples.max=max(n.samples)
    this.j=this.k=samp.type=ID=matrix(NA,N.session,n.samples.max)
    ID.marked=y=s=vector("list",N.session)
    n.M=n.UM=rep(NA,N.session)
    for(g in 1:N.session){
      this.j[g,1:n.samples[g]]=data[[g]]$this.j
      this.k[g,1:n.samples[g]]=data[[g]]$this.k
      samp.type[g,1:n.samples[g]]=data[[g]]$samp.type
      ID[g,1:n.samples[g]]=data[[g]]$ID
      n.M[g]=data[[g]]$n.M
      n.UM[g]=data[[g]]$n.UM
      y[[g]]=data[[g]]$y
      s[[g]]=data[[g]]$s
      ID.marked[[g]]=data[[g]]$ID.marked
    }
    
    if(any(tlocs>0)){
      locs=array(NA,dim=c(N.session,max(n.marked),max(tlocs),2))
      for(g in 1:N.session){
        locs[g,1:n.marked[g],1:tlocs[g],1:2]=data[[g]]$locs
      }
    }else{
      locs=NA
    }
    
    out<-list(this.j=this.j,this.k=this.k,samp.type=samp.type,ID.marked=ID.marked, #observed data
              n.marked=n.marked,locs=locs,n.M=n.M,n.UM=n.UM,
              y=y,s=s, ID=ID,N=N,#true data
              X=X,K=K,K1D=K1D,buff=buff,xlim=xlim,ylim=ylim)
    return(out)
  }