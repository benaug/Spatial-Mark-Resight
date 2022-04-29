e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.multisession=function(data,inits=NA,M1=NA,M2=NA,marktype="premarked",obstype="poisson"){
  N.session=nrow(data$this.j)
  n.samples=rowSums(!is.na(data$this.j))
  init.session=vector("list",N.session)
  
  #split inits by session
  inits.use=vector("list",N.session)
  parms=names(inits)
  for(g in 1:N.session){
    inits.use[[g]]=vector("list",length(parms))
    names(inits.use[[g]])=parms
    for(i in 1:length(parms)){
      inits.use[[g]][[i]]=inits[[i]][g]
    }
  }
  
  #initialize sessions one by one
  for(g in 1:N.session){
    if(all(is.na(data$locs))){
      locs.use=NA
    }else{
      tlocs.sess.max=max(rowSums(!is.na(data$locs[g,1:M1[g],,1])))
      locs.use=data$locs[g,1:M1[g],1:tlocs.sess.max,1:2]
    }
    data.use=list(this.j=data$this.j[g,1:n.samples[g]],this.k=data$this.k[g,1:n.samples[g]],samp.type=data$samp.type[g,1:n.samples[g]],
                  ID.marked=data$ID.marked[[g]],n.marked=data$n.marked[g],locs=locs.use,X=data$X[[g]],buff=data$buff[g],
                  K1D=data$K1D[[g]])
    init.session[[g]]=init.SMR(data.use,inits.use[[g]],M1=M1[g],M2=M2[g],marktype=marktype,obstype=obstype)
  }
  J=unlist(lapply(data$X,nrow))
  M.both=M1+M2
  maxM.both=max(M.both)
  s=array(NA,dim=c(N.session,maxM.both,2))
  z=matrix(NA,N.session,maxM.both)
  ID=matrix(NA,N.session,max(n.samples))
  y.full=array(NA,dim=c(N.session,maxM.both,max(J)))
  y.event=array(NA,dim=c(N.session,maxM.both,max(J),3))
  K1D=matrix(NA,N.session,max(J))
  n.fixed=rep(NA,N.session)
  samp.type=matrix(NA,N.session,max(n.samples))
  match=array(NA,dim=c(N.session,max(n.samples),maxM.both))
  if(!all(is.na(data$locs))){
    tel.inds=matrix(NA,N.session,dim(data$locs)[2],dim(data$locs)[3])
    n.locs.ind=matrix(NA,N.session,dim(data$locs)[2],dim(data$locs)[3])
    n.tel.inds=rep(NA,N.session)
  }else{
    tel.inds=n.locs.ind=n.tel.inds=NA
  }

  for(g in 1:N.session){
    s[g,1:M.both[g],]=init.session[[g]]$s
    z[g,1:M.both[g]]=init.session[[g]]$z
    ID[g,1:n.samples[g]]=init.session[[g]]$ID
    y.full[g,1:M.both[g],1:J[g]]=init.session[[g]]$y.full
    y.event[g,1:M.both[g],1:J[g],]=init.session[[g]]$y.event
    K1D[g,1:J[g]]=init.session[[g]]$K1D
    samp.type[g,1:n.samples[g]]=init.session[[g]]$samp.type
    match[g,1:n.samples[g],1:M.both[g]]=init.session[[g]]$match
    n.fixed[g]=init.session[[g]]$n.fixed
    if(!all(is.na(data$locs))){
      n.tel.inds[g]=sum(rowSums(!is.na(data$locs[g,,,1]))>0)
      tel.inds[g,1:n.tel.inds[g]]=init.session[[g]]$tel.inds
      n.locs.ind[g,1:n.tel.inds[g]]=init.session[[g]]$n.locs.ind
    }
  }
  #put X in ragged array
  X.new=array(NA,dim=c(N.session,max(J),2))
  for(g in 1:N.session){
    X.new[g,1:J[g],]=data$X[[g]]
  }


  return(list(s=s,z=z,ID=ID,y.full=y.full,y.event=y.event,K1D=K1D,J=J,X=X.new,
         n.samples=n.samples,n.fixed=n.fixed,samp.type=samp.type,this.j=data$this.j,match=match,
         xlim=data$xlim,ylim=data$ylim,locs=locs.use,tel.inds=tel.inds,n.locs.ind=n.locs.ind,n.tel.inds=n.tel.inds))

}