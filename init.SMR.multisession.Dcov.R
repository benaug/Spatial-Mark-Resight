e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.multisession.Dcov <- function(data,inits=NA,M1=NA,M2=NA,marktype="premarked",obstype="poisson"){
  N.session <- length(data)
  if(length(M1)!=N.session)stop("Must supply an M1 for each session.")
  if(length(M2)!=N.session)stop("Must supply an M2 for each session.")
  
  init.session <- vector("list",N.session)
  for(g in 1:N.session){
    inits.use <- inits #lam0 and sigma inits vary by session
    inits.use$lam0 <- inits.use$lam0[g]
    inits.use$sigma <- inits.use$sigma[g]
    init.session[[g]] <- init.SMR.Dcov(data[[g]],inits.use,M1=M1[g],M2=M2[g],marktype=marktype,obstype="poisson")
  }
  anyTelemetry <- FALSE
  locs.use <- vector("list",N.session)
  tlocs.sess.max <- rep(NA,N.session)
  for(g in 1:N.session){
    if(all(is.na(init.session[[g]]$locs))){
      locs.use <- NA
    }else{
      if(M1[g]>1){
        tlocs.sess.max[g] <- max(rowSums(!is.na(init.session[[g]]$locs[1:M1[g],,1])))
      }else{
        tlocs.sess.max[g] <- sum(!is.na(init.session[[g]]$locs[1:M1[g],,1]))
      }
      locs.use[[g]] <- init.session[[g]]$locs[1:M1[g],1:tlocs.sess.max[g],1:2]
      anyTelemetry <- TRUE
    }
  }
  
  if(anyTelemetry){
    n.tel.inds <- unlist(lapply(locs.use,nrow))
    n.locs.ind <- matrix(NA,N.session,max(n.tel.inds))
    tel.inds <- matrix(NA,N.session,max(n.tel.inds))
  }else{
    tel.inds <- n.locs.ind <- n.tel.inds <- NA
  }
 
  n.samples <- unlist(lapply(data,function(x){length(x$this.j)}))
  J <- unlist(lapply(data,function(x){nrow(x$X)}))
  n.marked <- unlist(lapply(data,function(x){x$n.marked}))
  M.both <- M1+M2
  maxM.both <- max(M.both)
  s <- array(NA,dim=c(N.session,maxM.both,2))
  z <- matrix(NA,N.session,maxM.both)
  ID <- matrix(NA,N.session,max(n.samples))
  y.full <- array(NA,dim=c(N.session,maxM.both,max(J)))
  y.event <- array(NA,dim=c(N.session,maxM.both,max(J),3))
  K1D <- matrix(NA,N.session,max(J))
  n.fixed <- rep(NA,N.session)
  samp.type <- matrix(NA,N.session,max(n.samples))
  match <- array(NA,dim=c(N.session,max(n.samples),maxM.both))
  
  n.cells <- unlist(lapply(data,function(x){x$n.cells}))
  n.cells.x <- unlist(lapply(data,function(x){x$n.cells.x}))
  n.cells.y <- unlist(lapply(data,function(x){x$n.cells.y}))
  n.cells.max <- max(n.cells)
  n.cells.x.max <- max(n.cells.x)
  n.cells.y.max <- max(n.cells.y)
  res <- unlist(lapply(data,function(x){x$res}))
  cellArea <- res^2
  xlim <- ylim <- matrix(NA,N.session,2)
  x.vals <- matrix(NA,N.session,n.cells.x.max)
  y.vals <- matrix(NA,N.session,n.cells.y.max)
  dSS <- array(NA,dim=c(N.session,n.cells.max,2))
  InSS <- array(0,dim=c(N.session,n.cells.max))
  D.cov <- array(NA,dim=c(N.session,n.cells.max))
  cells <- array(0,dim=c(N.session,n.cells.x.max,n.cells.y.max))
  
  for(g in 1:N.session){
    s[g,1:M.both[g],] <- init.session[[g]]$s
    z[g,1:M.both[g]] <- init.session[[g]]$z
    ID[g,1:n.samples[g]] <- init.session[[g]]$ID
    y.full[g,1:M.both[g],1:J[g]] <- init.session[[g]]$y.full
    y.event[g,1:M.both[g],1:J[g],] <- init.session[[g]]$y.event
    K1D[g,1:J[g]] <- init.session[[g]]$K1D
    samp.type[g,1:n.samples[g]] <- init.session[[g]]$samp.type
    match[g,1:n.samples[g],1:M.both[g]] <- init.session[[g]]$match
    n.fixed[g] <- init.session[[g]]$n.fixed
    xlim[g,] <- data[[g]]$xlim
    ylim[g,] <- data[[g]]$ylim
    K1D[g,1:J[g]] <- data[[g]]$K1D
    x.vals[g,1:n.cells.x[g]] <- data[[g]]$x.vals
    y.vals[g,1:n.cells.y[g]] <- data[[g]]$y.vals
    dSS[g,1:n.cells[g],] <- data[[g]]$dSS
    InSS[g,1:n.cells[g]] <- data[[g]]$InSS
    D.cov[g,1:n.cells[g]] <- data[[g]]$D.cov
    cells[g,1:n.cells.x[g],1:n.cells.y[g]] <- data[[g]]$cells
    if(anyTelemetry){
      # n.tel.inds[g] <- sum(rowSums(!is.na(data$locs[g,,,1]))>0)
      tel.inds[g,1:n.tel.inds[g]] <- init.session[[g]]$tel.inds
      n.locs.ind[g,1:n.tel.inds[g]] <- init.session[[g]]$n.locs.ind
    }
  }
  if(anyTelemetry){
    #reformat locs.use actually
    locs.use2 <- array(NA,dim=c(N.session,max(n.tel.inds,na.rm=TRUE),max(n.locs.ind,na.rm=TRUE),2))
    for(g in 1:N.session){
      locs.use2[g,1:nrow(locs.use[[g]]),1:ncol(locs.use[[g]]),] <- locs.use[[g]]
    }
    #remove unused telemetry dimensions if not all marked individuals telemetered 
    rem.idx <- which(colSums(is.na(n.locs.ind))==N.session)
    if(length(rem.idx)>0){
      n.locs.ind <- n.locs.ind[,-rem.idx]
      tel.inds <- tel.inds[,-rem.idx]
    }
  }else{
    locs.use2 <- NA
  }
  
  #put X in ragged array
  X.new <- array(NA,dim=c(N.session,max(J),2))
  for(g in 1:N.session){
    X.new[g,1:J[g],] <- data[[g]]$X
  }
  dummy.data <- matrix(0,N.session,maxM.both) #dummy data not used, doesn't really matter what the values are
  
  z.data <- matrix(NA,N.session,maxM.both)
  for(g in 1:N.session){
    z.data[g,1:n.marked[g]] <- 1
  }
  
  #format data for marginalized sampler, too
  y.mID <- array(NA,dim=c(N.session,max(M1),max(J)))
  y.mnoID <- y.um <- y.unk <- matrix(NA,N.session,max(J))
  for(g in 1:N.session){
    y.mID[g,1:M1[g],1:J[g]] <- y.event[g,1:M1[g],1:J[g],1]
    y.mnoID[g,1:J[g]] <- colSums((y.event[g,1:M1[g],1:J[g],2]))
    y.um[g,1:J[g]] <- colSums((y.event[g,(M1[g]+1):M.both[g],1:J[g],2]))
    y.unk[g,1:J[g]] <- colSums((y.event[g,1:M.both[g],1:J[g],3]))
  }
  
  return(list(y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,
              s.init=s,z.init=z,ID=ID,y.full=y.full,y.event=y.event,K1D=K1D,J=J,X=X.new,
              n.samples=n.samples,n.fixed=n.fixed,samp.type=samp.type,this.j=data$this.j,match=match,
              xlim=xlim,ylim=ylim,locs=locs.use2,tel.inds=tel.inds,n.locs.ind=n.locs.ind,n.tel.inds=n.tel.inds,
              res=res,cellArea=cellArea,x.vals=x.vals,xlim=xlim,ylim=ylim,
              y.vals=y.vals,dSS=dSS,InSS=InSS,cells=cells,n.cells=n.cells,n.cells.x=n.cells.x,
              n.cells.y=n.cells.y,D.cov=D.cov,dummy.data=dummy.data,z.data=z.data))

}