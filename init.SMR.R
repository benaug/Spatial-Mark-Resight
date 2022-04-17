e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR=function(data,inits=NA,M1=NA,M2=NA,marktype="premarked",obstype="poisson"){
  library(abind)
  #extract observed data
  this.j=data$this.j
  # this.k=data$this.k #not used in this 2D sampler
  samp.type=data$samp.type
  ID.marked=data$ID.marked
  n.marked=data$n.marked
  X<-as.matrix(data$X)
  J<-nrow(X)
  K<- data$K
  K1D=data$K1D
  buff<- data$buff
  M.both=M1+M2
  locs=data$locs
  
  xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  
  ##pull out initial values
  lam0=inits$lam0
  sigma<- inits$sigma
  
  n.samp1=sum(samp.type=="markedID")
  n.samp2=sum(samp.type=="markednoID")
  n.samp3=sum(samp.type=="unmarked")
  n.samp4=sum(samp.type=="unk")
  n.samples=length(this.j)
  
  useMarkednoID=FALSE
  if(n.samp2>0){
    useMarkednoID=TRUE
  }
  useUnk=FALSE
  if(n.samp4>0){
    useUnk=TRUE
  }
  
  #build y.marked
  y.marked=matrix(0,M1,J)
  for(l in 1:length(ID.marked)){
    y.marked[ID.marked[l],this.j[l]]=y.marked[ID.marked[l],this.j[l]]+1
  }
  
  #initialize unknown IDs
  G.true=matrix(c(rep(1,M1),rep(2,M2)),ncol=1) #individual marked and unmarked statuses
  ID=c(ID.marked,rep(NA,n.samples-length(ID.marked)))
  nextID=max(ID,na.rm=TRUE)+1
  
  y.true2D=apply(y.marked,c(1,2),sum)
  y.true2D=rbind(y.true2D,matrix(0,nrow=M2,ncol=J))
  if(M1<n.marked)stop("M1 must be larger than the number of marked individuals.")
  #marked noID first
  if(useMarkednoID){
    matches=1:M1
    if(marktype=="natural"){
      for(l in (n.samp1+1):(n.samp1+n.samp2)){
        sametrap=y.true2D[matches,this.j[l]]>0
        if(any(sametrap)){
          ID[l]=matches[which(sametrap)[1]]
          y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        }else{#must be new ID
          ID[l]=nextID
          y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
          nextID=nextID+1
        }
        if(nextID>M1)stop("Need to raise M1 to initialize data (marktype=natural)")
      }
    }else{ #premarked requires different approach because number of marks known
      for(l in (n.samp1+1):(n.samp1+n.samp2)){
        #chose match with closest other capture, if you match a captured ID
        assignedIDsthatmatch=which(ID%in%matches)
        if(length(assignedIDsthatmatch)>0){
          thesetraps=X[this.j[assignedIDsthatmatch],]
          if(length(assignedIDsthatmatch)>1){
            dists=sqrt((X[this.j[l],1]-thesetraps[,1])^2+(X[this.j[l],2]-thesetraps[,2])^2)
          }else{
            dists=sqrt((X[this.j[l],1]-thesetraps[1])^2+(X[this.j[l],2]-thesetraps[2])^2)
          }
          pick=which(dists==min(dists))[1]#pick first one if multiple to choose from
          ID[l]=ID[assignedIDsthatmatch[pick]]
          y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        }else{#if none of your matches were captured, pick first match
         ID[l]=matches[1]
         y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        }
      }
    }
  }
  #unmarked next...
  nextID=M1+1
  matches=(M1+1):M.both
  for(l in (n.samp1+n.samp2+1):(n.samp1+n.samp2+n.samp3)){
    #can you match an unmarked guy in same trap already assigned an ID?
    sametrap=y.true2D[matches,this.j[l]]>0
    if(any(sametrap)){
      ID[l]=matches[which(sametrap)[1]]
      y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
    }else{#must be new ID
      ID[l]=nextID
      y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
      nextID=nextID+1
    }
    if(nextID>M.both)stop("Need to raise M2 to initialize data.")
  }
  #then unknown...
  #keeping nextID where it is so we assign non matches to unmarked class
  if(useUnk){
    matches=1:M.both
    for(l in (n.samp1+n.samp2+n.samp3+1):(n.samp1+n.samp2+n.samp3+n.samp4)){
      #can you match an unmarked guy in same trap already assigned an ID?
      sametrap=y.true2D[matches,this.j[l]]>0
      if(any(sametrap)){
        ID[l]=matches[which(sametrap)[1]]
        y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
      }else{#must be new ID
        ID[l]=nextID
        y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        nextID=nextID+1
      }
    }
    if(nextID>M.both)stop("Need to raise M2 to initialize data.")
  }
  
  #initialize match
  match=matrix(FALSE,nrow=n.samples,ncol=M1+M2)
  #marked no ID samples can only match marked guys
  idx=which(samp.type=="markednoID")
  match[idx,1:M1]=TRUE
  #unmarked samples can only match unmarked guys
  idx=which(samp.type=="unmarked")
  match[idx,(M1+1):M.both]=TRUE
  #unk samps can match anyone
  idx=which(samp.type=="unk")
  match[idx,]=TRUE
  
  #initialize z
  z=1*(rowSums(y.true2D)>0)
  if(marktype=="premarked"){
    z[1:M1]=1
  }
  
  #intialize s
  s<- cbind(runif(M.both,xlim[1],xlim[2]), runif(M.both,ylim[1],ylim[2])) #assign random locations
  idx=which(rowSums(y.true2D)>0) #switch for those actually caught
  for(i in idx){
    trps<- matrix(X[y.true2D[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s[i,]<- trps
    }
  }
  
  #identify known ID (1), known status, unknown ID (2), unknown ID and status (3)
  G.type=rep(0,n.samples)
  G.type[samp.type=="markedID"]=1
  G.type[samp.type%in%c("markednoID","unmarked")]=2
  G.type[samp.type=="unk"]=3
  n.fixed=sum(samp.type=="markedID")
  
  if(!is.null(dim(data$locs))){
    max.locs=dim(locs)[2]
    tel.inds=which(rowSums(is.na(locs[,,1]))<max.locs)
    n.locs.ind=rowSums(!is.na(locs[,,1]))
    n.locs.ind=n.locs.ind[tel.inds]
    print("using telemetry to initialize telmetered s. Remove from data if not using in the model.")
    #update s starts for telemetry guys
    for(i in tel.inds){
      if(n.locs.ind[i]>1){
        s[i,]=colMeans(locs[i,1:n.locs.ind[i],])
      }else{
        s[i,]=locs[i,1,]
      }
      #make sure s is in state space
      if(s[i,1]<xlim[1]){
        s[i,1]=xlim[1]
      }
      if(s[i,1]>xlim[2]){
        s[i,1]=xlim[2]
      }
      if(s[i,2]<ylim[1]){
        s[i,2]=ylim[1]
      }
      if(s[i,2]>ylim[2]){
        s[i,2]=ylim[2]
      }
    }
  }else{
    tel.inds=NA
    n.locs.ind=NA
  }
  
  D=e2dist(s, X)
  lamd<- lam0*exp(-D*D/(2*sigma*sigma))
  #make 3D y.true
  y.true3D=array(0,dim=c(M.both,J,3))
  for(l in 1:n.samples){
    y.true3D[ID[l],this.j[l],G.type[l]]=y.true3D[ID[l],this.j[l],G.type[l]]+1
  }
  y.true2D=apply(y.true3D,c(1,2),sum)
  
  if(obstype=="poisson"){
    ll.y=dpois(y.true2D,K1D*lamd*z,log=TRUE)
  }else if(obstype=="negbin"){
    theta.d=inits$theta.d
    ll.y=y.true2D*0
    for(i in 1:M.both){
      if(z[i]==1){
        ll.y[i,]=dnbinom(y.true2D[i,],mu=lamd[i,],size=theta.d*K1D,log=TRUE)
      }
    }
  }else{
    stop("obstype not recognized")
  }
  
  if(!is.finite(sum(ll.y)))stop("Starting observation model likelihood not finite. Possible error in K1D (if supplied by user) or problem initializing data.")
  
  return(list(s=s,z=z,ID=ID,y.full=y.true2D,y.event=y.true3D,K1D=K1D,
         n.samples=n.samples,n.fixed=n.fixed,samp.type=G.type,this.j=this.j,match=match,
         xlim=xlim,ylim=ylim,locs=locs,tel.inds=tel.inds,n.locs.ind=n.locs.ind))

}