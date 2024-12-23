#multisession version of Poisson Dcov Marginal

library(nimble)
source("sim.SMR.multisession.Dcov.R")
source("sim.SMR.Dcov.R") #used by multisession simulator
source("NimbleModel SMR Multisession Poisson Dcov Marginal.R")
source("NimbleFunctions SMR Multisession Poisson Dcov Marginal.R")
source("init.SMR.multisession.Dcov.R")
source("init.SMR.Dcov.R") #used by multisession initializer
source("sSampler Multisession Poisson Dcov Marginal.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

####Simulate some data####
#Here, I'll simulate 3 populations with different n.marked, K, X, and state space areas
#sharing lam0, sigma so they can be shared during estimation
N.session <- 3
n.marked <- c(12,13,14)
lam0 <- rep(0.5,N.session)
sigma <- rep(0.5,N.session)
K <- c(10,11,9) #number of occasions
buff <- rep(2,N.session) #state space buffer

#theta is probability of observing each sample type for marked and unmarked individuals
theta.marked <- matrix(rep(c(0.75,0.15,0.1),N.session),nrow=N.session,byrow=TRUE) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked <- rep(0.75,N.session) #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)
marktype <- "premarked" #are individuals premarked, or naturally marked? This test script only handles premarked.
# marktype <- "natural"
obstype <- "poisson"
tlocs <- c(0,0,0) #number of telemetry locs/marked individual in each session. For "premarked"

#make an SCR trapping array. Making the trapping array size vary by session
X <- vector("list",N.session)
X[[1]] <- as.matrix(expand.grid(1:10,1:10))
X[[2]] <- as.matrix(expand.grid(1:9,1:9))
X[[3]] <- as.matrix(expand.grid(1:11,1:11))

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- ylim <- matrix(NA,N.session,2)
for(g in 1:N.session){
  xlim[g,] <- range(X[[g]][,1]) + c(-buff[g],buff[g])
  ylim[g,] <- range(X[[g]][,2]) + c(-buff[g],buff[g])
}

#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
for(g in 1:N.session){
  x.shift <- xlim[g,1]
  y.shift <- ylim[g,1]
  xlim[g,] <- xlim[g,] - x.shift
  ylim[g,] <- ylim[g,] - y.shift
  X[[g]][,1] <- X[[g]][,1]- x.shift
  X[[g]][,2] <- X[[g]][,2]- y.shift
}

res <- rep(0.25,N.session) #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- y.vals <- dSS <- cells <- vector("list",N.session)
n.cells <- n.cells.x <- n.cells.y <- rep(NA,N.session)
for(g in 1:N.session){
  x.vals[[g]] <- seq(xlim[g,1]+res[g]/2,xlim[g,2]-res[g]/2,res[g]) #x cell centroids
  y.vals[[g]] <- seq(ylim[g,1]+res[g]/2,ylim[g,2]-res[g]/2,res[g]) #y cell centroids
  dSS[[g]] <- as.matrix(cbind(expand.grid(x.vals[[g]],y.vals[[g]])))
  cells[[g]] <- matrix(1:nrow(dSS[[g]]),nrow=length(x.vals[[g]]),ncol=length(y.vals[[g]]))
  n.cells[g] <- nrow(dSS[[g]])
  n.cells.x[g] <- length(x.vals[[g]])
  n.cells.y[g] <- length(y.vals[[g]])
}

#create a density covariate - one for each session
library(geoR)
D.cov <- vector("list",N.session)
#need a simulated landscape with individuals living around traps to be captured
#these are pretty good
D.seeds <- c(13223,13216,13252)
for(g in 1:N.session){
  set.seed(D.seeds[g])
  D.cov.tmp <- grf(n.cells[g],grid=dSS[[g]],cov.pars=c(500,500),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
  D.cov.tmp <- as.numeric(scale(D.cov.tmp)) #scale
  par(mfrow=c(1,1),ask=FALSE)
  D.cov[[g]] <- D.cov.tmp
  image(x.vals[[g]],y.vals[[g]],matrix(D.cov[[g]],n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g," D.cov"),xlab="X",ylab="Y",col=cols1)
}

#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
InSS <- vector("list",N.session)
for(g in 1:N.session){
  dSS.tmp <- dSS[[g]] - res[g]/2 #convert back to grid locs
  InSS[[g]] <- rep(1,length(D.cov[[g]]))
  InSS[[g]][dSS.tmp[,1]<2&dSS.tmp[,2]<2] <- 0
  InSS[[g]][dSS.tmp[,1]<2&dSS.tmp[,2]>(ylim[g,2]-2)] <- 0
  InSS[[g]][dSS.tmp[,1]>(xlim[g,2]-2)&dSS.tmp[,2]<2] <- 0
  InSS[[g]][dSS.tmp[,1]>(xlim[g,2]-2)&dSS.tmp[,2]>(ylim[g,2]-2)] <- 0
  image(x.vals[[g]],y.vals[[g]],matrix(InSS[[g]],n.cells.x[g],n.cells.y[g]),main=paste("Session",g," Habitat"))
}

#Density covariates
D.beta0 <- rep(-1,N.session)
D.beta1 <- rep(0.5,N.session)
#what is implied expected N in state space?
for(g in 1:N.session){
  lambda.cell <- exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  print(sum(lambda.cell)) #expected N in state space
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell,n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g," Expected Density"),col=cols1)
  points(X[[g]],pch=4,cex=0.75)
}

#Simulate some data
#setting seed here because I am setting a seed to produce the D.cov and you will simulate the same
#data set over and over if you don't use different seeds here for each data set you simulate
set.seed(13435342) #change seed for new data set
data <- sim.SMR.multisession.Dcov(N.session=N.session,n.marked=n.marked,marktype=marktype,
             theta.marked=theta.marked,theta.unmarked=theta.unmarked,
             lam0=lam0,sigma=sigma,K=K,X=X,tlocs=tlocs,obstype=obstype,
             D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,res=res,xlim=xlim,ylim=ylim)

#What is the observed data?
g <- 1 #session to look at
head(t(data[[g]]$samp.type)) #vector of each samp.type of each detection. Must be in this order (I think).
table(data[[g]]$samp.type)
head(t(data[[g]]$this.j)) #trap of capture for each sample
head(t(data[[g]]$this.k)) #occasion of each capture for each sample (not used in this "2D" sampler)
head(data[[g]]$ID.marked) #true ID's for marked and identified samples ()
str(data[[g]]$locs) #possibly telemetry. N.session x n.marked x tlocs x 2 array (or ragged array if number of locs/ind and/or session differ). 

####Fit model in Nimble####
if(marktype=="natural"){
  M1 <- rep(75,N.session) #Augmentation level for marked if number marked unknown
}else{
  M1 <- n.marked #Set to n.marked if premarked.
}
M2 <- c(150,150,150) #Augmentation level for unmarked
#Monitor N.M and N, marked and total ind abundance to make sure N.M does not hit M1 (unless known)
#and N does not hit M1+M2 during sampling. If so, raise the offending M and run again.
M.both <- M1 + M2
#Need some inits to initialize data
#Use reasonable inits for lam0 and sigma since we check to make sure initial observation
#model likelihood is finite
#also use this function checks to make sure theta.marked and theta.unmarked inits are in
#the correct structure. 
inits <- list(lam0=rep(1,N.session),sigma=rep(1,N.session)) #initializing with 1 parameter per session, just set all to same value

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild <- init.SMR.multisession.Dcov(data,inits,M1=M1,M2=M2,marktype=marktype,obstype="poisson")

#inits for nimble
theta.unmarked.init <- matrix(c(0,0.5,0.5),N.session,3,byrow=TRUE)

N.init <- rowSums(nimbuild$z.init,na.rm=TRUE)
N.M.init <- N.UM.init <- rep(NA,N.session)
for(g in 1:N.session){
  N.UM.init[g] <- sum(nimbuild$z.init[g,(M1[g]+1):M.both[g]])
  N.M.init[g] <- sum(nimbuild$z.init[g,1:M1[g]],na.rm=TRUE)
}
(N.init-N.UM.init)==n.marked #should be n.marked[g] individuals in initialized data
N.M.init==n.marked #should be true for "premarked"
#init for total D0, but we split into marked and unmarked D0 to fit model
D0.init <- rowSums(nimbuild$z.init,na.rm=TRUE)/(rowSums(nimbuild$InSS)*nimbuild$cellArea)
#initializing assuming 50:50 split
D0.M.init <- D0.init*0.5
D0.UM.init <- D0.init*0.5

Niminits <- list(N=N.init,N.M=N.M.init,N.UM=N.UM.init,
                 z=nimbuild$z.init,N=rowSums(nimbuild$z.init),s=nimbuild$s.init,
                 lam0.fixed=0.5,sigma.fixed=0.5,D0.M=D0.M.init,D0.UM=D0.UM.init,
                 D.beta1=rep(0,N.session))

#constants for Nimble
J <- unlist(lapply(data,function(x){nrow(x$X)}))
#If you do not have telemetry use these. Make sure telemetry BUGS code is commented out.
constants <- list(N.session=N.session,M1=M1,M.both=M.both,J=J,K1D=nimbuild$K1D,
                xlim=nimbuild$xlim,ylim=nimbuild$ylim,D.cov=nimbuild$D.cov,res=nimbuild$res,
                cellArea=nimbuild$cellArea,n.cells=nimbuild$n.cells)

# Supply data to Nimble. marginalized data formatted in nimbuild object
Nimdata <- list(y.mID=nimbuild$y.mID, #marked with ID
                y.mnoID=nimbuild$y.mnoID, #marked without ID
                y.um=nimbuild$y.um, #unmarked
                y.unk=nimbuild$y.unk, #unk marked status
                z=nimbuild$z.data,X=nimbuild$X,cells=nimbuild$cells,
                dummy.data=nimbuild$dummy.data,InSS=nimbuild$InSS)

#If you have telemetry use these instead. Make sure to uncomment telemetry BUGS code.
# constants <- list(N.session=N.session,M1=M1,M.both=M.both,J=J,K1D=nimbuild$K1D,
#                   xlim=nimbuild$xlim,ylim=nimbuild$ylim,D.cov=nimbuild$D.cov,res=nimbuild$res,
#                   cellArea=nimbuild$cellArea,n.cells=nimbuild$n.cells,
#                   #telemetry stuff
#                   tel.inds=nimbuild$tel.inds,
#                   n.tel.inds=nimbuild$n.tel.inds,n.locs.ind=nimbuild$n.locs.ind)
# Nimdata <- list(y.mID=nimbuild$y.mID, #marked with ID
#                 y.mnoID=nimbuild$y.mnoID, #marked without ID
#                 y.um=nimbuild$y.um, #unmarked
#                 y.unk=nimbuild$y.unk, #unk marked status
#                 z=nimbuild$z.data,X=nimbuild$X,cells=nimbuild$cells,
#                 dummy.data=nimbuild$dummy.data,InSS=nimbuild$InSS,locs=nimbuild$locs)

# set parameters to monitor
parameters <- c('D0','D0.M','D0.UM','D.beta1','lambda.N.M','lambda.N.UM',
                'N.M','N.UM','N','lam0.fixed','sigma.fixed','theta.marked.fixed','theta.unmarked.fixed')
parameters2 <- c("lambda.cell","s.cell",'D0') #record D0 here for plotting
nt <- 1 #thinning rate for parameters
nt2 <- 10 #thinning rate for parameters2


# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#if you change the model file, make sure you make necessary changes here
config.nodes <- c('lam0.fixed','sigma.fixed','theta.marked.fixed','theta.unmarked.fixed[2:3]')
#this one for session-specific theta.marked and theta.unmarked
# config.nodes <- c('lam0.fixed','sigma.fixed','theta.marked',paste('theta.unmarked[1:',N.session,',2:3',']'))
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,monitors2=parameters2, thin2=nt2,nodes=config.nodes)

# how many z proposals per iteration per session for marked and unmarked individuals?
z.ups <- round(cbind(M1,M2)*0.25) #doing 25% of M1 and M2 here
z.ups #note, this is a matrix with one row per session. If number of marked inds is known (premarked scenario), 
#first column of z.ups ignored and marked z's not updated because known. updateMarked=FALSE then.

#do we need to update N.M?
if(marktype=="natural"){
  updateMarked <- TRUE
}else{
  updateMarked <- FALSE #if premarked, N.M is known
}
J <- nimbuild$J
for(g in 1:N.session){
  #nodes used for update, calcNodes + z nodes
  lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",1:",M.both[g],",1:",J[g],"]"))
  y.mID.nodes <- Rmodel$expandNodeNames(paste("y.mID[",g,",1:",M1[g],",1:",J[g],"]"))
  y.mnoID.nodes <- Rmodel$expandNodeNames(paste("y.mnoID[",g,",1:",J[g],"]"))
  y.um.nodes <- Rmodel$expandNodeNames(paste("y.um[",g,",1:",J[g],"]"))
  y.unk.nodes <- Rmodel$expandNodeNames(paste("y.unk[",g,",1:",J[g],"]"))
  bigLam.marked.nodes <- Rmodel$expandNodeNames(paste("bigLam.marked[",g,",1:",J[g],"]")) #only need this in calcNodes
  bigLam.unmarked.nodes <- Rmodel$expandNodeNames(paste("bigLam.unmarked[",g,",1:",J[g],"]")) #only need this in calcNodes
  lam.mnoID.nodes <- Rmodel$expandNodeNames(paste("lam.mnoID[",g,",1:",J[g],"]"))
  lam.um.nodes <- Rmodel$expandNodeNames(paste("lam.um[",g,",1:",J[g],"]"))
  lam.unk.nodes <- Rmodel$expandNodeNames(paste("lam.unk[",g,",1:",J[g],"]"))
  N.M.node <- Rmodel$expandNodeNames(paste("N.M[",g,"]"))
  N.UM.node <- Rmodel$expandNodeNames(paste("N.UM[",g,"]"))
  N.node <- Rmodel$expandNodeNames(paste("N[",g,"]"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",1:",M.both[g],"]"))
  calcNodes <- c(N.M.node,N.UM.node,N.node=N.node,lam.nodes,bigLam.marked.nodes,bigLam.unmarked.nodes,
                 lam.mnoID.nodes,lam.um.nodes,lam.unk.nodes,
                 y.mID.nodes,y.mnoID.nodes,y.um.nodes,y.unk.nodes)
  inds.detected <- which(rowSums(nimbuild$y.mID[g,,],na.rm=TRUE)>0) #inds with at least 1 mID detection
  conf$addSampler(target = paste("N.UM"),
                  type = 'zSampler',control = list(g=g,updateMarked=updateMarked,z.ups=z.ups[g,1:2],
                                                   J=J[g],M1=M1[g],M.both=M.both[g],inds.detected=inds.detected,
                                                   y.mID.nodes=y.mID.nodes,y.mnoID.nodes=y.mnoID.nodes,
                                                   y.um.nodes=y.um.nodes,y.unk.nodes=y.unk.nodes,
                                                   lam.nodes=lam.nodes,lam.mnoID.nodes=lam.mnoID.nodes,
                                                   lam.um.nodes=lam.um.nodes,lam.unk.nodes=lam.unk.nodes,
                                                   N.M.node=N.M.node,
                                                   N.UM.node=N.UM.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),
                  silent = TRUE)
}

#add sSampler
#if no telemetry,
for(g in 1:N.session){
  for(i in 1:M.both[g]){
    conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,g=g,J=J[g],res=nimbuild$res[g],n.cells=nimbuild$n.cells[g],
                                                   n.cells.x=nimbuild$n.cells.x[g],n.cells.y=nimbuild$n.cells.y[g],
                                                   xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],res=nimbuild$res[g],
                                                   M1=M1[g],n.locs.ind=0,scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}
#if telemetry
# for(g in 1:N.session){
#   for(i in 1:M.both[g]){
#     if(i %in% nimbuild$tel.inds[g,]){#inds with telemetry
#       conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
#                       type = 'sSampler',control=list(i=i,g=g,J=J[g],res=nimbuild$res[g],n.cells=nimbuild$n.cells[g],
#                                                      n.cells.x=nimbuild$n.cells.x[g],n.cells.y=nimbuild$n.cells.y[g],
#                                                      xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],res=nimbuild$res[g],
#                                                      n.locs.ind=nimbuild$n.locs.ind[g,i],M1=M1[g],scale=1),silent = TRUE)
#       #scale parameter here is just the starting scale. It will be tuned.
#     }else{ #inds with no telemetry
#       conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
#                       type = 'sSampler',control=list(i=i,g=g,J=J[g],res=nimbuild$res[g],n.cells=nimbuild$n.cells[g],
#                                                      n.cells.x=nimbuild$n.cells.x[g],n.cells.y=nimbuild$n.cells.y[g],
#                                                      xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],res=nimbuild$res[g],
#                                                      n.locs.ind=0,M1=M1[g],scale=1),silent = TRUE)
#     }
#   }
# }

#maybe replace independent lam0 and sigma samplers with block sampler better accommodating for posterior covariance
#should improve mixing and increase posterior effective sample size. AF_slice works better than block RW. 
#Need to not use this update or modify it when using lam0 or sigma covariates.
#This sampler is slower, so not worth it if data is not so sparse there is strong posterior correlation
#between lam0 and sigma.

#mixing seems generally good if you keep independent samplers and add block sampler (if data sparse enough that there is posterior correlation)
# conf$removeSampler(c("lam0.fixed","sigma.fixed"))
conf$addSampler(target = c("lam0.fixed","sigma.fixed"),type = 'RW_block',
                control = list(adaptive=TRUE),silent = TRUE)
#AF_slice pretty efficient for Dcov parameters. Block by session
for(g in 1:N.session){
  conf$addSampler(target = c(paste("D0.M[",g,"]"),paste("D0.UM[",g,"]"),paste("D.beta1[",g,"]")),
                  type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
}

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
#can ignore warning about pi.cell
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[250:nrow(mvSamples),]))

exp(D.beta0)
unlist(lapply(data,function(x){x$lambda.N}))#expected Ns
unlist(lapply(data,function(x){x$N})) #realized Ns


#Important! If N.UM[g] hits M2[g] during sampling, raise M2[g]. 
#For an unknown number of marked individuals, if N.M hits M1 during sampling, raise M1.

#Look at cell-level expected density estimates, compare to truth
mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx <- grep("lambda.cell",colnames(mvSamples2))
n.cells.max <- max(n.cells)
lambda.cell.idx <- matrix(lambda.cell.idx,N.session,n.cells.max)
D0.idx <- grep("D0",colnames(mvSamples2))
burnin2 <- 10 #consider nt2 thinning rate when setting burnin2

#compare expected D plot to truth
#image will show posterior means
n.iter.use <- burnin2:nrow(mvSamples2)
lambda.cell.post <- array(NA,dim=c(N.session,n.cells.max,length(n.iter.use)))
lambda.cell <- lambda.cell.ests <- array(NA,dim=c(N.session,n.cells.max))
lambda.cell.HPDs <- array(NA,dim=c(N.session,n.cells.max,2))
for(g in 1:N.session){
  lambda.cell[g,1:n.cells[g]] <- exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  lambda.cell.post[g,1:n.cells[g],] <- t(cellArea[g]*mvSamples2[n.iter.use,D0.idx[g]]*
                                           mvSamples2[n.iter.use,lambda.cell.idx[g,1:n.cells[g]]])
  lambda.cell.ests[g,1:n.cells[g]] <- rowMeans(lambda.cell.post[g,1:n.cells[g],])
  lambda.cell.HPDs[g,1:n.cells[g],] <- HPDinterval(mcmc(t(lambda.cell.post[g,1:n.cells[g],])))
  #remove nonhabitat (or not, comment out)
  lambda.cell[g,InSS[[g]]==0] <- NA
  lambda.cell.ests[g,InSS[[g]]==0] <- NA
}

par(mfrow=c(1,1),ask=FALSE)
for(g in 1:N.session){
  zlim <- range(c(lambda.cell[g,],lambda.cell.ests[g,]),na.rm=TRUE) #use same zlim for plots below
  #truth
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell[g,1:n.cells[g]],n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g,"True Expected Density"),zlim=zlim)
  points(X[[g]],pch=4)
  #estimate, posterior means
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell.ests[g,1:n.cells[g]],n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g,"Est Expected Density"),zlim=zlim)
  points(X[[g]],pch=4)
}
#cell ests and 95% HPDs vs. truth. 
#Need a lot of posterior samples for accurate 95% HPDs, if not, will look "jagged"
for(g in 1:N.session){
  idx <- order(lambda.cell[g,1:n.cells[g]])
  plot(lambda.cell.ests[g,1:n.cells[g]][idx]~lambda.cell[g,1:n.cells[g]][idx],type="l",lwd=2,
       main=paste("Session",g,"True vs. Estimated Density"),ylim=range(lambda.cell.HPDs[g,1:n.cells[g],]))
  lines(lambda.cell.HPDs[g,1:n.cells[g],1][idx]~lambda.cell[g,1:n.cells[g]][idx],lty=2)
  lines(lambda.cell.HPDs[g,1:n.cells[g],2][idx]~lambda.cell[g,1:n.cells[g]][idx],lty=2)
  abline(0,1,col="darkred",lwd=2) #1:1 expectation
}
