#This is an SMR data simulator and MCMC sampler that handles all sample types
#1) marked, known ID
#2) marked, unknown ID
#3) unmarked, unknown ID
#4) unknown marked status, unknown ID

#It handles both "premarked" scenarios where you know the number of marked individuals
#and "natural" where marks are from natural patterns on the animals so the number of marked 
#individuals is unknown. For "natural", consider using the random thinning model with partial ID covariates.

#This sampler also handles telemetry for marked individuals in the "premarked" scenario. Don't try
#to use telemetry with "natural" scenario. Not realistic and I'm not sure what my code will do!

#y[i,j,k] ~ Poisson(lam[i,j,k])
#y.event[i,j,k,1:3] ~ Multinomial(theta.marked[1:3],y[i,j,k]) for marked i
#y.event[i,j,k,1:3] ~ Multinomial(theta.unmarked[1:3],y[i,j,k]) for unmarked i

#event 1 is you know the ID (marked known ID samples)
#event 2 is you know the mark status, but not ID (marked, unknown ID or unmarked samples)
#event 3 is you don't know mark status or ID (unknown marked status samples)

#This version allows spatial density covariates and/or a habitat mask. Estimates were roughly unbiased in 2 simulation scenarios
#with nominal coverage using 20 marked individuals and a density covariate where high D areas are well covered by the trapping
#array. There was some bias in D.beta1 and D0 using a previous density covariate because many of the 20 marked individuals lived in the highest
#density areas that were farther away from the grid, leaving few to be detected. Still, the N estimates in this case were roughly
#unbiased. A good design is probably more important than when using
#regular SCR. Generally, I believe there will be more bias with fewer marked individuals and/or lower proportion of marked 
#individual samples being identified to individual (theta.marked[1] lower), and fewer samples being of known marked status.
#Also, more bias with natural marks vs. premarked because you then need to estimate the number of marked individuals.
#So, evaluating designs via simulation is a good idea.
#I think one challenge for D covs with latent ID models is that areas with highest D have the most home range overlap and 
#thus uncertainty in individual ID.

library(nimble)
source("sim.SMR.Dcov.R")
source("NimbleModel SMR Poisson Dcov DA2 Marginal.R")
source("NimbleFunctions SMR Poisson Dcov DA2 Marginal.R")
source("init.SMR.Dcov.R")
source("sSampler Poisson Dcov Marginal.R")
source("mask.check.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
n.marked <- 20 #for "premarked", this is the number marked. for "natural", this is the number of captured individuals that are identifiable
lam0 <- 0.5
sigma <- 0.5
K <- 10 #number of occasions
buff <- 2 #state space buffer
X <- expand.grid(3:11,3:11) #make a trapping array
#theta is probability of observing each sample type for marked and unmarked individuals
theta.marked <- c(0.75,0.15,0.1) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked <- 0.75 #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)
marktype <- "premarked" #are individuals premarked, or naturally marked?
# marktype <- "natural"
obstype <- "poisson"
tlocs <- 0 #number of telemetry locs/marked individual. For "premarked"

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- range(X[,1]) + c(-buff,buff)
ylim <- range(X[,2]) + c(-buff,buff)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim-x.shift
ylim <- ylim-y.shift
X[,1] <- X[,1]-x.shift
X[,2] <- X[,2]-y.shift

res <- 0.25 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#simulate a D.cov, higher cov.pars for large scale cov
#change seed to get new D.cov. trial and error to create one with good trapping array coverage
set.seed(152)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(500,500),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)
points(X,pch=4)

#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(1,length(D.cov))
InSS[dSS.tmp[,1]<2&dSS.tmp[,2]<2] <- 0
InSS[dSS.tmp[,1]<2&dSS.tmp[,2]>12] <- 0
InSS[dSS.tmp[,1]>12&dSS.tmp[,2]<2] <- 0
InSS[dSS.tmp[,1]>12&dSS.tmp[,2]>12] <- 0

image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y),main="Habitat")

#Density covariates
D.beta0 <- -0.5 #data simulator uses intercept for marked + unmarked
D.beta1 <- 0.5
#what is implied expected N in state space?
lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",col=cols1)
points(X,pch=4,cex=1,lwd=2)

#Simulate some data
#setting seed here because I am setting a seed to produce the D.cov and you will simulate the same
#data set over and over if you don't use different seeds here for each data set you simulate
set.seed(143532) #change seed for new data set
data <- sim.SMR.Dcov(D.beta0=D.beta0,D.beta1=D.beta1,res=res,
                D.cov=D.cov,InSS=InSS,n.marked=n.marked,marktype=marktype,
                theta.marked=theta.marked,theta.unmarked=theta.unmarked,
                lam0=lam0,sigma=sigma,K=K,X=X,xlim=xlim,ylim=ylim,tlocs=tlocs,
                obstype=obstype)
points(data$s,pch=16) #add activity centers


#What is the observed data?
head(data$samp.type) #vector of each samp.type of each detection. Must be in this order (I think).
table(data$samp.type)
head(data$this.j) #trap of capture for each sample
head(data$this.k) #occasion of each capture for each sample (not used in this "2D" sampler)
head(data$ID.marked) #true ID's for marked and identified samples ()
str(data$locs) #possibly telemetry. n.marked x tlocs x 2 array (or ragged array if number of locs/ind differ). 
#Rows are 1:n.marked individuals, columns are max telemetry points for a single
#individual, fill in NAs for inds with no telemetry and/or inds without max number of telemetry points.
#in latter case, order telemetry points first, then NAs

#function to test for errors in mask set up. 
mask.check(dSS=data$dSS,cells=data$cells,n.cells=data$n.cells,n.cells.x=data$n.cells.x,
           n.cells.y=data$n.cells.y,res=data$res,xlim=data$xlim,ylim=data$ylim,
           x.vals=data$x.vals,y.vals=data$y.vals)

####Fit model in Nimble####
if(marktype=="natural"){
  M1 <- 100 #Augmentation level for marked.
}else{
  M1 <- n.marked #Set to n.marked if premarked. psi1 will be estimated, but can be ignored.
}
M2 <- 250 #Augmentation level for unmarked
#Monitor N.M and N, marked and total ind abundance to make sure N.M does not hit M1 (unless known)
#and N does not hit M1+M2 during sampling. If so, raise the offending M and run again.
M.both <- M1 + M2
#Need some inits to initialize data
#Use reasonable inits for lam0 and sigma since we check to make sure initial observation
#model likelihood is finite
#also use this function checks to make sure theta.marked and theta.unmarked inits are in
#the correct structure. 
inits <- list(lam0=lam0,sigma=sigma)

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild <- init.SMR.Dcov(data,inits,M1=M1,M2=M2,marktype=marktype,obstype="poisson")
points(nimbuild$s,pch=16) #initialized activity centers

#inits for nimble
#init for total D0, but we split into marked and unmarked D0 to fit model
D0.init <- sum(nimbuild$z)/(sum(data$InSS)*data$res^2)
#initializing assuming 50:50 split
D0.M.init <- D0.init*0.5
D0.UM.init <- D0.init*0.5

#must initialize N.M and N.UM to be consistent with z. speeds converge to set consistent with lambda.N.M/UM
Niminits <- list(z=nimbuild$z,s=nimbuild$s,D0.M=D0.M.init,D0.UM=D0.UM.init,D.beta1=0,
                 N.M=sum(nimbuild$z[1:M1]),lambda.N.M=sum(nimbuild$z[1:M1]),
                 N.UM=sum(nimbuild$z[(M1+1):M.both]),lambda.N.UM=sum(nimbuild$z[(M1+1):M.both]),
                 theta.unmarked=c(0,0.5,0.5),lam0=inits$lam0,sigma=inits$sigma)


J <- nrow(data$X)
z.data <- c(rep(1,data$n.marked),rep(NA,M.both-data$n.marked))
dummy.data <- rep(0,M.both) #dummy data not used, doesn't really matter what the values are

#Use this if you do not have telemetry. Make sure telemetry commented out in model file
constants <- list(M1=M1,M.both=M.both,J=J,K1D=data$K1D,
                  D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
                  xlim=nimbuild$xlim,ylim=nimbuild$ylim,res=data$res)

Nimdata <- list(y.mID=nimbuild$y.event[1:M1,,1], #marked with ID
                y.mnoID=colSums((nimbuild$y.event[1:M1,,2])), #marked without ID
                y.um=colSums((nimbuild$y.event[(M1+1):M.both,,2])), #unmarked
                y.unk=colSums((nimbuild$y.event[1:M.both,,3])), #unk marked status
                z=z.data,X=as.matrix(X),
                dummy.data=dummy.data,cells=data$cells,InSS=data$InSS)

#If you have telemetry use these instead. Make sure to uncomment telemetry BUGS code.
# constants <- list(M1=M1,M.both=M.both,J=J,K1D=data$K1D,
#                   D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
#                   xlim=nimbuild$xlim,ylim=nimbuild$ylim,res=data$res,tel.inds=nimbuild$tel.inds,
#                   n.tel.inds=length(nimbuild$tel.inds),n.locs.ind=nimbuild$n.locs.ind)
# Nimdata <- list(y.mID=nimbuild$y.event[1:M1,,1], #marked with ID
#                 y.mnoID=colSums((nimbuild$y.event[1:M1,,2])), #marked without ID
#                 y.um=colSums((nimbuild$y.event[(M1+1):M.both,,2])), #unmarked
#                 y.unk=colSums((nimbuild$y.event[1:M.both,,3])), #unk marked status
#                 z=z.data,X=as.matrix(X),
#                 dummy.data=dummy.data,cells=data$cells,InSS=data$InSS,
#                 locs=data$locs)

# set parameters to monitor
parameters <- c('D0','D0.M','D0.UM','lambda.N.M','lambda.N.UM','lam0','sigma','theta.marked','theta.unmarked',
             'N.M','N.UM','N','D.beta1')
parameters2 <- c("lambda.cell","s.cell",'D0') #record D0 here for plotting

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#use block sampler below for 'D0.M','D0.UM','D.beta1'
config.nodes <- c('lam0','sigma','theta.marked','theta.unmarked[2:3]')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1,
                      monitors2=parameters2,thin2=10,nodes=config.nodes)


#add N/z updates
#do we need to update N.M?
if(marktype=="natural"){
  updateMarked <- TRUE
}else{
  updateMarked <- FALSE #if premarked, N.M is known
}

# how many z proposals per iteration per session for marked (if updated), unmarked?
z.ups <- round(c(M1*0.25,M2*0.25)) #doing 25% of M1 and M2 here
J <- nrow(data$X)

#nodes used for update, calcNodes + z nodes
lam.nodes <- Rmodel$expandNodeNames(paste("lam[1:",M.both,",1:",J,"]"))
y.mID.nodes <- Rmodel$expandNodeNames(paste("y.mID[1:",M1,",1:",J,"]"))
y.mnoID.nodes <- Rmodel$expandNodeNames(paste("y.mnoID[1:",J,"]"))
y.um.nodes <- Rmodel$expandNodeNames(paste("y.um[1:",J,"]"))
y.unk.nodes <- Rmodel$expandNodeNames(paste("y.unk[1:",J,"]"))
bigLam.marked.nodes <- Rmodel$expandNodeNames("bigLam.marked") #only need this in calcNodes
bigLam.unmarked.nodes <- Rmodel$expandNodeNames("bigLam.marked") #only need this in calcNodes
lam.mnoID.nodes <- Rmodel$expandNodeNames("lam.mnoID")
lam.um.nodes <- Rmodel$expandNodeNames("lam.um")
lam.unk.nodes <- Rmodel$expandNodeNames("lam.unk")
N.M.node <- Rmodel$expandNodeNames(paste("N.M"))
N.UM.node <- Rmodel$expandNodeNames(paste("N.UM"))
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M.both,"]"))
calcNodes <- c(N.M.node,N.UM.node,N.node=N.node,lam.nodes,bigLam.marked.nodes,bigLam.unmarked.nodes,
               lam.mnoID.nodes,lam.um.nodes,lam.unk.nodes,
               y.mID.nodes,y.mnoID.nodes,y.um.nodes,y.unk.nodes)
inds.detected <- which(rowSums(nimbuild$y.event[,,1])>0)
conf$addSampler(target = paste("N.UM"),
                type = 'zSampler',control = list(updateMarked=updateMarked,z.ups=z.ups[1:2],
                                                 J=J,M1=M1,M.both=M.both,inds.detected=inds.detected,
                                                 y.mID.nodes=y.mID.nodes,y.mnoID.nodes=y.mnoID.nodes,
                                                 y.um.nodes=y.um.nodes,y.unk.nodes=y.unk.nodes,
                                                 lam.nodes=lam.nodes,lam.mnoID.nodes=lam.mnoID.nodes,
                                                 lam.um.nodes=lam.um.nodes,lam.unk.nodes=lam.unk.nodes,
                                                 N.M.node=N.M.node,
                                                 N.UM.node=N.UM.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),
                silent = TRUE)

#add sSampler
#if no telemetry,
for(i in 1:M.both){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,J=J,res=data$res,n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                                 xlim=nimbuild$xlim,ylim=nimbuild$ylim,n.locs.ind=0,
                                                 M1=M1,scale=0.25),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}
#with telemetry,
# for(i in 1:M.both){
#   if(i %in% nimbuild$tel.inds){#inds with telemetry
#     conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
#                     type = 'sSampler',control=list(i=i,J=J,res=data$res,n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
#                                                    xlim=nimbuild$xlim,ylim=nimbuild$ylim,n.locs.ind=nimbuild$n.locs.ind[i],
#                                                    M1=M1,scale=0.25),silent = TRUE)
#   }else{ #inds with no telemetry
#     conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
#                     type = 'sSampler',control=list(i=i,J=J,res=data$res,n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
#                                                    xlim=nimbuild$xlim,ylim=nimbuild$ylim,n.locs.ind=0,
#                                                    M1=M1,scale=0.25),silent = TRUE)
#   }
# }

#can add block sampler if lam0, sigma, lambda.N.UM, and/or lambda.N.M (if N.M unknown) posteriors highly correlated
#RW_block faster, AFslice slower, but mixes better
conf$addSampler(target = c("lam0","sigma"),type = 'RW_block',
                control = list(adaptive=TRUE),silent = TRUE)
conf$addSampler(target = c("D0.M","D0.UM","D.beta1"),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time


library(coda)
mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[200:nrow(mvSamples),]))

exp(D.beta0)
data$N

cor(mvSamples[200:nrow(mvSamples),])

#Important! If N.UM hits M2 during sampling, raise M2. 
#For an unknown number of marked individuals, if N.M hits M1 during sampling, raise M1.

#plot density surface, etc.
mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx <- grep("lambda.cell",colnames(mvSamples2))
D0.idx <- grep("D0",colnames(mvSamples2))
burnin2 <- 10

#compare expected D plot to truth (for simulated data sets)
n.cells <- data$n.cells
lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
n.iter.use <- burnin2:nrow(mvSamples2)
lambda.cell.post <- t(cellArea*mvSamples2[n.iter.use,D0.idx]*mvSamples2[n.iter.use,lambda.cell.idx[1:n.cells]])
lambda.cell.ests <- rowMeans(lambda.cell.post[1:n.cells,])
lambda.cell.HPDs <- HPDinterval(mcmc(t(lambda.cell.post[1:n.cells,])))
#remove nonhabitat (or not, comment out)
lambda.cell[data$InSS==0] <- NA
lambda.cell.ests[data$InSS==0] <- NA

par(mfrow=c(1,1),ask=FALSE)
zlim <- range(c(lambda.cell,lambda.cell.ests),na.rm=TRUE) #use same zlim for plots below
#truth
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)
#estimate, posterior means
image(x.vals,y.vals,matrix(lambda.cell.ests,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)

#cell ests and 95% HPDs vs. truth. 
#Need a lot of posterior samples for accurate 95% HPDs, if not, will look "jagged"
idx <- order(lambda.cell)
plot(lambda.cell.ests[1:n.cells][idx]~lambda.cell[1:n.cells][idx],type="l",lwd=2,
     main="True vs. Estimated Density",ylim=range(lambda.cell.HPDs[1:n.cells,]))
lines(lambda.cell.HPDs[1:n.cells,1][idx]~lambda.cell[1:n.cells][idx],lty=2)
lines(lambda.cell.HPDs[1:n.cells,2][idx]~lambda.cell[1:n.cells][idx],lty=2)
abline(0,1,col="darkred",lwd=2) #1:1 expectation

