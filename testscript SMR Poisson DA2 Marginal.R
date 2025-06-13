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

library(nimble)
source("sim.SMR.R")
source("NimbleModel SMR Poisson DA2 Marginal.R")
source("NimbleFunctions SMR Poisson DA2 Marginal.R")
source("init.SMR.R")
source("sSampler Poisson Marginal.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
N <- 78
n.marked <- 12 #for "premarked", this is the number marked. for "natural", this is the number of captured individuals that are identifiable
lam0 <- 0.5
sigma <- 0.5
K <- 10 #number of occasions
buff <- 3 #state space buffer
X<- expand.grid(3:11,3:11) #make a trapping array
#theta is probability of observing each sample type for marked and unmarked individuals
theta.marked <- c(0.75,0.15,0.1) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked <- 0.75 #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)
marktype <- "premarked" #are individuals premarked, or naturally marked?
# marktype <- "natural"
obstype <- "poisson"
tlocs <- 0 #number of telemetry locs/marked individual. For "premarked"
data <- sim.SMR(N=N,n.marked=n.marked,marktype=marktype,
             theta.marked=theta.marked,theta.unmarked=theta.unmarked,
             lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,tlocs=tlocs,
             obstype=obstype)

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

####Fit model in Nimble####
if(marktype=="natural"){
  M1 <- 40 #Augmentation level for marked.
}else{
  M1 <- n.marked #Set to n.marked if premarked. psi1 will be estimated, but can be ignored.
}
M2 <- 125 #Augmentation level for unmarked
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
nimbuild <- init.SMR(data,inits,M1=M1,M2=M2,marktype=marktype,obstype="poisson")

#inits for nimble
#must initialize N.M and N.UM to be consistent with z. speeds converge to set consistent with lambda.N.M/UM
Niminits <- list(z=nimbuild$z,s=nimbuild$s,
                 N.M=sum(nimbuild$z[1:M1]),lambda.N.M=sum(nimbuild$z[1:M1]),
                 N.UM=sum(nimbuild$z[(M1+1):M.both]),lambda.N.UM=sum(nimbuild$z[(M1+1):M.both]),
                 theta.unmarked=c(0,0.5,0.5),lam0=inits$lam0,sigma=inits$sigma)

#constants for Nimble
J <- nrow(data$X)
constants <- list(M1=M1,M.both=M.both,J=J,K1D=data$K1D,xlim=data$xlim,ylim=data$ylim)

# Supply data to Nimble. Note, y.true and y.true.event are treated as completely latent (but known IDs enforced)
z.data <- c(rep(1,data$n.marked),rep(NA,M.both-data$n.marked))

Nimdata <- list(y.mID=nimbuild$y.event[1:M1,,1], #marked with ID
                y.mnoID=colSums((nimbuild$y.event[1:M1,,2])), #marked without ID
                y.um=colSums((nimbuild$y.event[(M1+1):M.both,,2])), #unmarked
                y.unk=colSums((nimbuild$y.event[1:M.both,,3])), #unk marked status
                z=z.data,X=as.matrix(X))
#all samples accounted for
sum(Nimdata$y.mID) + sum(Nimdata$y.mnoID) + sum(Nimdata$y.um) + sum(Nimdata$y.unk)
sum(nimbuild$y.event)

#If you have telemetry use these instead. Make sure to uncomment telemetry BUGS code.
# constants <- list(M1=M1,M2=M2,M.both=M.both,J=J,K=K,K1D=data$K1D,
#                   xlim=data$xlim,ylim=data$ylim,tel.inds=nimbuild$tel.inds,
#                   n.tel.inds=length(nimbuild$tel.inds),n.locs.ind=nimbuild$n.locs.ind)
# Nimdata <- list(y.mID=nimbuild$y.event[1:M1,,1], #marked with ID
#                 y.mnoID=rowSums((nimbuild$y.event[1:M1,,2])), #marked without ID
#                 y.um=rowSums((nimbuild$y.event[(M1+1):M.both,,2])), #unmarked
#                 y.unk=rowSums((nimbuild$y.event[1:M.both,,3])), #unk marked status
#                 z=z.data,X=as.matrix(X),locs=data$locs)

# set parameters to monitor
parameters <- c('lambda.N.M','lambda.N.UM','lam0','sigma','theta.marked','theta.unmarked',
             'N.M','N.UM','N')

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
config.nodes <- c('lambda.N.M','lambda.N.UM','lam0','sigma','theta.marked','theta.unmarked[2:3]')
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1,nodes=config.nodes)


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
# y.nodes <- Rmodel$expandNodeNames(paste("y.full[1:",M.both,",1:",J,"]"))
lam.nodes <- Rmodel$expandNodeNames(paste("lam[1:",M.both,",1:",J,"]"))
y.mID.nodes <- Rmodel$expandNodeNames(paste("y.mID[1:",M1,",1:",J,"]"))
y.mnoID.nodes <- Rmodel$expandNodeNames(paste("y.mnoID[1:",J,"]"))
y.um.nodes <- Rmodel$expandNodeNames(paste("y.um[1:",J,"]"))
y.unk.nodes <- Rmodel$expandNodeNames(paste("y.unk[1:",J,"]"))
bigLam.marked.nodes <- Rmodel$expandNodeNames("bigLam.marked") #only need this in calcNodes
bigLam.unmarked.nodes <- Rmodel$expandNodeNames("bigLam.unmarked") #only need this in calcNodes
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
for(i in 1:M.both){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,
                                                 M1=M1,scale=0.25),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

#can add block sampler if lam0, sigma, lambda.N.UM, and/or lambda.N.M (if N.M unknown) posteriors highly correlated
#RW_block faster, AFslice slower, but mixes better
# conf$addSampler(target = c("lam0","sigma"),type = 'RW_block',
#                 control = list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[250:nrow(mvSamples),]))


#Important! If N.UM hits M2 during sampling, raise M2. 
#For an unknown number of marked individuals, if N.M hits M1 during sampling, raise M1.