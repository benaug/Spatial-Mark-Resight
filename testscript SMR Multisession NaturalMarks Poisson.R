#Multisession SMR using an alternative approach to data augmentation

#For naturally marked populations. The difference is that you need to estimate 
#the number of marked individuals when it's not known. In this context, "marked"
#means "identifiable by natural marks". photos where you cannot tell if they are
#identifiable by marks or not are classed as "unknown marked status" and can match
#either marked or unmarked individuals. Consider using the random thinning model with
#categorical ID covariates (e.g. coat pattern, color morph) if these correlate with
#how identifiable individuals are.

#I have no idea how well this model works with a lot of unknown marked status
#samples and very sparse capture data. Can test via simulation.


#This testscript shows how to share lam0, sigma, and/or expected density across sessions
#can share theta.marked and theta.unmarked but not done here.

#See single session test script for more notes on basic SMR model structure
#Nimble sampler won't work as is if only 1 marked individual, can be fixed, email Ben

library(nimble)
source("sim.SMR.multisession.R")
source("sim.SMR.R") #used by multisession simulator
source("NimbleModel SMR Multisession NaturalMarks Poisson.R")
source("NimbleFunctions SMR Multisession NaturalMarks Poisson.R")
source("init.SMR.multisession.R")
source("init.SMR.R") #used by multisession initializer
source("sSampler Multisession.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
#Here, I'll simulate 3 populations with different n.marked, K, X, and state space areas
#sharing D, lam0, sigma so they can be shared during estimation
N.session <- 3
D <- rep(0.4,N.session) #expected density in units of sigma and X
n.marked <- c(12,13,14)
#data simulator for natural marks is sort of hacky.
#I simulate total N[g] for session g, then select n.marked[g]
#to be naturally marked. If N[g] is simulated less than n[g],
#you'll get an error.

lam0 <- rep(0.5,N.session)
sigma <- rep(0.5,N.session)
K <- c(5,6,7) #number of occasions
buff <- rep(2,N.session) #state space buffer
#make trapping arrays
X1 <- expand.grid(3:11,3:11)
X2 <- expand.grid(3:12,3:12)
X3 <- expand.grid(3:13,3:13)
X <- list(X1,X2,X3) #put in a list, one for each session

#See what expected N is for these expected D and state space areas
area <- getArea(X=X,buff=buff)
area #state space areas for each session resulting from X and buff
lambda <- D*area
lambda #expected N in each session

#theta is probability of observing each sample type for marked and unmarked individuals
theta.marked <- matrix(rep(c(0.75,0.15,0.1),N.session),nrow=N.session,byrow=TRUE) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked <- rep(0.75,N.session) #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)
marktype <- "natural" #are individuals premarked, or naturally marked? This test script only handles natural marks.
obstype <- "poisson"
data <- sim.SMR.multisession(N.session=N.session,lambda=lambda,n.marked=n.marked,marktype=marktype,
             theta.marked=theta.marked,theta.unmarked=theta.unmarked,
             lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,tlocs=rep(0,N.session),
             obstype=obstype)

#What is the observed data?
head(t(data$samp.type)) #vector of each samp.type of each detection. Must be in this order (I think).
table(data$samp.type)
head(t(data$this.j)) #trap of capture for each sample
head(t(data$this.k)) #occasion of each capture for each sample (not used in this "2D" sampler)
head(data$ID.marked) #true ID's for marked and identified samples ()



####Fit model in Nimble####
M1 <- rep(40,N.session) #Augmentation level for marked.
M2 <- c(115,125,135) #Augmentation level for unmarked
#Monitor N.M and N.UM, marked and unmarked ind abundance to make sure N.M does not hit M1
#and N.UM does not hit M1+M2 during sampling. If so, raise the offending M and run again.
M.both <- M1 + M2
#Need some inits to initialize data
#Use reasonable inits for lam0 and sigma since we check to make sure initial observation
#model likelihood is finite
#also use this function checks to make sure theta.marked and theta.unmarked inits are in
#the correct structure. 
inits <- list(lam0=lam0,sigma=sigma) #initializing with 1 parameter per session, just set all to same value

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild <- init.SMR.multisession(data,inits,M1=M1,M2=M2,marktype=marktype,obstype="poisson")

#inits for nimble
theta.unmarked.init <- matrix(c(0,0.5,0.5),N.session,3,byrow=TRUE)

N.init <- rowSums(nimbuild$z,na.rm=TRUE)
N.UM.init <- N.UM.init <- rep(NA,N.session)
for(g in 1:N.session){
  N.UM.init[g] <- sum(nimbuild$z[g,(M1[g]+1):M.both[g]])
}
N.M.init <- N.init - N.UM.init

Niminits <- list(N=N.init,N.UM=N.UM.init,N.M=N.M.init,
                 z=nimbuild$z,s=nimbuild$s,ID=nimbuild$ID,capcounts=apply(nimbuild$y.full,c(1,2),sum,na.rm=TRUE),
                 y.full=nimbuild$y.full,y.event=nimbuild$y.event,
                 theta.unmarked=theta.unmarked.init,lam0.fixed=0.5,sigma.fixed=0.5, #one init per lam0.fixed, sigma.fixed
                 D=0.5)


#constants for Nimble
J <- unlist(lapply(data$X,nrow))
constants <- list(N.session=N.session,M1=M1,M.both=M.both,J=J,K1D=nimbuild$K1D,n.samples=nimbuild$n.samples,
                xlim=data$xlim,ylim=data$ylim,area=area)

Nimdata <- list(y.full=array(NA,dim=c(N.session,max(M.both),max(J))),y.event=array(NA,c(N.session,max(M.both),max(J),3)),
              ID=matrix(NA,N.session,max(nimbuild$n.samples)),X=nimbuild$X,capcounts=matrix(NA,N.session,max(M.both)))


# set parameters to monitor
parameters <- c('lambda','lam0.fixed','sigma.fixed','theta.marked','theta.unmarked',
              'n.M','n.UM','N.UM','N.M','N',"D")
#other things we can monitor with separate thinning rate
parameters2 <- c("ID","s")

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
#much faster if we don't let nimble configure y and y.event nodes
#if you change the model file, make sure you make necessary changes here
config.nodes <- c('lam0.fixed','sigma.fixed','theta.marked',paste("theta.unmarked[1:",N.session,",2:3","]"),"D.M","D.UM")
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1,
                      monitors2=parameters2, thin2=10,
                      nodes=config.nodes,
                      useConjugacy = TRUE)

###Two *required* sampler replacement
##Here, we remove the default samplers for y.full and y.event, which are not correct
#and replace it with the custom "IDSampler"
# conf$removeSampler("y.full")
# conf$removeSampler("y.event")
for(g in 1:N.session){
  conf$addSampler(target = paste0("y.full[",g,",1:",M.both[g],",1:",J[g],"]"),
                  type = 'IDSampler',control = list(M1=M1[g],M2=M2[g],M.both=M.both[g],J=J[g],K1D=nimbuild$K1D[g,1:J[g]],
                                                    n.fixed=nimbuild$n.fixed[g],samp.type=nimbuild$samp.type[g,1:nimbuild$n.samples[g]],
                                                    n.samples=nimbuild$n.samples[g],
                                                    this.j=nimbuild$this.j[g,1:nimbuild$n.samples[g]],
                                                    match=nimbuild$match[g,1:nimbuild$n.samples[g],1:M.both[g]],
                                                    g=g),
                  silent = TRUE)
}

# how many z proposals per iteration per session for marked, unmarked?
z.ups <- rbind(round(M1*0.25),round(M2*0.25)) #doing 25% of M1 and M2 here
J <- nimbuild$J
# conf$removeSampler("N.M")
# conf$removeSampler("N.UM")
for(g in 1:N.session){
  #nodes used for update, calcNodes + z nodes
  y.nodes <- Rmodel$expandNodeNames(paste("y.full[",g,",","1:",M.both[g],",1:",J[g],"]"))
  lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",","1:",M.both[g],",1:",J[g],"]"))
  N.M.node <- Rmodel$expandNodeNames(paste("N.M[",g,"]"))
  N.UM.node <- Rmodel$expandNodeNames(paste("N.UM[",g,"]"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",","1:",M.both[g],"]"))
  calcNodes <- c(N.M.node,N.UM.node,y.nodes,lam.nodes)

  conf$addSampler(target = paste("N.UM[",g,"]"),
                  type = 'zSampler',control = list(z.ups=z.ups[1:2,g],J=J[g],M1=M1[g],M.both=M.both[g],g=g,
                                                   y.nodes=y.nodes,lam.nodes=lam.nodes,N.M.node=N.M.node,
                                                   N.UM.node=N.UM.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),
                  silent = TRUE)
}


###Two *optional* sampler replacements:

#replace default activity center sampler that updates x and y locations separately with a joint update
#a little more efficient. sSampler below only tunes s when z=1. Should better tune activity centers for 
#uncaptured individuals
conf$removeSampler("s")
for(g in 1:N.session){
  for(i in 1:M.both[g]){
    conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,g=g,xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}

#maybe replace independent lam0 and sigma samplers with block sampler better accommodating for posterior covariance
#should improve mixing and increase posterior effective sample size. AF_slice works better than block RW. 
#Need to not use this update or modify it when using lam0 or sigma covariates.
#This sampler is slower, so not worth it if data is not so sparse there is strong posterior correlation
#between lam0 and sigma.

#mixing seems generally good if you keep independent samplers and add block sampler (if data sparse enough that there is posterior correlation)
# conf$removeSampler(c("lam0.fixed","sigma.fixed"))
conf$addSampler(target = c("lam0.fixed","sigma.fixed"),type = 'RW_block',
control = list(adaptive=TRUE),silent = TRUE)

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
plot(mcmc(mvSamples[50:nrow(mvSamples),]))

data$N #realized Ns
lambda #expected N
data$n.M #true number of captured marked individuals
data$n.UM #true number of captured unmarked individuals

#Important! 
#If N.M[g] hits M1[g] during sampling, raise M1[g].
#If N.UM[g] hits M2[g] during sampling, raise M2[g]. 