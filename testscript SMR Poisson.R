#This is an SMR data simulator and MCMC sampler that handles all sample types
#1) marked, known ID
#2) marked, unknown ID
#3) unmarked, unknown ID
#4) unknown marked status, unknown ID
#It handles both "premarked" scenarios where you know the number of marked individuals
#and "natural" where marks are from natural patterns on the animals so the number of marked 
#individuals is unknown. For "premarked", consider using the random thinning model with partial ID covariates
#y[i,j,k] ~ Poisson(lam[i,j,k])
#y.event[i,j,k,1:3]~Multinomial(theta.marked[1:3],y[i,j,k]) for marked i
#y.event[i,j,k,1:3]~Multinomial(theta.unmarked[1:3],y[i,j,k]) for unmarked i
#event 1 is you know the ID
#event 2 is you know the mark status, but not ID
#event 3 is you don't know mark status or ID

library(nimble)
source("sim.SMR.R")
source("NimbleModel SMR Poisson.R")
source("NimbleFunctions SMR Poisson.R")
source("init.SMR.R")
source("sSampler.R")

#make sure to run this line!
nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)
nimbleOptions('MCMCjointlySamplePredictiveBranches') 

####Simulate some data####
N=78
n.marked=12
lam0=0.5
sigma=0.5
K=10 #number of occasions
buff=3 #state space buffer
X<- expand.grid(3:11,3:11) #make a trapping array
#theta is probability of observing each sample type for marked and unmarked individuals
theta.marked=c(0.75,0.15,0.1) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked=0.75 #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)
marktype="premarked" #are individuals premarked, or naturally marked?
# marktype="natural"
tlocs=0 #number of telemetry locs/marked individual. For "premarked"
data=sim.SMR(N=N,n.marked=n.marked,marktype=marktype,
             theta.marked=theta.marked,theta.unmarked=theta.unmarked,
             lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,tlocs=tlocs)

#What is the observed data?
head(data$samp.type) #vector of each samp.type of each detection. Must be in this order (I think).
table(data$samp.type)
head(data$this.j) #trap of capture for each sample
head(data$this.k) #occasion of each capture for each sample (not used in this "2D" sampler)
head(data$ID.marked) #true ID's for marked and identified samples ()
head(data$locs) #possibly telemetry

####Fit model in Nimble####
if(marktype=="natural"){
  M1=40 #Augmentation level for marked.
}else{
  M1=n.marked #Set to n.marked if premarked. psi1 will be estimated, but can be ignored.
}
M2=125 #Augmentation level for unmarked
#Monitor N.M and N.UM, marked and unmarked ind abundance to make sure N.M does not hit M1
#and N.UM does not hit M1+M2 during sampling. If so, raise the offending M and run again.
M.both=M1+M2
#Need some inits to initialize data
#Use reasonable inits for lam0 and sigma since we check to make sure initial observation
#model likelihood is finite
#also use this function checks to make sure theta.marked and theta.unmarked inits are in
#the correct structure. 
inits=list(lam0=lam0,sigma=sigma)

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild=init.SMR(data,inits,M1=M1,M2=M2,marktype=marktype)

#inits for nimble
#full inits. Nimble can initialize psi1 and psi2, but if sigma and lam0 initialized too far away
#from truth, it can stop adapting before convergence and mix very poorly.
Niminits <- list(z=nimbuild$z,s=nimbuild$s,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.full),
                 y.full=nimbuild$y.full,y.event=nimbuild$y.event,
                 theta.unmarked=c(0,0.5,0.5),lam0=inits$lam0,sigma=inits$sigma)

#constants for Nimble
J=nrow(data$X)
constants<-list(M1=M1,M2=M2,M.both=M.both,J=J,K=K,K1D=data$K1D,n.samples=nimbuild$n.samples,
                xlim=data$xlim,ylim=data$ylim)

# Supply data to Nimble. Note, y.true and y.true.event are treated as completely latent (but known IDs enforced)
z.data=c(rep(1,data$n.marked),rep(NA,M.both-data$n.marked))

Nimdata<-list(y.full=matrix(NA,nrow=M.both,ncol=J),y.event=array(NA,c(M.both,J,3)),
              ID=rep(NA,nimbuild$n.samples),z=z.data,X=as.matrix(X),capcounts=rep(NA,M.both))

# set parameters to monitor
parameters<-c('psi1','psi2','lam0','sigma','theta.marked','theta.unmarked',
              'n.M','n.UM','N.M','N.UM','N.tot')

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1, useConjugacy = TRUE)

###One *required* sampler replacement
##Here, we remove the default samplers for y.full and y.event, which are not correct
#and replace it with the custom "IDSampler"
conf$removeSampler("y.full")
conf$removeSampler("y.event")
conf$addSampler(target = paste0("y.full[1:",M.both,",1:",J,"]"),
                type = 'IDSampler',control = list(M1=M1,M2=M2,M.both=M.both,J=J,K1D=data$K1D,n.fixed=nimbuild$n.fixed,
                                                  samp.type=nimbuild$samp.type,n.samples=nimbuild$n.samples,
                                                  this.j=nimbuild$this.j,match=nimbuild$match),
                silent = TRUE)

###Two *optional* sampler replacements:

#replace default activity center sampler that updates x and y locations separately with a joint update
#a little more efficient. sSampler below only tunes s when z=1. Should better tune activity centers for 
#uncaptured individuals
conf$removeSampler(paste("s[1:",M.both,", 1:2]", sep=""))
for(i in 1:M.both){
  # conf$addSampler(target = paste("s[",i,", 1:2]", sep=""), #do not adapt covariance bc s's not deterministically linked to unmarked individuals
  #                 type = 'RW_block',control=list(adaptive=TRUE,adaptScaleOnly=TRUE,adaptInterval=500),silent = TRUE)
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,scale=0.25),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

#replace independent lam0 and sigma samplers with block sampler better accommodating for posterior covariance
#should improve mixing and increase posterior effective sample size. AF_slice works better than block RW. 
#Need to not use this update or modify it when using lam0 or sigma covariates.
conf$removeSampler(c("lam0","sigma"))
conf$addSampler(target = c("lam0","sigma"),type = 'AF_slice',
                control = list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2<-Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

data$n.M #true number of captured marked individuals
data$n.UM #true number of captured unmarked individuals

#Important! If N.UM hits M2 during sampling, raise M2. 
#For an unknown number of marked individuals, if N.M hits M1 during sampling, raise M1.
