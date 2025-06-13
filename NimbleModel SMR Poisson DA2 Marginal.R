NimModel <- nimbleCode({
  #detection function priors
  lam0 ~ dunif(0,15)
  sigma ~ dunif(0,10)
  
  #sample type observation model priors (Dirichlet)
  alpha.marked[1] <- 1
  alpha.marked[2] <- 1
  alpha.marked[3] <- 1
  alpha.unmarked[1] <- 1
  alpha.unmarked[2] <- 1
  theta.marked[1:3] ~ ddirch(alpha.marked[1:3])
  theta.unmarked[1] <- 0
  theta.unmarked[2:3] ~ ddirch(alpha.unmarked[1:2])
  
  lambda.N.M ~ dunif(0,1000) #Expected N, marked
  lambda.N.UM ~ dunif(0,1000) #Expected N, unmarked
  N.M ~ dpois(lambda.N.M) #realized N in state space
  N.UM ~ dpois(lambda.N.UM) #realized N in state space
  N <- N.M + N.UM #total realized N
  
  #Marked individuals first
  for(i in 1:M1) {
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.mID[i,1:J] ~ dPoissonVector(lam[i,1:J]*K1D[1:J]*theta.marked[1],z=z[i]) #marked and identified detections
  }#custom Metropolis-Hastings update for N.M/z[1:M1] 
  
  #Then unmarked individuals
  for(i in (M1+1):M.both){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
  }#custom Metropolis-Hastings update for N.UM/z[(M1+1):M.both]
  
  #Unidentified detections by type
  #1) marked with no ID detections
  bigLam.marked[1:J] <- GetbigLam(lam=lam[1:M1,1:J],z=z[1:M1])
  lam.mnoID[1:J] <- bigLam.marked[1:J]*K1D[1:J]*theta.marked[2]
  y.mnoID[1:J] ~ dPoissonVector(lam.mnoID[1:J],z=1) #plug in z=1 to reuse dPoissonVector
  
  #2) unmarked detections
  bigLam.unmarked[1:J] <- GetbigLam(lam=lam[(M1+1):M.both,1:J],z=z[(M1+1):M.both])
  lam.um[1:J] <- bigLam.unmarked[1:J]*K1D[1:J]*theta.unmarked[2]
  y.um[1:J] ~ dPoissonVector(lam.um[1:J],z=1) #plug in z=1 to reuse dPoissonVector
  
  #3) unknown marked status
  lam.unk[1:J] <- bigLam.marked[1:J]*K1D[1:J]*theta.marked[3] + bigLam.unmarked[1:J]*K1D[1:J]*theta.unmarked[3]
  y.unk[1:J] ~ dPoissonVector(lam.unk[1:J],z=1) #plug in z=1 to reuse dPoissonVector
  
  # #If you have telemetry
  # for(i in 1:n.tel.inds){
  #   for(m in 1:n.locs.ind[i]){
  #     locs[tel.inds[i],m,1] ~ dnorm(s[tel.inds[i],1],sd=sigma)
  #     locs[tel.inds[i],m,2] ~ dnorm(s[tel.inds[i],2],sd=sigma)
  #   }
  # }
})# end model