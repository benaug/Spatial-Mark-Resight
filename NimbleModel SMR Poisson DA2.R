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
    y.full[i,1:J] ~ dPoissonVector(lam[i,1:J]*K1D[1:J],z=z[i]) #vectorized obs mod
    #custom distribution that skips likelihood eval for the individuals currently with 0 captures.
    y.event[i,1:J,1:3] ~ dmulti2(y.full[i,1:J],prob=theta.marked[1:3],capcounts=capcounts[i])
  }
  #Then unmarked individuals
  for(i in (M1+1):M.both){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.full[i,1:J] ~ dPoissonVector(lam[i,1:J]*K1D[1:J],z=z[i]) #vectorized obs mod
    #custom distribution that skips likelihood eval for the individuals currently with 0 captures.
    y.event[i,1:J,2:3] ~ dmulti2(y.full[i,1:J],prob=theta.unmarked[2:3],capcounts=capcounts[i])
  }
  
  # #If you have telemetry
  # for(i in 1:n.tel.inds){
  #   for(m in 1:n.locs.ind[i]){
  #     locs[tel.inds[i],m,1] ~ dnorm(s[tel.inds[i],1],sd=sigma)
  #     locs[tel.inds[i],m,2] ~ dnorm(s[tel.inds[i],2],sd=sigma)
  #   }
  # }
  
  #calculate number of marked and unmarked inds captured and abundance
  capcounts[1:M.both] <- Getcapcounts(y.full=y.full[1:M.both,1:J])
  n.M <- Getncap(capcounts=capcounts[1:M1],ID=ID[1:n.samples])
  n.UM <- Getncap(capcounts=capcounts[(M1+1):M.both],ID=ID[1:n.samples])
})# end model
