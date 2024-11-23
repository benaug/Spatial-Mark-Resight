NimModel <- nimbleCode({
  #detection function priors
  lam0 ~ dunif(0,15)
  theta.d ~ dunif(0,25) #careful with this prior. Too much prior mass near 0 gives very strong prior weight to high overdispersion
  sigma ~ dunif(0,10)
  #data augmentation priors for marked (1) and unmarked (2) individuals
  psi1 ~ dunif(0,1)
  psi2 ~ dunif(0,1)
  #sample type observation model priors (Dirichlet)
  alpha.marked[1] <- 1
  alpha.marked[2] <- 1
  alpha.marked[3] <- 1
  alpha.unmarked[1] <- 1
  alpha.unmarked[2] <- 1
  theta.marked[1:3] ~ ddirch(alpha.marked[1:3])
  theta.unmarked[1] <- 0
  theta.unmarked[2:3] ~ ddirch(alpha.unmarked[1:2])

  #likelihoods (except for s/z priors)
  #Marked individuals first
  for(i in 1:M1) {
    z[i] ~ dbern(psi1)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    p[i,1:J] <- theta.d/(theta.d+lam[i,1:J])
    y.full[i,1:J] ~ dNBVector(p=p[i,1:J],theta.d=theta.d*K1D[1:J],z=z[i]) #vectorized obs mod. trap op: sum of N NB RVs is NB with theta.d=N*theta.d
    #custom distribution that skips likelihood eval for the individuals currently with 0 captures.
    y.event[i,1:J,1:3] ~ dmulti2(y.full[i,1:J],prob=theta.marked[1:3],capcounts=capcounts[i])
  }

  #Then unmarked individuals
  for(i in (M1+1):M.both){
    z[i] ~ dbern(psi2)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    p[i,1:J] <- theta.d/(theta.d+lam[i,1:J])
    y.full[i,1:J] ~ dNBVector(p=p[i,1:J],theta.d=theta.d*K1D[1:J],z=z[i]) #vectorized obs mod. trap op: sum of N NB RVs is NB with theta.d=N*theta.d
    #custom distribution that skips likelihood eval for the individuals currently with 0 captures.
    y.event[i,1:J,2:3] ~ dmulti2(y.full[i,1:J],prob=theta.unmarked[2:3],capcounts=capcounts[i])
  }
  
  #If you have telemetry
  for(i in 1:n.tel.inds){
    for(m in 1:n.locs.ind[i]){
      locs[tel.inds[i],m,1] ~ dnorm(s[tel.inds[i],1],sd=sigma)
      locs[tel.inds[i],m,2] ~ dnorm(s[tel.inds[i],2],sd=sigma)
    }
  }
  
  #calculate number of marked and unmarked inds captured and abundance
  capcounts[1:M.both] <- Getcapcounts(y.full=y.full[1:M.both,1:J])
  n.M <- Getncap(capcounts=capcounts[1:M1],ID=ID[1:n.samples])
  n.UM <- Getncap(capcounts=capcounts[(M1+1):M.both],ID=ID[1:n.samples])
  N.M <- sum(z[1:M1])
  N.UM <- sum(z[(M1+1):M.both])
  N.tot<- N.M + N.UM
})# end model
