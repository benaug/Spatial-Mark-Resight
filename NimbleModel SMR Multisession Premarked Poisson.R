NimModel <- nimbleCode({
  #detection function priors - shared across sessions
  lam0.fixed~dunif(0,15)
  sigma.fixed~dunif(0,10)
  #Expected density for marked + unmarked individuals
  D ~ dunif(0,10) #Expected density
  for(g in 1:N.session){
    #plug in shared df parameter for each session. Must use lam0[g] and sigma[g] here for custom update.
    #alternatively, can be estimated separately or with random effects.
    lam0[g] <- lam0.fixed
    sigma[g] <- sigma.fixed
    lambda[g] <- D*area[g] #expected N
    N[g] ~ dpois(lambda[g]) #realized N
    N.UM[g] <- N[g] - M1[g]
    
    #sample type observation model priors (Dirichlet)
    #If sharing across sessions, model as fixed above and plug theta.marked and theta.unmarked
    #into g indices here
    alpha.marked[g,1] <- 1
    alpha.marked[g,2] <- 1
    alpha.marked[g,3] <- 1
    alpha.unmarked[g,1] <- 1
    alpha.unmarked[g,2] <- 1
    theta.marked[g,1:3]~ddirch(alpha.marked[g,1:3])
    theta.unmarked[g,1] <- 0
    theta.unmarked[g,2:3]~ddirch(alpha.unmarked[g,1:2])
    
    #likelihoods (mostly)
    #Marked individuals first
    for(i in 1:M1[g]) {
      s[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      lam[g,i,1:J[g]] <- GetDetectionRate(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma[g], lam0=lam0[g], z=z[g,i])
      y.full[g,i,1:J[g]] ~ dPoissonVector(lam[g,i,1:J[g]]*K1D[g,1:J[g]],z=z[g,i]) #vectorized obs mod
      #custom distribution that skips likelihood eval for the individuals currently with 0 captures.
      y.event[g,i,1:J[g],1:3] ~ dmulti2(y.full[g,i,1:J[g]],prob=theta.marked[g,1:3],capcounts=capcounts[g,i])
    }
    
    #Then unmarked individuals
    for(i in (M1[g]+1):M.both[g]){
      s[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      lam[g,i,1:J[g]] <- GetDetectionRate(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma[g], lam0=lam0[g], z=z[g,i])
      y.full[g,i,1:J[g]] ~ dPoissonVector(lam[g,i,1:J[g]]*K1D[g,1:J[g]],z=z[g,i]) #vectorized obs mod
      #custom distribution that skips likelihood eval for the individuals currently with 0 captures.
      y.event[g,i,1:J[g],2:3] ~ dmulti2(y.full[g,i,1:J[g]],prob=theta.unmarked[g,2:3],capcounts=capcounts[g,i])
    }
    
    # #If you have telemetry
    # for(i in 1:n.tel.inds[g]){
    #   for(m in 1:n.locs.ind[g,i]){
    #     locs[g,tel.inds[g,i],m,1]~dnorm(s[g,tel.inds[g,i],1],sd=sigma[g])
    #     locs[g,tel.inds[g,i],m,2]~dnorm(s[g,tel.inds[g,i],2],sd=sigma[g])
    #   }
    # }
    
    #calculate number of marked and unmarked inds captured and abundance
    capcounts[g,1:M.both[g]] <- Getcapcounts(y.full=y.full[g,1:M.both[g],1:J[g]])
    n.M[g] <- Getncap(capcounts=capcounts[g,1:M1[g]],ID=ID[g,1:n.samples[g]])
    n.UM[g] <- Getncap(capcounts=capcounts[g,(M1[g]+1):M.both[g]],ID=ID[g,1:n.samples[g]])
  }
})# end model
