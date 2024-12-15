NimModel <- nimbleCode({
  #Density covariates, marked and unmarked have separate intercepts
  D0.M ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  D0.UM ~ dunif(0,100)
  D0 <- D0.M + D0.UM #total D0
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  
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

  #Density model
  D.intercept.M <- D0.M*cellArea
  D.intercept.UM <- D0.UM*cellArea
  # D.intercept <- exp(D.beta0)*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells])
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda.N.M <- D.intercept.M*pi.denom #Expected N, marked
  lambda.N.UM <- D.intercept.UM*pi.denom #Expected N, unmarked
  N.M ~ dpois(lambda.N.M) #realized marked N in state space
  N.UM ~ dpois(lambda.N.UM) #realized unmarked N in state space
  N <- N.M + N.UM #total realized N
  
  #Marked individuals first
  for(i in 1:M1) {
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]],InSS=InSS[s.cell[i]])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.mID[i,1:J] ~ dPoissonVector(lam[i,1:J]*K1D[1:J]*theta.marked[1],z=z[i]) #marked and identified detections
  }#custom Metropolis-Hastings update for N.M/z[1:M1] 
  
  #Then unmarked individuals
  for(i in (M1+1):M.both){
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]],InSS=InSS[s.cell[i]])
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