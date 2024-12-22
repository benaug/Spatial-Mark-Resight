NimModel <- nimbleCode({
  #detection function priors - shared across sessions
  lam0.fixed ~ dunif(0,15)
  sigma.fixed ~ dunif(0,10)
  
  #sample type observation probabilities, shared over sessions.
  alpha.marked[1] <- 1
  alpha.marked[2] <- 1
  alpha.marked[3] <- 1
  alpha.unmarked[1] <- 1
  alpha.unmarked[2] <- 1
  theta.marked.fixed[1:3] ~ ddirch(alpha.marked[1:3])
  theta.unmarked.fixed[1] <- 0
  theta.unmarked.fixed[2:3] ~ ddirch(alpha.unmarked[1:2])
  
  for(g in 1:N.session){
    #not sharing Dcov parameters over sessions, but can do that
    D0.M[g] ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
    D0.UM[g] ~ dunif(0,100)
    D0[g] <- D0.M[g] + D0.UM[g] #total D0
    # D.beta0[g] ~ dnorm(0,sd=10)
    D.beta1[g] ~ dnorm(0,sd=10)
    #Density model
    D.intercept.M[g] <- D0.M[g]*cellArea[g]
    D.intercept.UM[g] <- D0.UM[g]*cellArea[g]
    lambda.cell[g,1:n.cells[g]] <- InSS[g,1:n.cells[g]]*exp(D.beta1[g]*D.cov[g,1:n.cells[g]])
    pi.cell[g,1:n.cells[g]] <- lambda.cell[g,1:n.cells[g]]/pi.denom[g] #expected proportion of total N in cell c
    pi.denom[g] <- sum(lambda.cell[g,1:n.cells[g]])
    lambda.N.M[g] <- D.intercept.M[g]*pi.denom[g] #Expected N, marked
    lambda.N.UM[g] <- D.intercept.UM[g]*pi.denom[g] #Expected N, unmarked
    N.M[g] ~ dpois(lambda.N.M[g]) #realized marked N in state space
    N.UM[g] ~ dpois(lambda.N.UM[g]) #realized unmarked N in state space
    N[g] <- N.M[g] + N.UM[g] #total realized N
    
    #plug in shared df parameter for each session. Must use lam0[g] and sigma[g] here for custom update.
    #alternatively, can be estimated separately or with random effects.
    lam0[g] <- lam0.fixed
    sigma[g] <- sigma.fixed
    
    #sample type observation model priors (Dirichlet)
    #If not shared across sessions, use this
    #into g indices here
    # alpha.marked[g,1] <- 1
    # alpha.marked[g,2] <- 1
    # alpha.marked[g,3] <- 1
    # alpha.unmarked[g,1] <- 1
    # alpha.unmarked[g,2] <- 1
    # theta.marked[g,1:3] ~ ddirch(alpha.marked[g,1:3])
    # theta.unmarked[g,1] <- 0
    # theta.unmarked[g,2:3] ~ ddirch(alpha.unmarked[g,1:2])
    #if shared over sessions use this
    theta.marked[g,1:3] <- theta.marked.fixed[1:3]
    theta.unmarked[g,1:3] <- theta.unmarked.fixed[1:3]
    
    #Marked individuals first
    for(i in 1:M1[g]){
      #dunif() here implies uniform distribution within a grid cell
      #also tells nimble s's are in continuous space, not discrete
      s[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      #get cell s_i lives in using look-up table
      s.cell[g,i] <- cells[g,trunc(s[g,i,1]/res[g])+1,trunc(s[g,i,2]/res[g])+1]
      #categorical likelihood for this cell, equivalent to zero's trick
      #also disallowing s's in non-habitat
      dummy.data[g,i] ~ dCell(pi.cell[g,s.cell[g,i]],InSS=InSS[g,s.cell[g,i]])
      lam[g,i,1:J[g]] <- GetDetectionRate(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma[g], lam0=lam0[g], z=z[g,i])
      y.mID[g,i,1:J[g]] ~ dPoissonVector(lam[g,i,1:J[g]]*K1D[g,1:J[g]]*theta.marked[g,1],z=z[g,i]) #marked and identified detections
    }
    
    #Then unmarked individuals
    for(i in (M1[g]+1):M.both[g]){
      #dunif() here implies uniform distribution within a grid cell
      #also tells nimble s's are in continuous space, not discrete
      s[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      #get cell s_i lives in using look-up table
      s.cell[g,i] <- cells[g,trunc(s[g,i,1]/res[g])+1,trunc(s[g,i,2]/res[g])+1]
      #categorical likelihood for this cell, equivalent to zero's trick
      #also disallowing s's in non-habitat
      dummy.data[g,i] ~ dCell(pi.cell[g,s.cell[g,i]],InSS=InSS[g,s.cell[g,i]])
      lam[g,i,1:J[g]] <- GetDetectionRate(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma[g], lam0=lam0[g], z=z[g,i])
    }
    
    #Unidentified detections by type
    #1) marked with no ID detections
    bigLam.marked[g,1:J[g]] <- GetbigLam(lam=lam[g,1:M1[g],1:J[g]],z=z[g,1:M1[g]])
    lam.mnoID[g,1:J[g]] <- bigLam.marked[g,1:J[g]]*K1D[g,1:J[g]]*theta.marked[g,2]
    y.mnoID[g,1:J[g]] ~ dPoissonVector(lam.mnoID[g,1:J[g]],z=1) #plug in z=1 to reuse dPoissonVector
    
    #2) unmarked detections
    bigLam.unmarked[g,1:J[g]] <- GetbigLam(lam=lam[g,(M1[g]+1):M.both[g],1:J[g]],z=z[g,(M1[g]+1):M.both[g]])
    lam.um[g,1:J[g]] <- bigLam.unmarked[g,1:J[g]]*K1D[g,1:J[g]]*theta.unmarked[g,2]
    y.um[g,1:J[g]] ~ dPoissonVector(lam.um[g,1:J[g]],z=1) #plug in z=1 to reuse dPoissonVector
    
    #3) unknown marked status
    lam.unk[g,1:J[g]] <- bigLam.marked[g,1:J[g]]*K1D[g,1:J[g]]*theta.marked[g,3] + 
      bigLam.unmarked[g,1:J[g]]*K1D[g,1:J[g]]*theta.unmarked[g,3]
    y.unk[g,1:J[g]] ~ dPoissonVector(lam.unk[g,1:J[g]],z=1) #plug in z=1 to reuse dPoissonVector
    
    #If you have telemetry
    # for(i in 1:n.tel.inds[g]){
    #   for(m in 1:n.locs.ind[g,i]){
    #     locs[g,tel.inds[g,i],m,1] ~ dnorm(s[g,tel.inds[g,i],1],sd=sigma[g])
    #     locs[g,tel.inds[g,i],m,2] ~ dnorm(s[g,tel.inds[g,i],2],sd=sigma[g])
    #   }
    # }
  }
})# end model