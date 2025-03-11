# Function to calculate detection rate, but skip when z=0
GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     ans <- lam0*exp(-d2/(2*sigma^2))
     return(ans)
    }
  }
)
#Vectorized observation model that also prevents z from being turned off if an unmarked ind currently has samples.
#also skips likelihood eval when z=0
dNBVector <- nimbleFunction(
  run = function(x = double(1), p = double(1), theta.d = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dnbinom(x, p = p, size = theta.d, log = TRUE))
      return(logProb)
    }
  }
)

#dummy random vector generator to make nimble happy
rNBVector <- nimbleFunction(
  run = function(n = integer(0),p = double(1),theta.d = double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(p)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)
#custom multinomial distribution to skip calcs when an ind has 0 samples (most of them!).
dmulti2 <- nimbleFunction(
  run = function(x = double(2), size = double(1), prob = double(1), capcounts = double(0),
                 log = integer(0)) {
    returnType(double(0))
    levels <- nimDim(prob)[1]
    J <- nimDim(size)[1]
    if(capcounts==0){
      return(0)
    }else{
      logProb <- 0
      for(j in 1:J){
        if(size[j]>0){
          logProb <- logProb + dmulti(x[j,1:levels], size=size[j], prob=prob, log = TRUE)
        }
      }
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rmulti2 <- nimbleFunction(
  run = function(n=integer(0),size = double(1), prob = double(1), capcounts = double(0)) {
    returnType(double(2))
    J <- nimDim(size)[1]
    out <- matrix(J,3,value=0)
    return(out)
  }
)
#make dummy random vector generator to make nimble happy
Getcapcounts <- nimbleFunction(
  run = function(y.full=double(2)){
    returnType(double(1))
    M.both <- nimDim(y.full)[1]
    J <- nimDim(y.full)[2]
    capcounts <- numeric(M.both, value = 0)
    for(i in 1:M.both){
      capcounts[i] <- sum(y.full[i,1:J])
    }
    return(capcounts)
  }
)

#calculates how many samples each individual is currently allocated.
Getncap <- nimbleFunction(
  run = function(capcounts=double(1),ID=double(1)){ #don't need ID, but nimble requires is it used in a function 
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    nstate <- numeric(M, value = 0)
    for(i in 1:M){
      if(capcounts[i]>0){
        nstate[i] <- 1
      }
    }
    n.cap <- sum(nstate)
    return(n.cap)
  }
)

#------------------------------------------------------------------
# Customer sampler to update latent IDs, y.full, and y.event
#------------------------------------------------------------------
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M1 <- control$M1
    M2 <- control$M2
    M.both <- control$M.both
    J <- control$J
    K1D <- control$K1D
    n.fixed <- control$n.fixed
    samp.type <- control$samp.type
    match <- control$match
    n.samples <- control$n.samples
    this.j <- control$this.j
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model$z
    y.full <- model$y.full
    y.event <- model$y.event

    #precalculate log likelihoods
    ll.y <- matrix(0,nrow=M.both,ncol=J)
    ll.y.event <- matrix(0,nrow=M.both,ncol=J)
    for(i in 1:M.both){
      if(z[i]==1){
        for(j in 1:J) {
          #if theta is a function of individual or trap covariates, need to fix the line below, e.g., size=model$theta[i,j]*model$K1D[j]
          ll.y[i,j] <-  dnbinom(y.full[i,j],size=model$theta.d[1]*model$K1D[j],prob=model$p[i,j], log = TRUE)
        }
      }
    }
    for(i in 1:M1){
      if(z[i]==1){
        for(j in 1:J) {
          if(y.full[i,j]>0){
            ll.y.event[i,j] <- dmulti(y.event[i,j,1:3],y.full[i,j],model$theta.marked,log=TRUE)
          }
        }
      }
    }
    for(i in (M1+1):M.both){
      if(z[i]==1){
        for(j in 1:J) {
          if(y.full[i,j]>0){
            ll.y.event[i,j] <- dmulti(y.event[i,j,1:3],y.full[i,j],model$theta.unmarked,log=TRUE)
          }
        }
      }
    }
    ll.y.cand <- ll.y
    ll.y.event.cand <- ll.y.event
    ID.curr <- model$ID
    
    ###update IDs
    ID.cand <- ID.curr
    y.full.cand <- y.full
    y.event.cand <- y.event
    for(l in (n.fixed+1):n.samples){#for all samples without known IDs
      propprobs <- model$lam[1:M.both,this.j[l]]
      for(i in 1:M.both){ #zero out nonmatches and z=0
        if(!match[l,i] | z[i]==0){
          propprobs[i] <- 0
        }
      }
      denom <- sum(propprobs) #abort if propprobs sum to 0. No matches anywhere nearby.
      if(denom>0){
        propprobs <- propprobs/denom
        ID.cand[l] <- rcat(1,prob=propprobs)
        if(ID.cand[l]!=ID.curr[l]){
          swapped <- c(ID.curr[l],ID.cand[l])
          #new sample proposal probabilities
          forprob <- propprobs[swapped[2]]
          backprob <- propprobs[swapped[1]]
          #new y.full's - move sample from ID to ID.cand
          y.event.cand[ID.curr[l],this.j[l],samp.type[l]] <- y.event[ID.curr[l],this.j[l],samp.type[l]]-1
          y.event.cand[ID.cand[l],this.j[l],samp.type[l]] <- y.event[ID.cand[l],this.j[l],samp.type[l]]+1
          y.full.cand[ID.curr[l],this.j[l]] <- y.full[ID.curr[l],this.j[l]]-1
          y.full.cand[ID.cand[l],this.j[l]] <- y.full[ID.cand[l],this.j[l]]+1
          #if theta is a function of individual or trap covariates, need to fix these 2 lines below, e.g., size=model$theta[swapped[1],this.j[l]]*model$K1D[this.j[l]]
          ll.y.cand[swapped[1],this.j[l]] <- dnbinom(y.full.cand[swapped[1],this.j[l]],size=model$theta.d[1]*model$K1D[this.j[l]],prob=model$p[swapped[1],this.j[l]],log=TRUE)
          ll.y.cand[swapped[2],this.j[l]] <- dnbinom(y.full.cand[swapped[2],this.j[l]],size=model$theta.d[1]*model$K1D[this.j[l]],prob=model$p[swapped[2],this.j[l]],log=TRUE)
          #old ID theta likelihood
          if(swapped[1]<=M1){#marked guy
            if(y.full.cand[swapped[1],this.j[l]]==0){
              ll.y.event.cand[swapped[1],this.j[l]] <- 0
            }else{
              ll.y.event.cand[swapped[1],this.j[l]] <- dmulti(y.event.cand[swapped[1],this.j[l],1:3],
                                                           y.full.cand[swapped[1],this.j[l]],model$theta.marked,log=TRUE)
            }
          }else{#unmarked guy
            if(y.full.cand[swapped[1],this.j[l]]==0){
              ll.y.event.cand[swapped[1],this.j[l]] <- 0
            }else{
              ll.y.event.cand[swapped[1],this.j[l]] <- dmulti(y.event.cand[swapped[1],this.j[l],1:3],
                                                           y.full.cand[swapped[1],this.j[l]],model$theta.unmarked,log=TRUE)
            }
          }
          #new ID theta likelihood
          if(swapped[2]<=M1){#marked guy
            if(y.full.cand[swapped[2],this.j[l]]==0){
              ll.y.event.cand[swapped[2],this.j[l]] <- 0
            }else{
              ll.y.event.cand[swapped[2],this.j[l]] <- dmulti(y.event.cand[swapped[2],this.j[l],1:3],
                                                           y.full.cand[swapped[2],this.j[l]],model$theta.marked,log=TRUE)
            }
          }else{#unmarked guy
            if(y.full.cand[swapped[2],this.j[l]]==0){
              ll.y.event.cand[swapped[2],this.j[l]] <- 0
            }else{
              ll.y.event.cand[swapped[2],this.j[l]] <- dmulti(y.event.cand[swapped[2],this.j[l],1:3],
                                                           y.full.cand[swapped[2],this.j[l]],model$theta.unmarked,log=TRUE)
            }
          }
          #select sample to move proposal probabilities
          #P(select a sample of this type for this ID)*P(select this j|sample of this type and this ID)
          #n.samples cancels out in MH ratio. Wrong number of samples (includes known ID that are not updated), but doesn't matter.
          #Including for clarity
          focalprob <- (sum(ID.curr==swapped[1]&samp.type==samp.type[l])/n.samples)*
            y.event[swapped[1],this.j[l],samp.type[l]]/sum(y.event[swapped[1],1:J,samp.type[l]])
          focalbackprob <- (sum(ID.cand==swapped[2]&samp.type==samp.type[l])/n.samples)*
            y.event.cand[swapped[2],this.j[l],samp.type[l]]/sum(y.event.cand[swapped[2],1:J,samp.type[l]])
          
          #sum log likelihoods and do MH step
          lp_initial <- sum(ll.y[swapped,this.j[l]])+sum(ll.y.event[swapped,this.j[l]])
          lp_proposed <- sum(ll.y.cand[swapped,this.j[l]])+sum(ll.y.event.cand[swapped,this.j[l]])
          log_MH_ratio <- (lp_proposed+log(backprob)+log(focalbackprob)) - (lp_initial+log(forprob)+log(focalprob))
          accept <- decide(log_MH_ratio)
          if(accept){
            y.event[swapped[1],this.j[l],samp.type[l]] <- y.event.cand[swapped[1],this.j[l],samp.type[l]]
            y.event[swapped[2],this.j[l],samp.type[l]] <- y.event.cand[swapped[2],this.j[l],samp.type[l]]
            y.full[swapped[1],this.j[l]] <- y.full.cand[swapped[1],this.j[l]]
            y.full[swapped[2],this.j[l]] <- y.full.cand[swapped[2],this.j[l]]
            ll.y[swapped[1],this.j[l]] <- ll.y.cand[swapped[1],this.j[l]]
            ll.y[swapped[2],this.j[l]] <- ll.y.cand[swapped[2],this.j[l]]
            ll.y.event[swapped[1],this.j[l]] <- ll.y.event.cand[swapped[1],this.j[l]]
            ll.y.event[swapped[2],this.j[l]] <- ll.y.event.cand[swapped[2],this.j[l]]
            ID.curr[l] <- ID.cand[l]
          }else{
            #set these back
            y.event.cand[swapped[1],this.j[l],samp.type[l]] <- y.event[swapped[1],this.j[l],samp.type[l]]
            y.event.cand[swapped[2],this.j[l],samp.type[l]] <- y.event[swapped[2],this.j[l],samp.type[l]]
            y.full.cand[swapped[1],this.j[l]] <- y.full[swapped[1],this.j[l]]
            y.full.cand[swapped[2],this.j[l]] <- y.full[swapped[2],this.j[l]]
            ID.cand[l] <- ID.curr[l]
          }
        }
      }
    }
    #put everything back into the model$stuff
    model$y.full <<- y.full
    model$y.event <<- y.event
    # model$capcounts <<- capcounts
    model$ID <<- ID.curr
    model.lp.proposed <- model$calculate(calcNodes) #update logprob
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)
