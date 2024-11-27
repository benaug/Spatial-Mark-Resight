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
dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lambda = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dpois(x, lambda = lambda, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0),lambda = double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(lambda)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)
GetbigLam <- nimbleFunction(
  run = function(lam = double(2), z = double(1)){ 
    returnType(double(1))
    M <- nimDim(lam)[1]
    J <- nimDim(lam)[2]
    bigLam <- rep(0,J)
    for(i in 1:M){
      if(z[i]==1){
        bigLam <- bigLam + lam[i,]
      }
    }
    return(bigLam)
  }
)

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    J <- control$J
    M1 <- control$M1
    M.both <- control$M.both
    z.ups <- control$z.ups
    inds.detected <- inds.detected
    y.mID.nodes <- control$y.mID.nodes
    y.mnoID.nodes <- control$y.mnoID.nodes
    y.um.nodes <- control$y.um.nodes
    y.unk.nodes <- control$y.unk.nodes
    lam.nodes <- control$lam.nodes
    lam.mnoID.nodes <- control$lam.mnoID.nodes
    lam.um.nodes <- control$lam.um.nodes
    lam.unk.nodes <- control$lam.unk.nodes
    N.M.node <- control$N.M.node
    N.UM.node <- control$N.UM.node
    z.nodes <- control$z.nodes
    updateMarked <- control$updateMarked
    calcNodes <- control$calcNodes
  },
  run = function(){
    if(updateMarked){
      bigLam.marked.initial <- model$bigLam.marked
      ###1)  Do marked individuals first###
      for(up in 1:z.ups[1]){ #how many updates per iteration?
        # propose to add/subtract 1
        updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
        if(updown==0){#subtract
          reject <- FALSE #we auto reject if you select a detected individual

          #find all z's currently on *excluding unmarked individuals*
          z.on <- which(model$z[1:M1]==1)

          n.z.on <- length(z.on)
          if(n.z.on>0){ #skip if no marked z's to turn off, otherwise nimble will crash
            pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
            pick <- z.on[pick]

            #prereject turning off detected individuals
            if(any(pick==inds.detected)){ #is this individual detected?
              reject <- TRUE
            }
            if(!reject){
              #get initial logprobs for N and y
              lp.initial.N <- model$getLogProb(N.M.node)
              lp.initial.y.mID <- model$getLogProb(y.mID.nodes[pick])
              lp.initial.y.mnoID <- model$getLogProb(y.mnoID.nodes)
              lp.initial.y.unk <- model$getLogProb(y.unk.nodes)

              #propose new N/z
              model$N.M[1] <<-  model$N.M[1] - 1
              model$N[1] <<-  model$N[1] - 1
              model$z[pick] <<- 0
              
              #turn off
              bigLam.marked.proposed <- bigLam.marked.initial - model$lam[pick,] #subtract these out before calculate
              model$calculate(lam.nodes[pick])
              model$bigLam.marked <<- bigLam.marked.proposed
              model$calculate(lam.mnoID.nodes)
              model$calculate(lam.unk.nodes)
              
              #get proposed logprobs for N and y
              lp.proposed.N <- model$calculate(N.M.node)
              lp.proposed.y.mID <- model$calculate(y.mID.nodes[pick]) #will always be 0
              lp.proposed.y.mnoID <- model$calculate(y.mnoID.nodes)
              lp.proposed.y.unk <- model$calculate(y.unk.nodes)

              #MH step
              log_MH_ratio <- (lp.proposed.N + lp.proposed.y.mID + lp.proposed.y.mnoID + lp.proposed.y.unk) -
                (lp.initial.N + lp.initial.y.mID + lp.initial.y.mnoID + lp.initial.y.unk)
              accept <- decide(log_MH_ratio)
              if(accept) {
                mvSaved["N",1][1] <<- model[["N"]]
                mvSaved["N.M",1][1] <<- model[["N.M"]]
                mvSaved["z",1][pick] <<- model[["z"]][pick]
                mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
                mvSaved["bigLam.marked",1][1:J] <<- model[["bigLam.marked"]][1:J]
                mvSaved["lam.mnoID",1][1:J] <<- model[["lam.mnoID"]][1:J]
                mvSaved["lam.unk",1][1:J] <<- model[["lam.unk"]][1:J]
                bigLam.marked.initial <- bigLam.marked.proposed
              }else{
                model[["N"]] <<- mvSaved["N",1][1]
                model[["N.M"]] <<- mvSaved["N.M",1][1]
                model[["z"]][pick] <<- mvSaved["z",1][pick]
                model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
                model[["bigLam.marked"]][1:J] <<- mvSaved["bigLam.marked",1][1:J]
                model[["lam.mnoID"]][1:J] <<- mvSaved["lam.mnoID",1][1:J]
                model[["lam.unk"]][1:J] <<- mvSaved["lam.unk",1][1:J]
                model$calculate(y.mID.nodes[pick])
                model$calculate(y.mnoID.nodes)
                model$calculate(y.unk.nodes)
                model$calculate(N.M.node)
              }
            }
          }
        }else{ #add
          if(model$N.M[1] < M1){ #cannot update if z maxed out. Need to raise M1
            #find all z's currently off. marked guys only
            z.off <- which(model$z[1:M1]==0)

            n.z.off <- length(z.off)
            pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
            pick <- z.off[pick]

            #get initial logprobs for N and y
            lp.initial.N <- model$getLogProb(N.M.node)
            lp.initial.y.mID <- model$getLogProb(y.mID.nodes[pick])
            lp.initial.y.mnoID <- model$getLogProb(y.mnoID.nodes)
            lp.initial.y.unk <- model$getLogProb(y.unk.nodes)

            #propose new N/z
            model$N.M[1] <<-  model$N.M[1] + 1
            model$N[1] <<-  model$N[1] + 1
            model$z[pick] <<- 1

            #turn on
            model$calculate(lam.nodes[pick])
            bigLam.marked.proposed <- bigLam.marked.initial + model$lam[pick,] #add these in after calculate
            model$bigLam.marked <<- bigLam.marked.proposed
            model$calculate(lam.mnoID.nodes)
            model$calculate(lam.unk.nodes)

            #get proposed logprobs for N and y
            lp.proposed.N <- model$calculate(N.M.node)
            lp.proposed.y.mID <- model$calculate(y.mID.nodes[pick])
            lp.proposed.y.mnoID <- model$calculate(y.mnoID.nodes)
            lp.proposed.y.unk <- model$calculate(y.unk.nodes)

            #MH step
            log_MH_ratio <- (lp.proposed.N + lp.proposed.y.mID + lp.proposed.y.mnoID + lp.proposed.y.unk) -
              (lp.initial.N + lp.initial.y.mID + lp.initial.y.mnoID + lp.initial.y.unk)
            accept <- decide(log_MH_ratio)
            if(accept){
              mvSaved["N",1][1] <<- model[["N"]]
              mvSaved["N.M",1][1] <<- model[["N.M"]]
              mvSaved["z",1][pick] <<- model[["z"]][pick]
              mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
              mvSaved["bigLam.marked",1][1:J] <<- model[["bigLam.marked"]][1:J]
              mvSaved["lam.mnoID",1][1:J] <<- model[["lam.mnoID"]][1:J]
              mvSaved["lam.unk",1][1:J] <<- model[["lam.unk"]][1:J]
              bigLam.marked.initial <- bigLam.marked.proposed
            }else{
              model[["N"]] <<- mvSaved["N",1][1]
              model[["N.M"]] <<- mvSaved["N.M",1][1]
              model[["z"]][pick] <<- mvSaved["z",1][pick]
              model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
              model[["bigLam.marked"]][1:J] <<- mvSaved["bigLam.marked",1][1:J]
              model[["lam.mnoID"]][1:J] <<- mvSaved["lam.mnoID",1][1:J]
              model[["lam.unk"]][1:J] <<- mvSaved["lam.unk",1][1:J]
              model$calculate(y.mID.nodes[pick])
              model$calculate(y.mnoID.nodes)
              model$calculate(y.unk.nodes)
              model$calculate(N.M.node)
            }
          }
        }
      }
    }
    ###2)  Do unmarked individuals second###
    bigLam.unmarked.initial <- model$bigLam.unmarked
    for(up in 1:z.ups[2]){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      if(updown==0){#subtract
        reject <- FALSE #we auto reject if you select a detected individual

        #find all z's currently on *excluding marked individuals*
        z.on2 <- which(model$z[(M1+1):M.both]==1) + M1 #cannot reuse objects created from R functions, (z.on2, instead of z.on)

        n.z.on <- length(z.on2)
        if(n.z.on>0){ #skip if no unmarked z's to turn off, otherwise nimble will crash
          pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
          pick <- z.on2[pick]

          #prereject turning off all unmarked individuals
          if(model$N.UM[1]==1){ #is this the last unmarked individual?
            reject <- TRUE
          }
          if(!reject){
            #get initial logprobs for N and y
            lp.initial.N <- model$getLogProb(N.UM.node)
            lp.initial.y.um <- model$getLogProb(y.um.nodes)
            lp.initial.y.unk <- model$getLogProb(y.unk.nodes)

            #propose new N/z
            model$N.UM[1] <<-  model$N.UM[1] - 1
            model$N[1] <<-  model$N[1] - 1
            model$z[pick] <<- 0

            #turn off
            bigLam.unmarked.proposed <- bigLam.unmarked.initial - model$lam[pick,] #subtract these out before calculate
            model$calculate(lam.nodes[pick])
            model$bigLam.unmarked <<- bigLam.unmarked.proposed
            model$calculate(lam.um.nodes)
            model$calculate(lam.unk.nodes)

            #get proposed logprobs for N and y
            lp.proposed.N <- model$calculate(N.UM.node)
            lp.proposed.y.um <- model$calculate(y.um.nodes)
            lp.proposed.y.unk <- model$calculate(y.unk.nodes)

            #MH step
            log_MH_ratio <- (lp.proposed.N + lp.proposed.y.um + lp.proposed.y.unk) -
              (lp.initial.N + lp.initial.y.um + lp.initial.y.unk)
            accept <- decide(log_MH_ratio)
            if(accept) {
              mvSaved["N",1][1] <<- model[["N"]]
              mvSaved["N.UM",1][1] <<- model[["N.UM"]]
              mvSaved["z",1][pick] <<- model[["z"]][pick]
              mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
              mvSaved["bigLam.unmarked",1][1:J] <<- model[["bigLam.unmarked"]][1:J]
              mvSaved["lam.um",1][1:J] <<- model[["lam.um"]][1:J]
              mvSaved["lam.unk",1][1:J] <<- model[["lam.unk"]][1:J]
              bigLam.unmarked.initial <- bigLam.unmarked.proposed
            }else{
              model[["N"]] <<- mvSaved["N",1][1]
              model[["N.UM"]] <<- mvSaved["N.UM",1][1]
              model[["z"]][pick] <<- mvSaved["z",1][pick]
              model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
              model[["bigLam.unmarked"]][1:J] <<- mvSaved["bigLam.unmarked",1][1:J]
              model[["lam.um"]][1:J] <<- mvSaved["lam.um",1][1:J]
              model[["lam.unk"]][1:J] <<- mvSaved["lam.unk",1][1:J]
              model$calculate(y.um.nodes)
              model$calculate(y.unk.nodes)
              model$calculate(N.UM.node)
            }
          }
        }
      }else{#add
        if(model$N.UM[1] < M.both-M1 ){ #cannot update if z maxed out. Need to raise M2

          #find all z's currently off. Marked inds excluded here bc always on.
          z.off2 <- which(model$z[(M1+1):M.both]==0) + M1

          n.z.off <- length(z.off2)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off2[pick]

          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.UM.node)
          lp.initial.y.um <- model$getLogProb(y.um.nodes)
          lp.initial.y.unk <- model$getLogProb(y.unk.nodes)

          #propose new N/z
          model$N.UM[1] <<-  model$N.UM[1] + 1
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1

          #turn on
          model$calculate(lam.nodes[pick])
          bigLam.unmarked.proposed <- bigLam.unmarked.initial + model$lam[pick,] #add these in after calculate
          model$bigLam.unmarked <<- bigLam.unmarked.proposed
          model$calculate(lam.um.nodes)
          model$calculate(lam.unk.nodes)

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.UM.node)
          lp.proposed.y.um <- model$calculate(y.um.nodes)
          lp.proposed.y.unk <- model$calculate(y.unk.nodes)

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y.um + lp.proposed.y.unk) -
            (lp.initial.N + lp.initial.y.um + lp.initial.y.unk)
          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["N.UM",1][1] <<- model[["N.UM"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
            mvSaved["bigLam.unmarked",1][1:J] <<- model[["bigLam.unmarked"]][1:J]
            mvSaved["lam.um",1][1:J] <<- model[["lam.um"]][1:J]
            mvSaved["lam.unk",1][1:J] <<- model[["lam.unk"]][1:J]
            bigLam.unmarked.initial <- bigLam.unmarked.proposed
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["N.UM"]] <<- mvSaved["N.UM",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
            model[["bigLam.unmarked"]][1:J] <<- mvSaved["bigLam.unmarked",1][1:J]
            model[["lam.um"]][1:J] <<- mvSaved["lam.um",1][1:J]
            model[["lam.unk"]][1:J] <<- mvSaved["lam.unk",1][1:J]
            model$calculate(y.um.nodes)
            model$calculate(y.unk.nodes)
            model$calculate(N.UM.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)