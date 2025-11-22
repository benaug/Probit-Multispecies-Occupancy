#this version has conjugate samplers for w and RW updates for z/w
NimModel <- nimbleCode({
  #priors
  for(i in 1:S){
    #for RW updates
    B[i] ~ dnorm(0,sd=10) # expanded intercept. too diffuse, mean psi can get stuck at 0 or 1
    #can use this with intercept only model to keep psi between approx 0.001 an 0.999 so it doesn't get stuck at 0 or 1
    # B[i] ~ dunif(-3.1,3.1)
    #for conjugate updates (currenly only works for intercept only)
    # B[i] ~ dflat() #assumption for B conjugate sampler
    beta0.derived[i] <- B[i]/sqrt(Sigma[i,i]) #identifiable intercept
    p[i] ~ dbeta(1, 1) #detection probability
  }
  #Huang–Wand prior
  nu <- 2  # for uniform marginal on r_{j,k} Dont change this!
  for(i in 1:S){
    st[i] <- 10  # large for uninformative half-t on sqrt(sigma_{j,j})
    a[i] ~ dgamma(shape=0.5,rate =1/st[i]^2)  #half-t on sqrt(σ_{j,j})
    A[i,i] <- a[i] # A = diag(a1, ..., aS)
  }
  Rscale[1:S,1:S] <- inverse(2*nu * A[1:S,1:S])
  Omega[1:S,1:S] ~ dwish(R=Rscale[1:S,1:S],df=nu+S-1) #precision matrix
  for(j in 1:J){
    w[1:S,j] ~ dmnorm(B[1:S],prec=Omega[1:S,1:S])
  }
  #Derive Correlation matrix
  Sigma[1:S,1:S] <- inverse(Omega[1:S,1:S]) #covariance matrix
  for(i in 1:S){
    for(k in 1:S){
      R[i,k] <- Sigma[i,k]/(sqrt(Sigma[i,i])*sqrt(Sigma[k,k]))
    }
  }
  #occupancy states and detection process
  for(i in 1:S){
    for(j in 1:J){
      z[i,j] <- step(w[i,j]) #0 if w < 0, 1 if w > 0
      y[i,j] ~ dbinom(p[i]*z[i,j],K[j])
    }
  }
})

#using this for simulation
uppertri_mult_diag <- nimbleFunction(
  run = function(mat = double(2),vec = double(1)) {
    returnType(double(2))
    p <- length(vec)
    out <- matrix(nrow = p,ncol = p,init = FALSE)
    for(i in 1:p)
      out[ ,i] <- mat[ ,i] * vec[i]
    return(out)
  })

#conjugate sampler for a=diag(A). Stolen from Dorazio et al.
aConjugateSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    # Defined stuff
    S <- control$S
    calcNodes <- control$calcNodes
  },
  run = function(){
    for(i in 1:S){
      shape <- (model$nu[1]+S)/2
      rate <- 1/model$st[i]^2+model$nu[1]*model$Omega[i,i]
      model$a[i] <<- rgamma(1,shape=shape,rate=rate) #full conditionals
    }
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#conjugate sampler B, intercept only
BConjugateSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    # Defined stuff
    S <- control$S
    J <- control$J
    calcNodes <- model$getDependencies(target)
  },
  run = function(){
    B.mean <- rep(NA,S)
    for(i in 1:S){
      B.mean[i] <- mean(model$w[i,])
    }
    #can't figure out how to get MVN RNG to work, using noncentered parameterization instead
    chol.sigma <- chol(model$Sigma/J) 
    B.prop <- B.mean + chol.sigma %*% rnorm(S,0,1)
    z <- rnorm(S,0,1)
    model$B <<- (B.mean + chol.sigma %*% z)[,1]
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#nimble function to call rtmvnorm, used in WConjugateSampler1 below
rtmvnormNim <- nimbleRcall(function(n=integer(0),mean=double(1),sigma=double(2),lower=double(1),upper=double(1),algorithm=character(0),
                                    start.value=double(1),burn.in.samples=double(0)){},
                           Rfun = 'rtmvnorm',
                           returnType = double(1),
                           where = .GlobalEnv)

#Conjugate sampler for w[,j]|z[,j]
WConjugateSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    # Defined stuff
    S <- control$S
    j <- control$j
    K <- control$K
    burn.in.samples <- control$burn.in.samples
    #since we never change z states in this sampler, we don't need downstream nodes for detection
    calcNodes <- model$expandNodeNames(paste0("w[1:", S,",",j,"]"))
  },
  run = function(){
    #update w[,j]|z[,j] from full conditional
    wStart <- model$w[,j]
    wLower <- rep(-Inf,S)
    wUpper <- rep(Inf,S)
    for(i in 1:S){
      if(model$z[i,j]==1){
        wLower[i] <- 0
      }else{
        wUpper[i] <- 0
      }
    }
    W.prop <- rtmvnormNim(1, mean=model$B, sigma=model$Sigma, lower=wLower, upper=wUpper,
                          algorithm='gibbs',start.value=wStart,
                          burn.in.samples=burn.in.samples)
    if(all(!is.nan(W.prop))){  # check if sampling of truncated multivariate normal distribution failed)
      model$w[,j] <<- W.prop
    }else{
      #keep current w, pretend this never happened
      print("W prop failed")
    }
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#nimble function to call rtmvnormJ below, used in WConjugateSampler2 below
#apologies for poor file names, something was going wrong with 2 functions that start with "rtmvnormNim"
#probably something to do with nimble recognizing functions starting with "r" as RNGs.
allJrtmvnormNim <- nimbleRcall(function(n=integer(0),mean=double(1),sigma=double(2),lower=double(2),upper=double(2),algorithm=character(0),
                                    start.value=double(2),burn.in.samples=double(0)){},
                           Rfun = 'allJrtmvnorm',
                           returnType = double(2),
                           where = .GlobalEnv)
# 
#R function to call rtmvnorm J times and return output
allJrtmvnorm <- function(n=integer(0),mean=double(2),sigma=double(2),lower=double(2),upper=double(2),algorithm=character(0),
                      start.value=double(2),burn.in.samples=double(0)){
  S <- nrow(start.value)
  J <- ncol(start.value)
  W.prop <- matrix(0,S,J)
  for(j in 1:J){
    W.prop[,j] <- rtmvnorm(1, mean=mean, sigma=sigma, lower=lower[,j], upper=upper[,j],
                           algorithm='gibbs',start.value=start.value[,j],
                           burn.in.samples=burn.in.samples)
  }
  return(W.prop)
}

#Conjugate sampler for w[,j]|z[,j] - all J together
WConjugateSampler2 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    # Defined stuff
    S <- control$S
    J <- control$J
    K <- control$K
    burn.in.samples <- control$burn.in.samples
    #since we never change z states in this sampler, we don't need downstream nodes for detection
    calcNodes <- model$expandNodeNames(paste0("w[1:", S,",1:",J,"]"))
  },
  run = function(){
    #update w[,j]|z[,j] from full conditional for all J in one function and 1 R call
    #get w upper and lower for all J
    wLower <- matrix(-Inf,S,J)
    wUpper <- matrix(Inf,S,J)
    for(j in 1:J){
      for(i in 1:S){
        if(model$z[i,j]==1){
          wLower[i,j] <- 0
        }else{
          wUpper[i,j] <- 0
        }
      }
    }
    #get w.prop from R for all j
    w.prop <- allJrtmvnormNim(n=1, mean=model$B, sigma=model$Sigma, lower=wLower, upper=wUpper,
                          algorithm='gibbs',start.value=model$w,#start from current values
                          burn.in.samples=burn.in.samples)
    for(j in 1:J){
      if(all(!is.nan(w.prop[,j]))){ # check if sampling of truncated multivariate normal distribution failed
        model$w[,j] <<- w.prop[,j]
      }else{
        #keep current w, pretend this never happened
        print("w prop failed")
      }
    }
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)


