#this version samples unidentified betas
NimModel <- nimbleCode({
  #priors
  for(i in 1:S){
    #for RW updates
    # B[i] ~ dnorm(0,sd=10) # expanded intercept. too diffuse, mean psi can get stuck at 0 or 1
    #can use this with intercept only model to keep psi between approx 0.001 an 0.999 so it doesn't get stuck at 0 or 1
    # B[i] ~ dunif(-3.1,3.1)
    #for conjugate updates (currenly only works for intercept only)
    B[i] ~ dflat() #assumption for B conjugate sampler
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
  setup = function(model, mvSaved, target, control) {
    # Defined stuff
    S <- control$S
    calcNodes <- control$calcNodes
  },
  run = function(){
    shape <- (model$nu[1]+S)/2
    for(i in 1:S){
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