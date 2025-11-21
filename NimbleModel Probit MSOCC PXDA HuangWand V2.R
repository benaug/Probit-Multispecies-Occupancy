#this version has conjugate samplers for w and RW updates for z/w
NimModel <- nimbleCode({
  #priors
  for(i in 1:S){
    B[i] ~ dnorm(0,sd=10) # expanded intercept
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

#nimble function to call rtmvnorm, used in WConjugateSampler
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