#Without PXDA
NimModel <- nimbleCode({
  #priors
  for(i in 1:S){
    beta0[i] ~ dnorm(0,sd=10)#unidentifiable intercept
    beta0.derived[i] <- beta0[i]/sds[i] #beta0/sds, identifiable intercept
    p[i] ~ dbeta(1,1) #detection probabilities
  }
  
  #Wishart prior
  Sigma[1:S,1:S] ~ dinvwish(S=Ident[1:S,1:S], df=S+2) # weakly informative
  
  for(j in 1:J){
    mu[1:S,j] ~ dmnorm(beta0[1:S],cov=Sigma[1:S,1:S])
  }
  #get sds and R
  for(i in 1:S){
    sds[i] <- sqrt(Sigma[i,i])
  }
  for(i in 1:S){
    for(k in 1:S){
      R[i,k] <- Sigma[i,k]/(sds[i]*sds[k])
    }
  }
  #occupancy states and detection process
  for(i in 1:S){
    for(j in 1:J){
      z[i,j] <- step(mu[i,j]) #0 if mu < 0, 1 if mu > 0
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