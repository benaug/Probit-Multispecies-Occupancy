#This version uses the Huang-Wand hierarchical prior for the precision matrix 
#V2 uses: 
#1) full conditional w[,j] updates conditioned on z[,j] from Dorazio. Custom updates use R package tmvtnorm
#and it is very slow. I provide code to generate them one j at a time or all J at once and the latter is
#very marginally faster (at least with 100 sites) and should use less RAM

#2) to allow z states to change, this version uses a simple RW update for w[i,j] in cases where y[i,j]=0.
#This was just the easiest thing to do. I will change it to use a proposal that always swaps the z state,
#which will improve efficiency further.

#W proposals can fail and if so, it will print an error each time. I
#saw it fail sometimes using Dorazio's z updates that lead to poor starting
#values for w simulation when a z has been updated and the current w is inconsistent,
#but I haven't seen it so far when not using that z update where current values of w
#are always consistent with current z

library(mvtnorm)
library(tmvtnorm)
library(truncnorm)
library(nimble)
library(coda)
source("NimbleModel Probit MSOCC PXDA HuangWand V2.R")

#Data dimensions (definitions below may be specific to eDNA interpretation)
S <- 5 #species
J <- 250 #sites. set lower than the other test scripts due to longer run times
K <- rep(5,J) #detection occasions per site

set.seed(1123)

#parameter values
beta0 <- rnorm(S,0,1)
sds <- rep(2,S)
Ustar <- rlkj_corr_cholesky(n=1,eta=1,p=S) #simulate upper triangular Cholesky factor of correlation matrix
R <- t(Ustar) %*% Ustar #compute correlation matrix
D <- diag(sds) 
Sigma <- D%*%R%*%D #compute covariance matrix
p <- rep(0.25,S) #fixed detection probabilities
#uncertainty/mixing will get worse as cumulative p is lower.
1-(1-p)^K[1] #cumulative capture probs for each species for 1st trap or all traps if K is the same

#correlated site occupancy
z <- matrix(0,S,J)
eps <- w <- matrix(0,S,J)
for(j in 1:J){
  w[,j] <- rmvnorm(1,mean=beta0,sigma=Sigma)[1,]
  z[,j] <- 1*(w[,j]>0)
}

#these need to be informative, not too close to 0 or 1
rowMeans(z)

#detection process
y <- array(0,dim=c(S,J))
for(i in 1:S){
  for(j in 1:J){
    y[i,j] <- rbinom(1,K[j],z[i,j]*p[i])
  }
}

#initialize
z.init <- 1*(y>0)
#Following Dorazio using glm to set B inits. modify if you add covariates
B.init <- rep(NA,S)
z.obs <- 1*(y>0)
for(s in 1:S){
  zvec = z.obs[s,]
  fit = glm(zvec ~ 1, family=binomial(link='probit'))
  B.init[s] = fit$coefficients
}

R.init <- diag(rep(1,S))
D.init = diag(S)

w.init <- matrix(0,S,J)
for(i in 1:S){
  for(j in 1:J) {
    if(z.init[i,j]==1){
      w.init[i,j] <- rtruncnorm(1,mean=B.init[i],sd=D.init[i,i],a=0)
    }else{
      w.init[i,j] <- rtruncnorm(1,mean=B.init[i],sd=D.init[i,i],b=0)
    }
  }                   
}
#more Dorazio initializatoion code
Sigma.init <- D.init %*% R.init %*% D.init
Omega.init <- chol2inv(chol(Sigma.init))
Eps <- t(w.init - B.init)
crossProdMat <- crossprod(Eps, Eps)
nu.Sigma <- 2
scale.Sigma <- 10
a.init <- rgamma(S, shape=(nu.Sigma+S)/2, rate=(1/scale.Sigma^2) + nu.Sigma*diag(Omega.init))
Omega.init <- rWishart(1, df=nu.Sigma+S-1 + J, Sigma=chol2inv(chol(2*nu.Sigma*diag(a.init) + crossProdMat)))[,, 1]
Sigma.init <- chol2inv(chol(Omega.init))


Niminits <- list(B=B.init,z=z.init,w=w.init,a=a.init,Omega=Omega.init)
constants <- list(S=S,J=J,K=K)

z.data <- z.init
z.data[z.data==0] <- NA
A.data <- matrix(0,S,S)
diag(A.data) <- NA

Nimdata <- list(y=y,z=z.data,A=A.data)

#set parameters to monitor
parameters <- c('beta0.derived',"p")
parameters2 <- c('R',"w")

#Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters,thin=2,
                      monitors2=parameters2,thin2=2)

conf$removeSampler("w") #remove block RW sampler for w that mixes poorly
#independent w samplers that can switch z states for y[i,j]=0
for(i in 1:S){
  for(j in 1:J){
    if(y[i,j]==0){
      conf$addSampler(target = paste0("w[", i,",",j,"]"), type = "RW",
                      control = list(S=S,i=i,j=j,K=K[j]))
    }
  }
}

#full conditional w updates conditioned on z states using rtmvnorm in R that uses gibbs sampling to generate them.
burn.in.samples <- 500 #can change the number of burn in samples for rtmvnorm() proposals. no idea what is necessary,
#but 500 may be overkill. seem to converge right away.
# #this version does 1 j at a time with J R calls
# for(j in 1:J){
#   conf$addSampler(target = paste0("w[1:", S,",",j,"]"), type = "WConjugateSampler",
#                   control = list(S=S,j=j,K=K[j],burn.in.samples=burn.in.samples))
# }
#this version does all J at in one update with 1 R call
#very marginally faster, should use slightly less RAM
conf$addSampler(target = paste0("w[1:", S,",1:",J,"]"), type = "WConjugateSampler2",
                  control = list(S=S,J=J,K=K,burn.in.samples=burn.in.samples))


#remove a samplers, replace with custom conjugate samplers that nimble didn't recognize
conf$removeSampler("a")
calcNodes <- Rmodel$getDependencies(paste0("a[1:",S,"]"))
calcNodes <- c(calcNodes,Rmodel$expandNodeNames("Sigma"),
               Rmodel$expandNodeNames("R"),
               Rmodel$expandNodeNames("w"))
conf$addSampler(target = paste0("a[1:",S,"]"), type = "aConjugateSampler",
                control = list(S=S,calcNodes=calcNodes))

#remove B samplers, replace with custom conjugate samplers that nimble didn't recognize
#mixes similarly to RW updates, but are faster to compute
conf$removeSampler("B")
conf$addSampler(target = paste0("B[1:",S,"]"), type = "BConjugateSampler",
                control = list(S=S,J=J))

#Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
#can ignore warnings about NAs in A
Cmcmc$run(2000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)

burnin <- 250

#plot derived betas, detection parameters
plot(coda::mcmc(mvSamples[-c(1:burnin),])) #might take a while for all parameters to converge
round(beta0/sds,2) #targets for beta0.derived
p

#plot mean occupancy
rowMeans(z) #average occupancy targets
idx.bd <- grep("beta0.derived",colnames(mvSamples))
plot(coda::mcmc(pnorm(mvSamples[-c(1:burnin),idx.bd])))


mvSamples2 <- as.matrix(Cmcmc$mvSamples2)
idx.R <- grep("R",colnames(mvSamples2))
idx.w <- grep("w",colnames(mvSamples2))
burnin2 <- 250

#can look at R and w posteriors
plot(coda::mcmc(mvSamples2[-c(1:burnin2),idx.R[-seq(1,S*S,S+1)]])) #removing diagonals that are all 1
#w's that are constrained >0 by observed occupancy states
plot(coda::mcmc(mvSamples2[-c(1:burnin2),idx.w[which(y>0)]]))
#unconstrainted w's
plot(coda::mcmc(mvSamples2[-c(1:burnin2),idx.w[which(y==0)]]))
#z states (still says w on the plots, ignore)
plot(coda::mcmc(1*(mvSamples2[-c(1:burnin2),idx.w[which(y==0)]]>0)))

#get point and interval estimates for correlations, compare to truth
ests2 <- matrix(colMeans(mvSamples2[-c(1:burnin2),idx.R]),S,S)
HPDs2 <- array(HPDinterval(mcmc(mvSamples2[-c(1:burnin2),idx.R])),dim=c(S,S,2))

par(mfrow=c(1,1),ask=FALSE)
plot(NA,xlim=c(-1,1),ylim=c(-1,1),pch=16,xlab="True Species Correlations",ylab="Estimated Species Correlations")
for(s1 in 1:S){
  if(s1<S){
    for(s2 in (s1+1):S){
      points(ests2[s1,s2]~R[s1,s2],pch=16)
      lines(x=c(R[s1,s2],R[s1,s2]),y=c(HPDs2[s1,s2,1],HPDs2[s1,s2,2]))
    }
  }
}
abline(0,1)

#coverage of R upper triangular matrix
covers <- HPDs2[,,1]<R&HPDs2[,,2]>R
mean(covers[upper.tri(covers)]) 
