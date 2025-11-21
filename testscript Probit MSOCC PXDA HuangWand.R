#This version uses the Huang-Wand hierarchical prior for the precision matrix 

library(mvtnorm)
library(nimble)
library(coda)
source("NimbleModel Probit MSOCC PXDA HuangWand.R")

#Data dimensions (definitions below may be specific to eDNA interpretation)
S <- 10 #species
J <- 500 #sites
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
#uncertainty/mixing will get worse as cumulative p is lower. Dorazio notes bias in large correlations when cumulative p is low
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
w.init <- matrix(0, S, J)
w.init[z.init==0] <- -abs(rnorm(sum(z.init==0),0,1))
w.init[z.init==1] <- abs(rnorm(sum(z.init==1),0,1))

#Following Dorazio using glm to set B inits. modify if you add covariates
#if you don't provide decent B inits, they are drawn from prior, which is currently
#very diffuse and can take a long time to converge

B.init <- rep(NA,S)
z.obs <- 1*(y>0)
for(s in 1:S){
  zvec = z.obs[s,]
  fit = glm(zvec ~ 1, family=binomial(link='probit'))
  B.init[s] = fit$coefficients
}

Niminits <- list(B=B.init,z=z.init,w=w.init)
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

#remove a samplers, replace with custom conjugate samplers that nimble didn't recognize
conf$removeSampler("a")
calcNodes <- Rmodel$getDependencies(paste0("a[1:",S,"]"))
calcNodes <- c(calcNodes,Rmodel$expandNodeNames("Sigma"),Rmodel$expandNodeNames("R"))
conf$addSampler(target = paste0("a[1:",S,"]"), type = "aConjugateSampler",
                control = list(S=S,calcNodes=calcNodes))

#could replace default RW samplers for B with slice samplers, not sure which is more efficient
# conf$removeSamplers("B")
# for(i in 1:S){
#   conf$addSampler(target = paste("B[",i,"]"),
#                   type = 'slice',control = list(adaptive=TRUE),silent = TRUE)
# }

#Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
#can ignore warnings about NAs in A
Cmcmc$run(5000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)

burnin <- 1000

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
burnin2 <- 1000

#can look at R and w posteriors
plot(coda::mcmc(mvSamples2[-c(1:burnin2),idx.R[-seq(1,S*S,S+1)]])) #removing diagonals that are all 1
#w's that are constrained >0 by observed occupancy states
plot(coda::mcmc(mvSamples2[-c(1:burnin2),idx.w[which(y>0)]]))
#unconstrainted w's
plot(coda::mcmc(mvSamples2[-c(1:burnin2),idx.w[which(y==0)]]))

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