#This version doesn't use PX. Same as Tobler et al.

library(mvtnorm)
library(nimble)
library(coda)
source("NimbleModel Probit MSOCC DA Wishart.R")

#Data dimensions
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
#mixing will get worse as cumulative p is lower
1-(1-p)^K[1] #cumulative capture probs for each species for 1st trap or all traps if K is the same

#correlated site occupancy
z <- matrix(0,S,J)
eps <- mu <- matrix(0,S,J)
for(j in 1:J){
  mu[,j] <- rmvnorm(1,mean=beta0,sigma=Sigma)[1,]
  z[,j] <- 1*(mu[,j]>0)
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
mu.init <- matrix(0, S, J)
mu.init[z.init==0] <- -abs(rnorm(sum(z.init==0),0,1))
mu.init[z.init==1] <- abs(rnorm(sum(z.init==1),0,1))

beta0.init <- rnorm(S,0,0.5)
sds.init <- runif(S, 0.5, 3)

Niminits <- list(beta0=beta0.init,z=z.init,mu=mu.init,sds=sds.init)
constants <- list(S=S,J=J,K=K)

z.data <- z.init
z.data[z.data==0] <- NA

Nimdata <- list(y=y,z=z.data,Ident=diag(S))

#set parameters to monitor
parameters <- c('beta0.derived',"p","sds")
parameters2 <- c('R',"mu")

#Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters,thin=2,
                      monitors2=parameters2,thin2=2)

#slice samplers may be better than RW for unidentified betas? Didn't compare, either works.
conf$removeSamplers("beta0")
for(i in 1:S){
  conf$addSampler(target = paste("beta0[",i,"]"),
                  type = 'slice',control = list(adaptive=TRUE),silent = TRUE)
}

#Build and compile
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)

burnin <- 1000

#plot derived betas
plot(coda::mcmc(mvSamples[-c(1:burnin),]))
beta0/sds #targets for beta0.derived
p

#plot mean occupancy
rowMeans(z) #average occupancy targets
idx.bd <- grep("beta0.derived",colnames(mvSamples))
plot(coda::mcmc(pnorm(mvSamples[-c(1:burnin),idx.bd])))


mvSamples2 <- as.matrix(Cmcmc$mvSamples2)
idx.R <- grep("R",colnames(mvSamples2))
idx.mu <- grep("mu",colnames(mvSamples2))
burnin2 <- 1000

#can look at R and mu posteriors
plot(coda::mcmc(mvSamples2[-c(1:burnin2),idx.R[-seq(1,S*S,S+1)]])) #removing diagonals that are all 1
#mu's that are constrained >0 by observed occupancy states
plot(coda::mcmc(mvSamples2[-c(1:burnin2),idx.mu[which(y>0)]]))
#unconstrainted mu's
plot(coda::mcmc(mvSamples2[-c(1:burnin2),idx.mu[which(y==0)]]))

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
