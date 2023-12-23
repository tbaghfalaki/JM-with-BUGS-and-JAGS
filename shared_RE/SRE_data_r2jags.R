rm(list = ls())
library(R2jags)

load(file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated_data/long.data_sre_1.RData")
load(file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated_data/surv.data_sre_1.RData")

# Number of patients and number of longitudinal observations per patient
n <- length(surv.data$id)
M <- table(long.data$id)
# Survival and censoring times
st <- surv.data$survtime
death <- surv.data$death # death=1 means observed

# Longitudinal information in matrix format
time <- Y <- X1 <- X2 <- matrix(NA, n, max(M))
for (i in 1:n) {
  time[i, 1:M[i]] <- long.data$obstime[long.data$id == i]
  Y[i, 1:M[i]] <- long.data$Y1[long.data$id == i]
  X1[i, 1:M[i]] <- long.data$x1[long.data$id == i]
  X2[i, 1:M[i]] <- long.data$x2[long.data$id == i]
}

X1 <- X1[, 1]
X2 <- X2[, 1]

W <- model.matrix(~ surv.data$w1+surv.data$w2) # Fixed effects

X <- array(1, dim = c(n, max(M), 4)) # Fixed effects
X[, , 2] <- time
X[, , 3] <- X1
X[, , 4] <- X2

Z <- array(1, dim = c(n, max(M), 2)) # Random effects
Z[, , 2] <- time
########  BUGS code  ########
sink("model_file")
cat("model{
  for(i in 1:n){
    #Longitudinal sub-model
    for(j in 1:M[i]){
      Y[i,j]~dnorm(mu[i,j],tau)
      mu[i,j]<-inprod(beta[],X[i,j,])+inprod(b[i,],Z[i,j,])
    }
    t[i] ~ dweib(nu,mut[i]) 
    log(mut[i])<-inprod(alpha[],W[i,])+inprod(b[i,],gamma[])
    
    is.censored[i]~dinterval(t[i],c[i])
    
    b[i,1:Nb]~dmnorm(mub[],Omega[,])
  }
  #Prior distributions
  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  for(l in 1:2){ 
    gamma[l]~dnorm(0,0.001)
  }
  nu~dgamma(0.01,0.01)
  tau~dgamma(0.01,0.01)
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  #Determinetic nodes
  lambda0<-exp(alpha[1])
  sigma<-1/tau
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])
}", fill = TRUE)
sink()
#### Running JAGS
is.censored=1-death
t=c=rep(0,n)
for(i in 1:n){
  if(death[i]==0)((t[i]=NA) & (c[i]=st[i]))
  if(death[i]==1)((t[i]=st[i]) & (c[i]=max(st)))
}


d.jags <- list(
  n = n, M = M, W = W, Y = Y, X = X, Z = Z,
  mub = rep(0, 2), V = diag(1, 2), Nb = 2,
  Nbeta = dim(X)[3], Nalpha = ncol(W), is.censored=is.censored,c=c,t=t
)


i.jags <- function() {
  list(
    alpha = rnorm(ncol(W)), gamma = rnorm(2),
    beta = rnorm(dim(X)[3]), tau = runif(1), Omega = diag(runif(2))
  )
}
parameters <- c("alpha", "gamma", "beta", "sigma", "Sigma", "nu")


sim <- jags(
  data = d.jags, inits = i.jags, parameters,
  n.iter = 15000, model.file =  "model_file"
)

print(sim)

##### Results
print(sim)
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 15000 iterations (first 7500 discarded), n.thin = 7
n.sims = 3213 iterations saved
mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
Sigma[1,1]    0.819   0.089    0.656    0.757    0.815    0.876    1.007 1.001  3200
Sigma[2,1]    0.562   0.090    0.387    0.501    0.562    0.622    0.735 1.008   310
Sigma[1,2]    0.562   0.090    0.387    0.501    0.562    0.622    0.735 1.008   310
Sigma[2,2]    0.839   0.156    0.573    0.730    0.822    0.933    1.186 1.043    53
alpha[1]      1.179   0.161    0.910    1.070    1.163    1.271    1.577 1.124    23
alpha[2]     -1.143   0.134   -1.466   -1.217   -1.131   -1.054   -0.919 1.124    24
alpha[3]     -0.589   0.135   -0.863   -0.675   -0.589   -0.498   -0.327 1.016   150
beta[1]      -0.529   0.081   -0.686   -0.583   -0.528   -0.475   -0.369 1.002  3200
beta[2]       0.625   0.122    0.380    0.543    0.626    0.709    0.861 1.004   670
beta[3]       0.621   0.104    0.417    0.548    0.620    0.691    0.821 1.002  1700
beta[4]       0.436   0.051    0.337    0.402    0.436    0.471    0.537 1.002  1800
gamma[1]      1.151   0.369    0.577    0.902    1.092    1.321    2.066 1.228    14
gamma[2]     -1.155   0.414   -2.193   -1.356   -1.090   -0.864   -0.542 1.371    10
nu            1.111   0.101    0.950    1.040    1.101    1.164    1.361 1.156    19
sigma         1.033   0.043    0.951    1.003    1.032    1.061    1.119 1.004   650
deviance   4921.173  70.447 4765.664 4879.351 4926.158 4969.726 5047.103 1.099    28

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 2302.0 and DIC = 7223.2
DIC is an estimate of expected predictive error (lower deviance is better).
> 