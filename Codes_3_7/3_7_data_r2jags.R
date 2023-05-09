rm(list = ls())
library(R2jags)

load(file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated data/long.data_1.RData")
load(file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated data/surv.data_1.RData")

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

W <- model.matrix(~ surv.data$w1) # Fixed effects

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
  n.iter = 1000, model.file =  "model_file"
)

print(sim)

##### Results
print(sim)
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 1000 iterations (first 500 discarded)
n.sims = 1500 iterations saved
mu.vect sd.vect      2.5%       25%       50%       75%     97.5%  Rhat n.eff
Sigma[1,1]     1.055   0.080     0.909     1.002     1.049     1.106     1.224 1.001  1500
Sigma[2,1]     0.531   0.068     0.401     0.487     0.530     0.574     0.672 1.035    63
Sigma[1,2]     0.531   0.068     0.401     0.487     0.530     0.574     0.672 1.035    63
Sigma[2,2]     0.897   0.099     0.729     0.824     0.890     0.958     1.110 1.028    83
alpha[1]      -1.153   0.087    -1.333    -1.205    -1.151    -1.097    -0.984 1.037    75
alpha[2]       1.116   0.086     0.954     1.055     1.112     1.176     1.291 1.012   170
beta[1]        0.461   0.100     0.259     0.395     0.463     0.528     0.645 1.252    12
beta[2]       -0.556   0.093    -0.705    -0.627    -0.558    -0.499    -0.343 2.004     5
beta[3]        1.047   0.105     0.851     0.978     1.040     1.120     1.248 1.369     9
beta[4]        0.965   0.063     0.855     0.914     0.962     1.015     1.071 1.699     6
gamma[1]      -0.403   0.130    -0.678    -0.482    -0.406    -0.326    -0.141 1.243    13
gamma[2]       0.151   0.185    -0.193     0.028     0.145     0.268     0.544 1.339    10
nu             0.880   0.049     0.777     0.848     0.881     0.916     0.971 1.062    42
sigma          1.002   0.019     0.964     0.988     1.001     1.015     1.041 1.004   450
deviance   19194.182  53.521 19086.294 19160.559 19193.764 19232.081 19295.876 1.028    98

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1404.7 and DIC = 20598.9
DIC is an estimate of expected predictive error (lower deviance is better).
