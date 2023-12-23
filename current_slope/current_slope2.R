rm(list = ls())
library(R2jags)
library(nnet)
library(statmod)
library(splines2)

load(file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated_data/long.data_2.RData")
load(file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated_data/surv.data_2.RData")

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


# Gauss-Legendre quadrature (15 points)
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes # Nodes
wk <- glq$weights # Weights
K <- length(xk) # K-points


########  BUGS code  ########
sink("model_file")
cat("model{
  for(i in 1:n){
    # Longitudinal sub-model
    for(j in 1:M[i]){
      Y[i,j]~dnorm(mu[i,j],tau)
      mu[i,j]<-inprod(beta[],X[i,j,])+inprod(b[i,],Z[i,j,])
    }
    Alpha0[i]<- gamma[1]*(beta[1]+beta[3]*X1[i]+beta[4]*X2[i]+b[i,1])+gamma[2]*(beta[2]+b[i,2])
    Alpha1[i]<- gamma[1]*(beta[2]+b[i,2])

    # Survival sub-model
    # Hazard function
    haz[i]<- kappa*pow(st[i],kappa-1)*exp(inprod(alpha[],W[i,])+Alpha0[i]+Alpha1[i]*st[i])
    # Cumulative hazard function
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature 
      xk1[i,j]<-(xk[j]+1)/2*st[i] 
      wk1[i,j]<- wk[j]*st[i]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[i,j]<- kappa*pow(xk1[i,j],kappa-1)*exp(inprod(alpha[],W[i,])+Alpha0[i]+Alpha1[i]*xk1[i,j])
    }
    
    # Log-survival function
    logSurv[i]<- -inprod(wk1[i,],chaz[i,])
    # Definition of the survival log-likelihood using zeros trick
    phi[i]<-100000-death[i]*log(haz[i])-logSurv[i]
    zeros[i]~dpois(phi[i])
    # Random effects
    b[i,1:Nb]~dmnorm(mub[],Omega[,])
  }
  # Prior distributions
  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
}

 for(l in 1:2){ 
  gamma[l]~dnorm(0,0.001)
 }
  kappa~dgamma(0.01,0.01)
  tau~dgamma(0.01,0.01)
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  # Determinetic nodes
  sigma<-1/tau
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])

}", fill = TRUE)
sink()
#### Running JAGS #######

d.jags <- list(
  n = n, M = M, st = st, W = W, Y = Y, X = X, Z = Z, K = K, X1 = X1, X2 = X2,
  death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
  Nbeta = dim(X)[3], Nalpha = dim(W)[2], xk = xk, wk = wk
)

i.jags <- function() {
  list(
    gamma = rnorm(2), alpha = rnorm(dim(W)[2]),
    beta = rnorm(dim(X)[3]), tau = runif(1), Omega = diag(runif(2)), kappa = 1
  )
}

parameters <- c("alpha", "gamma", "beta", "sigma", "Sigma", "kappa")

sim <- jags(
  data = d.jags, inits = i.jags, parameters,
  n.iter = 5000, model.file = "model_file"
)

print(sim)

### Results
print(sim)
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.026   0.098         0.851         0.956         1.020         1.091         1.238 1.009   240
Sigma[2,1]         0.445   0.095         0.250         0.384         0.448         0.509         0.624 1.007   640
Sigma[1,2]         0.445   0.095         0.250         0.384         0.448         0.509         0.624 1.007   640
Sigma[2,2]         1.209   0.216         0.870         1.055         1.177         1.335         1.704 1.014   170
alpha[1]          -1.590   0.191        -1.998        -1.705        -1.583        -1.459        -1.242 1.001  3000
alpha[2]           0.973   0.144         0.701         0.876         0.967         1.069         1.263 1.002  3800
alpha[3]          -1.081   0.260        -1.608        -1.251        -1.075        -0.902        -0.599 1.004   560
beta[1]           -0.413   0.090        -0.578        -0.476        -0.415        -0.351        -0.235 1.082    29
beta[2]            0.563   0.105         0.343         0.495         0.565         0.634         0.765 1.035    81
beta[3]            0.528   0.057         0.414         0.490         0.530         0.567         0.639 1.016   170
beta[4]            0.444   0.110         0.233         0.368         0.439         0.521         0.655 1.030    71
gamma[1]          -0.290   0.154        -0.593        -0.391        -0.294        -0.187         0.007 1.008  1500
gamma[2]          -0.348   0.390        -1.176        -0.602        -0.326        -0.082         0.368 1.010  1600
kappa              1.752   0.171         1.444         1.635         1.740         1.854         2.119 1.002  1800
sigma              0.962   0.033         0.897         0.939         0.960         0.984         1.027 1.006   400
deviance   100007206.809  55.570 100007086.724 100007171.457 100007209.570 100007244.318 100007312.487 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1525.0 and DIC = 100008731.8
DIC is an estimate of expected predictive error (lower deviance is better).
