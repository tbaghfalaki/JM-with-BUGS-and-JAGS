rm(list = ls())
library(R2jags)
library(nnet)
library(statmod)
library(splines2)

load(file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated data/long.data_2.RData")
load(file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated data/surv.data_2.RData")

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
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.136   0.123         0.915         1.051         1.131         1.211         1.398 1.028    76
Sigma[2,1]         0.290   0.158         0.050         0.154         0.282         0.408         0.612 2.011     5
Sigma[1,2]         0.290   0.158         0.050         0.154         0.282         0.408         0.612 2.011     5
Sigma[2,2]         0.815   0.477         0.131         0.263         0.903         1.206         1.591 5.329     3
alpha[1]           4.171   1.281         2.371         3.131         3.969         5.069         7.271 2.101     5
alpha[2]          -0.238   0.359        -1.091        -0.450        -0.237        -0.031         0.514 1.100    49
beta[1]            0.443   0.141         0.180         0.344         0.437         0.540         0.728 1.029    75
beta[2]           -0.919   0.452        -1.600        -1.284        -1.032        -0.459        -0.124 3.877     3
beta[3]            1.073   0.113         0.843         0.998         1.076         1.150         1.289 1.037   150
beta[4]            1.147   0.120         0.913         1.066         1.146         1.228         1.380 1.042    56
gamma[1]          -0.746   0.300        -1.533        -0.896        -0.711        -0.535        -0.280 1.369    10
gamma[2]          -3.486   1.645        -6.667        -4.887        -3.641        -1.952        -0.662 4.557     3
kappa              4.725   1.454         2.303         3.633         4.254         6.301         7.322 3.057     4
sigma              1.019   0.067         0.894         0.973         1.016         1.062         1.158 1.018   120
deviance   100001679.425 274.425 100001211.331 100001417.926 100001744.000 100001858.953 100002243.944 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 11840.8 and DIC = 100013520.2
DIC is an estimate of expected predictive error (lower deviance is better).

