rm(list = ls())
library(R2jags)
library(nnet)
library(statmod)
library(splines2)

load(file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated_data/long.data_1.RData")
load(file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated_data/surv.data_1.RData")

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
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect          2.5%           25%          50%          75%         97.5%  Rhat n.eff
Sigma[1,1]  1.03500e+00   0.104         0.849         0.963  1.03000e+00  1.10000e+00         1.258 1.010   240
Sigma[2,1]  4.92000e-01   0.111         0.295         0.416  4.83000e-01  5.61000e-01         0.735 1.029    75
Sigma[1,2]  4.92000e-01   0.111         0.295         0.416  4.83000e-01  5.61000e-01         0.735 1.029    75
Sigma[2,2]  7.39000e-01   0.155         0.483         0.627  7.26000e-01  8.34000e-01         1.082 1.093    30
alpha[1]    7.99000e-01   0.156         0.485         0.699  7.94000e-01  9.07000e-01         1.102 1.036    92
alpha[2]    5.58000e-01   0.068         0.421         0.514  5.59000e-01  6.04000e-01         0.690 1.005   420
alpha[3]   -5.58000e-01   0.118        -0.793        -0.636 -5.55000e-01 -4.79000e-01        -0.328 1.007   310
beta[1]    -8.24000e-01   0.095        -1.013        -0.890 -8.21000e-01 -7.58000e-01        -0.641 1.046    50
beta[2]     5.80000e-01   0.146         0.296         0.468  5.98000e-01  6.87000e-01         0.846 1.226    14
beta[3]     9.58000e-01   0.120         0.735         0.875  9.51000e-01  1.04100e+00         1.197 1.039    57
beta[4]     4.16000e-01   0.059         0.300         0.373  4.16000e-01  4.56000e-01         0.530 1.013   170
gamma[1]   -5.57000e-01   0.078        -0.721        -0.607 -5.54000e-01 -5.04000e-01        -0.416 1.038    76
gamma[2]   -2.85000e-01   0.226        -0.714        -0.442 -2.80000e-01 -1.32000e-01         0.146 1.072    43
kappa       8.92000e-01   0.046         0.807         0.861  8.91000e-01  9.21000e-01         0.986 1.031    75
sigma       1.01100e+00   0.038         0.937         0.985  1.01000e+00  1.03700e+00         1.089 1.011   200
deviance    1.00006e+08  51.230 100005865.494 100005937.842  1.00006e+08  1.00006e+08 100006064.092 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1288.3 and DIC = 100007259.7
DIC is an estimate of expected predictive error (lower deviance is better).
