rm(list = ls())
library(R2jags)
library(statmod)

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

########  Gauss-Legendre quadrature (15 points)  ########
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
    Alpha0[i]<- gamma*(beta[1]+beta[3]*X1[i]+beta[4]*X2[i]+b[i,1])
    Alpha1[i]<- gamma*(beta[2]+b[i,2])
    # Survival sub-model
    # Hazard function 
    haz[i]<- nu*pow(st[i],nu-1)*exp(inprod(alpha[],W[i,])+
                                            (Alpha0[i]+Alpha1[i]*st[i]))
    #Cumulative hazard function 
    
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature 
      xk1[i,j]<-(xk[j]+1)/2*st[i] 
      wk1[i,j]<- wk[j]*st[i]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[i,j]<- nu*pow(xk1[i,j],nu-1)*exp(inprod(alpha[],W[i,])+
                                                    (Alpha0[i]+Alpha1[i]*xk1[i,j]))
    }
    
    #Log-survival function with Gauss-Kronrod/Legendre requadrature
    logSurv[i]<- -inprod(wk1[i,],chaz[i,])
    #Definition of the survival log-likelihood using zeros trick
    phi[i]<-100000-death[i]*log(haz[i])-logSurv[i]
    zeros[i]~dpois(phi[i])
    #Random effects
    b[i,1:Nb]~dmnorm(mub[],Omega[,])
  }
  #Prior distributions
  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  gamma~dnorm(0,0.001)
  nu~dgamma(0.01,0.01)
  tau~dgamma(0.01,0.01)
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  #Determinetic nodes
  sigma<-1/tau
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])
}", fill = TRUE)
sink()
#### Running JAGS
d.jags <- list(
  n = n, M = M, st = st, W = W, Y = Y, X = X, Z = Z, X1 = X1, X2 = X2,
  death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
  Nbeta = dim(X)[3], Nalpha = ncol(W),xk=xk,wk=wk,K=K
)

i.jags <- function() {
  list(
    alpha = rnorm(ncol(W)), gamma = rnorm(1),
    beta = rnorm(dim(X)[3]), tau = runif(1), Omega = diag(runif(2))
  )
}
parameters <- c("alpha", "gamma", "beta", "sigma", "Sigma", "nu")

sim <- jags(
  data = d.jags, inits = i.jags, parameters,
  n.iter = 2000, model.file = "model_file"
)

print(sim)


##### Results ##### 
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 2000 iterations (first 1000 discarded)
n.sims = 3000 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.002   0.093         0.827         0.940         0.998         1.060         1.198 1.019   110
Sigma[2,1]         0.475   0.089         0.302         0.415         0.476         0.534         0.658 1.089    29
Sigma[1,2]         0.475   0.089         0.302         0.415         0.476         0.534         0.658 1.089    29
Sigma[2,2]         1.115   0.186         0.771         0.988         1.101         1.227         1.517 1.024   140
alpha[1]          -1.666   0.154        -1.991        -1.763        -1.664        -1.564        -1.381 1.005  1100
alpha[2]           0.953   0.137         0.689         0.864         0.952         1.040         1.235 1.002  1700
alpha[3]          -1.030   0.256        -1.562        -1.196        -1.024        -0.866        -0.540 1.001  3000
beta[1]           -0.432   0.088        -0.607        -0.489        -0.433        -0.372        -0.263 1.008   500
beta[2]            0.582   0.081         0.432         0.525         0.580         0.639         0.740 1.004   600
beta[3]            0.530   0.055         0.410         0.496         0.531         0.568         0.630 1.042    85
beta[4]            0.451   0.109         0.242         0.375         0.458         0.523         0.661 1.006  1900
gamma             -0.413   0.071        -0.556        -0.461        -0.410        -0.364        -0.281 1.003   740
nu                 1.723   0.148         1.450         1.622         1.718         1.818         2.040 1.001  3000
sigma              0.965   0.034         0.900         0.941         0.964         0.987         1.034 1.003   860
deviance   100007233.312  48.283 100007140.250 100007201.101 100007232.610 100007265.061 100007329.575 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1160.0 and DIC = 100008393.3
DIC is an estimate of expected predictive error (lower deviance is better).
