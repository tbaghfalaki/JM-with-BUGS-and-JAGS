rm(list = ls())
library(R2jags)
library(statmod)

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

W <- model.matrix(~ surv.data$w) # Fixed effects

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
pInference for Bugs model at "model_file", fit using jags,
3 chains, each with 2000 iterations (first 1000 discarded)
n.sims = 3000 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.026   0.093         0.849         0.963         1.023         1.086         1.216 1.005   490
Sigma[2,1]         0.444   0.086         0.278         0.388         0.444         0.500         0.614 1.004   690
Sigma[1,2]         0.444   0.086         0.278         0.388         0.444         0.500         0.614 1.004   690
Sigma[2,2]         0.973   0.166         0.684         0.854         0.964         1.079         1.323 1.042    52
alpha[1]          -2.182   0.161        -2.507        -2.289        -2.173        -2.071        -1.892 1.021   110
alpha[2]           1.003   0.132         0.744         0.916         1.001         1.087         1.267 1.007   290
beta[1]           -0.538   0.075        -0.684        -0.589        -0.539        -0.487        -0.385 1.079    32
beta[2]            0.584   0.080         0.425         0.534         0.587         0.636         0.739 1.050    59
beta[3]            0.467   0.055         0.366         0.430         0.463         0.502         0.582 1.121    22
beta[4]            0.730   0.111         0.501         0.655         0.731         0.807         0.934 1.073    67
gamma             -0.633   0.088        -0.810        -0.691        -0.633        -0.572        -0.466 1.019   130
nu                 1.715   0.142         1.452         1.618         1.711         1.808         2.002 1.005   490
sigma              0.959   0.033         0.897         0.936         0.958         0.981         1.028 1.012   180
deviance   100007190.802  49.031 100007099.453 100007156.552 100007190.393 100007223.004 100007287.437 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1164.0 and DIC = 100008354.8
DIC is an estimate of expected predictive error (lower deviance is better).
