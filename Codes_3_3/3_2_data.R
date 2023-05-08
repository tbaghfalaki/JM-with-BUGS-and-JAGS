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

W <- model.matrix(~ surv.data$w1) # Fixed effects

X <- array(1, dim = c(n, max(M), 4)) # Fixed effects
X[, , 2] <- time
X[, , 3] <- X1
X[, , 4] <- X2

Z <- array(1, dim = c(n, max(M), 2)) # Random effects
Z[, , 2] <- time


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
  lambda0<-exp(alpha[1])
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
parameters <- c("alpha", "gamma", "beta", "sigma", "Sigma", "lambda0", "nu")

sim <- jags(
  data = d.jags, inits = i.jags, parameters,
  n.iter = 5000, model.file = "model_file"
)

print(sim)


##### Results ##### 
print(sim)
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.114   0.125         0.883         1.030         1.108         1.193         1.374 1.002  1400
Sigma[2,1]         0.308   0.211        -0.065         0.172         0.307         0.432         0.779 1.082    31
Sigma[1,2]         0.308   0.211        -0.065         0.172         0.307         0.432         0.779 1.082    31
Sigma[2,2]         0.971   0.487         0.285         0.624         0.818         1.285         2.149 1.127    23
alpha[1]           2.904   0.164         2.595         2.787         2.902         3.018         3.218 1.009   250
alpha[2]           0.041   0.124        -0.197        -0.045         0.041         0.122         0.283 1.002  1800
beta[1]            0.440   0.139         0.144         0.350         0.446         0.536         0.699 1.032    94
beta[2]           -0.329   0.289        -0.902        -0.525        -0.325        -0.127         0.223 1.094    26
beta[3]            1.057   0.107         0.849         0.985         1.054         1.128         1.278 1.056    50
beta[4]            1.181   0.125         0.936         1.097         1.182         1.266         1.422 1.004  1300
gamma             -0.569   0.058        -0.687        -0.607        -0.568        -0.529        -0.459 1.006   400
lambda0           18.490   3.073        13.401        16.236        18.212        20.456        24.976 1.009   250
nu                 1.908   0.076         1.760         1.856         1.907         1.959         2.056 1.003   950
sigma              0.983   0.066         0.855         0.938         0.982         1.027         1.115 1.006   390
deviance   100002414.091  49.401 100002316.426 100002381.472 100002414.206 100002447.242 100002509.580 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1199.2 and DIC = 100003613.3
DIC is an estimate of expected predictive error (lower deviance is better).
