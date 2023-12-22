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
    # Longitudinal sub-model
    for(j in 1:M[i]){
      Y[i,j]~dnorm(mu[i,j],tau)
      mu[i,j]<-inprod(beta[],X[i,j,])+inprod(b[i,],Z[i,j,])
    }
    Alpha0[i]<- gamma*(beta[1]+beta[3]*X1[i]+beta[4]*X2[i]+b[i,1])
    Alpha1[i]<- gamma*(beta[2]+b[i,2])
    # Survival sub-model
    # Hazard function
    haz[i]<- exp(inprod(alpha[],W[i,])+
                   (Alpha0[i]+Alpha1[i]*st[i]))
    # Cumulative hazard function
    chaz[i]<- exp(inprod(alpha[],W[i,])+Alpha0[i])*(exp(st[i]*Alpha1[i])-1)/Alpha1[i]
    # Log-survival function
    logSurv[i]<- -chaz[i]
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
  gamma~dnorm(0,0.001)
  kappa~dgamma(0.01,0.01)
  tau~dgamma(0.01,0.01)
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  # Determinetic nodes
  lambda0<-exp(alpha[1])
  sigma<-1/tau
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])
}", fill = TRUE)
sink()
#### Running JAGS
d.jags <- list(
  n = n, M = M, st = st, W = W, Y = Y, X = X, Z = Z, X1 = X1, X2 = X2,
  death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
  Nbeta = dim(X)[3], Nalpha = ncol(W)
)

i.jags <- function() {
  list(
    alpha = rnorm(ncol(W)), gamma = rnorm(1),
    beta = rnorm(dim(X)[3]), tau = runif(1), Omega = diag(runif(2))
  )
}
parameters <- c("alpha", "gamma", "beta", "sigma", "Sigma", "lambda0")

sim <- jags(
  data = d.jags, inits = i.jags, parameters,
  n.iter = 5000, model.file = "model_file"
)

print(sim)
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.068   0.103         0.873         0.996         1.066         1.136         1.278 1.009   250
Sigma[2,1]         0.529   0.114         0.312         0.453         0.525         0.601         0.762 1.013   410
Sigma[1,2]         0.529   0.114         0.312         0.453         0.525         0.601         0.762 1.013   410
Sigma[2,2]         1.330   0.190         1.003         1.202         1.313         1.443         1.743 1.026    94
alpha[1]           0.153   0.065         0.026         0.110         0.154         0.197         0.281 1.001  2600
alpha[2]          -0.856   0.068        -0.989        -0.902        -0.856        -0.810        -0.723 1.011   200
beta[1]           -0.493   0.100        -0.682        -0.559        -0.496        -0.430        -0.293 1.027   340
beta[2]            0.553   0.108         0.352         0.474         0.546         0.628         0.771 1.014   250
beta[3]            0.472   0.122         0.215         0.396         0.476         0.554         0.696 1.028   180
beta[4]            0.455   0.056         0.343         0.416         0.455         0.495         0.563 1.005  1100
gamma             -0.947   0.058        -1.064        -0.987        -0.945        -0.906        -0.834 1.006   350
lambda0            1.168   0.076         1.027         1.116         1.166         1.218         1.325 1.001  2600
sigma              0.992   0.032         0.931         0.969         0.992         1.014         1.057 1.001  3800
deviance   100007476.564  45.956 100007385.557 100007445.697 100007476.564 100007508.032 100007564.937 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1051.7 and DIC = 100008528.3
DIC is an estimate of expected predictive error (lower deviance is better).
> 