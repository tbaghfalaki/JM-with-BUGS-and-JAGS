rm(list = ls())
library(R2jags)
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


##### Results ##### 
> print(sim)
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.028   0.100         0.844         0.959         1.024         1.093         1.239 1.006   760
Sigma[2,1]         0.444   0.094         0.259         0.379         0.441         0.509         0.627 1.007   390
Sigma[1,2]         0.444   0.094         0.259         0.379         0.441         0.509         0.627 1.007   390
Sigma[2,2]         1.001   0.157         0.717         0.893         0.991         1.099         1.334 1.042    53
alpha[1]          -2.037   0.149        -2.347        -2.133        -2.033        -1.938        -1.752 1.004   550
alpha[2]           0.819   0.120         0.585         0.738         0.819         0.899         1.056 1.001  3800
beta[1]           -0.518   0.073        -0.665        -0.565        -0.519        -0.469        -0.375 1.025   140
beta[2]            0.560   0.091         0.380         0.499         0.559         0.621         0.735 1.074    32
beta[3]            0.442   0.055         0.325         0.408         0.443         0.478         0.549 1.053    44
beta[4]            0.710   0.115         0.475         0.633         0.710         0.793         0.919 1.070    40
gamma             -0.554   0.092        -0.733        -0.615        -0.553        -0.490        -0.372 1.006   390
lambda0            0.132   0.020         0.096         0.118         0.131         0.144         0.173 1.004   550
sigma              0.958   0.033         0.895         0.936         0.957         0.980         1.024 1.003   700
deviance   100007221.781  48.638 100007131.817 100007188.790 100007220.720 100007254.233 100007319.222 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1171.7 and DIC = 100008393.5
DIC is an estimate of expected predictive error (lower deviance is better).
