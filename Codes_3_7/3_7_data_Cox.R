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
    
    # Survival sub-model
    # Hazard function 
    haz[i]<- nu*pow(st[i],nu-1)*exp(inprod(alpha[],W[i,])+inprod(b[i,],gamma[]))
    #Cumulative hazard function 
    chaz[i]<- nu*pow(st[i],nu)*exp(inprod(alpha[],W[i,])+inprod(b[i,],gamma[]))
    #Log-survival function with Gauss-Kronrod/Legendre requadrature
    logSurv[i]<- -chaz[i]
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
    for(l in 1:2){ 
    gamma[l]~dnorm(0,0.001)
  }
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
  n = n, M = M, st = st, W = W, Y = Y, X = X, Z = Z,
  death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
  Nbeta = dim(X)[3], Nalpha = ncol(W)
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
  n.iter = 5000, model.file = "model_file"
)

print(sim)

##### Results
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect         2.5%          25%          50%           75%         97.5%  Rhat n.eff
Sigma[1,1]  1.04100e+00   0.086  8.83000e-01  9.82000e-01  1.03600e+00         1.096         1.223 1.001  2200
Sigma[2,1]  5.44000e-01   0.065  4.22000e-01  5.01000e-01  5.42000e-01         0.586         0.677 1.001  2200
Sigma[1,2]  5.44000e-01   0.065  4.22000e-01  5.01000e-01  5.42000e-01         0.586         0.677 1.001  2200
Sigma[2,2]  8.69000e-01   0.085  7.18000e-01  8.08000e-01  8.64000e-01         0.926         1.052 1.017   120
alpha[1]    8.23000e+00   1.912  5.12800e+00  6.64200e+00  8.00500e+00         9.734        11.973 2.573     4
alpha[2]    4.99000e-01   0.066  3.74000e-01  4.54000e-01  4.98000e-01         0.544         0.633 1.001  3800
beta[1]     4.95000e-01   0.080  3.34000e-01  4.42000e-01  4.97000e-01         0.549         0.642 1.064    41
beta[2]    -5.73000e-01   0.064 -7.06000e-01 -6.13000e-01 -5.73000e-01        -0.530        -0.447 1.126    21
beta[3]     1.01600e+00   0.103  8.12000e-01  9.57000e-01  1.01400e+00         1.074         1.244 1.033    66
beta[4]     1.01300e+00   0.050  9.23000e-01  9.79000e-01  1.01000e+00         1.047         1.121 1.016   160
gamma[1]   -1.13000e-01   0.089 -2.96000e-01 -1.70000e-01 -1.14000e-01        -0.054         0.059 1.001  3000
gamma[2]    1.09000e-01   0.101 -8.50000e-02  4.10000e-02  1.09000e-01         0.177         0.306 1.001  3800
nu          0.00000e+00   0.001  0.00000e+00  0.00000e+00  0.00000e+00         0.001         0.003 2.707     4
sigma       1.00400e+00   0.019  9.68000e-01  9.91000e-01  1.00400e+00         1.016         1.043 1.003   810
deviance    1.00019e+08  43.897  1.00019e+08  1.00019e+08  1.00019e+08 100019068.705 100019125.757 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 958.3 and DIC = 100019997.5
DIC is an estimate of expected predictive error (lower deviance is better).

