rm(list = ls())
library(R2jags)

load(file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated_data/long.data_sre_1.RData")
load(file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated_data/surv.data_sre_1.RData")

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
    chaz[i]<- pow(st[i],nu)*exp(inprod(alpha[],W[i,])+inprod(b[i,],gamma[]))
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
  nu~dgamma(0.1,0.1)
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
mu.vect sd.vect          2.5%          25%          50%           75%         97.5%  Rhat n.eff
Sigma[1,1]  8.20000e-01   0.091         0.656  7.59000e-01  8.13000e-01         0.874         1.011 1.018   120
Sigma[2,1]  5.65000e-01   0.083         0.406  5.08000e-01  5.64000e-01         0.620         0.730 1.038    61
Sigma[1,2]  5.65000e-01   0.083         0.406  5.08000e-01  5.64000e-01         0.620         0.730 1.038    61
Sigma[2,2]  8.20000e-01   0.134         0.586  7.30000e-01  8.08000e-01         0.898         1.127 1.030    73
alpha[1]    1.22900e+00   0.170         0.945  1.11200e+00  1.21300e+00         1.328         1.644 1.042    68
alpha[2]   -1.17400e+00   0.135        -1.471 -1.25600e+00 -1.16200e+00        -1.078        -0.949 1.044    56
alpha[3]   -6.14000e-01   0.146        -0.929 -7.06000e-01 -6.08000e-01        -0.512        -0.353 1.049    57
beta[1]    -5.26000e-01   0.080        -0.687 -5.78000e-01 -5.25000e-01        -0.471        -0.375 1.007  1100
beta[2]     6.11000e-01   0.106         0.401  5.40000e-01  6.13000e-01         0.686         0.809 1.063    39
beta[3]     6.21000e-01   0.108         0.408  5.51000e-01  6.19000e-01         0.690         0.844 1.006   360
beta[4]     4.33000e-01   0.051         0.333  4.00000e-01  4.32000e-01         0.467         0.534 1.016   160
gamma[1]    1.22200e+00   0.337         0.660  1.00100e+00  1.18300e+00         1.375         2.023 1.181    16
gamma[2]   -1.24100e+00   0.343        -2.008 -1.43100e+00 -1.19800e+00        -1.004        -0.641 1.187    17
nu          1.13800e+00   0.099         0.972  1.06700e+00  1.13200e+00         1.198         1.365 1.034    76
sigma       1.03400e+00   0.042         0.954  1.00600e+00  1.03300e+00         1.063         1.119 1.012   170
deviance    1.00005e+08  76.771 100004845.394  1.00005e+08  1.00005e+08 100005063.339 100005144.203 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 2870.0 and DIC = 100007877.0
DIC is an estimate of expected predictive error (lower deviance is better).
