rm(list = ls())
library(R2jags)
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


##### Results ##### 
print(sim)
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect          2.5%          25%          50%          75%         97.5%  Rhat n.eff
Sigma[1,1]  1.02700e+00   0.101         0.844  9.56000e-01  1.02300e+00  1.09200e+00         1.235 1.018   120
Sigma[2,1]  4.70000e-01   0.111         0.282  3.92000e-01  4.62000e-01  5.33000e-01         0.733 1.114    24
Sigma[1,2]  4.70000e-01   0.111         0.282  3.92000e-01  4.62000e-01  5.33000e-01         0.733 1.114    24
Sigma[2,2]  7.07000e-01   0.149         0.450  6.00000e-01  6.97000e-01  7.95000e-01         1.029 1.027   140
alpha[1]    7.16000e-01   0.077         0.564  6.63000e-01  7.16000e-01  7.67000e-01         0.868 1.005   460
alpha[2]    5.99000e-01   0.063         0.479  5.57000e-01  6.00000e-01  6.41000e-01         0.722 1.001  3700
alpha[3]   -5.76000e-01   0.113        -0.803 -6.54000e-01 -5.73000e-01 -4.98000e-01        -0.361 1.001  3800
beta[1]    -8.23000e-01   0.084        -0.997 -8.74000e-01 -8.21000e-01 -7.68000e-01        -0.659 1.053    43
beta[2]     6.95000e-01   0.118         0.454  6.15000e-01  7.00000e-01  7.77000e-01         0.905 1.129    24
beta[3]     9.24000e-01   0.107         0.712  8.53000e-01  9.22000e-01  9.95000e-01         1.132 1.049    48
beta[4]     4.05000e-01   0.058         0.293  3.65000e-01  4.06000e-01  4.47000e-01         0.514 1.020   120
gamma      -7.07000e-01   0.048        -0.804 -7.38000e-01 -7.07000e-01 -6.75000e-01        -0.617 1.005   500
lambda0     2.05200e+00   0.159         1.758  1.94100e+00  2.04700e+00  2.15400e+00         2.382 1.005   490
sigma       1.01400e+00   0.039         0.943  9.88000e-01  1.01300e+00  1.04000e+00         1.094 1.006   380
deviance    1.00006e+08  43.427 100005916.850  1.00006e+08  1.00006e+08  1.00006e+08 100006086.777 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 936.6 and DIC = 100006935.7
DIC is an estimate of expected predictive error (lower deviance is better).
> 