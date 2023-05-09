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
  time[i, 1:M[i]] <- long.data$obstime[long.data$id==i]
  Y[i, 1:M[i]] <- long.data$Y1[long.data$id==i]
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
  death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0,n),
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
  n.iter = 5000, model.file =  "model_file"
)

print(sim)


### Results
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.050   0.081         0.899         0.996         1.045         1.102         1.221 1.002  1500
Sigma[2,1]         0.548   0.066         0.422         0.503         0.546         0.591         0.680 1.012   170
Sigma[1,2]         0.548   0.066         0.422         0.503         0.546         0.591         0.680 1.012   170
Sigma[2,2]         0.863   0.085         0.707         0.804         0.859         0.916         1.042 1.011   190
alpha[1]          -0.983   0.130        -1.242        -1.067        -0.982        -0.894        -0.725 1.009   240
alpha[2]           1.175   0.087         1.001         1.116         1.175         1.234         1.342 1.001  3800
beta[1]            0.466   0.068         0.335         0.418         0.464         0.514         0.597 1.029    73
beta[2]           -0.606   0.067        -0.747        -0.648        -0.605        -0.560        -0.472 1.074    42
beta[3]            1.016   0.090         0.828         0.953         1.022         1.079         1.185 1.051    51
beta[4]            1.024   0.050         0.922         0.989         1.024         1.059         1.119 1.015   140
gamma[1]          -0.187   0.051        -0.285        -0.221        -0.187        -0.153        -0.089 1.011   220
gamma[2]           0.073   0.144        -0.210        -0.023         0.070         0.167         0.366 1.014   190
kappa              0.826   0.047         0.735         0.795         0.827         0.858         0.919 1.001  3800
sigma              1.004   0.019         0.967         0.991         1.004         1.016         1.043 1.005   510
deviance   100019470.165  45.063 100019377.276 100019440.860 100019469.458 100019500.334 100019557.372 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1004.1 and DIC = 100020474.3
DIC is an estimate of expected predictive error (lower deviance is better).
> 