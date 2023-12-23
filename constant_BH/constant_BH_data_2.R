rm(list = ls())
library(R2jags)
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
print(sim)
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.014   0.099         0.825         0.947         1.010         1.077         1.218 1.002  2400
Sigma[2,1]         0.455   0.098         0.262         0.389         0.456         0.521         0.643 1.006   390
Sigma[1,2]         0.455   0.098         0.262         0.389         0.456         0.521         0.643 1.006   390
Sigma[2,2]         1.169   0.179         0.845         1.043         1.155         1.282         1.553 1.046    50
alpha[1]          -1.557   0.148        -1.861        -1.655        -1.553        -1.452        -1.289 1.002  1500
alpha[2]           0.818   0.125         0.576         0.734         0.817         0.900         1.067 1.002  1300
alpha[3]          -0.971   0.241        -1.454        -1.133        -0.962        -0.809        -0.507 1.002  2200
beta[1]           -0.420   0.089        -0.589        -0.483        -0.420        -0.355        -0.245 1.043    54
beta[2]            0.565   0.084         0.399         0.509         0.566         0.622         0.729 1.003   690
beta[3]            0.525   0.056         0.420         0.487         0.524         0.562         0.639 1.011   190
beta[4]            0.451   0.121         0.225         0.359         0.453         0.543         0.670 1.050    51
gamma             -0.370   0.077        -0.522        -0.421        -0.371        -0.318        -0.219 1.002  1500
lambda0            0.213   0.031         0.155         0.191         0.212         0.234         0.276 1.002  1500
sigma              0.963   0.033         0.898         0.940         0.961         0.985         1.030 1.004   600
deviance   100007257.098  48.535 100007164.035 100007224.149 100007256.511 100007289.407 100007352.837 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1173.7 and DIC = 100008430.8
DIC is an estimate of expected predictive error (lower deviance is better).
> 