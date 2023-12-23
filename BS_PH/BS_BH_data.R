rm(list = ls())
library(R2jags)
library(nnet)
library(statmod)
library(splines2)

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
W=W[,-1]
X <- array(1, dim = c(n, max(M), 4)) # Fixed effects
X[, , 2] <- time
X[, , 3] <- X1
X[, , 4] <- X2

Z <- array(1, dim = c(n, max(M), 2)) # Random effects
Z[, , 2] <- time

# B-spline for survival time
knots <- quantile(st, prob = seq(0.3, 0.7, length = 2))
st_bs <- bSpline(st, knots = knots, degree = 3, intercept = TRUE)

# B-spline for Gauss-Legendre quadrature
# Gauss-Legendre quadrature (15 points)
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes # Nodes
wk <- glq$weights # Weights
K <- length(xk) # K-points


XK_sc <- matrix(0, K, n)
for (i in 1:n) {
  XK_sc[, i] <- (xk + 1) / 2 * st[i]
}


XK_t <- array(0, c(K, dim(st_bs)[2], n))
for (i in 1:n) {
  knots <- quantile(XK_sc[, i], prob = seq(0.3, .8, length = 2)) # B-spline
  
  XK_t[, , i] <- bSpline(XK_sc[, i], knots = knots, degree = 3, intercept = TRUE)
}
dim(XK_t)

#### Psi1 and Psi2
K1 <- matrix(0, (dim(st_bs)[2] - 1), dim(st_bs)[2])
for (i in 1:(dim(st_bs)[2] - 1)) {
  K1[i, i] <- 1
  K1[i, i + 1] <- -1
}
Psi1 <- t(K1) %*% K1


K2 <- matrix(0, (dim(st_bs)[2] - 2), dim(st_bs)[2])
for (i in 1:(dim(st_bs)[2] - 2)) {
  K2[i, i] <- 1
  K2[i, i + 1] <- -2
  K2[i, i + 2] <- 1
}
Psi2 <- t(K2) %*% K2
########  BUGS code  ########
sink("model_file")
cat("model{
  for(i in 1:n){
    # Longitudinal sub-model
    for(j in 1:M[i]){
      Y[i,j]~dnorm(mu[i,j],tau)
      mu[i,j]<-inprod(beta[],X[i,j,])+inprod(b[i,],Z[i,j,])
    }
    Alpha0[i]<- inprod(alpha[],W[i,])+gamma*(beta[1]+beta[3]*X1[i]+beta[4]*X2[i]+b[i,1])
    Alpha1[i]<- gamma*(beta[2]+b[i,2])

    # Survival sub-model
    # Hazard function
    haz[i]<- exp(inprod(st_bs[i,],gamma_bs[])+Alpha0[i]+Alpha1[i]*st[i])
    # Cumulative hazard function
    for(j in 1:K){
    # Scaling Gauss Kronrod (see wikipedia)
    x_scaled[i,j]<-(xk[j]+1)/2*st[i]
    w_scaled[i,j]<- wk[j]*st[i]/2
    chaz[i,j]<-exp(inprod(XK_t[j,,i],gamma_bs[])+Alpha0[i]+Alpha1[i]*x_scaled[i,j])
}
    # Log-survival function
    logSurv[i]<- -inprod(w_scaled[i,],chaz[i,])
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

gamma_bs[1:Nbs]~dmnorm(mu0[], Sigmasp[,])

for(k1 in 1:Nbs){
for(k2 in 1:Nbs){
Sigmasp[k1,k2]<- taubs*Psi1[k1,k2]
}}
taubs~dgamma(1,0.005)

  gamma~dnorm(0,0.001)
  tau~dgamma(0.01,0.01)
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  # Determinetic nodes
  sigma<-1/tau
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])

}", fill = TRUE)
sink()
#### Running JAGS #######

d.jags <- list(
  n = n, M = M, st = st, W = W, Y = Y, X = X, Z = Z, K = K, X1 = X1, X2 = X2,
  death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
  Nbeta = dim(X)[3], Nalpha = dim(W)[2], st_bs = st_bs, xk = xk, wk = wk,
  XK_t = XK_t, Nbs = dim(st_bs)[2], Psi1 = Psi1, mu0 = rep(0, dim(st_bs)[2])
)

i.jags <- function() {
  list(
    gamma = rnorm(1), alpha = rnorm(dim(W)[2]),
    beta = rnorm(dim(X)[3]), tau = runif(1), Omega = diag(runif(2)), gamma_bs = rnorm(dim(st_bs)[2]), taubs1 = 1
  )
}

parameters <- c("alpha", "gamma", "beta", "sigma", "Sigma", "gamma_bs")

sim <- jags(
  data = d.jags, inits = i.jags, parameters,
  n.iter = 2000, model.file = "model_file"
)

print(sim)

> print(sim)
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 2000 iterations (first 1000 discarded)
n.sims = 3000 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%        97.5%  Rhat n.eff
Sigma[1,1]          1.044   0.104         0.860         0.972         1.037         1.111  1.26600e+00 1.053    43
Sigma[2,1]          0.569   0.097         0.361         0.511         0.569         0.629  7.61000e-01 1.114    30
Sigma[1,2]          0.569   0.097         0.361         0.511         0.569         0.629  7.61000e-01 1.114    30
Sigma[2,2]          0.799   0.139         0.559         0.704         0.792         0.884  1.11000e+00 1.112    24
alpha[1]            0.611   0.064         0.484         0.568         0.614         0.654  7.38000e-01 1.001  2200
alpha[2]           -0.572   0.106        -0.787        -0.644        -0.571        -0.498 -3.65000e-01 1.036    62
beta[1]            -0.808   0.081        -0.962        -0.864        -0.810        -0.750 -6.50000e-01 1.056    45
beta[2]             0.576   0.079         0.420         0.524         0.575         0.629  7.36000e-01 1.053    62
beta[3]             0.920   0.113         0.682         0.848         0.925         0.996  1.12600e+00 1.057    93
beta[4]             0.395   0.054         0.295         0.358         0.393         0.431  5.07000e-01 1.005  1300
gamma              -0.720   0.046        -0.814        -0.751        -0.719        -0.687 -6.35000e-01 1.017   150
gamma_bs[1]         1.037   0.106         0.846         0.956         1.036         1.105  1.25400e+00 1.093    29
gamma_bs[2]         1.003   0.083         0.859         0.940         1.005         1.059  1.17800e+00 1.064    42
gamma_bs[3]         0.800   0.076         0.644         0.743         0.799         0.858  9.42000e-01 1.073    37
gamma_bs[4]         0.426   0.085         0.246         0.365         0.428         0.484  5.86000e-01 1.095    28
gamma_bs[5]         0.213   0.094         0.041         0.144         0.207         0.271  4.18000e-01 1.074    33
gamma_bs[6]         0.189   0.110        -0.017         0.104         0.189         0.263  4.33000e-01 1.019   140
sigma               1.013   0.038         0.942         0.987         1.012         1.038  1.09000e+00 1.006   360
deviance    100005876.652  45.217 100005786.176 100005846.529 100005876.366 100005906.089  1.00006e+08 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1016.7 and DIC = 100006893.4
DIC is an estimate of expected predictive error (lower deviance is better).
> 