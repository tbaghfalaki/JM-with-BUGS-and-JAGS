rm(list = ls())
library(R2jags)
library(nnet)
library(statmod)
library(splines2)

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

# B-spline for survival time
knots <- quantile(st, prob = seq(0.3, 0.8, length = 2))
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
taubs<-1/taubs1
taubs1~dgamma(1,0.005)

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
  n.iter = 5000, model.file = "model_file"
)

print(sim)

Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]          1.042   0.084         0.887         0.985         1.040         1.094         1.215 1.001  3800
Sigma[2,1]          0.559   0.065         0.438         0.514         0.557         0.603         0.690 1.005   450
Sigma[1,2]          0.559   0.065         0.438         0.514         0.557         0.603         0.690 1.005   450
Sigma[2,2]          0.865   0.091         0.700         0.800         0.861         0.925         1.057 1.009   250
alpha[1]           -1.448   1.015        -3.129        -2.319        -1.474        -0.610         0.498 1.542     7
alpha[2]            1.359   0.084         1.197         1.302         1.359         1.416         1.522 1.002  1600
beta[1]             0.529   0.074         0.385         0.480         0.532         0.581         0.668 1.043    58
beta[2]            -0.649   0.052        -0.749        -0.686        -0.649        -0.615        -0.545 1.019   120
beta[3]             0.993   0.095         0.801         0.928         0.996         1.058         1.172 1.041    57
beta[4]             1.035   0.046         0.944         1.003         1.035         1.069         1.119 1.103    25
gamma              -0.372   0.046        -0.463        -0.405        -0.373        -0.341        -0.282 1.002  1200
gamma_bs[1]         1.931   1.034        -0.129         1.093         2.010         2.828         3.665 1.508     7
gamma_bs[2]         0.848   0.991        -0.924         0.042         0.809         1.720         2.558 1.483     8
gamma_bs[3]         2.406   1.204         0.209         1.551         2.270         3.353         4.626 1.590     7
gamma_bs[4]        -2.804   1.139        -4.853        -3.646        -2.855        -1.946        -0.619 1.209    14
gamma_bs[5]        -6.613   1.352        -9.621        -7.480        -6.423        -5.666        -4.160 1.081    41
gamma_bs[6]        -7.769   2.697       -13.569        -9.270        -7.272        -5.863        -3.187 1.074    62
sigma               1.006   0.019         0.971         0.994         1.007         1.019         1.044 1.002  1600
deviance    100019131.193  43.337 100019050.981 100019101.346 100019128.481 100019158.608 100019221.155 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 936.9 and DIC = 100020068.0
DIC is an estimate of expected predictive error (lower deviance is better).

