rm(list = ls())
library(PermAlgo)
library(mvtnorm)
library(MASS)
library(R2jags)
library(abind)
library(statmod)
library(splines2)

# number of indivduals
n <- 500
# real values
Beta <- c(0.5, -0.2, 0.5, -0.5)
# residual error
sigmae <- sqrt(0.5)
# association parameter
gamma <- -0.5
# covariance matrix of random effects
rho <- 0.5
Sigma <- matrix(c(1, rho, rho, 1), 2, 2)

gapLongi <- 0.2 # gap between longi measurements
gap <- 0.02 # used to generate a lot of time points because the permutation
# algorithm chooses among those time points to define survival times
followup <- 3 # follow-up time
mestime <- seq(0, followup, gap) # measurement times
timesLongi <- mestime[which(round(mestime - round(mestime / gapLongi, 0) * gapLongi, 6) == 0)] # visit times
time <- rep(mestime, n) # time column
nmesindiv <- followup / gap + 1 # max. number of individual measurements
nmesy <- nmesindiv * n # max. total number of longi measurements
# the number is reduced because of censoring by the terminal event
idY <- rep(1:n, each = nmesindiv) # individual id
###############

b <- rmvnorm(n, rep(0, 2), Sigma)
b0 <- rep(b[, 1], each = nmesindiv) # random intercept Y1
b1 <- rep(b[, 2], each = nmesindiv) # random slope Y1

x1 <- rnorm(n, 1, 0.5)
x2 <- rbinom(n, 1, 0.5)
X1 <- rep(x1, each = nmesindiv) # continuous covariate
X2 <- rep(x2, each = nmesindiv) # binary covariate

# linear predictors
linPredY1 <- (Beta[1] + b0) + (Beta[2] + b1) * time + Beta[3] * X1 + Beta[4] * X2

# outcome Y1
Y1 <- rnorm(nmesy, linPredY1, sigmae)

# Permutation algorithm to generate survival times that depends on the linear predictors
DatTmp <- permalgorithm(n, nmesindiv,
                        Xmat = linPredY1,
                        eventRandom = round(rexp(n, 0.003) + 1, 0), # ~40% death
                        censorRandom = runif(n, 1, nmesindiv), # uniform random censoring
                        XmatNames = c("linPredY1"), # association
                        betas = gamma
) # association parameters

# extract last line for each Id (= death/censoring time)
DatTmp2 <- DatTmp[c(which(diff(DatTmp[, "Id"]) == 1), dim(DatTmp)[1]), c("Id", "Event", "Stop")]
DatTmp2$deathTimes <- mestime[DatTmp2$Stop + 1] # deathtimes
survDat <- DatTmp2[, c("Id", "deathTimes", "Event")]
DatTmp$time <- mestime[DatTmp$Start + 1] # measurements times of the biomarker
DatTmp$Uid <- paste(DatTmp$Id, DatTmp$time) # unique identifier to match covariates and observed biomarker values
longDat3 <- merge(DatTmp[, c("Uid", "Id", "time")], cbind("Uid" = paste(idY, time), X1, X2, Y1), by = c("Uid"))
longDat <- sapply(longDat3[longDat3$time %in% timesLongi, -1], as.numeric)
longDat <- as.data.frame(longDat[order(longDat[, "Id"], longDat[, "time"]), ])
summary(survDat) # survival dataset
summary(longDat) # longitudinal dataset
names(survDat)
names(longDat)
median(table(longDat$Id))
min(table(longDat$Id))
max(table(longDat$Id))
sum(survDat$Event) / length(survDat$Event)
#######################################
# Number of patients and number of longitudinal observations per patient
n
M <- table(longDat$Id)
max(M)
# Survival and censoring times
st <- survDat$deathTimes
death <- survDat$Event # death=1 means observed

# Longitudinal information in matrix format
time1 <- Y1n <- matrix(NA, n, max(M))
for (i in 1:n) {
  time1[i, 1:M[i]] <- longDat$time[longDat$Id == i]
  Y1n[i, 1:M[i]] <- longDat$Y1[longDat$Id == i]
}

treat <- x2 # Reference = placebo
W <- model.matrix(~x2) # Fixed effects
W[, 1] <- x1
X <- array(1, dim = c(n, max(M), 4)) # Fixed effects
X[, , 2] <- time1
X[, , 3] <- x1
X[, , 4] <- x2

ZL <- array(1, dim = c(n, max(M), 2)) # Random effects
ZL[, , 2] <- time1
############
# Gauss-Legendre quadrature (15 points)
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes # Nodes
wk <- glq$weights # Weights
K <- length(xk) # K-points

knots <- quantile(st, prob = seq(0.1, .9, .15)) # B-spline
st_bs <- bSpline(st, knots = knots, degree = 4, intercept = TRUE)

XK_sc <- matrix(0, K, n)
for (i in 1:n) {
  XK_sc[, i] <- (xk + 1) / 2 * st[i]
}


XK_t <- array(0, c(K, dim(st_bs)[2], n))
for (i in 1:n) {
  knots <- quantile(XK_sc[, i], prob = seq(0.1, .9, .15)) # B-spline
  
  XK_t[, , i] <- bSpline(XK_sc[, i], knots = knots, degree = 4, intercept = TRUE)
}
dim(XK_t)

K1 <- matrix(0, (dim(st_bs)[2] - 1), dim(st_bs)[2])
for (i in 1:(dim(st_bs)[2] - 1)) {
  K1[i, i] <- 1
  K1[i, i + 1] <- -1
}
Psi1 <- t(K1) %*% K1
######
K2 <- matrix(0, (dim(st_bs)[2] - 2), dim(st_bs)[2])
for (i in 1:(dim(st_bs)[2] - 2)) {
  K2[i, i] <- 1
  K2[i, i + 1] <- -2
  K2[i, i + 2] <- 1
}
Psi2 <- t(K2) %*% K2
#######
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
  #Prior distributions
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
  #Determinetic nodes
  sigma<-1/tau
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])
}", fill = TRUE)
sink()
#######
d.jags <- list(
  n = n, M = M, st = st, W = W, Y = Y1n, X = X, Z = ZL, K = K, X1 = X1, X2 = X2,
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
mu.vect sd.vect         2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]           0.984   0.075  8.49000e-01         0.932         0.980         1.032         1.141 1.003   860
Sigma[2,1]           0.479   0.064  3.58000e-01         0.436         0.479         0.520         0.609 1.003   970
Sigma[1,2]           0.479   0.064  3.58000e-01         0.436         0.479         0.520         0.609 1.003   970
Sigma[2,2]           1.021   0.101  8.32000e-01         0.949         1.016         1.087         1.231 1.002  1500
alpha[1]             0.142   0.218 -2.78000e-01        -0.009         0.140         0.291         0.566 1.012   680
alpha[2]             0.339   0.230 -1.07000e-01         0.182         0.339         0.492         0.802 1.001  3800
beta[1]              0.736   0.122  4.87000e-01         0.654         0.741         0.819         0.962 1.065    81
beta[2]             -0.204   0.053 -3.09000e-01        -0.240        -0.203        -0.168        -0.099 1.018   190
beta[3]              0.153   0.090 -2.10000e-02         0.089         0.152         0.214         0.339 1.016   140
beta[4]             -0.293   0.103 -5.16000e-01        -0.356        -0.293        -0.223        -0.097 1.015   140
gamma               -0.684   0.103 -8.89000e-01        -0.752        -0.683        -0.615        -0.488 1.004  1300
gamma_bs[1]         -1.425   0.638 -2.79300e+00        -1.827        -1.372        -0.988        -0.238 1.033    73
gamma_bs[2]         -0.585   0.585 -1.69400e+00        -0.970        -0.608        -0.163         0.578 1.029    96
gamma_bs[3]         -1.792   0.684 -3.19400e+00        -2.193        -1.777        -1.335        -0.538 1.052   130
gamma_bs[4]         -1.717   0.667 -2.99100e+00        -2.181        -1.701        -1.242        -0.424 1.035   100
gamma_bs[5]         -2.378   0.669 -3.78000e+00        -2.782        -2.356        -1.920        -1.174 1.046   110
gamma_bs[6]         -2.659   0.651 -4.05100e+00        -3.072        -2.626        -2.200        -1.497 1.051    58
gamma_bs[7]         -2.074   0.685 -3.28400e+00        -2.539        -2.102        -1.661        -0.722 1.043   640
gamma_bs[8]         -3.661   0.773 -5.25600e+00        -4.163        -3.668        -3.113        -2.276 1.024   850
gamma_bs[9]         -4.475   1.044 -6.56000e+00        -5.105        -4.383        -3.812        -2.566 1.065    63
gamma_bs[10]        -5.937   1.631 -9.94000e+00        -6.831        -5.728        -4.856        -3.223 1.172    25
gamma_bs[11]        -6.872   2.740 -1.50050e+01        -7.636        -6.394        -5.112        -3.300 1.362    12
sigma                0.508   0.014  4.82000e-01         0.498         0.508         0.517         0.535 1.001  3800
deviance     100008139.406  47.526  1.00008e+08 100008106.251 100008139.786 100008172.524 100008229.882 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1128.2 and DIC = 100009267.6
DIC is an estimate of expected predictive error (lower deviance is better).

