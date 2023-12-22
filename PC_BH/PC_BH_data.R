rm(list = ls())
library(R2jags)
library(nnet)

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
    #Longitudinal sub-model
    for(j in 1:M[i]){
      Y[i,j]~dnorm(mu[i,j],tau)
      mu[i,j]<-inprod(beta[],X[i,j,])+inprod(b[i,],Z[i,j,])
    }
    Alpha0[i]<- inprod(alpha[],W[i,])+gamma*(beta[1]+beta[3]*X1[i]+beta[4]*X2[i]+b[i,1])
    Alpha1[i]<- gamma*(beta[2]+b[i,2])
    #Survival sub-model
    #Hazard function
    haz[i]<- ((h[1]*step(s[1]-st[i]))+
             (h[2]*step(st[i]-s[1])*step(s[2]-st[i]))+
             (h[3]*step(st[i]-s[2])*step(s[3]-st[i]))+
             (h[4]*step(st[i]-s[3])*step(s[4]-st[i]))+
             (h[5]*step(st[i]-s[4])))*exp(Alpha0[i]+Alpha1[i]*st[i])
  
      #Cumulative hazard function
    chaz1[i]<- h[1]*(exp(Alpha1[i]*st[i])-1)*step(s[1]-st[i])+
      (h[1]*(exp(Alpha1[i]*s[1])-1)+h[2]*(exp(Alpha1[i]*st[i])-exp(Alpha1[i]*s[1])))*(step(st[i]-s[1])*step(s[2]-st[i]))+
      (h[1]*(exp(Alpha1[i]*s[1])-1)+h[2]*(exp(Alpha1[i]*s[2])-exp(Alpha1[i]*s[1]))+h[3]*(exp(Alpha1[i]*st[i])-exp(Alpha1[i]*s[2])))*(step(st[i]-s[2])*step(s[3]-st[i]))+
      (h[1]*(exp(Alpha1[i]*s[1])-1)+h[2]*(exp(Alpha1[i]*s[2])-exp(Alpha1[i]*s[1]))+h[3]*(exp(Alpha1[i]*s[3])-exp(Alpha1[i]*s[2]))+h[4]*(exp(Alpha1[i]*st[i])-exp(Alpha1[i]*s[3])))*(step(st[i]-s[3])*step(s[4]-st[i]))+
      (h[1]*(exp(Alpha1[i]*s[1])-1)+h[2]*(exp(Alpha1[i]*s[2])-exp(Alpha1[i]*s[1]))+h[3]*(exp(Alpha1[i]*s[3])-exp(Alpha1[i]*s[2]))+
      h[4]*(exp(Alpha1[i]*s[4])-exp(Alpha1[i]*s[3]))+h[5]*(exp(Alpha1[i]*st[i])-exp(Alpha1[i]*s[4])))*step(st[i]-s[4])

    chaz[i]<-exp(Alpha0[i])*chaz1[i]/Alpha1[i]
    #Log-survival function
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

  for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

  gamma~dnorm(0,0.001)
  tau~dgamma(0.01,0.01)
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  #Determinetic nodes
  sigma<-1/tau
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])
}", fill = TRUE)
sink()
#### Running JAGS
s <- quantile(st, seq(0.2, .8, length = 4))

delta <- rep(0, n)
for (i in 1:n) {
  if (st[i] <= s[1]) (delta[i] <- 1)
  if (st[i] > s[1] & st[i] <= s[2]) (delta[i] <- 2)
  if (st[i] > s[2] & st[i] <= s[3]) (delta[i] <- 3)
  if (st[i] > s[3] & st[i] <= s[4]) (delta[i] <- 4)
  if (st[i] > s[4]) (delta[i] <- 5)
}
table(delta)
Delta <- class.ind(delta)

d.jags <- list(
  n = n, M = M, st = st, W = W, Y = Y, X = X, Z = Z, X1 = X1, X2 = X2,
  death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
  Nbeta = dim(X)[3], Nalpha = ncol(W), delta = Delta, s = s, J = length(s) + 1
)

i.jags <- function() {
  list(
    alpha = rnorm(ncol(W)), gamma = rnorm(1),
    beta = rnorm(dim(X)[3]), tau = runif(1), Omega = diag(runif(2))
  )
}
parameters <- c("alpha", "gamma", "beta", "sigma", "Sigma", "h")

sim <- jags(
  data = d.jags, inits = i.jags, parameters,
  n.iter = 2000, model.file = "model_file"
)

print(sim)

##### Results #####
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 2000 iterations (first 1000 discarded)
n.sims = 3000 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.037   0.105         0.843         0.964         1.035         1.106         1.250 1.050    53
Sigma[2,1]         0.453   0.094         0.276         0.387         0.451         0.516         0.647 1.038    65
Sigma[1,2]         0.453   0.094         0.276         0.387         0.451         0.516         0.647 1.038    65
Sigma[2,2]         0.944   0.154         0.673         0.840         0.930         1.038         1.287 1.036    94
alpha[1]          -1.668   1.163        -3.923        -2.333        -1.589        -1.021         0.729 1.197    20
alpha[2]           0.995   0.136         0.734         0.903         0.994         1.085         1.265 1.005   550
beta[1]           -0.528   0.081        -0.678        -0.582        -0.534        -0.478        -0.355 1.019   180
beta[2]            0.594   0.080         0.418         0.548         0.600         0.645         0.734 1.091    41
beta[3]            0.453   0.056         0.353         0.411         0.450         0.491         0.566 1.139    19
beta[4]            0.713   0.104         0.489         0.645         0.721         0.788         0.896 1.029   220
gamma             -0.636   0.091        -0.816        -0.698        -0.635        -0.575        -0.459 1.033    75
h[1]               0.351   0.463         0.015         0.090         0.166         0.363         1.783 1.172    23
h[2]               1.618   2.102         0.073         0.428         0.768         1.643         7.824 1.184    22
h[3]               1.320   1.827         0.055         0.320         0.614         1.272         7.030 1.189    21
h[4]               1.543   2.032         0.064         0.399         0.720         1.645         7.894 1.180    21
h[5]               2.873   3.847         0.124         0.709         1.416         3.037        14.640 1.153    25
sigma              0.961   0.034         0.899         0.938         0.959         0.983         1.029 1.003   730
deviance   100007198.312  49.606 100007102.273 100007164.492 100007198.233 100007230.671 100007301.156 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1219.4 and DIC = 100008417.7
DIC is an estimate of expected predictive error (lower deviance is better).
