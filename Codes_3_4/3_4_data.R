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
    #Longitudinal sub-model
    for(j in 1:M[i]){
      Y[i,j]~dnorm(mu[i,j],tau)
      mu[i,j]<-inprod(beta[],X[i,j,])+inprod(b[i,],Z[i,j,])
    }
    Alpha0[i]<- inprod(alpha[],W[i,])+gamma*(beta[1]+beta[3]*X1[i]+beta[4]*X2[i]+b[i,1])
    Alpha1[i]<- gamma*(beta[2]+b[i,2])
    #Survival sub-model
    #Hazard function
    haz[i]<- (h[1]*equals(delta[i,1],1)+h[2]*equals(delta[i,2],1)+h[3]*equals(delta[i,3],1)+
                h[4]*equals(delta[i,4],1)+h[5]*equals(delta[i,5],1))*exp(Alpha0[i]+Alpha1[i]*st[i])
    #Cumulative hazard function
    chaz1[i]<- h[1]*(exp(Alpha1[i]*st[i])-1)*equals(delta[i,1],1)+
      (h[1]*(exp(Alpha1[i]*s[1])-1)+h[2]*(exp(Alpha1[i]*st[i])-exp(Alpha1[i]*s[1])))*equals(delta[i,2],1)+
      (h[1]*(exp(Alpha1[i]*s[1])-1)+h[2]*(exp(Alpha1[i]*s[2])-exp(Alpha1[i]*s[1]))+h[3]*(exp(Alpha1[i]*st[i])-exp(Alpha1[i]*s[2])))*equals(delta[i,3],1)+
      (h[1]*(exp(Alpha1[i]*s[1])-1)+h[2]*(exp(Alpha1[i]*s[2])-exp(Alpha1[i]*s[1]))+h[3]*(exp(Alpha1[i]*s[3])-exp(Alpha1[i]*s[2]))+h[4]*(exp(Alpha1[i]*st[i])-exp(Alpha1[i]*s[3])))*equals(delta[i,4],1)+
      (h[1]*(exp(Alpha1[i]*s[1])-1)+h[2]*(exp(Alpha1[i]*s[2])-exp(Alpha1[i]*s[1]))+h[3]*(exp(Alpha1[i]*s[3])-exp(Alpha1[i]*s[2]))+h[4]*(exp(Alpha1[i]*s[4])-exp(Alpha1[i]*s[3]))+h[5]*(exp(Alpha1[i]*st[i])-exp(Alpha1[i]*s[4])))*equals(delta[i,5],1)

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
  n.iter = 5000, model.file = "model_file"
)

print(sim)

##### Results #####
> print(sim)
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect          2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.102   0.125         0.877         1.014         1.093         1.181         1.366 1.028    77
Sigma[2,1]         0.423   0.242        -0.028         0.246         0.424         0.592         0.898 1.063    38
Sigma[1,2]         0.423   0.242        -0.028         0.246         0.424         0.592         0.898 1.063    38
Sigma[2,2]         1.117   0.516         0.286         0.738         1.055         1.428         2.343 1.525     8
alpha[1]           3.180   1.969        -0.435         1.310         3.637         4.822         5.877 4.890     3
alpha[2]           0.026   0.114        -0.203        -0.048         0.027         0.101         0.246 1.010   210
beta[1]            0.457   0.149         0.174         0.355         0.451         0.558         0.746 1.040   230
beta[2]           -0.523   0.305        -1.151        -0.722        -0.519        -0.316         0.049 1.115    23
beta[3]            1.079   0.118         0.848         0.997         1.080         1.160         1.308 1.010   210
beta[4]            1.190   0.122         0.956         1.102         1.189         1.278         1.420 1.005   420
gamma             -0.393   0.051        -0.496        -0.426        -0.393        -0.359        -0.292 1.025    86
h[1]               0.825   1.500         0.010         0.028         0.095         0.949         5.802 4.821     3
h[2]               1.290   2.314         0.016         0.045         0.150         1.500         8.645 4.853     3
h[3]               1.546   2.803         0.018         0.053         0.176         1.799        10.147 4.824     3
h[4]               1.619   2.964         0.019         0.056         0.183         1.851        11.303 4.792     3
h[5]               2.398   4.369         0.028         0.082         0.268         2.721        17.104 4.758     3
sigma              0.962   0.065         0.843         0.917         0.958         1.005         1.095 1.015   140
deviance   100002582.609  48.073 100002488.415 100002550.538 100002581.500 100002615.309 100002679.683 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1113.3 and DIC = 100003695.9
DIC is an estimate of expected predictive error (lower deviance is better).
