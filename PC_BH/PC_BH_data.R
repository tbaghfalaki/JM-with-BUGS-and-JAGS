rm(list = ls())
library(R2jags)
library(nnet)

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
W <- W[,-1]

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



d.jags <- list(
  n = n, M = M, st = st, W = W, Y = Y, X = X, Z = Z, X1 = X1, X2 = X2,
  death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
  Nbeta = dim(X)[3], Nalpha = ncol(W),  s = s, J = length(s) + 1
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
Sigma[1,1]         1.016   0.095         0.845         0.951         1.011         1.077         1.209 1.029    78
Sigma[2,1]         0.455   0.086         0.293         0.395         0.456         0.515         0.625 1.018   120
Sigma[1,2]         0.455   0.086         0.293         0.395         0.456         0.515         0.625 1.018   120
Sigma[2,2]         1.144   0.170         0.862         1.021         1.127         1.248         1.520 1.042    60
alpha[1]           0.940   0.138         0.679         0.842         0.940         1.032         1.208 1.002  1500
alpha[2]          -1.007   0.251        -1.506        -1.167        -1.008        -0.837        -0.525 1.002  1300
beta[1]           -0.434   0.083        -0.585        -0.497        -0.432        -0.372        -0.276 1.076    36
beta[2]            0.563   0.079         0.394         0.516         0.567         0.616         0.706 1.140    20
beta[3]            0.541   0.057         0.435         0.500         0.538         0.579         0.658 1.178    16
beta[4]            0.463   0.095         0.286         0.392         0.464         0.536         0.631 1.024    91
gamma             -0.429   0.074        -0.578        -0.479        -0.428        -0.378        -0.286 1.006   350
h[1]               0.069   0.019         0.036         0.055         0.067         0.080         0.112 1.004   540
h[2]               0.235   0.053         0.144         0.198         0.229         0.267         0.350 1.003   780
h[3]               0.170   0.050         0.086         0.134         0.164         0.200         0.280 1.006   360
h[4]               0.240   0.085         0.099         0.178         0.231         0.292         0.427 1.001  3000
h[5]               0.896   0.278         0.442         0.694         0.866         1.077         1.503 1.002  1500
sigma              0.964   0.032         0.904         0.941         0.963         0.984         1.027 1.012   170
deviance   100007229.959  44.821 100007144.017 100007199.292 100007228.407 100007261.053 100007319.071 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 973.4 and DIC = 100008203.4
DIC is an estimate of expected predictive error (lower deviance is better).
