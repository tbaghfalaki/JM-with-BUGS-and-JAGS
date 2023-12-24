rm(list = ls())
library(R2jags)
load("/Users/taban/Desktop/Taban/joint modeling bugs/CR/long.data.RData")
load("/Users/taban/Desktop/Taban/joint modeling bugs/CR/CR.data.RData")
#######################################
# Number of patients and number of longitudinal observations per patient
n <- length(CR.data$id)
M <- table(long.data$id)
max(M)
# Survival and censoring times
Time <- CR.data$survtime
CR <- CR.data$CR

# Longitudinal information in matrix format
time1 <- Y1n <- X1n <- X2n <- matrix(NA, n, max(M))
for (i in 1:n) {
  time1[i, 1:M[i]] <- long.data$obstime[long.data$id == i]
  Y1n[i, 1:M[i]] <- long.data$Y1[long.data$id == i]
  X1n[i, 1:M[i]] <- long.data$x1[long.data$id == i]
  X2n[i, 1:M[i]] <- long.data$x2[long.data$id == i]
}

W <- model.matrix(~ CR.data$w1 + CR.data$w2) # Fixed effects

XL <- array(1, dim = c(n, max(M), 4)) # Fixed effects
XL[, , 2] <- time1
XL[, , 3] <- X1n
XL[, , 4] <- X2n
ZL <- array(1, dim = c(n, max(M), 2)) # Random effects
ZL[, , 2] <- time1
########  BUGS code  ########
sink("model_file")
cat("model{
  for(i in 1:n){
    #Longitudinalobservations
    for(j in 1:M[i]){
      Y1[i,j]~dnorm(mu[i,j],tau)
      mu[i,j]<-inprod(betaL1[],XL1[i,j,])+inprod(b[i,1:2],ZL1[i,j,])
         }
    #Survival and censoring times
    # 1th cause
    Alpha01[i]<- inprod(alpha1[],W[i,])+gamma1*(betaL1[1]+betaL1[3]*x1[i]+betaL1[4]*x2[i]+b[i,1])
    Alpha11[i]<- gamma1*(betaL1[2]+b[i,2])
    haz1[i]<- exp(Alpha01[i]+Alpha11[i]*Time[i])
    chaz1[i]<- exp(Alpha01[i])/Alpha11[i]*(exp(Alpha11[i]*Time[i])-1)
    logSurv1[i]<- -chaz1[i]
    
    # 2th cause
    Alpha02[i]<- inprod(alpha2[],W[i,])+gamma2*(betaL1[1]+betaL1[3]*x1[i]+betaL1[4]*x2[i]+b[i,1])
    Alpha12[i]<- gamma2*(betaL1[2]+b[i,2])
    haz2[i]<- exp(Alpha02[i]+Alpha12[i]*Time[i])
    chaz2[i]<- exp(Alpha02[i])/Alpha12[i]*(exp(Alpha12[i]*Time[i])-1)

    logSurv2[i]<- -chaz2[i]   
    #Definition of the survival log-likelihood using zeros trick
    phi[i]<-100000-equals(CR[i],1)*log(haz1[i])-equals(CR[i],2)*log(haz2[i])-logSurv1[i]-logSurv2[i]
    zeros[i]~dpois(phi[i])
    #Random effects
    b[i,1:Nb]~dmnorm(mub[],Omega[,])
  }
  #Prior distributions
  for(l in 1:NbetasL){
    betaL1[l]~dnorm(0,0.001)
  }
  for(k in 1:Nalpha){ 
    alpha1[k]~dnorm(0,0.001)
    alpha2[k]~dnorm(0,0.001)
  }
    gamma1~dnorm(0,0.001)
    gamma2~dnorm(0,0.001)
    
    tau~dgamma(0.01,0.01)
    sigma<-1/tau
  
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])

  lambda01<-exp(alpha1[1])
  lambda02<-exp(alpha2[1])
  
}", fill = TRUE)
sink()
#### Running JAGS
i.jags <- function() {
  list(
    alpha1 = rnorm(dim(W)[2]), alpha2 = rnorm(dim(W)[2]), gamma1 = rnorm(1), gamma2 = rnorm(1),
    betaL1 = rnorm(dim(XL)[3]),
    tau = 1, Omega = diag(runif(2))
  )
}
d.jags <- list(
  n = n, M = M, Time = Time, W = W, Y1 = Y1n,
  XL1 = XL, ZL1 = ZL,
  CR = CR, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
  NbetasL = dim(XL)[3], Nalpha = dim(W)[2], x1 = CR.data$x1,
  x2 = CR.data$x2
)
parameters <- c(
  "alpha1", "alpha2", "gamma1", "gamma2", "betaL1",
  "sigma", "Sigma", "lambda01", "lambda02"
)

library(R2jags)
sim123 <- jags(
  data = d.jags, inits = i.jags, parameters,
  n.iter = 2000, model.file = "model_file"
)
print(sim123)
### Results
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 2000 iterations (first 1000 discarded)
n.sims = 3000 iterations saved
mu.vect sd.vect         2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.016   0.081  8.72000e-01         0.958         1.012         1.067         1.185 1.014   150
Sigma[2,1]         0.082   0.086 -9.70000e-02         0.027         0.086         0.139         0.249 1.110    24
Sigma[1,2]         0.082   0.086 -9.70000e-02         0.027         0.086         0.139         0.249 1.110    24
Sigma[2,2]         0.780   0.121  5.63000e-01         0.695         0.773         0.858         1.040 1.011   490
alpha1[1]         -0.003   0.152 -3.07000e-01        -0.100         0.002         0.104         0.295 1.022   110
alpha1[2]          0.959   0.156  6.62000e-01         0.856         0.953         1.057         1.297 1.009   230
alpha1[3]         -0.905   0.078 -1.06300e+00        -0.959        -0.903        -0.854        -0.749 1.001  3000
alpha2[1]         -0.040   0.188 -4.09000e-01        -0.172        -0.033         0.088         0.311 1.003   910
alpha2[2]         -1.187   0.244 -1.68400e+00        -1.348        -1.183        -1.023        -0.716 1.001  3000
alpha2[3]          1.005   0.132  7.49000e-01         0.918         1.006         1.095         1.269 1.004   570
betaL1[1]          0.837   0.074  7.06000e-01         0.781         0.838         0.884         0.988 1.052    49
betaL1[2]          1.228   0.088  1.05700e+00         1.168         1.226         1.290         1.404 1.092    27
betaL1[3]          1.102   0.108  9.02000e-01         1.023         1.091         1.189         1.298 1.143    21
betaL1[4]          0.989   0.061  8.82000e-01         0.945         0.984         1.029         1.117 1.121    23
gamma1            -0.466   0.046 -5.55000e-01        -0.497        -0.467        -0.435        -0.377 1.026    82
gamma2            -0.543   0.069 -6.76000e-01        -0.591        -0.542        -0.496        -0.405 1.002  1100
lambda01           1.009   0.153  7.36000e-01         0.905         1.002         1.109         1.343 1.022   110
lambda02           0.977   0.184  6.64000e-01         0.842         0.968         1.092         1.365 1.003   910
sigma              0.475   0.015  4.47000e-01         0.465         0.475         0.485         0.504 1.005   470
deviance   100006110.029  43.811  1.00006e+08 100006079.006 100006111.122 100006140.598 100006196.627 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 945.7 and DIC = 100007055.8
DIC is an estimate of expected predictive error (lower deviance is better).
