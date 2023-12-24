rm(list = ls())
library(PermAlgo) # Permutation algorithm to generate survival times dependent on time-varying covariates
library(mvtnorm)
library(MASS)
library(R2jags)

nsample <- 500 # number of indivduals

# Y1 (continuous)
Beta1 <- c(1, -1, 0.5, 0.5)
Beta2 <- c(-0.5, 0.5, 0.5, 0.5)
Beta3 <- c(-0.5, 0.5, 0.5, 0.5)
Beta4 <- c(-0.5, 0.5, 0.5, 0.5)
Beta5 <- c(-0.5, 0.5, 0.5, 0.5)
alpha <- 1
sigma <- sqrt(0.5) # residual error
# Survival
gamma <- c(-.5, 0.5, -0.2, 0.5, -0.2)
gapLongi <- 0.2 # gap between longi measurements
gap <- 0.02 # used to generate a lot of time points because the permutation
# algorithm chooses among those time points to define survival times
followup <- 2 # follow-up time
mestime <- seq(0, followup, gap) # measurement times
timesLongi <- mestime[which(round(mestime - round(mestime / gapLongi, 0) * gapLongi, 6) == 0)] # visit times
time <- rep(mestime, nsample) # time column
nmesindiv <- followup / gap + 1 # max. number of individual measurements
nmesy <- nmesindiv * nsample # max. total number of longi measurements
# the number is reduced because of censoring by the terminal event
idY <- rep(1:nsample, each = nmesindiv) # individual id
S <- seq(0, 3.5, length = 4)

# random effects variance and covariance matrix
rho <- .5
rho1 <- 0.1
Sigma <- Matrix::bdiag(
  matrix(c(1, rho, rho, 1), 2, 2), matrix(c(1, rho, rho, 1), 2, 2), matrix(c(1, rho, rho, 1), 2, 2),
  matrix(c(1, rho, rho, 1), 2, 2), matrix(c(1, rho, rho, 1), 2, 2)
)
Sigma[Sigma == 0] <- rho1
Sigma_r <- Sigma <- matrix(Sigma, 10, 10)
############### real
nmark <- 5
u <- rmvnorm(nsample, rep(0, 2 * nmark), Sigma)
b1_int <- rep(u[, 1], each = nmesindiv) # random intercept Y1
b1_slo <- rep(u[, 2], each = nmesindiv) # random slope Y1
b2_int <- rep(u[, 3], each = nmesindiv) # random intercept Y2
b2_slo <- rep(u[, 4], each = nmesindiv) # random slope Y2
b3_int <- rep(u[, 5], each = nmesindiv) # random intercept Y3
b3_slo <- rep(u[, 6], each = nmesindiv) # random intercept Y3
b4_int <- rep(u[, 7], each = nmesindiv) # random intercept Y3
b4_slo <- rep(u[, 8], each = nmesindiv) # random intercept Y3
b5_int <- rep(u[, 9], each = nmesindiv) # random intercept Y3
b5_slo <- rep(u[, 10], each = nmesindiv) # random intercept Y3

x1 <- rnorm(nsample, 1, 0.5)
X1 <- rep(x1, each = nmesindiv) # continuous covariate
x2 <- rbinom(nsample, 1, 0.5)
X2 <- rep(x2, each = nmesindiv) # binary covariate

# linear predictors
linPredY1 <- (Beta1[1] + b1_int) + (Beta1[2] + b1_slo) * time + Beta1[3] * X1 + Beta1[4] * X2
linPredY2 <- (Beta2[1] + b2_int) + (Beta2[2] + b2_slo) * time + Beta2[3] * X1 + Beta2[4] * X2
linPredY3 <- (Beta3[1] + b3_int) + (Beta3[2] + b3_slo) * time + Beta3[3] * X1 + Beta3[4] * X2
linPredY4 <- (Beta4[1] + b4_int) + (Beta4[2] + b4_slo) * time + Beta4[3] * X1 + Beta4[4] * X2
linPredY5 <- (Beta5[1] + b5_int) + (Beta5[2] + b5_slo) * time + Beta5[3] * X1 + Beta5[4] * X2

# binary outcome Y1
Y1 <- rbinom(nmesy, 1, exp(linPredY1) / (1 + exp(linPredY1)))
# continuous outcomes
Y2 <- rnorm(nmesy, linPredY2, sigma)
Y3 <- rnorm(nmesy, linPredY3, sigma)
Y4 <- rnorm(nmesy, linPredY4, sigma)
Y5 <- rnorm(nmesy, linPredY5, sigma)
# Permutation algorithm to generate survival times that depends on the linear predictors
DatTmp <- permalgorithm(nsample, nmesindiv,
  Xmat = matrix(c(linPredY1, linPredY2, linPredY3, linPredY4, linPredY5, X1), nrow = nsample * nmesindiv),
  eventRandom = round(rexp(nsample, 0.003) + 1, 0), # ~40% death
  censorRandom = runif(nsample, 1, nmesindiv), # uniform random censoring
  XmatNames = c("linPredY1", "linPredY2", "linPredY3", "linPredY4", "linPredY5", "X1"), # association
  betas = c(gamma, alpha)
) # association parameters

# extract last line for each Id (= death/censoring time)
DatTmp2 <- DatTmp[c(which(diff(DatTmp[, "Id"]) == 1), dim(DatTmp)[1]), c("Id", "Event", "Stop")]
DatTmp2$deathTimes <- mestime[DatTmp2$Stop + 1] # deathtimes
survDat <- DatTmp2[, c("Id", "deathTimes", "Event")]
DatTmp$time <- mestime[DatTmp$Start + 1] # measurements times of the biomarker
DatTmp$Uid <- paste(DatTmp$Id, DatTmp$time) # unique identifier to match covariates and observed biomarker values
longDat3 <- merge(DatTmp[, c("Uid", "Id", "time")], cbind("Uid" = paste(idY, time), X1, X2, Y1, Y2, Y3, Y4, Y5), by = c("Uid"))
longDat <- sapply(longDat3[longDat3$time %in% timesLongi, -1], as.numeric)
longDat <- as.data.frame(longDat[order(longDat[, "Id"], longDat[, "time"]), ])
summary(survDat) # survival dataset
summary(longDat) # longitudinal dataset
head(longDat)
names(survDat)
names(longDat)
median(table(longDat$Id))
min(table(longDat$Id))
max(table(longDat$Id))
sum(survDat$Event) / length(survDat$Event)
#######################################
# Number of patients and number of longitudinal observations per patient
n <- nsample
M <- table(longDat$Id)
max(M)
# Survival and censoring times
Time <- survDat$deathTimes
death <- survDat$Event # death=1 means observed
survDat$x1 <- x1
# Longitudinal information in matrix format
time1 <- Y1n <- Y2n <- Y3n <- Y4n <- Y5n <- matrix(NA, n, max(M))
for (i in 1:n) {
  time1[i, 1:M[i]] <- longDat$time[longDat$Id == i]
  Y1n[i, 1:M[i]] <- longDat$Y1[longDat$Id == i]
  Y2n[i, 1:M[i]] <- longDat$Y2[longDat$Id == i]
  Y3n[i, 1:M[i]] <- longDat$Y3[longDat$Id == i]
  Y4n[i, 1:M[i]] <- longDat$Y4[longDat$Id == i]
  Y5n[i, 1:M[i]] <- longDat$Y5[longDat$Id == i]
}

treat <- rep(1, n) # Reference = placebo
W <- model.matrix(~x1) # Fixed effects
XL <- array(1, dim = c(n, max(M), 4)) # Fixed effects
XL[, , 2] <- time1
XL[, , 3] <- x1
XL[, , 4] <- x2
ZL <- array(1, dim = c(n, max(M), 2)) # Random effects
ZL[, , 2] <- time1
########  BUGS code  ########
sink("model_file")
cat("model{
  for(i in 1:n){
    #Longitudinalobservations
    for(j in 1:M[i]){
      Y1[i,j]~dbern(p1[i,j])
      logit(p1[i,j])<- mu1[i,j]
      mu1[i,j]<-inprod(betaL1[],XL1[i,j,])+inprod(b[i,1:2],ZL1[i,j,])
      Y2[i,j]~dnorm(mu2[i,j],tau[1])
      mu2[i,j]<-inprod(betaL2[],XL2[i,j,])+inprod(b[i,3:4],ZL2[i,j,])
      Y3[i,j]~dnorm(mu3[i,j],tau[2])
      mu3[i,j]<-inprod(betaL3[],XL3[i,j,])+inprod(b[i,5:6],ZL3[i,j,])
      Y4[i,j]~dnorm(mu4[i,j],tau[3])
      mu4[i,j]<-inprod(betaL4[],XL4[i,j,])+inprod(b[i,7:8],ZL4[i,j,])
      Y5[i,j]~dnorm(mu5[i,j],tau[4])
      mu5[i,j]<-inprod(betaL5[],XL5[i,j,])+inprod(b[i,9:10],ZL5[i,j,])
      
      
    }
    #Survival and censoring times
    #Hazard function
    Alpha0[i]<- inprod(alpha[],W[i,])+gamma[1]*(betaL1[1]+b[i,1])+gamma[2]*(betaL2[1]+b[i,3])+
      gamma[3]*(betaL3[1]+b[i,5])+gamma[4]*(betaL4[1]+b[i,7])+gamma[5]*(betaL5[1]+b[i,9])
    Alpha1[i]<- gamma[1]*(betaL1[2]+b[i,2])+gamma[2]*(betaL2[2]+b[i,4])+gamma[3]*(betaL3[2]+b[i,6])+
    gamma[4]*(betaL4[2]+b[i,8])+gamma[5]*(betaL5[2]+b[i,10])

    haz[i] <- exp(Alpha0[i]+Alpha1[i]*Time[i])
    chaz[i]<- exp(Alpha0[i])/Alpha1[i]*(exp(Alpha1[i]*Time[i])-1)
    
    #Log-survival function log(S)=-H(t) 
    logSurv[i]<- -chaz[i]
    #Definition of the survival log-likelihood using zeros trick
    phi[i]<-100000-death[i]*log(haz[i])-logSurv[i]
    zeros[i]~dpois(phi[i])
    #Random effects
    b[i,1:Nb]~dmnorm(mub[],Omega[,])
  }
  #Prior distributions
  for(l in 1:NbetasL){
    betaL1[l]~dnorm(0,0.001)
    betaL2[l]~dnorm(0,0.001)
    betaL3[l]~dnorm(0,0.001)
    betaL4[l]~dnorm(0,0.001)
    betaL5[l]~dnorm(0,0.001)
  }
  for(k in 1:Nalpha){ 
  alpha[k]~dnorm(0,0.001)
  }

  for(k in 1:nmark){
  gamma[k]~dnorm(0,0.001)
 }

 for(k in 1:(nmark-1)){ 
 tau[k]~dgamma(0.01,0.01)
 sigma[k]<-1/tau[k]
 }

  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])
  lambda0<-exp(alpha[1])
}", fill = TRUE)
sink()
#### Running JAGS
#### Running JAGS
i.jags <- function() {
  list(
    alpha = rnorm(dim(W)[2]), gamma = rnorm(nmark),
    betaL1 = rnorm(dim(XL)[3]), betaL2 = rnorm(dim(XL)[3]), betaL3 = rnorm(dim(XL)[3]),
    betaL4 = rnorm(dim(XL)[3]), betaL5 = rnorm(dim(XL)[3]),
    tau = rep(1,(nmark-1)), Omega = diag(runif(10))
  )
}
d.jags <- list(
  n = nsample / 2, M = M, Time = Time, W = W, Y1 = Y1n, Y2 = Y2n, Y3 = Y3n, Y4 = Y4n, Y5 = Y5n,
  XL1 = XL, ZL1 = ZL, XL2 = XL, ZL2 = ZL, XL3 = XL, ZL3 = ZL, XL4 = XL, ZL4 = ZL, XL5 = XL, ZL5 = ZL,
  death = death, mub = rep(0, 2 * nmark), V = diag(1, 2 * nmark), Nb = 2 * nmark, zeros = rep(0, n),
  NbetasL = dim(XL)[3], Nalpha = dim(W)[2],nmark=nmark
)
parameters <- c(
  "alpha", "gamma", "betaL1", "betaL2", "betaL3", 
  "betaL4", "betaL5","lambda0",
  "sigma"#, "Sigma"
)

sim123 <- jags(
  data = d.jags, inits = i.jags, parameters,
  n.iter = 2000, model.file = "model_file"
)

print(sim123)


##### Results ##### 
> print(sim123)
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 2000 iterations (first 1000 discarded)
n.sims = 3000 iterations saved
mu.vect sd.vect         2.5%          25%          50%          75%        97.5%  Rhat n.eff
alpha[1]        -3.766   0.563       -5.009       -4.120       -3.711       -3.391       -2.789 1.008  1400
alpha[2]         1.656   0.406        0.910        1.375        1.637        1.929        2.475 1.025   180
betaL1[1]        1.085   0.237        0.635        0.918        1.082        1.238        1.553 1.149    18
betaL1[2]       -1.102   0.170       -1.413       -1.220       -1.105       -0.992       -0.739 1.040    61
betaL1[3]        0.500   0.178        0.150        0.375        0.503        0.620        0.846 1.057    48
betaL1[4]        0.473   0.176        0.124        0.360        0.478        0.592        0.806 1.007   410
betaL2[1]       -0.489   0.126       -0.779       -0.564       -0.479       -0.398       -0.280 1.139    23
betaL2[2]        0.541   0.077        0.389        0.489        0.538        0.595        0.694 1.097    26
betaL2[3]        0.578   0.113        0.375        0.500        0.575        0.649        0.822 1.481     8
betaL2[4]        0.439   0.083        0.290        0.382        0.434        0.493        0.610 1.026   110
betaL3[1]       -0.174   0.143       -0.441       -0.261       -0.179       -0.097        0.159 1.130   160
betaL3[2]        0.604   0.071        0.466        0.559        0.600        0.650        0.748 1.129    21
betaL3[3]        0.336   0.120        0.121        0.255        0.321        0.407        0.587 1.095    45
betaL3[4]        0.283   0.134        0.031        0.188        0.278        0.372        0.566 1.397     9
betaL4[1]       -0.541   0.151       -0.776       -0.647       -0.569       -0.443       -0.208 1.392     9
betaL4[2]        0.530   0.086        0.372        0.472        0.527        0.584        0.718 1.070   100
betaL4[3]        0.564   0.095        0.391        0.495        0.558        0.632        0.750 1.111    25
betaL4[4]        0.498   0.121        0.244        0.426        0.505        0.581        0.713 1.301    11
betaL5[1]       -0.583   0.105       -0.802       -0.646       -0.582       -0.519       -0.367 1.087    45
betaL5[2]        0.439   0.096        0.279        0.375        0.427        0.489        0.681 1.075    70
betaL5[3]        0.614   0.128        0.392        0.511        0.610        0.708        0.856 1.013   190
betaL5[4]        0.369   0.095        0.180        0.307        0.373        0.436        0.550 1.142    26
gamma[1]        -0.087   0.189       -0.458       -0.218       -0.089        0.043        0.277 1.129    21
gamma[2]         0.462   0.119        0.232        0.380        0.463        0.542        0.702 1.009   260
gamma[3]        -0.343   0.127       -0.594       -0.425       -0.341       -0.259       -0.099 1.008   330
gamma[4]         0.484   0.127        0.236        0.400        0.482        0.569        0.731 1.009   590
gamma[5]        -0.043   0.122       -0.284       -0.123       -0.045        0.042        0.188 1.027    82
lambda0          0.027   0.015        0.007        0.016        0.024        0.034        0.061 1.008  1400
sigma[1]         0.544   0.027        0.494        0.526        0.543        0.561        0.600 1.015   160
sigma[2]         0.519   0.025        0.471        0.502        0.518        0.535        0.570 1.068    34
sigma[3]         0.544   0.028        0.493        0.524        0.543        0.563        0.603 1.087    28
sigma[4]         0.516   0.026        0.469        0.497        0.514        0.533        0.569 1.023   110
deviance  50011779.417  98.256 50011600.054 50011713.199 50011775.516 50011843.690 50011983.340 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 4388.0 and DIC = 50016167.4
DIC is an estimate of expected predictive error (lower deviance is better).
> 
