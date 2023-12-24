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
########  Gauss-Legendre quadrature (15 points)  ########
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes # Nodes
wk <- glq$weights # Weights
K <- length(xk) # K-points
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

     haz[i]<- nu*pow(Time[i],nu-1)*exp(Alpha0[i]+Alpha1[i]*Time[i])
    
    
    #Cumulative hazard function 
    
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature 
      xk1[i,j]<-(xk[j]+1)/2*Time[i] 
      wk1[i,j]<- wk[j]*Time[i]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[i,j]<- nu*pow(xk1[i,j],nu-1)*exp(Alpha0[i]+Alpha1[i]*xk1[i,j])
    }
    
    #Log-survival function with Gauss-Kronrod/Legendre requadrature
    logSurv[i]<- -inprod(wk1[i,],chaz[i,])
 
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
  nu~dgamma(.1,.1)
}", fill = TRUE)
sink()
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
  NbetasL = dim(XL)[3], Nalpha = dim(W)[2],xk=xk,wk=wk,K=K,nmark=nmark
)
parameters <- c(
  "alpha", "gamma", "betaL1", "betaL2", "betaL3", "nu",
  "betaL4", "betaL5",
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
alpha[1]        -3.682   1.513       -5.581       -4.569       -3.985       -3.332        0.706 1.484     9
alpha[2]         1.592   0.537        0.515        1.240        1.621        1.945        2.664 1.115    24
betaL1[1]        0.883   0.225        0.472        0.730        0.882        1.027        1.350 1.078    33
betaL1[2]       -1.037   0.182       -1.404       -1.159       -1.032       -0.920       -0.676 1.099    26
betaL1[3]        0.634   0.199        0.242        0.504        0.641        0.769        1.020 1.015   210
betaL1[4]        0.374   0.186        0.005        0.253        0.376        0.500        0.734 1.018   470
betaL2[1]       -0.523   0.143       -0.758       -0.633       -0.546       -0.411       -0.241 1.708     6
betaL2[2]        0.506   0.087        0.349        0.443        0.504        0.563        0.693 1.020   430
betaL2[3]        0.572   0.128        0.334        0.472        0.573        0.674        0.797 2.294     4
betaL2[4]        0.379   0.143        0.138        0.267        0.372        0.474        0.692 1.166    16
betaL3[1]       -0.157   0.121       -0.442       -0.212       -0.154       -0.087        0.082 1.150    27
betaL3[2]        0.387   0.077        0.242        0.332        0.387        0.441        0.537 1.045    95
betaL3[3]        0.329   0.105        0.140        0.253        0.317        0.410        0.538 1.040    58
betaL3[4]        0.318   0.110        0.104        0.240        0.319        0.396        0.530 1.079    44
betaL4[1]       -0.486   0.127       -0.740       -0.571       -0.482       -0.401       -0.249 1.145    20
betaL4[2]        0.250   0.077        0.107        0.197        0.246        0.299        0.408 1.054    62
betaL4[3]        0.549   0.151        0.279        0.435        0.543        0.658        0.838 1.204    16
betaL4[4]        0.355   0.106        0.160        0.278        0.355        0.429        0.561 1.060    39
betaL5[1]       -0.669   0.119       -0.860       -0.760       -0.687       -0.578       -0.426 1.630     7
betaL5[2]        0.623   0.084        0.453        0.570        0.626        0.680        0.783 1.511     7
betaL5[3]        0.528   0.092        0.321        0.473        0.534        0.592        0.696 1.108    34
betaL5[4]        0.506   0.118        0.289        0.418        0.510        0.587        0.739 1.150    19
gamma[1]        -1.027   0.571       -2.220       -1.427       -0.926       -0.633       -0.017 2.258     5
gamma[2]         0.831   0.261        0.384        0.646        0.799        1.001        1.392 1.496     8
gamma[3]        -0.234   0.177       -0.591       -0.344       -0.234       -0.116        0.105 1.042    61
gamma[4]         1.019   0.246        0.554        0.852        1.009        1.187        1.521 1.199    15
gamma[5]        -0.349   0.196       -0.780       -0.470       -0.338       -0.214        0.007 1.084    44
nu               0.654   0.387        0.001        0.412        0.692        0.925        1.354 1.887     6
sigma[1]         0.570   0.026        0.522        0.552        0.570        0.588        0.623 1.007   290
sigma[2]         0.526   0.025        0.479        0.508        0.526        0.543        0.577 1.047    50
sigma[3]         0.544   0.025        0.497        0.527        0.544        0.561        0.593 1.072    33
sigma[4]         0.533   0.023        0.489        0.517        0.533        0.548        0.581 1.010   220
deviance  50012856.372  70.643 50012710.501 50012813.488 50012863.668 50012904.354 50012978.340 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 2330.4 and DIC = 50015186.8
DIC is an estimate of expected predictive error (lower deviance is better).

