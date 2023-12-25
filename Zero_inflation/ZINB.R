rm(list = ls())
library(mvtnorm)
library(countreg)
library(R2jags)
library(statmod)

n <- 1000
m <- 9
Nobs <- n * m
x1 <- rnorm(n)
x2 <- rbinom(n, 1, .5)
D <- matrix(c(1, .5, .5, 1), 2, 2)
b <- rmvnorm(n, c(0, 0), D)

y <- z <- muy <- muz <- matrix(0, n, m)
Beta1 <- c(2, -1, 1, -1)
Beta2 <- c(1, 1, -1)
follow_up <- 1
t <- seq(from = 0, to = follow_up, length = m)

for (i in 1:n) {
  for (j in 1:m) {
    muy[i, j] <- Beta1[1] + Beta1[2] * t[j] + Beta1[3] * x1[i] + Beta1[4] * x2[i] + b[i, 1]
    muz[i, j] <- Beta2[1] + Beta2[2] * t[j] + Beta2[3] * x1[i] + b[i, 2]

    z[i, j] <- rbinom(1, 1, exp(muz[i, j]) / (1 + exp(muz[i, j])))
  }
}
theta <- 3

for (i in 1:n) {
  for (j in 1:m) {
    if (z[i, j] == 1) (y[i, j] <- 0)
    if (z[i, j] == 0) (y[i, j] <- rztnbinom(1, mu = exp(muy[i, j]), theta = theta))
  }
}
table(z)
gamma=c(-1,1)
alpha <- c(0.5, -0.5)
mut <- surt <- rep(0, n)
omega <- 1
for (i in 1:n) {
  mut[i] <- (alpha[1] + alpha[2] * x1[i] + gamma[1] * b[i, 1] + gamma[2] * b[i, 2])
}

for (i in 1:n) {
  surt[i] <- rweibull(1, omega, exp(-mut[i] / omega))
}
is.censored <- rep(0, n)

scale <- 2
shape <- 1
td <- rweibull(n, shape, scale)

for (i in 1:n) {
  if (surt[i] > td[i]) ((surt[i] <- td[i]) & (is.censored[i] <- 1))
}
table(is.censored) / n
M <- c()
for (i in 1:n) {
  M[i] <- 0
  for (j in 1:m) {
    if ((surt[i] > t[j])) (M[i] <- M[i] + 1)
  }
  if (M[i] < m) ((y[i, (M[i] + 1):m] <- NA) & (z[i, (M[i] + 1):m] <- NA))
}

 colnames(y)<-c("y1","y2","y3","y4","y5","y6","y7","y8","y9")
y


Data <- cbind(y, surt, is.censored, x1, x2)
head(Data)
#################################
library(R2jags)
n <- length(y[, 1])
m <- length(y[1, ])

U0 <- c(0, 0)
R <- structure(.Data = c(1, 0, 0, 1), .Dim = c(2, 2))

c <- rep(0, n)
for (i in 1:n) {
  if (is.censored[i] == 1) ((c[i] <- surt[i]) & (surt[i] <- NA))
  if (is.censored[i] == 0) ((c[i] <- 1000))
}
cbind(surt, c, is.censored)
z <- matrix(0, n, m)
z[y == 0] <- 1

########  BUGS code  ########
sink("model_file")
cat("model{
  K<-1000
  for (i in 1:n) {
    for (j in 1:M[i]) {
      zeros[i,j]~dpois(phi[i,j])
      phi[i,j]<-  - ll[i,j]+K				
      ll[i,j]<-z[i,j]*log(pi[i,j]) +
        (1-z[i,j])*(log(1-pi[i,j])+loggam(r+y[i,j])-loggam(r)-loggam(y[i,j]+1)+
                      r*log(r/(r+lambda[i, j]))+
                      y[i,j]*log(lambda[i, j]/(lambda[i, j]+r))-log(1-pow(r/(r+lambda[i, j]),r)))
      
      log(lambda[i, j])<-beta1[1]+beta1[2]*t[j]+beta1[3]*x1[i]+beta1[4]*x2[i]+
        U[i,1]
      logit(pi[i, j])<-beta2[1]+beta2[2]*t[j]+beta2[3]*x1[i]+U[i,2]
    }
    surt[i] ~ dweib(omega,mut[i])    
    is.censored[i]~dinterval(surt[i],c[i])
    log(mut[i])<-alpha[1]+alpha[2]*x1[i]+gamma[1]*U[i, 1]+gamma[2]*U[i, 2]
    U[i,1:2] ~ dmnorm(U0[],tau[,])   
  }  
  r~ dgamma(0.1,0.1) 
  omega ~ dgamma(0.1,0.1) 
  sigma[1:2,1:2]<-inverse(tau[,])
  
  #priors
  tau[1:2,1:2] ~ dwish(R[,], 2)
  for(k in 1:Nbeta1){ 
    beta1[k]~dnorm(0, 0.0001)}
  for(k in 1:Nbeta2){ 
    beta2[k]~dnorm(0, 0.0001)}
  for(k in 1:Nalpha){ 
    alpha[k]~dnorm(0, 0.0001)}
  for(k in 1:Ngamma){ 
    gamma[k]~dnorm(0, 0.0001)}
}", fill = TRUE)
sink()
####################### MULTI
i.jags <- function() {
  list(
    beta1 = rep(0, 4), beta2 = rep(0, 3), alpha = rep(0, 2), gamma = c(1, 1), omega = 1, r = 1,
    U = matrix(0, n, 2)
  )
}
d.jags <- list(
  n = n, t = t, x1 = x1, x2 = x2, R = R, M = M,Nbeta1=4,Nbeta2=3,Nalpha=2,Ngamma=2,
  y = y, z = z, surt = surt, c = c, U0 = U0, is.censored = is.censored,
  zeros = matrix(0, n, m)
)
parameters <- c("beta1", "beta2", "alpha", "gamma", "sigma", "omega", "r")

# file.show(model.file)
sim123 <- jags(
  data = d.jags, inits = i.jags, parameters,
  n.iter = 3000, model.file = "model_file"
)

print(sim123)

### Results
n=500
Inference for Bugs model at "/Users/taban/Desktop/bugs_paper/ZINB/zinb_bugs.R", fit using jags,
3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
n.sims = 3750 iterations saved
mu.vect sd.vect        2.5%         25%         50%         75%       97.5%  Rhat n.eff
alpha[1]         0.461   0.097       0.286       0.393       0.457       0.523       0.664 1.043    53
alpha[2]        -0.430   0.081      -0.602      -0.482      -0.427      -0.374      -0.276 1.008   370
beta1[1]         1.896   0.130       1.639       1.808       1.895       1.982       2.157 1.047    61
beta1[2]        -0.953   0.125      -1.202      -1.035      -0.952      -0.875      -0.706 1.002  1100
beta1[3]         1.058   0.092       0.868       1.000       1.058       1.117       1.239 1.039   110
beta1[4]        -1.010   0.130      -1.259      -1.101      -1.010      -0.920      -0.748 1.027    84
beta2[1]         1.030   0.105       0.825       0.959       1.030       1.101       1.238 1.003   700
beta2[2]         0.901   0.212       0.490       0.756       0.895       1.044       1.320 1.021   120
beta2[3]        -1.014   0.091      -1.196      -1.076      -1.010      -0.950      -0.840 1.007   360
gamma[1]        -0.918   0.194      -1.309      -1.058      -0.913      -0.773      -0.577 1.013   210
gamma[2]         0.881   0.214       0.540       0.718       0.855       1.018       1.344 1.060    39
omega            1.020   0.081       0.878       0.962       1.014       1.073       1.193 1.043    52
r                3.413   0.387       2.685       3.138       3.396       3.673       4.214 1.004   520
sigma[1,1]       0.984   0.135       0.747       0.889       0.973       1.066       1.279 1.026    89
sigma[2,1]       0.385   0.144       0.103       0.291       0.386       0.478       0.677 1.013   990
sigma[1,2]       0.385   0.144       0.103       0.291       0.386       0.478       0.677 1.013   990
sigma[2,2]       1.053   0.221       0.687       0.900       1.027       1.186       1.552 1.012   170
deviance   3929377.228  55.047 3929267.658 3929339.711 3929376.879 3929414.081 3929484.957 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 1468.3 and DIC = 3930845.5
DIC is an estimate of expected predictive error (lower deviance is better).

n=1000
print(sim123)
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 3000 iterations (first 1500 discarded)
n.sims = 4500 iterations saved
mu.vect sd.vect        2.5%         25%         50%         75%       97.5%  Rhat n.eff
alpha[1]         0.401   0.072       0.270       0.351       0.399       0.444       0.568 1.089    28
alpha[2]        -0.474   0.062      -0.595      -0.516      -0.472      -0.432      -0.354 1.024    91
beta1[1]         2.090   0.083       1.906       2.038       2.092       2.147       2.243 1.048    49
beta1[2]        -0.962   0.092      -1.137      -1.024      -0.965      -0.901      -0.778 1.007   720
beta1[3]         0.852   0.058       0.733       0.815       0.854       0.893       0.958 1.010   210
beta1[4]        -1.050   0.090      -1.221      -1.110      -1.052      -0.989      -0.874 1.011   190
beta2[1]         0.914   0.072       0.773       0.867       0.914       0.963       1.057 1.005   460
beta2[2]         1.175   0.148       0.886       1.074       1.172       1.274       1.457 1.023    93
beta2[3]        -1.040   0.070      -1.183      -1.086      -1.038      -0.992      -0.911 1.063    43
gamma[1]        -1.018   0.163      -1.388      -1.103      -1.001      -0.907      -0.742 1.178    17
gamma[2]         1.200   0.232       0.810       1.026       1.156       1.376       1.714 1.146    20
omega            0.993   0.060       0.890       0.954       0.989       1.025       1.136 1.110    26
r                3.199   0.282       2.681       3.003       3.184       3.382       3.789 1.015   140
sigma[1,1]       1.083   0.097       0.904       1.016       1.081       1.146       1.284 1.005   520
sigma[2,1]       0.613   0.098       0.426       0.547       0.607       0.675       0.822 1.010   240
sigma[1,2]       0.613   0.098       0.426       0.547       0.607       0.675       0.822 1.010   240
sigma[2,2]       1.010   0.181       0.697       0.875       1.001       1.132       1.386 1.100    27
deviance   7822690.522  81.337 7822521.918 7822640.122 7822693.014 7822746.982 7822844.318 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 3008.6 and DIC = 7825699.1
DIC is an estimate of expected predictive error (lower deviance is better).
