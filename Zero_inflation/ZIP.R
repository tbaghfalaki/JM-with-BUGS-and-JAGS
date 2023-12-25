rm(list = ls())
library(mvtnorm)
library(extraDistr)
library(R2jags)
library(statmod)
set.seed(9)
n <- 1000
m <- 9
Nobs <- n * m
x1 <- rnorm(n)
x2 <- rbinom(n, 1, .5)
D <- matrix(c(1, .5, .5, 1), 2, 2)
b <- rmvnorm(n, c(0, 0), D)

y <- z <- muy <- muz <- matrix(0, n, m)
Beta1 <- c(-1, -1, 1, -1)
Beta2 <- c(-1, 1, -1) 
follow_up <- 1
t <- seq(from = 0, to = follow_up, length = m)

for (i in 1:n) {
  for (j in 1:m) {
    muy[i, j] <- Beta1[1] + Beta1[2] * t[j] + Beta1[3] * x1[i] + Beta1[4] * x2[i] + b[i, 1]
    muz[i, j] <- Beta2[1] + Beta2[2] * t[j] + Beta2[3] * x1[i] + b[i, 2]

    z[i, j] <- rbinom(1, 1, exp(muz[i, j]) / (1 + exp(muz[i, j])))
  }
}

for (i in 1:n) {
  for (j in 1:m) {
    if (z[i, j] == 1) (y[i, j] <- 0)
    if (z[i, j] == 0) (y[i, j] <- rtpois(1, exp(muy[i, j]), a = 0, b = Inf))
  }
}
table(z)
gamma <- c(-1, 1)
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

scale <- 3
shape <- 1
td <- rweibull(n, shape, scale)
mean(td)
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
colnames(y) <- c("y1", "y2", "y3", "y4", "y5", "y6", "y7", "y8", "y9")
y

table(M)/n
mean(M)

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
   
        ll[i,j]<-z[i,j]*log(pi[i,j]) + (1-z[i,j])*(log(1-pi[i,j])+y[i,j]*log(lambda[i,j])- 
                lambda[i,j] - loggam(y[i,j]+1)-log(1-exp(-lambda[i,j])))

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
  n = n, t = t, x1 = x1, x2 = x2, R = R, M = M, Nbeta1 = 4, Nbeta2 = 3, Nalpha = 2, Ngamma = 2,
  y = y, z = z, surt = surt, c = c, U0 = U0, is.censored = is.censored,
  zeros = matrix(0, n, m)
)
parameters <- c("beta1", "beta2", "alpha", "gamma", "sigma", "omega")

# file.show(model.file)
sim123 <- jags(
  data = d.jags, inits = i.jags, parameters,
  n.iter = 10000, model.file = "model_file"
)

print(sim123)

### Results
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 10000 iterations (first 5000 discarded), n.thin = 5
n.sims = 3000 iterations saved
mu.vect sd.vect        2.5%         25%         50%         75%       97.5%  Rhat n.eff
alpha[1]         0.464   0.069       0.340       0.418       0.459       0.507       0.611 1.012   180
alpha[2]        -0.536   0.060      -0.654      -0.575      -0.535      -0.493      -0.423 1.008   280
beta1[1]        -1.110   0.107      -1.322      -1.181      -1.109      -1.037      -0.907 1.016   130
beta1[2]        -1.109   0.102      -1.309      -1.179      -1.111      -1.040      -0.903 1.004   590
beta1[3]         1.039   0.081       0.889       0.983       1.037       1.092       1.203 1.022   120
beta1[4]        -1.045   0.125      -1.300      -1.127      -1.042      -0.961      -0.802 1.001  3000
beta2[1]        -0.920   0.066      -1.044      -0.964      -0.922      -0.874      -0.788 1.001  3000
beta2[2]         0.596   0.150       0.297       0.495       0.597       0.700       0.875 1.004   560
beta2[3]        -0.973   0.061      -1.093      -1.013      -0.973      -0.932      -0.856 1.002  1300
gamma[1]        -0.802   0.131      -1.072      -0.889      -0.798      -0.706      -0.572 1.023   110
gamma[2]         0.750   0.141       0.501       0.650       0.741       0.839       1.067 1.035    65
omega            0.945   0.052       0.853       0.909       0.941       0.979       1.059 1.014   180
sigma[1,1]       1.178   0.147       0.923       1.075       1.165       1.270       1.510 1.006   350
sigma[2,1]       0.361   0.102       0.170       0.292       0.358       0.428       0.575 1.006   360
sigma[1,2]       0.361   0.102       0.170       0.292       0.358       0.428       0.575 1.006   360
sigma[2,2]       0.936   0.133       0.698       0.842       0.931       1.024       1.204 1.004   840
deviance   8137425.047  81.413 8137253.032 8137372.884 8137428.709 8137481.707 8137570.399 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 3296.3 and DIC = 8140721.4
DIC is an estimate of expected predictive error (lower deviance is better).
