rm(list = ls()) 
library(R2OpenBUGS)

load(file = "C:/Users/BitBin/Downloads/archive (1)/long.data_sre_1.RData")
load(file = "C:/Users/BitBin/Downloads/archive (1)/surv.data_sre_1.RData")

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
    st_n[i] ~ dweib(nu,mut[i])I(st_cen[i],)
    log(mut[i])<-inprod(alpha[],W[i,])+inprod(b[i,],gamma[])
    
    b[i,1:Nb]~dmnorm(mub[],Omega[,])
  }
  #Prior distributions
  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  for(l in 1:2){ 
    gamma[l]~dnorm(0,0.001)
  }
  nu~dgamma(0.01,0.01)
  tau~dgamma(0.01,0.01)
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  #Determinetic nodes
  lambda0<-exp(alpha[1])
  sigma<-1/tau
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])
}", fill = TRUE)
sink()
#### Running R2Openbugs
st_n=st_cen=rep(0,n)
for(i in 1:n){
  if(death[i]==0)((st_n[i]=NA) & (st_cen[i]=st[i]))
  if(death[i]==1)((st_n[i]=st[i]) & (st_cen[i]=0))
}

mub = rep(0, 2)
V = diag(1, 2)
Nb = 2
Nbeta = dim(X)[3]
Nalpha = ncol(W)
M=as.numeric(M)
data <- list("n",  "M", "W", "Y", "X","Z",
             "mub" , "V" , "Nb",
             "Nbeta" , "Nalpha" ,"st_n","st_cen")



inits <- function() {
  list(
    alpha = rnorm(ncol(W)), gamma = rnorm(2),nu=1,
    beta = rnorm(dim(X)[3]), tau = runif(1), Omega = diag(runif(2))
  )
}


parameters <- c("alpha","beta", "gamma", "beta", "sigma", "Sigma", "nu")


sim_openbugs <- bugs(data, inits, parameters, model.file="model_file",
                     n.chains=2, n.iter=2000,debug=TRUE)
print(sim_openbugs)


Inference for Bugs model at "model_file", 
Current: 2 chains, each with 2000 iterations (first 1000 discarded)
Cumulative: n.sims = 2000 iterations saved
mean   sd   2.5%    25%    50%    75%  97.5% Rhat n.eff
alpha[1]      1.2  0.1    0.9    1.1    1.2    1.3    1.5  1.1    75
alpha[2]     -1.2  0.1   -1.4   -1.2   -1.2   -1.1   -0.9  1.1   100
alpha[3]     -0.6  0.1   -0.9   -0.7   -0.6   -0.5   -0.4  1.0  1300
beta[1]      -0.5  0.1   -0.7   -0.6   -0.5   -0.5   -0.4  1.1    32
beta[2]       0.6  0.1    0.3    0.5    0.6    0.7    0.9  1.1   160
beta[3]       0.6  0.1    0.4    0.6    0.6    0.7    0.8  1.1    26
beta[4]       0.4  0.1    0.3    0.4    0.4    0.5    0.5  1.0    64
gamma[1]      1.2  0.4    0.5    0.9    1.2    1.4    2.0  1.1    55
gamma[2]     -1.2  0.5   -2.3   -1.5   -1.2   -0.9   -0.5  1.2    19
beta[1]      -0.5  0.1   -0.7   -0.6   -0.5   -0.5   -0.4  1.1    32
beta[2]       0.6  0.1    0.3    0.5    0.6    0.7    0.9  1.1   160
beta[3]       0.6  0.1    0.4    0.6    0.6    0.7    0.8  1.1    26
beta[4]       0.4  0.1    0.3    0.4    0.4    0.5    0.5  1.0    64
sigma         1.0  0.0    0.9    1.0    1.0    1.1    1.1  1.0   190
Sigma[1,1]    0.8  0.1    0.7    0.8    0.8    0.9    1.0  1.0   360
Sigma[1,2]    0.5  0.1    0.3    0.5    0.6    0.6    0.7  1.0   510
Sigma[2,1]    0.5  0.1    0.3    0.5    0.6    0.6    0.7  1.0   510
Sigma[2,2]    0.8  0.2    0.6    0.7    0.8    0.9    1.3  1.0    45
nu            1.1  0.1    1.0    1.1    1.1    1.2    1.3  1.1    35
deviance   5003.9 73.3 4851.0 4954.0 5006.0 5051.0 5144.0  1.1    21

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = Dbar-Dhat)
pD = 509.5 and DIC = 5513.0
DIC is an estimate of expected predictive error (lower deviance is better).
> 