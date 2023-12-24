rm(list=ls())
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
time1 =Y1n=X1n=X2n=matrix(NA, n, max(M))
for(i in 1:n) {
  time1[i, 1:M[i]] <- long.data$obstime[long.data$id==i]
  Y1n[i, 1:M[i]] <- long.data$Y1[long.data$id==i]
  X1n[i, 1:M[i]] <- long.data$x1[long.data$id==i]
  X2n[i, 1:M[i]] <- long.data$x2[long.data$id==i]
}

W <- model.matrix(~ CR.data$w1 + CR.data$w2) # Fixed effects

XL <- array(1, dim = c(n, max(M), 4)) # Fixed effects
XL[, , 2] <- time1; 
XL[, , 3] <- X1n; 
XL[, , 4] <- X2n; 
ZL <- array(1, dim = c(n, max(M), 2)) # Random effects
ZL[, , 2] <- time1
########  Gauss-Legendre quadrature (15 points)  ########
glq <- statmod::gauss.quad(15, kind = "legendre")
xk <- glq$nodes # Nodes
wk <- glq$weights # Weights
K <- length(xk) # K-points
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
    haz1[i]<- nu1*pow(Time[i],nu1-1)*exp(Alpha01[i]+Alpha11[i]*Time[i])
    for(j in 1:K){
      xk1[i,j]<-(xk[j]+1)/2*Time[i] 
      wk1[i,j]<- wk[j]*Time[i]/2
      chaz1[i,j]<- nu1*pow(xk1[i,j],nu1-1)*exp(Alpha01[i]+Alpha11[i]*xk1[i,j])
      chaz2[i,j]<- nu2*pow(xk1[i,j],nu2-1)*exp(Alpha02[i]+Alpha12[i]*xk1[i,j])
    }
    logSurv1[i]<- -inprod(wk1[i,],chaz1[i,])
    # 2th cause
    Alpha02[i]<- inprod(alpha2[],W[i,])+gamma2*(betaL1[1]+betaL1[3]*x1[i]+betaL1[4]*x2[i]+b[i,1])
    Alpha12[i]<- gamma2*(betaL1[2]+b[i,2])
    haz2[i]<- nu2*pow(Time[i],nu2-1)*exp(Alpha02[i]+Alpha12[i]*Time[i])
    logSurv2[i]<- -inprod(wk1[i,],chaz2[i,])    
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
  nu1~dgamma(.1,.1)
  nu2~dgamma(.1,.1)
  
}", fill = TRUE)
sink()
#### Running JAGS
i.jags <- function() {
  list(alpha1 = rnorm(dim(W)[2]),alpha2 = rnorm(dim(W)[2]),  gamma1= rnorm(1), gamma2= rnorm(1), 
       betaL1 = rnorm(dim(XL)[3]),
       tau = 1, Omega = diag(runif(2)))
}
d.jags <- list(n = n, M = M, Time = Time, W = W, Y1 = Y1n, 
               XL1 =XL, ZL1 = ZL,
               CR = CR, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0,n),
               NbetasL = dim(XL)[3],Nalpha=dim(W)[2],xk=xk,wk=wk,K=K,x1=CR.data$x1,
               x2=CR.data$x2)
parameters <- c("alpha1","alpha2", "gamma1", "gamma2",  "betaL1",
                "sigma",  "Sigma")
#file.show(model.file)
library(R2jags)
sim123 <- jags(data=d.jags, inits=i.jags, parameters,
               n.iter=2000, model.file="model_file")
print(sim123)
### Results
Inference for Bugs model at "model_file", fit using jags,
3 chains, each with 2000 iterations (first 1000 discarded)
n.sims = 3000 iterations saved
mu.vect sd.vect         2.5%           25%           50%           75%         97.5%  Rhat n.eff
Sigma[1,1]         1.015   0.079  8.71000e-01         0.957         1.012         1.068         1.176 1.005   410
Sigma[2,1]         0.088   0.092 -9.70000e-02         0.026         0.089         0.152         0.255 1.055    44
Sigma[1,2]         0.088   0.092 -9.70000e-02         0.026         0.089         0.152         0.255 1.055    44
Sigma[2,2]         0.744   0.138  5.26000e-01         0.646         0.727         0.827         1.071 1.055    41
alpha1[1]          0.065   0.157 -2.48000e-01        -0.039         0.064         0.173         0.377 1.014   150
alpha1[2]          0.965   0.153  6.66000e-01         0.859         0.965         1.066         1.268 1.005   500
alpha1[3]         -0.924   0.082 -1.08700e+00        -0.981        -0.924        -0.870        -0.761 1.002  1500
alpha2[1]         -0.028   0.219 -4.57000e-01        -0.171        -0.029         0.125         0.388 1.001  3000
alpha2[2]         -1.182   0.245 -1.68300e+00        -1.343        -1.173        -1.012        -0.716 1.002  1200
alpha2[3]          1.024   0.128  7.77000e-01         0.942         1.018         1.108         1.276 1.002  1100
betaL1[1]          0.848   0.082  7.04000e-01         0.785         0.840         0.912         1.004 1.198    15
betaL1[2]          1.204   0.088  1.01500e+00         1.155         1.204         1.256         1.396 1.033    88
betaL1[3]          1.116   0.111  9.02000e-01         1.032         1.125         1.196         1.325 1.267    12
betaL1[4]          0.987   0.046  8.92000e-01         0.956         0.987         1.020         1.076 1.031    85
gamma1            -0.487   0.052 -5.83000e-01        -0.522        -0.488        -0.451        -0.382 1.003   780
gamma2            -0.554   0.076 -6.93000e-01        -0.607        -0.556        -0.504        -0.402 1.007   290
sigma              0.476   0.015  4.46000e-01         0.466         0.476         0.486         0.508 1.006   400
deviance   100006118.537  45.020  1.00006e+08 100006087.159 100006116.553 100006147.897 100006208.303 1.000     1

For each parameter, n.eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor (at convergence, Rhat=1).

DIC info (using the rule, pD = var(deviance)/2)
pD = 999.2 and DIC = 100007117.8
DIC is an estimate of expected predictive error (lower deviance is better).