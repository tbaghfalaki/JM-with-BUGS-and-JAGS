rm(list = ls())

# load packages
library(joineRML)
library(R2jags)

# load data
data(pbc2)
attach(pbc2)
names(pbc2)
############# data_preparation
n <- length(unique(pbc2$id))
Surv1 <- Stat <- Age <- Drug <- id <- matrix(NA, n, 20)
ni <- c()
for (i in 1:n) {
  id[i, 1] <- i
  ni[i] <- length(pbc2$id[pbc2$id == i])
  Surv1[i, 1:ni[i]] <- pbc2$years[pbc2$id == i]
  Stat[i, 1:ni[i]] <- pbc2$status2[pbc2$id == i]
  Age[i, 1:ni[i]] <- pbc2$age[pbc2$id == i]
  Drug[i, 1:ni[i]] <- pbc2$drug[pbc2$id == i]
}
years <- Surv1[, 1]
status2 <- Stat[, 1]
age <- Age[, 1]
drug <- Drug[, 1]
id <- id[, 1]
long.data <- data.frame(id, years, status2, age, drug)

Time <- long.data$years
death <- long.data$status2 # death=1 means observed
table(death)
Age <- (long.data$age - mean(long.data$age)) / sd(long.data$age)
treat <- long.data$drug - 1 # Reference = placebo
######
Albumin1 <- pbc2$albumin
Year1 <- pbc2$year
ID <- pbc2$id
Alkaline1 <- log(pbc2$alkaline)
Sgot1 <- log(pbc2$SGOT)
Platelets1 <- log(pbc2$platelets)
Prothrombin1 <- log(pbc2$prothrombin)
SerBilir1 <- log(pbc2$serBilir)

ni <- M <- table(ID)
n <- dim(table(ID))
Albumin <- Alkaline <- SGOT1 <- Platelets <- Prothrombin <- SerBilir <- Year <- matrix(NA, n, max(ni))
NUM <- ID
unqNUM <- unique(ID)

for (i in 1:n) {
  Albumin[i, 1:ni[i]] <- Albumin1[NUM == unqNUM[i]]
  Alkaline[i, 1:ni[i]] <- Alkaline1[NUM == unqNUM[i]]
  SGOT1[i, 1:ni[i]] <- Sgot1[NUM == unqNUM[i]]
  Platelets[i, 1:ni[i]] <- Platelets1[NUM == unqNUM[i]]
  Prothrombin[i, 1:ni[i]] <- Prothrombin1[NUM == unqNUM[i]]
  SerBilir[i, 1:ni[i]] <- SerBilir1[NUM == unqNUM[i]]
  Year[i, 1:ni[i]] <- Year1[NUM == unqNUM[i]]
}
#######################################
W <- matrix(1, n, 3) # Fixed effects
W[, 2] <- treat
W[, 3] <- Age

XL <- array(1, dim = c(n, max(M), 2)) # Fixed effects
XL[, , 2] <- Year
ZL <- array(1, dim = c(n, max(M), 2)) # Random effects
ZL[, , 2] <- Year
########  Gauss-Legendre quadrature (15 points)  ########
glq <- statmod::gauss.quad(15, kind = "legendre")
xk <- glq$nodes # Nodes
wk <- glq$weights # Weights
K <- length(xk) # K-points
############ JM
sink("model_file")
cat("model{
  for(i in 1:n){
    #Longitudinalobservations
    for(j in 1:M[i]){
      Y1[i,j]~dnorm(mu1[i,j],tau[1])
      mu1[i,j]<-inprod(beta1[],XL1[i,j,])+inprod(b[i,1:2],ZL1[i,j,])
      Y2[i,j]~dnorm(mu2[i,j],tau[2])
      mu2[i,j]<-inprod(beta2[],XL1[i,j,])+inprod(b[i,3:4],ZL1[i,j,])
      Y3[i,j]~dnorm(mu3[i,j],tau[3])
      mu3[i,j]<-inprod(beta3[],XL1[i,j,])+inprod(b[i,5:6],ZL1[i,j,])
      Y4[i,j]~dnorm(mu4[i,j],tau[4])
      mu4[i,j]<-inprod(beta4[],XL1[i,j,])+inprod(b[i,7:8],ZL1[i,j,])
      Y5[i,j]~dnorm(mu5[i,j],tau[5])
      mu5[i,j]<-inprod(beta5[],XL1[i,j,])+inprod(b[i,9:10],ZL1[i,j,])
    }
    #Survival and censoring times
    #Hazard function
    Alpha0[i]<- inprod(alpha,W[i,])+gamma[1]*(beta1[1]+b[i,1])+
      gamma[2]*(beta2[1]+b[i,3])+
      gamma[3]*(beta3[1]+b[i,5])+
      gamma[4]*(beta4[1]+b[i,7])+
      gamma[5]*(beta5[1]+b[i,9])
    
    Alpha1[i]<- gamma[1]*(beta1[2]+b[i,2])+
      gamma[2]*(beta2[2]+b[i,4])+
      gamma[3]*(beta3[2]+b[i,6])+
      gamma[4]*(beta4[2]+b[i,8])+
      gamma[5]*(beta5[2]+b[i,10])
    
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
    phi[i]<-1000-death[i]*log(haz[i])-logSurv[i]
    zeros[i]~dpois(phi[i])
    #Random effects
    b[i,1:Nb]~dmnorm(mub[],Omega[,])
  }
  #Prior distributions
  for(l in 1:NbetasL){
    beta1[l]~dnorm(0,0.001)
    beta2[l]~dnorm(0,0.001)
    beta3[l]~dnorm(0,0.001)
    beta4[l]~dnorm(0,0.001)
    beta5[l]~dnorm(0,0.001)
  }
  
  for(l in 1:5){ 
    gamma[l]~dnorm(0,0.001)
    tau[l]~dgamma(0.01,0.01)
    sigma[l]<-1/tau[l]
  }  
  
  for(l in 1:Nbetas){
    alpha[l]~dnorm(0,0.001)
  }
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])
  nu~dgamma(0.01,0.01)
}", fill = TRUE)
sink()

# Run ‘JAGS’ from R using R2jags package
# data as list 
d.jags <- list(
  n = n, M = M, Time = Time, Y1 = Albumin, Y2 = Alkaline,
  Y3 = SGOT1, Y4 = Platelets, Y5 = SerBilir,
  XL1 = XL, ZL1 = ZL, W = W,
  death = death, mub = rep(0, 10), V = diag(1, 10), Nb = 10, zeros = rep(0, n),
  NbetasL = dim(XL)[3], Nbetas = dim(W)[2], xk = xk, wk = wk, K = K
)
# Initial value 
i.jags <- function() {
  list(gamma = rnorm(5),beta1 = rnorm(dim(XL)[3]),beta2 = rnorm(dim(XL)[3]), 
       beta3 = rnorm(dim(XL)[3]),beta4 = rnorm(dim(XL)[3]),
       beta5 = rnorm(dim(XL)[3]), 
       alpha = rnorm(dim(W)[2]), tau = rep(1, 5), Omega = diag(runif(10))
  )
}
# Parameters of interest
parameters <- c(
  "beta1", "beta2", "beta3", "beta4", "beta5", "alpha",
  "gamma", "nu", "sigma", "Sigma"
)

# The main fumction for rung BUGS code
jm <- jags(
  data = d.jags, inits = i.jags, parameters,
  n.chains=3, n.iter = 100000, 
  n.burnin=80000,
  n.thin=10,model.file = "model_file"
)
print(jm)
#save(jm, file="JMpbc_newfeb1.RData")

## Some MCMC plots
library(mcmcplots)
mcmc_output=as.mcmc(jm)
traplot(mcmc_output, parms = c("beta1","beta2","beta3","beta4","beta5","alpha","gamma"))#,auto.layout =FALSE)
denplot(mcmc_output, parms = c("beta1","beta2","beta3","beta4","beta5","alpha","gamma"))


BB=jm$BUGSoutput$summary[101:125,c(1:3,7:8)]


mcmcplots::caterplot(mcmc_output, parms = c("gamma"), denstrip=TRUE, 
                     col="blue", pch=NA,
                     labels.loc=expression(paste(gamma[1],gamma[2],gamma[3],gamma[4],gamma[5])))
#mcmcplots::caterplot(mcmc_output, parms = c("gamma"), col="white", pch=NA, add=TRUE)


