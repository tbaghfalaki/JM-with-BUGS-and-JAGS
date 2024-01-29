rm(list = ls())
library(joineRML)
library(R2jags)
data(pbc2)
attach(pbc2)
names(pbc2)
head(pbc2)
CR=c()
CR[pbc2$status=="alive"]=0
CR[pbc2$status=="transplanted"]=1
CR[pbc2$status=="dead"]=2

############# data_preparation
n <- length(unique(pbc2$id))
Surv1 <- Stat <- Age <- Drug <- id <- matrix(NA, n, 20)
ni <- c()
for (i in 1:n) {
  id[i, 1] <- i
  ni[i] <- length(pbc2$id[pbc2$id == i])
  Surv1[i, 1:ni[i]] <- pbc2$years[pbc2$id == i]
  Stat[i, 1:ni[i]] <- CR[pbc2$id == i]
  Age[i, 1:ni[i]] <- pbc2$age[pbc2$id == i]
  Drug[i, 1:ni[i]] <- pbc2$drug[pbc2$id == i]
}
nmark=5
years <- Surv1[, 1]
status2 <- Stat[, 1]
age <- Age[, 1]
drug <- Drug[, 1]
id <- id[, 1]
long.data <- data.frame(id, years, status2, age, drug)
Time <- long.data$years
CR <- long.data$status2 
table(CR)/n
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
XS <- matrix(1, n, 3) # Fixed effects
XS[, 2] <- treat
XS[, 3] <- Age

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
      mu1[i,j]<-inprod(betaL1[],XL1[i,j,])+inprod(b[i,1:2],ZL1[i,j,])
      Y2[i,j]~dnorm(mu2[i,j],tau[2])
      mu2[i,j]<-inprod(betaL2[],XL1[i,j,])+inprod(b[i,3:4],ZL1[i,j,])
      Y3[i,j]~dnorm(mu3[i,j],tau[3])
      mu3[i,j]<-inprod(betaL3[],XL1[i,j,])+inprod(b[i,5:6],ZL1[i,j,])
      Y4[i,j]~dnorm(mu4[i,j],tau[4])
      mu4[i,j]<-inprod(betaL4[],XL1[i,j,])+inprod(b[i,7:8],ZL1[i,j,])
      Y5[i,j]~dnorm(mu5[i,j],tau[5])
      mu5[i,j]<-inprod(betaL5[],XL1[i,j,])+inprod(b[i,9:10],ZL1[i,j,])
    }
    
    #Survival and censoring times
    # 1th cause
    Alpha01[i]<- inprod(alpha1[],XS[i,])+gamma1[1]*(betaL1[1]+b[i,1])+
      gamma1[2]*(betaL2[1]+b[i,3])+
      gamma1[3]*(betaL3[1]+b[i,5])+
      gamma1[4]*(betaL4[1]+b[i,7])+
      gamma1[5]*(betaL5[1]+b[i,9])
    
    Alpha11[i]<- gamma1[1]*(betaL1[2]+b[i,2])+
      gamma1[2]*(betaL2[2]+b[i,4])+
      gamma1[3]*(betaL3[2]+b[i,6])+
      gamma1[4]*(betaL4[2]+b[i,8])+
      gamma1[5]*(betaL5[2]+b[i,10])
    
    
    haz1[i]<- nu1*pow(Time[i],nu1-1)*exp(Alpha01[i]+Alpha11[i]*Time[i])
    #Cumulative hazard function 
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature 
      xk1[i,j]<-(xk[j]+1)/2*Time[i] 
      wk1[i,j]<- wk[j]*Time[i]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz1[i,j]<- nu1*pow(xk1[i,j],nu1-1)*exp(Alpha01[i]+Alpha11[i]*xk1[i,j])
    }
    
    logSurv1[i]<- -inprod(wk1[i,],chaz1[i,])
    
    
    # 2th cause
    Alpha02[i]<- inprod(alpha2[],XS[i,])+gamma2[1]*(betaL1[1]+b[i,1])+
      gamma2[2]*(betaL2[1]+b[i,3])+
      gamma2[3]*(betaL3[1]+b[i,5])+
      gamma2[4]*(betaL4[1]+b[i,7])+
      gamma2[5]*(betaL5[1]+b[i,9])
    
    Alpha12[i]<- gamma2[1]*(betaL1[2]+b[i,2])+
      gamma2[2]*(betaL2[2]+b[i,4])+
      gamma2[3]*(betaL3[2]+b[i,6])+
      gamma2[4]*(betaL4[2]+b[i,8])+
      gamma2[5]*(betaL5[2]+b[i,10])
    
    
    haz2[i]<- nu2*pow(Time[i],nu2-1)*exp(Alpha02[i]+Alpha12[i]*Time[i])
    #Cumulative hazard function 
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature 
      xk2[i,j]<-(xk[j]+1)/2*Time[i] 
      wk2[i,j]<- wk[j]*Time[i]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz2[i,j]<- nu2*pow(xk2[i,j],nu2-1)*exp(Alpha02[i]+Alpha12[i]*xk2[i,j])
    }
    
    logSurv2[i]<- -inprod(wk2[i,],chaz2[i,])
    
    
    #Definition of the survival log-likelihood using zeros trick
    phi[i]<-100000-equals(CR[i],1)*log(haz1[i])-equals(CR[i],2)*log(haz2[i])-logSurv1[i]-logSurv2[i]
    zeros[i]~dpois(phi[i])
    #Random effects
    b[i,1:Nb]~dmnorm(mub[],Omega[,])
  }
  #Prior distributions
  for(l in 1:NbetaL){
    betaL1[l]~dnorm(0,0.001)
    betaL2[l]~dnorm(0,0.001)
    betaL3[l]~dnorm(0,0.001)
    betaL4[l]~dnorm(0,0.001)
    betaL5[l]~dnorm(0,0.001)
  }
  
  for(l in 1:nmark){ 
    gamma1[l]~dnorm(0,0.001)
    gamma2[l]~dnorm(0,0.001)
    tau[l]~dgamma(0.01,0.01)
    sigma[l]<-1/tau[l]
  }  
  
  for(l in 1:Nalpha){
    alpha1[l]~dnorm(0,0.001)
    alpha2[l]~dnorm(0,0.001)
    
  }
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])
  nu1~dgamma(0.01,0.01)
  nu2~dgamma(0.01,0.01)
  
}", fill = TRUE)
sink()

d.jags <- list(
  n = n, M = M, Time = Time, Y1 = Albumin, Y2 = Alkaline,
  Y3 = SGOT1, Y4 = Platelets, Y5 = SerBilir, nmark = nmark,
  XL1 = XL, ZL1 = ZL, XS = XS,
  CR = CR, mub = rep(0, 10), V = diag(1, 10), Nb = 10, zeros = rep(0, n),
  NbetaL = dim(XL)[3], Nalpha = dim(XS)[2], xk = xk, wk = wk, K = K
)

i.jags <- function() {
  list(
    gamma1 = rnorm(nmark),    gamma2 = rnorm(nmark),

    betaL1 = rnorm(dim(XL)[3]),
    betaL2 = rnorm(dim(XL)[3]),
    betaL3 = rnorm(dim(XL)[3]),
    betaL4 = rnorm(dim(XL)[3]),
    betaL5 = rnorm(dim(XL)[3]),
    alpha1 = rnorm(dim(XS)[2]),
    alpha2 = rnorm(dim(XS)[2]), tau = rep(1, 5), Omega = diag(runif(10))
  )
}

parameters <- c(
  "betaL1", "betaL2", "betaL3", "betaL4", "betaL5", "alpha1", "alpha2",
  "gamma1", "gamma2", "nu1","nu2", "sigma", "Sigma"
)
jm <- jags(
  data = d.jags, inits = i.jags, parameters, n.chains = 2,
  n.iter = 70000, model.file = "model_file"
)
print(jm)

save(jm, file="/Users/taban/Desktop/Taban/joint modeling bugs/PBC-CR/CRpbc.RData")
#jm1 <- update(jm, n.iter = 20000, n.burnin = 10000)


BB=jm$BUGSoutput$summary[101:134,c(1:3,7:8)]
xtable::xtable(BB, type = "latex", digits = 3, file = "/Users/taban/Desktop/fff/filename2.tex")


mcmc_output=as.mcmc(jm)
par(mfrow = c(1,2)) 
mcmcplots::caterplot(mcmc_output, parms = c("gamma1"),main="transplanted")
mcmcplots::caterplot(mcmc_output, parms = c("gamma2"),main="dead")



par(mfrow = c(1,2)) 
mcmcplots::caterplot(mcmc_output, parms = c("gamma1"), denstrip=TRUE, 
                     col="blue", pch=NA,main="222")
mcmcplots::caterplot(mcmc_output, parms = c("gamma1"), col="white", pch=NA, add=TRUE)

mcmcplots::caterplot(mcmc_output, parms = c("gamma2"), denstrip=TRUE, col="blue", pch=NA)
mcmcplots::caterplot(mcmc_output, parms = c("gamma2"), col="white", pch=NA, add=TRUE)

