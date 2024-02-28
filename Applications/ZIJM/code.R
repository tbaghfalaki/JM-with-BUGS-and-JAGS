rm(list=ls())
# load packages
library(mvtnorm)
library(pec)
library(timeROC)
library(survival)

# load data
setwd("/Users/taban/Desktop/Taban/joint modeling bugs/Applications/Zero_inflation JM")
Long=read.table("Long.txt",header=TRUE)
Surv=read.table("surv.txt",header=TRUE)
attach(Long)
attach(Surv)

## preparing data
ID=unique(SubjectID)
m=table(SubjectID)
n=dim(Surv)[1]
y=obstime=LibrarySize1=HH=matrix(NA,n,max(m)) 
for(i in 1:n){
  y[i,1:length(Prevotella[SubjectID==ID[i]])]=Prevotella[SubjectID==ID[i]]
  obstime[i,1:length(GWColl[SubjectID==ID[i]])]=GWColl[SubjectID==ID[i]]
  LibrarySize1[i,1:length(LibrarySize[SubjectID==ID[i]])]=LibrarySize[SubjectID==ID[i]]
  HH[i,1:length(History_of_preterm_delivery[SubjectID==ID[i]])]=History_of_preterm_delivery[SubjectID==ID[i]]
}

z=matrix(0,n,max(m))
z[y==0]=1
R=diag(2)*10;U0=c(0,0)
non_hispanic_whites=as.numeric(non_hispanic_whites)
Preeclampsias=as.numeric(Preeclampsias)
high_incomes=as.numeric(high_incomes)
History=HH[,1]
zeros=matrix(0,n,max(m))

########  BUGS code  ########
sink("model_file_ZINB")
cat("model{
  K<-1000
  for (i in 1:n) {

    for (j in 1:m[i]) {
      zeros[i,j]~dpois(phi[i,j])
      phi[i,j]<-  -ll[i,j]+K				
      ll[i,j]<-z[i,j]*log(pi[i,j]) +
        (1-z[i,j])*(
          log(1-pi[i,j])+loggam(r+y[i,j])-loggam(r)-loggam(y[i,j]+1)+r*log(r/(r+lambda[i, j]))+y[i,j]*log(lambda[i, j]/(lambda[i, j]+r))-log(1-pow(r/(r+lambda[i, j]),r))   )
      
      
      log(lambda[i, j])<-beta1[1]+beta1[2]*obstime[i,j]+beta1[3]*Preeclampsias[i]+log(LibrarySize1[i,j])+U[i,1]
      logit(pi[i, j])<-beta2[1]+beta2[2]*obstime[i,j]+beta2[3]*Preeclampsias[i]+U[i,2]
      
    }
    surt[i] ~ dweib(p,mut[i])     
    is.censored[i]~dinterval(surt[i],surt.cen[i])
    
    log(mut[i])<-beta3[1]+beta3[2]*Preeclampsias[i]+beta3[3]*non_hispanic_whites[i]+
      beta3[4]*high_incomes[i]+beta3[5]*History[i]+r1*U[i, 1]+r2*U[i, 2]
    
    U[i,1:2] ~ dmnorm(U0[],tau[,])   
  }  
  r~ dgamma(.1,.1) 
  p~ dgamma(.1,.1) 
  sigma[1:2,1:2]<-inverse(tau[,])
  
  #priors
  tau[1:2,1:2] ~ dwish(R[,], 2)
  for(j in 1:5){ 
    beta3[j]~dnorm(0,.0001)}
  for(j in 1:3){
    beta1[j]~dnorm(0,.0001)
    beta2[j]~dnorm(0,.0001)
    
  }
  r1~dnorm(0, 0.0001)
  r2~dnorm(0, 0.0001)
}
", fill = TRUE)
sink()
# creation of jags data  
surt=GWDels
status=rep(1,n)
surt.cen=rep(0,n)
for(i in 1:n){
  if(status[i]==1)(surt.cen[i]=800) 
  if(status[i]==0)((surt.cen[i]=surt[i]) & (surt[i]=NA))
}
time00 = Sys.time()
zeros=matrix(0,n,max(m))
is.censored=1-status
data <- list (n=40,obstime=obstime, Preeclampsias=Preeclampsias,non_hispanic_whites=non_hispanic_whites,zeros=zeros,
              high_incomes=high_incomes,m=m,History=History,R=R,is.censored=is.censored,LibrarySize1=LibrarySize1,
              y=y,z=z,surt=surt,surt.cen=surt.cen,U0=U0)
U=matrix(0,40,2)
inits <- function(){
  list(beta1=rep(0,3),beta2=rep(0,3),beta3=rep(0,5),r1=1,r2=-1,p=2,r=1,
       U=U)
  list(beta1=rep(0,3),beta2=rep(0,3),beta3=rep(0.1,5),r1=1,r2=-1,p=2,r=1,
       U=U)
}
parameters <- c("beta1","beta2","beta3","r1","r2","tau","sigma","r","p")

sim <- R2jags::jags(data, inits, parameters,n.chains=2,n.thin = 30,
                    n.iter=200000,model.file="model_file_ZINB")
save(sim,file="simzinb.RData")

