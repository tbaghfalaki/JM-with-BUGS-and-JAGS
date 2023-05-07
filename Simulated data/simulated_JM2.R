rm(list=ls())
library(PermAlgo)
library(mvtnorm)


nsample=500 # number of indivduals
# real values of parameters
Beta=c(0.5, -0.5, 1,1)
sigma=1 
gamma=-0.5 
rho=.5 
Sigma=matrix(c(1,rho,rho,1),2,2)


# gap between longitudinal measurements 
gapLongi=0.2 
gap=0.1 
# observed times for longitudinal marker
followup=2 
t <- seq(from = 0, to = followup, by = gap)
timesLongi=t[which(round(t-round(t/gapLongi,0)*gapLongi,6)==0)] # visit times
time=rep(t, nsample) 
# max. number of individual measurements
nmesindiv=followup/gap+1 
# max. total number of longitudinal measurements
nmesy= nmesindiv*nsample 
#  individual id for longitudinal
idY<-rep(1:nsample, each=nmesindiv) 

###############

b <- rmvnorm(nsample, rep(0, 2), Sigma)
b1 <- rep(b[,1], each=nmesindiv) # random intercept Y1
b2 <- rep(b[,2], each=nmesindiv) # random slope Y1

x1=rnorm(nsample,1, 0.5)
x2=rbinom(nsample,1, 0.5)
X1=rep(x1 , each=nmesindiv)  
X2=rep(x2, each=nmesindiv)  
# linear predictor
eta <- (Beta[1]+b1) + (Beta[2]+b2)*time + Beta[3]*X1 + Beta[4]*X2

Y1 <- rnorm(nmesy, eta, sigma)
# Permutation algorithm to generate survival times 
Data_permu <- permalgorithm(nsample,nmesindiv,Xmat=eta,
                            eventRandom = round(rweibull(nsample,1,3)+1,0), 
                            censorRandom=runif(nsample,1,nmesindiv),  
                            XmatNames=c("eta"),  
                            betas=gamma)  

# extract last line for each id (= death/censoring time)
Data_permu2=Data_permu[c(which(diff(Data_permu[,"Id"])==1), dim(Data_permu)[1]), c("Id","Event","Stop")]
Data_permu2$survtime <- t[Data_permu2$Stop+1] # survtime
surv.data <- Data_permu2[, c("Id", "survtime", "Event")]
surv.data$x1<-x1
Data_permu$time <- t[Data_permu$Start+1] # measurements times of the biomarker
Data_permu$Uid <- paste(Data_permu$Id, Data_permu$time) # unique identifier to match covariates and observed biomarker values
long.data3 <- merge(Data_permu[,c("Uid", "Id", "time")], cbind("Uid"=paste(idY,time), X1, X2, Y1), by=c("Uid"))
long.data <- sapply(long.data3[long.data3$time%in%timesLongi,-1], as.numeric)
long.data <- as.data.frame(long.data[order(long.data[, "Id"], long.data[, "time"]),])
colnames(long.data)<-c("id","obstime","x1","x2","Y1")
colnames(surv.data)<-c("id","survtime","death","w1")

head(long.data)
head(surv.data)


save(long.data, file = "/Users/taban/Desktop/Taban/joint modeling bugs/new_codes/long.data_2.RData")
save(surv.data, file = "/Users/taban/Desktop/Taban/joint modeling bugs/new_codes/surv.data_2.RData")

