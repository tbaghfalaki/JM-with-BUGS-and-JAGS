  rm(list = ls())
  library(PermAlgo)
  library(mvtnorm)
  library(survival)
  
  nsample <- 500 # number of indivduals
  # real values of parameters
  gamma <- -0.5
  Beta <- c(-0.5, 0.5, 0.5, 0.5)
  sigma <- 1
  alpha <- c(1,-1)
  rho <- 0.5
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  
  set.seed(12345)
  # gap between longitudinal measurements
  gapLongi <- 0.2
  gap <- 0.1
  # observed times for longitudinal marker
  followup <- 2
  t <- seq(from = 0, to = followup, by = gap)
  timesLongi <- t[which(round(t - round(t / gapLongi, 0) * gapLongi, 6) == 0)] # visit times
  time <- rep(t, nsample)
  # max. number of individual measurements
  nmesindiv <- followup / gap + 1
  # max. total number of longitudinal measurements
  nmesy <- nmesindiv * nsample
  #  individual id for longitudinal
  idY <- rep(1:nsample, each = nmesindiv)
  
  ###############
  
  b <- rmvnorm(nsample, rep(0, 2), Sigma)
  b1 <- rep(b[, 1], each = nmesindiv) # random intercept Y1
  b2 <- rep(b[, 2], each = nmesindiv) # random slope Y1
  
  x1 <- rnorm(nsample, 0, 1)
  x2 <- rbinom(nsample, 1, 0.6)
  X1 <- rep(x1, each = nmesindiv)
  X2 <- rep(x2, each = nmesindiv)
  w1 <- rnorm(nsample, 0, 1)
  W1<- rep(w1, each = nmesindiv)
  w2 <- rbinom(nsample, 1, 0.5)
  W2<- rep(w2, each = nmesindiv)
  
  # linear predictor
  eta <- (Beta[1] + b1) + (Beta[2] + b2) * time + Beta[3] * X1 + Beta[4] * X2
  
  Y1 <- rnorm(nmesy, eta, sigma)
  # Permutation algorithm to generate survival times
  Data_permu <- permalgorithm(nsample, nmesindiv,
    Xmat = cbind(eta,W1,W2),
    eventRandom = round(rexp(nsample, 0.02) + 1, 0),
    censorRandom = runif(nsample, 1, nmesindiv),
    XmatNames = c("eta","W1","W2"),
    betas = c(gamma,alpha)
  )

# extract last line for each id (= death/censoring time)
Data_permu2 <- Data_permu[c(which(diff(Data_permu[, "Id"]) == 1), dim(Data_permu)[1]), c("Id", "Event", "Stop")]
Data_permu2$survtime <- t[Data_permu2$Stop + 1] # survtime
surv.data <- Data_permu2[, c("Id", "survtime", "Event")]
surv.data$w1 <- w1
surv.data$w2 <- w2
Data_permu$time <- t[Data_permu$Start + 1] # measurements times of the biomarker
Data_permu$Uid <- paste(Data_permu$Id, Data_permu$time) # unique identifier to match covariates and observed biomarker values
long.data3 <- merge(Data_permu[, c("Uid", "Id", "time")], cbind("Uid" = paste(idY, time), X1, X2, Y1), by = c("Uid"))
long.data <- sapply(long.data3[long.data3$time %in% timesLongi, -1], as.numeric)
long.data <- as.data.frame(long.data[order(long.data[, "Id"], long.data[, "time"]), ])
colnames(long.data) <- c("id", "obstime", "x1", "x2", "Y1")
colnames(surv.data) <- c("id", "survtime", "death", "w1", "w2")

head(long.data)
head(surv.data)


save(long.data, file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated_data/long.data_2.RData")
save(surv.data, file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated_data/surv.data_2.RData")



library(ggplot2)
p1=ggplot(data = long.data, aes(x = obstime, y = Y1, group = id)) +  
  geom_line()+xlab("Time")+stat_smooth(aes(group = 1), method = "lm")+
  ylab("Longitudinal measurements")+ggtitle("(a)")


library(ggfortify)
km_fit <- survfit(Surv(survtime, death) ~ 1, data=surv.data)
p2=autoplot(km_fit,ylab="Survival Probability",xlab="Time")+ggtitle("(b)")


library(gridExtra)
library(grid)
grid.arrange(p1, p2,
             ncol=2, nrow = 1)




