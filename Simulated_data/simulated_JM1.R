rm(list = ls())
library(MASS)
library(survival)

# Random number generation
set.seed(12345)
# number of indivduals
nsample <- 500
# real values of parameters
mu <- rep(0, 2)
rho <- 0.5
D <- matrix(c(1, rho, rho, 1), 2, 2)
gamma <- -0.5
Beta <- c(-0.5, 0.5, 0.5, 0.5)
sigma <- 1
alpha <- c(1,-1)
# observed times for longitudinal marker
t <- seq(from = 0, to = 2.0, by = 0.2)

p <- 1
C <- rep(NA, nsample)
x1 <- x2 <- w1 <- w2 <- rep(NA, nsample)
death <- rep(0, nsample) # 1=observed, 0=censored
Y1 <- matrix(NA, ncol = length(t), nrow = nsample)
time <- rep(NA, nsample)
id <- rep(NA, nsample)
ct <- rep(NA, nsample)
mat <- matrix(NA, ncol = 7, nrow = nsample * length(t))
evtim <- td <- rep(0,nsample)
u <- matrix(0, nsample, 2)
A0 <- A1 <- c()
for (i in 1:nsample) {
  x1[i] <- rbinom(1, 1, 0.6)
  x2[i] <- rnorm(1)
  w1[i] <- rnorm(1)
  w2[i] <- rbinom(1, 1, 0.5)
  
  u[i, ] <- mvrnorm(1, mu, D)
  Alpha1 <- gamma * (Beta[2] + u[i, 2])
  Alpha0 <-  alpha[1] * w1[i] +alpha[2] * w2[i] + gamma * (Beta[1] + Beta[3] * x1[i] + Beta[4] * x2[i] + u[i, 1])
  
  
  td[i] <- rexp(1, 2)
  
  for (j in 1:length(t)) {
    Y1[i, j] <- Beta[1] + Beta[2] * t[j] + Beta[3] * x1[i] + Beta[4] * x2[i] + u[i, 1] + u[i, 2] * t[j] + rnorm(1, 0, sigma)
  }
  
  lambda0 <- 1
  temp <- (Alpha1 * (-log(runif(1, min = 0, max = 1))) / lambda0 + exp(Alpha0))
  if(temp>0){
    evtim[i] <- (log(temp) - Alpha0) / Alpha1
    time[i] <- min(evtim[i], td[i])
    if (time[i] >= 2) {
      time[i]=2
    } else {
      death[i] <- 1
    }
  }else{
    time[i]=2 
    death[i] <- 0
  }
  
  for (j in 1:length(t)) {
    if (time[i] < t[j]) {
      Y1[i, j] <- NA
    }
  }
  ct[i] <- length(Y1[i, ][is.na(Y1[i, ]) == FALSE])
  id[i] <- i
  jk <- 1
  for (jj in p:(ct[i] + p - 1)) {
    mat[jj, 1] <- t[jk]
    mat[jj, 2] <- i
    mat[jj, 3] <- x1[i]
    mat[jj, 4] <- x2[i]
    mat[jj, 5] <- death[i]
    mat[jj, 6] <- time[i]
    mat[jj, 7] <- Y1[i, jk]
    
    jk <- jk + 1
  }
  p <- p + ct[i]
}
# 1 - sum(death) / nsample # censoring rate
colnames(mat) <- c("obstime", "id", "x1", "x2", "death", "survtime", "Y1")
####### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
surv.data <- data.frame(id)
surv.data$id <- id
surv.data$w1 <- w1
surv.data$w2 <- w2
surv.data$survtime <- time
surv.data$death <- death
long.data <- data.frame(mat)
long.data <- na.omit(long.data)
head(long.data)
head(surv.data)

save(long.data, file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated_data/long.data_1.RData")
save(surv.data, file = "/Users/taban/Desktop/Taban/joint modeling bugs/Simulated_data/surv.data_1.RData")


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








