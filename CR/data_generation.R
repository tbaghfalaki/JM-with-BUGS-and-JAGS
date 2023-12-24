rm(list = ls())
library(mvtnorm)
library(MASS)
library(survival)
library(Rsolnp)
library(Matrix)
set.seed(9)
nsample <- 500
follow_up <- 1
t <- seq(from = 0, to = follow_up, 0.1)
mu <- rep(0, 2)
rho <- 0.1
Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
gamma1 <- -0.5
gamma2 <- -0.5
alpha1 <- c(1, -1)
alpha2 <- c(-1, 1)
Beta <- rep(1, 4)
sigma <- sqrt(.5)

p <- 1
C <- rep(NA, nsample)
w1 <- w2 <- x1 <- x2 <- rep(NA, nsample)
CR <- rep(0, nsample) # 1,2=observed, 0=censored
Y1 <- matrix(NA, ncol = length(t), nrow = nsample)
time <- rep(NA, nsample)
id <- rep(NA, nsample)
ct <- rep(NA, nsample)
mat <- matrix(NA, ncol = 7, nrow = nsample * length(t))

evtim1 <- evtim2 <- td <- c()
u <- matrix(0, nsample, 2)
A10 <- A20 <- A11 <- A21 <- c()
for (i in 1:nsample) {
  u[i, ] <- mvrnorm(1, mu, Sigma)
  x1[i] <- rbinom(1, 1, 0.6)
  x2[i] <- rnorm(1, 0, 1)
  w1[i] <- rbinom(1, 1, 0.5)
  w2[i] <- rnorm(1, 0, 1)
  Alpha10 <- alpha1[1] * w1[i] + alpha1[2] * w2[i] + gamma1 * (Beta[1] + Beta[3] * x1[i] + Beta[4] * x2[i] + u[i, 1])
  Alpha11 <- gamma1 * (Beta[2] + u[i, 2])
  
  
  Alpha20 <- alpha2[1] * w1[i] + alpha2[2] * w2[i] + gamma2 * (Beta[1] + Beta[3] * x1[i] + Beta[4] * x2[i] + u[i, 1])
  Alpha21 <- gamma2 * (Beta[2] + u[i, 2])
  
  
  scale <- 1
  shape <- 0.5
  td[i] <- rweibull(1, shape, scale)
  for (j in 1:length(t)) {
    Y1[i, j] <- Beta[1] + Beta[2] * t[j] + Beta[3] * x1[i] + Beta[4] * x2[i] + u[i, 1] + u[i, 2] * t[j] + rnorm(1, 0, sigma)
  }
  
  
  
  
  
  lambda0 <- 1
  temp1 <- (Alpha11 * (-log(runif(1, min = 0, max = 1))) / lambda0 + exp(Alpha10))
  if(temp1> 0){
    evtim1[i] <- (log(temp1) - Alpha10) / Alpha11
  }else{
    evtim1[i] <- 10000
    CR[i] <- 0
    td[i]=follow_up
    
  }
  
  
  temp2 <- (Alpha21 * (-log(runif(1, min = 0, max = 1))) / lambda0 + exp(Alpha20))
  if(temp2> 0){
    evtim2[i] <- (log(temp2) - Alpha20) / Alpha21
  }else{
    evtim2[i] <- 10000
    CR[i] <- 0
    td[i]=follow_up
  }
  
  time[i] <- min(evtim1[i],evtim2[i], td[i],2)
  if (time[i] == follow_up | time[i] ==  td[i]) {
    CR[i] <- 0
  } 
  
  if (time[i] ==evtim1[i]) {
    CR[i] <- 1
  }
  
  if (time[i] == evtim2[i]) {
    CR[i] <- 2
    
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
    mat[jj, 5] <- CR[i]
    mat[jj, 6] <- time[i]
    mat[jj, 7] <- Y1[i, jk]
    jk <- jk + 1
  }
  p <- p + ct[i]
}
colnames(mat) <- c("obstime", "id", "x1", "x2", "CR", "survtime", "Y1")
summary(time)
table(CR) / nsample

####### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
long.data <- data.frame(mat)
long.data <- na.omit(long.data)
head(long.data)
####### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
CR.data <- data.frame(1:nsample, w1, w2, x1, x2, time, CR)
colnames(CR.data) <- c("id", "w1", "w2", "x1", "x2", "survtime", "CR")
head(CR.data)
#####################
save(long.data, file = "/Users/taban/Desktop/Taban/joint modeling bugs/CR/long.data.RData")
save(CR.data, file = "/Users/taban/Desktop/Taban/joint modeling bugs/CR/CR.data.RData")

table(CR) / nsample



library(ggplot2)
p1 <- ggplot(data = long.data, aes(x = obstime, y = Y1, group = id)) +
  geom_line() +
  xlab("Time") +
  stat_smooth(aes(group = 1), method = "lm") +
  ylab("Longitudinal measurements") +
  ggtitle("(a)") +
  theme(plot.title = element_text(size = 16, face = "bold"))



require(cmprsk)
library(survminer)
library(cowplot)
library(ggsci)
print(fit <- cmprsk::cuminc(ftime = CR.data$survtime, fstatus = CR.data$CR, group = ""))

p2 <- ggcompetingrisks(fit, conf.int = TRUE, multiple_panels = FALSE, palette = "npg", subset = "", legend.title = "") + xlab("Time") +
  ylab("CIF") + ggtitle("(b)") +
  theme_cowplot() + scale_fill_jco()


library(gridExtra)
library(grid)
grid.arrange(p1, p2,
             ncol = 2, nrow = 1
)

