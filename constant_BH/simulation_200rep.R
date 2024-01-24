rm(list = ls())
library(MASS)
library(survival)
library(R2jags)

# Random number generation
# number of indivduals

resultsss <- list()
NNNN <- 200
for (ij in 1:NNNN) {
  set.seed(ij)
  nsample <- 500
  # real values of parameters
  mu <- rep(0, 2)
  rho <- 0.5
  D <- matrix(c(1, rho, rho, 1), 2, 2)
  gamma <- -0.5
  Beta <- c(-0.5, 0.5, 0.5, 0.5)
  sigma <- 1
  alpha <- c(1, -1)
  Realpara <- c(
    Beta, alpha, gamma,
    (sigma)^2, D[1, ], D[2, 2]
  )

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
  evtim <- td <- rep(0, nsample)
  u <- matrix(0, nsample, 2)
  A0 <- A1 <- c()
  for (i in 1:nsample) {
    x1[i] <- rbinom(1, 1, 0.6)
    x2[i] <- rnorm(1)
    w1[i] <- rnorm(1)
    w2[i] <- rbinom(1, 1, 0.5)

    u[i, ] <- mvrnorm(1, mu, D)
    Alpha1 <- gamma * (Beta[2] + u[i, 2])
    Alpha0 <- alpha[1] * w1[i] + alpha[2] * w2[i] + gamma * (Beta[1] + Beta[3] * x1[i] + Beta[4] * x2[i] + u[i, 1])

    for (j in 1:length(t)) {
      Y1[i, j] <- Beta[1] + Beta[2] * t[j] + Beta[3] * x1[i] + Beta[4] * x2[i] + u[i, 1] + u[i, 2] * t[j] + rnorm(1, 0, sigma)
    }

    lambda0 <- 1
    temp <- (Alpha1 * (-log(runif(1, min = 0, max = 1))) / lambda0 + exp(Alpha0))
    if (temp > 0) {
      time[i] <- (log(temp) - Alpha0) / Alpha1
      if (time[i] >= 2) {
        time[i] <- 2
        death[i] <- 0
      } else {
        death[i] <- 1
      }
    } else {
      time[i] <- 2
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
  sum(death) / nsample
  # Number of patients and number of longitudinal observations per patient
  n <- length(surv.data$id)
  M <- table(long.data$id)
  mean(M)
  median(M)

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

  W <- model.matrix(~ surv.data$w1 + surv.data$w2) # Fixed effects

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
    # Longitudinal sub-model
    for(j in 1:M[i]){
      Y[i,j]~dnorm(mu[i,j],tau)
      mu[i,j]<-inprod(beta[],X[i,j,])+inprod(b[i,],Z[i,j,])
    }
    Alpha0[i]<- gamma*(beta[1]+beta[3]*X1[i]+beta[4]*X2[i]+b[i,1])
    Alpha1[i]<- gamma*(beta[2]+b[i,2])
    # Survival sub-model
    # Hazard function
    haz[i]<- exp(inprod(alpha[],W[i,])+
                   (Alpha0[i]+Alpha1[i]*st[i]))
    # Cumulative hazard function
    chaz[i]<- exp(inprod(alpha[],W[i,])+Alpha0[i])*(exp(st[i]*Alpha1[i])-1)/Alpha1[i]
    # Log-survival function
    logSurv[i]<- -chaz[i]
    # Definition of the survival log-likelihood using zeros trick
    phi[i]<-100000-death[i]*log(haz[i])-logSurv[i]
    zeros[i]~dpois(phi[i])
    # Random effects
    b[i,1:Nb]~dmnorm(mub[],Omega[,])
  }
  # Prior distributions
  for(l in 1:Nbeta){
    beta[l]~dnorm(0,0.001)
  }
  for(l in 1:Nalpha){
    alpha[l]~dnorm(0,0.001)
  }
  gamma~dnorm(0,0.001)
  kappa~dgamma(0.01,0.01)
  tau~dgamma(0.01,0.01)
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  # Determinetic nodes
  lambda0<-exp(alpha[1])
  sigma<-1/tau
  Sigma[1:Nb,1:Nb]<-inverse(Omega[,])
}", fill = TRUE)
  sink()
  #### Running JAGS
  d.jags <- list(
    n = n, M = M, st = st, W = W, Y = Y, X = X, Z = Z, X1 = X1, X2 = X2,
    death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
    Nbeta = dim(X)[3], Nalpha = ncol(W)
  )

  i.jags <- function() {
    list(
      alpha = rnorm(ncol(W)), gamma = rnorm(1),
      beta = rnorm(dim(X)[3]), tau = runif(1), Omega = diag(runif(2))
    )
  }
  parameters <- c("alpha", "gamma", "beta", "sigma", "Sigma", "lambda0")

  sim1 <- jags(
    data = d.jags, inits = i.jags, parameters, n.chains = 1,
    n.iter = 3000, model.file = "model_file"
  )

  print(sim1)

  ################
  Realpara <- c(
    Beta, alpha, gamma,
    (sigma)^2, D[1, ], D[2, 2]
  )
  Est <- c(
    sim1$BUGSoutput$mean$beta, sim1$BUGSoutput$mean$alpha[-1], sim1$BUGSoutput$mean$gamma,
    sim1$BUGSoutput$mean$sigma,
    sim1$BUGSoutput$mean$Sigma[1, ], sim1$BUGSoutput$mean$Sigma[2, 2]
  )
  Est_l <- c(
    apply(sim1$BUGSoutput$sims.list$beta, 2, quantile, 0.025),
    apply(sim1$BUGSoutput$sims.list$alpha, 2, quantile, 0.025)[-1],
    apply(sim1$BUGSoutput$sims.list$gamma, 2, quantile, 0.025),
    apply(sim1$BUGSoutput$sims.list$sigma, 2, quantile, 0.025),
    apply(sim1$BUGSoutput$sims.list$Sigma[, , 1], 2, quantile, 0.025),
    apply(sim1$BUGSoutput$sims.list$Sigma[, , 2], 2, quantile, 0.025)[2]
  )

  Est_u <- c(
    apply(sim1$BUGSoutput$sims.list$beta, 2, quantile, 0.975),
    apply(sim1$BUGSoutput$sims.list$alpha, 2, quantile, 0.975)[-1],
    apply(sim1$BUGSoutput$sims.list$gamma, 2, quantile, 0.975),
    apply(sim1$BUGSoutput$sims.list$sigma, 2, quantile, 0.975),
    apply(sim1$BUGSoutput$sims.list$Sigma[, , 1], 2, quantile, 0.975),
    apply(sim1$BUGSoutput$sims.list$Sigma[, , 2], 2, quantile, 0.975)[2]
  )


  sdd <- c(
    sim1$BUGSoutput$sd$beta,
    sim1$BUGSoutput$sd$alpha[-1], sim1$BUGSoutput$sd$gamma, sim1$BUGSoutput$sd$sigma,
    sim1$BUGSoutput$sd$Sigma[1, ], sim1$BUGSoutput$sd$Sigma[2, 2]
  )
  CR <- rep(0, length(Realpara))

  for (kk in 1:length(Realpara)) {
    if (Realpara[kk] > Est_l[kk] & Realpara[kk] < Est_u[kk]) (CR[kk] <- 1)
  }

  cbind(Realpara, Est, Est_l, Est_u, CR)

  resultsss[[ij]] <- list(CR = CR, Est = Est, Est_l = Est_l, Est_u = Est_u, sdd = sdd)

  print(ij)
}

MSE <- cr <- est <- sd <- matrix(NA, NNNN, length(Realpara))
for (j in 1:NNNN) {
  MSE[j, ] <- (resultsss[[j]]$Est - Realpara)^2
  cr[j, ] <- resultsss[[j]]$CR
  est[j, ] <- resultsss[[j]]$Est
  sd[j, ] <- resultsss[[j]]$sdd
}



BB <- cbind(Realpara, apply(est, 2, mean), apply(est, 2, sd), (apply(est, 2, mean) - Realpara) / Realpara, sqrt(apply(MSE, 2, mean)), apply(sd, 2, mean), apply(cr, 2, mean))
colnames(BB) <- c("Real", "Est", "se", "Bias", "RMSE", "ese", "CR")

xtable::xtable(BB, type = "latex", digits = 3, file = "/Users/taban/Desktop/fff/filename2.tex")




Real        Est         se         Bias       RMSE        ese    CR
[1,] -0.5 -0.4946052 0.09527520 -0.010789632 0.09518971 0.08200852 0.900
[2,]  0.5  0.5058958 0.07597355  0.011791649 0.07601238 0.06825725 0.915
[3,]  0.5  0.4945095 0.12143364 -0.010981086 0.12125405 0.10503914 0.895
[4,]  0.5  0.4982298 0.05510085 -0.003540340 0.05499142 0.05282170 0.910
[5,]  1.0  1.0073511 0.06777871  0.007351087 0.06800752 0.06605791 0.920
[6,] -1.0 -1.0028293 0.11921530  0.002829285 0.11895054 0.11987429 0.940
[7,] -0.5 -0.5020777 0.04564233  0.004155448 0.04557547 0.04497996 0.940
[8,]  1.0  1.0034742 0.03034837  0.003474183 0.03047110 0.02902721 0.930
[9,]  1.0  1.0011179 0.10012805  0.001117928 0.09988367 0.09270015 0.925
[10,]  0.5  0.4987655 0.08542366 -0.002468959 0.08521878 0.08071672 0.945
[11,]  1.0  1.0147057 0.12087473  0.014705740 0.12146566 0.11709410 0.935
