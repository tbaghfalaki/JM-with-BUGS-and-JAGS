### This code was written by Antoine Barbieri, and I have included it here with his permission. ###
# load packages
library(survival)
library(ggfortify)
library(ggplot2)
library(gridExtra)
library(grid)

# load data
data("aids", package = "JMbayes2")
data("aids.id", package = "JMbayes2")
data.id <- aids.id
data.long <- aids

# survival fit
km_fit <- survfit(Surv(data.id$Time, data.id$death) ~ 1)
p2 = autoplot(km_fit, ylab="Survival Probability", xlab="Time", ylim = c(0,1)) + ggtitle("(b)")

# longitudinal plot
p1 = ggplot(data = data.long, aes(x = obstime, y = CD4, group = patient)) +
  geom_line() +
  xlab("Time") +
  geom_abline(slope = -0.188, intercept = 5.656, lty = 1, lwd = 1.2, color="blue") +
  geom_abline(slope = -0.132, intercept = 10.413, lty = 2, lwd = 1.2, color="blue") +
  ylab("Longitudinal measurements") +
  ggtitle("(a)")

# both plots
grid.arrange(p1, p2, ncol=2, nrow = 1)



# Definition of models from formulas
#----------------------------------


# load data
data("aids", package = "JMbayes2")
data("aids.id", package = "JMbayes2")
data.id <- aids.id
data.long <- aids


# creation of jags data list starting with longitudinal part
jags.data <- list(
  y = data.long$CD4,
  X = cbind(1, data.long$obstime),
  Z = cbind(1, data.long$obstime),
  ncX = 2,
  ncZ = 2,
  id_row = as.vector(c(1, 
                       1+cumsum(tapply(as.integer(data.long$patient),
                                       as.integer(data.long$patient),
                                       length))))
)


# Working data for survival integral approximation using Gauss-Kronrod
sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
        0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
        -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
        0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329)
wk15 <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
          0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
          0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
          0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)
# update jags data with zeros trick
zeros <- numeric(length(data.id$Time))
# update jags data with the survival data
jags.data <- c(jags.data, 
               list(n = length(data.id$Time),
                    st = data.id$Time,
                    delta = data.id$death,
                    W2 = cbind(1, 
                               data.id$drug, 
                               data.id$gender, 
                               data.id$prevOI, 
                               data.id$AZT),
                    ncW2 = 5,
                    K = length(wk15), 
                    wk = wk15, 
                    xk = sk,
                    zeros = zeros))


# Define the partially latent variable U for class membership
u <- data.id$death
u[which(u==0)] <- NA
# zero-tail constraint
u[which(data.id$death==0 & 
          data.id$Time>max(data.id$Time[which(data.id$death==1)]))] <- 0
# update jags data with the incidence data
jags.data <- c(jags.data, 
               list(W1 = cbind(1, 
                               data.id$drug, 
                               data.id$gender, 
                               data.id$prevOI, 
                               data.id$AZT),
                    ncW1 = 5,
                    u = u))


################################################################################
# Priors 
################################################################################


# Choice of precision
precision <- 0.01
# prior parameters for longitudinal submodel
priorMean.beta <- as.numeric(rep(0, jags.data$ncX))
priorTau.beta <- diag(rep(precision, length(priorMean.beta)))
mu0 <- rep(0, jags.data$ncZ)
priorR.prec.mat <- diag(rep(precision, jags.data$ncZ))
priorK.prec.mat <- jags.data$ncZ
# prior parameters for incidence submodel
priorMean.xi <- as.numeric(rep(0, jags.data$ncW1))
priorTau.xi <- diag(rep(precision, length(priorMean.xi)))
# prior parameters for latency submodel
priorMean.alpha <- as.numeric(rep(0, jags.data$ncW2))
priorTau.alpha <- diag(rep(precision, length(priorMean.alpha)))

# update jags data with hyper parameters
jags.data <- c(jags.data, 
               list(precision = precision, 
                    priorMean.beta = priorMean.beta,
                    priorTau.beta = priorTau.beta,
                    priorMean.xi = priorMean.xi,
                    priorTau.xi = priorTau.xi,
                    priorMean.alpha = priorMean.alpha,
                    priorTau.alpha = priorTau.alpha,
                    mu0 = mu0,
                    priorR.prec.mat = priorR.prec.mat,
                    priorK.prec.mat = priorK.prec.mat))


jags.model <- "
model{
    # Likelihood
    for (i in 1:n) {
      # Individual contribution to the likelihood
      # Longitudinal submodel
      # Distribution of Y
      for(j in id_row[i]:(id_row[i+1]-1)){
        y[j] ~ dnorm(mu[j], prec.tau2[index[i]])
        mu[j] <- inprod(beta[index[i], 1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncZ], Z[j, 1:ncZ])
      }
      # Distribution of random effects
      b[i, 1:ncZ] ~ dmnorm(mu0, prec.mat[index[i], , ])
      # Incidence submodel
      logit(pi[i]) <- inprod(xi[1:ncW1], W1[i, 1:ncW1])
      u[i] ~ dbern(pi[i])
      # define the index for class specific parameters
      index[i] <- 1*equals(u[i],1) + 2*equals(u[i],0)
      # Latency submodel
      A0[i] <- gamma * ( beta[index[i], 1] + b[i, 1] )
      A1[i] <- gamma * ( beta[index[i], 2] + b[i, 2] )
      # Hazard function
      eta_baseline[i] <- inprod(alpha[1:ncW2], W2[i, 1:ncW2])
      haz[i]<- nu * pow( st[i], nu-1) * exp(eta_baseline[i] + A0[i] + A1[i] * st[i] )
      # Cumulative hazard function
      for (k in 1:K) {
        # Scaling Gauss-Kronrod / Legendre quadrature
        xk1[i, k] <- (xk[k]+1)*st[i]/2
        wk1[i, k] <- wk[k]*st[i]/2
        # Hazard function at Gauss-Kronrod / Legendre nodes
        chaz[i, k] <- nu*pow(xk1[i, k], nu-1) * exp(eta_baseline[i] + A0[i] + A1[i]*xk1[i, k] )
      }
      # Log-survival function with Gauss-Kronrod / Legendrer requadrature
      logSurv[i] <- -inprod( wk1[i, 1:K], chaz[i, 1:K] )
      # Definition of the survival log-likelihood using zeros trick
      logL[i] <- u[i]*log(pi[i]) + u[i]*delta[i]*log(haz[i]) + u[i]*logSurv[i] + (1-delta[i])*(1-u[i])*log(1-pi[i])
      phi[i] <- 100000 - logL[i]
      zeros[i] ~ dpois(phi[i])
    }
    # contribution to prior
    # Longitudinal part
    prec.tau2[1] ~ dgamma(0.01, 0.01)
    prec.tau2[2] ~ dgamma(0.01, 0.01)
    beta[1,1:ncX] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])         
    beta[2,1:ncX] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])          
    prec.mat[1, 1:ncZ, 1:ncZ] ~ dwish(priorR.prec.mat[, ], priorK.prec.mat)  
    prec.mat[2, 1:ncZ, 1:ncZ] ~ dwish(priorR.prec.mat[, ], priorK.prec.mat)  
    # Incidence priors
    xi[1:ncW1] ~ dmnorm(priorMean.xi[], priorTau.xi[, ])
    # Latency/survival priors
    alpha[1:ncW2] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
    nu ~ dgamma(0.01, 0.01)
    gamma ~ dnorm(0, precision)
    # determinetic nodes
    sigma2[1] <- 1/prec.tau2[1]
    sigma2[2] <- 1/prec.tau2[2]
    covariance.b[1, 1:ncZ, 1:ncZ] <- inverse(prec.mat[1, 1:ncZ, 1:ncZ])
    covariance.b[2, 1:ncZ, 1:ncZ] <- inverse(prec.mat[2, 1:ncZ, 1:ncZ])
  }
"



# parameters to save in the sampling step
parms_to_save <- c("beta", "sigma2", "xi", "alpha", "gamma", "nu", "covariance.b", "b")
# using jagsUI
out_jags = jagsUI::jags(data = jags.data,
                        parameters.to.save = parms_to_save,
                        model.file = textConnection(jags.model),
                        n.chains = 3,
                        parallel = FALSE,
                        n.adapt = NULL,
                        n.iter = 200000,
                        n.burnin = 50000,
                        n.thin = 1,
                        DIC = F)

jagsUI::traceplot(out_jags, 
                  parameters = c("beta", "xi", "alpha", "nu", "gamma"), 
                  layout = c(4,4))

jagsUI::traceplot(out_jags, parameters = c("covariance.b", "sigma2"), 
                  layout = c(5,2))