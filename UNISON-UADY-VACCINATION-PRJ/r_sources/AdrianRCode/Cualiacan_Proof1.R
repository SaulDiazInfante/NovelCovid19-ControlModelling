library(MHadaptive)
library(MASS)
library(mvtnorm)
library(deSolve)
library(invgamma)
library(coda)

Data <- read.csv("InformationData.csv", header = TRUE)[1:23, c(7, 3)]
set.seed(69)
#Ode system
OdeSystem <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    delta_e <- 1/5.1
    mu_s <- 1/8
    alpha_a <- 1/8
    ##
    N_star <- S + E + Is + Ia + R
    ##
    dS <- -(beta_a * Ia + beta_s * Is) * (S / N_star)
    dE <- (beta_a * Ia + beta_s * Is) * (S / N_star) - delta_e * E
    dIs <- rho * delta_e * E - (alpha_s + mu_s) * Is
    dIa <- (1 - rho) * delta_e * E - alpha_a * Ia
    dR <- alpha_a * Ia + alpha_s * Is
    dD <- mu_s * Is
    dCIs <- rho * delta_e * E
    list(c(dS, dE, dIs, dIa, dR, dD, dCIs)) })
}
#### priors ####
Prior <- function(pars){
  beta_a <- pars[1]
  beta_s <- pars[2]
  rho <- pars[3]
  alpha_s <- pars[4]
  tau <- pars[5]
  #
  logpitheta <- dunif(beta_a, min = 1.0, max = exp(2), log = TRUE) +
    dunif(beta_s, min = 1.0, max = exp(2), log = TRUE) +
    dunif(rho, min = 0, max = 1, log = TRUE) +
    dunif(alpha_s, min = exp(0.05), max = exp(0.25), log = TRUE) +
    dunif(tau, min = exp(0.2), max = exp(0.6), log = TRUE)
    # dinvgamma(tau, shape = 2, scale = 1/0.1, log = TRUE)
  return(logpitheta)
}

log_likelihood <- function(pars, data, X_Initial){
  parameters <- c(beta_a = pars[1],
                  beta_s = pars[2],
                  rho = pars[3],
                  alpha_s = pars[4])
  tau <- pars[5]
  Disc_time <- seq(0, 22, 1)
  TimeSolve <- seq(0, 22, 0.01)
  Output_Ode <- ode(X_Initial,
                    TimeSolve,
                    OdeSystem,
                    parameters,
                    method = "lsoda")
  Cum_D_Day  <- Output_Ode[, 7][Output_Ode[, 1] %in% Disc_time]
  Cum_Is_Day <- Output_Ode[, 8][Output_Ode[, 1] %in% Disc_time]
  ####
  Log_Temp_Is <- sum(log(gamma(data[,2] + tau)) - log(gamma(tau))) +
                     tau * sum(log(tau) - log(tau + Cum_Is_Day)) +
                     sum(data[, 2] * (log(Cum_Is_Day) - log(tau + Cum_Is_Day)))
  TT <- is.nan(Log_Temp_Is)
  if (TT == TRUE) {
    print(paste(sum(log(tau + Cum_Is_Day)), 1))
    print(paste(sum(log(Cum_Is_Day)), 1))
    print(Cum_Is_Day)
    print(c(pars[1], pars[2], pars[3], pars[4]))
    break
  }
  return(Log_Temp_Is)
}

FuncLogAproxPosterior <- function(pars, data) {
  N_Total <- 858638
  E0   <- 0
  Is0  <- data[1,2]
  Ia0  <- 0
  R0   <- 0
  D0   <- data[1,1]
  CIs0 <- data[1,2]
  S0   <- N_Total - Is0 - E0 - Ia0

  X_Init <- c(S = S0,
              E = E0,
              Is = Is0,
              Ia = Ia0,
              R = R0,
              D = D0,
              CIs = CIs0)

  Logpitheta <- Prior(pars)
  ValLogVero <- log_likelihood(pars, Data, X_Init)

  LogAproxPosterior <- ValLogVero + Logpitheta
  return(LogAproxPosterior)
}

MatrizCovInst <- matrix(c(0.0001, 0, 0, 0, 0,
                          0, 0.001, 0, 0, 0,
                          0, 0, 0.0001, 0, 0,
                          0, 0, 0, 0.0001, 0,
                          0, 0, 0, 0, 0.001), 5)

mcmc_r <- Metro_Hastings(li_func = FuncLogAproxPosterior,
                         pars = c(2, 2, 1, 1.25, 1.25),
                         prop_sigma = MatrizCovInst,
                         iterations = 1000,
                         burn_in = 900,
                         par_names = c('beta_a','beta_s','rho','alpha_s','tau'),
                         data = Data)

posterior <- mcmc_r$trace
par(mfrow = c(3, 2))
plot(posterior[, 1],type = "l")
plot(posterior[, 2],type = "l")
plot(posterior[, 3],type = "l")
plot(posterior[, 4],type = "l")
plot(posterior[, 5],type = "l")

par(mfrow = c(3, 2))
hist(posterior[, 1])
hist(posterior[, 2])
hist(posterior[, 3])
hist(posterior[, 4])
hist(posterior[, 5])

mcmc_r <- mcmc_thin(mcmc_r)

mcmc_r <- Metro_Hastings(li_func = FuncLogAproxPosterior,
                         pars = c(2, 2, 1, exp(0.25), exp(0.6)),
                         prop_sigma = mcmc_r$prop_sigma,
                         iterations = 10000,
                         burn_in = 1000,
                         par_names = c('beta_a','beta_s','rho','alpha_s','tau'),
                         data = Data)
####
mcmc_r <- mcmc_thin(mcmc_r)
plotMH(mcmc_r)
BCI(mcmc_r)

#### Post analysis ####
coda::autocorr.plot(as.mcmc(posterior[1,]))
coda::autocorr.diag(as.mcmc(posterior[1,]))
coda::effectiveSize(as.mcmc(posterior[1,]))
