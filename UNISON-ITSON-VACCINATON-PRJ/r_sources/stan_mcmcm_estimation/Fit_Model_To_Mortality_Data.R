# Title :     Estimation Program
# Objective : Parameter estimation before mitigation measures were implemented,
#             that is, from February 19 to October 31, 2020
# Created by: Biomathematics Sonora Team
# Created on: 25/01/2021
#
library(deSolve)
library(dplyr)
library(rstan)
library(outbreaks)
library(bayesplot)
library(data.table)
library(ggpubr)
library(lubridate)

#
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#
#
covid19_CDMX_MX_data <- fread("CDMX_Data.csv",
                              select = c("Date",
                                         "Deaths_PD",
                                         "Cumulative_Deaths"))
covid19_CDMX_MX_data <- data.frame(covid19_CDMX_MX_data)
onset <- covid19_CDMX_MX_data %>% select(Date)
#
#
cum_cases   <- covid19_CDMX_MX_data %>% select(Deaths_PD)
cum_cases   <- unlist(cum_cases, use.names = FALSE)
N           <- nrow(onset) # Number of days observed throughout the outbreak
pop         <- 26446435    # Population
sample_time <- 1:(N + 13)
#
#### Modify data into a form suitable for Stan
covid19_data <- list(n_obs = N,
                     n_theta = 5,
                     n_difeq = 6,
                     n_pop = pop,
                     y = cum_cases,
                     t0 = 0,
                     ts = sample_time)
#
#### Specify parameters to monitor
parameters <- c("y_hat", "y_init", "theta", "R_0", "phi_inv")
#
#### Model 1 - Poisson model ####
m1 <- stan_model("vaccination_bneg.stan")
#
#### mcmcm parameters ####
n_chains  <- 5
n_warmups <- 500
n_iter    <- 100500
n_thin    <- 50
set.seed(1234)
#
#### Set initial values:
ini_1 <- function(){
  list(theta = c(rnorm(1, mean = 0.5, sd = 0.1),
                 rnorm(1, mean = 0.5, sd = 0.1),
                 runif(1, 0.05, 0.25),
                 runif(1, 0.05, 1),
                 runif(1, 0.05, 1)
  ),
  E0  = runif(1, 1, 15),
  Ia0 = runif(1, 1, 15),
  Is0 = runif(1, 1, 5)
  )
}
#
#### Sampling
nuts_fit_1 <- sampling(m1,
                       data = covid19_data,
                       pars = parameters,
                       init = ini_1,
                       chains = n_chains,
                       warmup = n_warmups,
                       iter = n_iter,
                       thin = n_thin,
                       seed = 13219,
                       control = list(adapt_delta = 0.95, max_treedepth = 15))
#
#### Summary Posterior Distributions
nuts_fit_1_summary <- summary(nuts_fit_1,
                              pars = c("lp__",
                                       "theta[1]",
                                       "theta[2]",
                                       "theta[3]",
                                       "theta[4]",
                                       "theta[5]",
                                       "phi_inv",
                                       "R_0",
                                       "y_init[2]",
                                       "y_init[3]",
                                       "y_init[4]"))$summary
print(nuts_fit_1_summary, scientific = FALSE, digits = 4)
#
#### Post analysis
posts_1 <-  rstan::extract(nuts_fit_1)
mod1_diagnostics  <- rstan::get_sampler_params(nuts_fit_1)
#
#### Check for divergent transitions
rstan::check_divergences(nuts_fit_1)
posterior_1 <- as.array(nuts_fit_1)
color_scheme_set("viridis")
#
#### Markov chain traceplots
mcmc_trace(posterior_1, pars = "lp__")
mcmc_trace(posterior_1, pars = c("theta[1]",
                                 "theta[2]",
                                 "theta[3]",
                                 "theta[4]",
                                 "theta[5]",
                                 "y_init[2]",
                                 "y_init[3]",
                                 "y_init[4]"))
mcmc_trace(posterior_1, pars = "R_0")
#
#### Univariate and bivariate marginal posterior distributions
pairs(nuts_fit_1,
      pars = c("lp__",
               "theta[1]",
               "theta[2]",
               "theta[3]",
               "theta[4]",
               "theta[5]",
               "R_0",
               "y_init[2]",
               "y_init[3]",
               "y_init[4]"
      ),
      labels = c(
          "lp_", 
          "beta_a",
          "beta_s",
          "rho",
          "eps_1",
          "eps_2",
          "R_0",
          "E0",
          "Is0", 
          "Ia0"),
      cex.labels = 1.5,
      font.labels = 9,
      condition = "accept_stat__")
#
#### Kernel density estimates of each Markov chain separately, overlaid
mcmc_dens_overlay(posterior_1,
                  pars = c("theta[1]", "theta[2]", "y_init[1]"))
#
#### Central posterior uncertainty intervals
mcmc_intervals(posterior_1 ,
               pars = c("theta[1]",
                        "theta[2]",
                        "theta[3]"
               ))

mcmc_intervals(posterior_1 , 
               pars = c("y_init[2]",
                        "y_init[4]",
                        "R_0"
               ))
#

### Eliminate string that does not converge
Temp_New_1 <- as.matrix(nuts_fit_1, pars = c("y_init", "theta", "R_0"))
Temp_New_2 <- R0_TT_1[2001:10000,]
###

write.table(Temp_New_2[,7:12], file = "theta_values_v1.csv", sep = ",", 
            col.names = c("beta_a", "beta_s", "rho", "eps_1", "eps_2", "R_0"),
            row.names = FALSE)

write.table(Temp_New_2[,1:6], file = "IC_v1.csv", sep = ",", 
            col.names = c("S0", "E0", "Is0", "Ia0", "R0", "D0"),
            row.names = FALSE)


theta_val    <- read.csv("theta_values_v1.csv", head = TRUE)
IC_val       <- read.csv("IC_v1.csv", head = TRUE)
Data_CDMXEST <- read.csv("CDMX_DataAll.csv", head = TRUE)

## Solve ODE system
ode1 <- function(t, x1, par1){ 
  S   = x1[1]
  E   = x1[2]
  Is  = x1[3]
  Ia  = x1[4]
  R   = x1[5]
  D   = x1[6]
  CIs = x1[7]
  ##  
  beta_a  = par1[1]
  beta_s  = par1[2]
  delta_e = par1[3]
  rho     = par1[4]
  alpha_s = par1[5]
  alpha_a = par1[6]
  mu_h    = par1[7]
  alpha_r = par1[8]
  eps_1   = par1[9]
  eps_2   = par1[10]
  ## Infection force
  N_star = S + E + Is + Ia + R
  FOI_1  = (beta_a * Ia + beta_s * Is) * S  / N_star
  #
  if (t <= 33){
    dS   = mu_h * N_star - FOI_1 - mu_h * S + alpha_r * R
    dE   = FOI_1 - (delta_e + mu_h) * E
    dIs  = rho * delta_e * E - (alpha_s + mu_h) * Is
    dIa  = (1 - rho) * delta_e * E - (alpha_a + mu_h) * Ia
    dR   = alpha_a * Ia + 0.89 * alpha_s * Is - (alpha_r + mu_h) * R
    dD   = 0.11 * alpha_s * Is
    dCIs = rho * delta_e * E
  }
  else if (33 < t && t <= 71)
  {
    dS   = mu_h * N_star - eps_1 * FOI_1 - mu_h * S + alpha_r * R
    dE   = eps_1 * FOI_1 - (delta_e + mu_h) * E
    dIs  = rho * delta_e * E - (alpha_s + mu_h) * Is
    dIa  = (1 - rho) * delta_e * E - (alpha_a + mu_h) * Ia
    dR   = alpha_a * Ia + 0.89 * alpha_s * Is - (alpha_r + mu_h) * R
    dD   = 0.11 * alpha_s * Is
    dCIs = rho * delta_e * E
  }
  else 
  {
    dS   = mu_h * N_star - eps_2 * eps_1 * FOI_1 - mu_h * S + alpha_r * R
    dE   = eps_2 * eps_1 * FOI_1 - (delta_e + mu_h) * E
    dIs  = rho * delta_e * E - (alpha_s + mu_h) * Is
    dIa  = (1 - rho) * delta_e * E - (alpha_a + mu_h) * Ia
    dR   = alpha_a * Ia + 0.89 * alpha_s * Is - (alpha_r + mu_h) * R
    dD   = 0.11 * alpha_s * Is
    dCIs = rho * delta_e * E
  }

  return(list(c(dS, dE, dIs, dIa, dR, dD, dCIs)))
}

N  = 26446435

di = as.Date("2020-02-19")
df = as.Date("2020-10-31")

Date1 = seq(di, df, 1)

nd1 = length(Date1)

disc_time_1 = seq(0, nd1-1, 1)

t1 = seq(0, nd1-1, 0.1)

vect <- seq(1, 8000,1)
MIs  <- matrix(data = 0, nrow = length(disc_time_1), ncol = length(vect))
MD   <- matrix(data = 0, nrow = length(disc_time_1), ncol = length(vect))
MIs1 <- matrix(data = 0, nrow = length(disc_time_1), ncol = length(vect))
MD1  <- matrix(data = 0, nrow = length(disc_time_1), ncol = length(vect))
MSt  <- matrix(data = 0, nrow = length(vect), ncol = 7)

j <- 1

#### Solve ODE system for each parameter set
for(i in vect){
  par1 = c(beta_a = theta_val[i,1], beta_s = theta_val[i,2], delta_e = 1/5.1,
           rho = theta_val[i,3], alpha_s = 1/10.81, alpha_a = 1/5.97,
           mu_h = 0.000039139, alpha_r = 1/365, eps_1 = theta_val[i,4], eps_2 = theta_val[i,5])
  
  z1   = c(IC_val[i,1], IC_val[i,2], IC_val[i,3], IC_val[i,4], IC_val[i,5], IC_val[i,6], IC_val[i,3])
  
  Y1   = as.matrix(ode(func = ode1, y = z1, times = t1, parms = par1[1:10], method = "rk4"))[,1:8]
  
  # Cumulative per day
  CD1  = Y1[,7][Y1[,1] %in% disc_time_1]
  CIs1 = Y1[,8][Y1[,1] %in% disc_time_1]
  
  # Incidence per day
  NCD1  = (CD1[2:length(disc_time_1)]  - CD1[1:(length(disc_time_1) - 1)])
  NCIs1 = (CIs1[2:length(disc_time_1)] - CIs1[1:(length(disc_time_1) - 1)])
  
  MSt[j,] <- Y1[length(t1),2:8]
  MIs[,j] <- CIs1
  MD[,j]  <- CD1
  
  MIs1[,j] <- c(Y1[1,8], NCIs1)  
  MD1[,j]  <- c(Y1[1,7], NCD1)  
  
  j <- j + 1
}

M1 <- apply(MD1, 1, quantile, probs = c(0.0005, 0.5, 0.9995))
M2 <- apply(MD, 1, quantile, probs = c(0.0005, 0.5, 0.9995))

## Figure
nd1   = length(Date1)

di2 = as.Date("2020-02-19")
df2 = as.Date("2020-10-31")

Date2 = seq(di2, df2, 1)
nd2   = length(Date2)

plIs1 <- ggplot() +
  geom_point(aes(x = Date2, y = Data_CDMXEST[1:nd2,8]), size = 1,
             color = "red", alpha = 0.75) +  
  geom_line(aes(x = Date1, y = M1[1,1:nd1], color = "Quantile 2.5%"), size = 0.75) +
  geom_line(aes(x = Date1, y = M1[2,1:nd1], color = "Quantile 50%"), size = 0.75) +
  geom_line(aes(x = Date1, y = M1[3,1:nd1], color = "Quantile 97.5%"), size = 0.75) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8,face="bold"),
        legend.position=c(0.2, 0.8), legend.title = element_blank(),
        legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"))  +
  labs(x = "", y = "Is Per Day", color = "Variable")+
  scale_color_manual(name="Variable",
                     breaks = c("Quantile 2.5%", "Quantile 50%", "Quantile 97.5%"),
                     values = c("Quantile 50%" = "blue", "Quantile 2.5%" = "grey75",
                                "Quantile 97.5%" = "grey55")) + 
  ggtitle("From February 19 to March 23")

plD1 <- ggplot() +
  geom_point(aes(x = Date2, y = Data_CDMXEST[1:nd2,9]), size = 1,
             color = "red", alpha = 0.75) +  
  geom_line(aes(x = Date1, y = M2[1,1:nd1], color = "Quantile 2.5%"), size = 0.75) +
  geom_line(aes(x = Date1, y = M2[2,1:nd1], color = "Quantile 50%"), size = 0.75) +
  geom_line(aes(x = Date1, y = M2[3,1:nd1], color = "Quantile 97.5%"), size = 0.75) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8,face="bold"),
        legend.position=c(0.2, 0.8), legend.title = element_blank(),
        legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"))  +
  labs(x = "", y = "Is Per Day", color = "Variable")+
  scale_color_manual(name="Variable",
                     breaks = c("Quantile 2.5%", "Quantile 50%", "Quantile 97.5%"),
                     values = c("Quantile 50%" = "blue", "Quantile 2.5%" = "grey75",
                                "Quantile 97.5%" = "grey55")) + 
  ggtitle("From February 19 to March 23")

cairo_ps(filename = "Figu_CIs_v1.eps",
         width = 9, height = 4.5, pointsize = 8,
         fallback_resolution = 1080)
print(plIs1)
dev.off()