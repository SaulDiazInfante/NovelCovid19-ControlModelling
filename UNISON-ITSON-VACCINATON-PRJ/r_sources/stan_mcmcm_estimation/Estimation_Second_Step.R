# Title :     
# Objective : Parameter estimation after mitigation measures were implemented,
#             that is, from March 23 to April 23
# Created by: Biomathematics Sonoran Team
# Created on: 09/03/2020
#
library(deSolve)
library(dplyr)
library(rstan)
library(outbreaks)
library(bayesplot)
library(data.table)
library(deSolve)
library(ggpubr)
library(lubridate)

#
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#
#
covid19_CDMX_MX_data <- fread("Data_CDMX_Estado_0309_v2.csv",
                              select = c("Date",
                                         "Onset_PD",
                                         "Cumulative_Onset"))
covid19_CDMX_MX_data <- data.frame(covid19_CDMX_MX_data)
onset <- covid19_CDMX_MX_data %>% select(Date)
#
#
cum_cases   <- covid19_CDMX_MX_data %>% select(Onset_PD)
cum_cases   <- unlist(cum_cases, use.names = FALSE)
N           <- nrow(onset) # Number of days observed throughout the outbreak
pop         <- 26446435    # Population
sample_time <- 1:N
#
#### Modify data into a form suitable for Stan
covid19_data <- list(n_obs = N,
                     n_theta = 1,
                     n_difeq = 7,
                     n_pop = pop,
                     y = cum_cases,
                     t0 = 0,
                     ts = sample_time)
#
#### Specify parameters to monitor
parameters <- c("y_hat", "y_init", "theta",  "R_0")
#
#### Model 1 - Poisson model ####
m1 <- stan_model("vaccination_v2.stan")
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
  list(theta = array(runif(1, 0.25, 0.75)))
}
#
#### Sampling
#time.start_nuts1 <- Sys.time()
nuts_fit_1 <- sampling(m1,
                       data = covid19_data,
                       pars = parameters,
                       init = ini_1,
                       chains = n_chains,
                       warmup = n_warmups,
                       iter = n_iter,
                       thin = n_thin,
                       seed = 13219,
                       control = list(adapt_delta = 0.99, max_treedepth = 15))
#time.end_nuts1 <- Sys.time()
#duration_nuts1 <- time.end_nuts1 - time.start_nuts1
#save.image(file = "Data2304_v1_70.RData")
#
#### Summary Posterior Distributions
nuts_fit_1_summary <- summary(nuts_fit_1,
                              pars = c("lp__",
                                       "theta[1]",
                                       "R_0"))$summary
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
mcmc_trace(posterior_1, pars = "theta[1]")
mcmc_trace(posterior_1, pars = "R_0")
#
#### Univariate and bivariate marginal posterior distributions
pairs(nuts_fit_1,
      pars = c("theta[1]",
               "R_0"
      ),
      labels = c("epsilon", "R0"),
      cex.labels = 1.5,
      font.labels = 9,
      condition = "accept_stat__")
#
#### Kernel density estimates of each Markov chain separately, overlaid
mcmc_dens_overlay(posterior_1,
                  pars = c("theta[1]", "theta[2]", "y_init[1]"))
#
#Central posterior uncertainty intervals
mcmc_intervals(posterior_1 ,
               pars = "theta[1]")

mcmc_intervals(posterior_1 , 
               pars = "R_0")
#
#### Save Estimate parameters and initial conditions
Mat_theta <- posts_1$theta
Mat_IC    <- posts_1$y_init

write.table(Mat_theta, file = "theta_values_v2.csv", sep = ",", 
            col.names = c("epsilon"),
            row.names = FALSE)

write.table(Mat_IC, file = "IC_v2.csv", sep = ",", 
            col.names = c("S0", "E0", "Is0", "Ia0", "R0", "D0", "CIs0"),
            row.names = FALSE)

## Solve ODE system
theta_val <- read.csv("theta_values_v2.csv", head = TRUE)
IC_val    <- read.csv("IC_v2.csv", head = TRUE)
Data_CDMXEST <- read.csv("Data_CDMX_Estado_0309_v2.csv", head = TRUE)

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
  delta_r = par1[8]
  epsilon = par1[9]
  ## Infection force
  N_star = S + E + Is + Ia + R
  FOI_1  = (beta_a * Ia + beta_s * Is) * S  / N_star
  
  dS   = mu_h * N_star - epsilon * FOI_1 - mu_h * S + delta_r * R
  dE   = epsilon * FOI_1 - (delta_e + mu_h) * E
  dIs  = rho * delta_e * E - (alpha_s + mu_h) * Is
  dIa  = (1 - rho) * delta_e * E - (alpha_a + mu_h) * Ia
  dR   = alpha_a * Ia + alpha_s * Is - (delta_r + mu_h) * R
  dD   = 0
  dCIs = rho * delta_e * E
  
  return(list(c(dS, dE, dIs, dIa, dR, dD, dCIs)))
}

N  = 26446435

di = as.Date("2020-03-23")
df = as.Date("2020-04-23")

Date1 = seq(di, df, 1)

nd1 = length(Date1)

disc_time_1 = seq(0, nd1-1, 1)

t1 = seq(0, nd1-1, 0.1)

vect <- seq(1,10000,1)
MIs  <- matrix(data = 0, nrow = length(disc_time_1), ncol = length(vect))

j <- 1

#### Solve ODE system for each parameter set
for(i in vect){
  par1 = c(beta_a = 0.6435, beta_s = 0.9322, delta_e = 1/5.1,
           rho = 0.1227, alpha_s = 1/10.81, alpha_a = 1/5.97,
           mu_h = 0.000039139, delta_r = 1/365, epsilon = theta_val[i,1])
  
  z1   = c(IC_val[i,1], IC_val[i,2], IC_val[i,3], IC_val[i,4], IC_val[i,5], IC_val[i,6],
           IC_val[i,7])
  
  Y1   = as.matrix(ode(func = ode1, y = z1, times = t1, parms = par1[1:9], method = "rk4"))[,1:8]
  
  # Cumulative per day
  CIs1 = Y1[,8][Y1[,1] %in% disc_time_1]
  
  # Incidence per day
  NCIs1 = (CIs1[2:length(disc_time_1)]  - CIs1[1:(length(disc_time_1) - 1)])
  
  MIs[,j] <- c(Y1[1,8],NCIs1)
  j <- j + 1
}

M1 <- apply(MIs, 1, quantile, probs = c(0.025, 0.5, 0.975))

## Figure
nd1   = length(Date1)

plIs1 <- ggplot() +
  geom_point(aes(x = Date1, y = Data_CDMXEST[1:nd1,2]), size = 1,
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
  ggtitle("From March 23 to April 23")

cairo_ps(filename = "Figu_CIs_v2.eps",
         width = 9, height = 4.5, pointsize = 8,
         fallback_resolution = 1080)
print(pl_Is1)
dev.off()