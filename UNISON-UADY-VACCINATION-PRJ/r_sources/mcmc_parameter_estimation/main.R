library(deSolve)
library(dplyr)
library(rstan)
library(gridExtra)
library(outbreaks)
library(bayesplot)
library(data.table)
library(knitr)
library(kableExtra)
library(bayesplot)
source("load_data.R")
source("init_seir.R")
source("mcmc_stain_summary.R")
source("mcmc_post_analysis.R")
source("divergence_plots.R")
source("model_fit.R")
#
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#
data_cdmx <- load_data()
onset <- data_cdmx[[1]]
cum_cases <- data_cdmx[[2]]
cum_cases <- unlist(cum_cases, use.names = FALSE)
N <- nrow(onset)
pop <- 8918653 # Population
sample_time <- 1:N
#
#### SIR rstan model ####
mdl <- stan_model("seirh.stan")
#
#### mcmc parameters ####
n_chains <- 10
n_warmups <- 100
n_iter <- 50500
n_thin <- 50
set.seed(972198)
covid19_data <- list(n_obs = N,
                 n_theta = 10,
                 n_difeq = 8,
                 n_pop = pop,
                 y = cum_cases,
                 t0 = 0,
                 ts = sample_time)
parameters <- c("y_hat", "y_init", "theta", "R_0")
time.start_nuts <- Sys.time()
nuts_fit <-
        sampling(mdl,
            data = covid19_data,
            pars = parameters,
            init = init_seirh,
            chains = n_chains,
            warmup = n_warmups,
            iter = n_iter,
            thin = n_thin,
            control = list(adapt_delta = 0.85))
time.end_nuts <- Sys.time()
duration_nuts <- time.end_nuts - time.start_nuts
parameters = c("theta", "R_0", "y_init")
nuts_fit_summary <- summary(nuts_fit, pars = parameters)$summary
print(nuts_fit_summary,
    scientific = FALSE,
    digits = 4)
    sub_path_1 <- "/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19"
    sub_path_2 <- "NovelCovid19-ControlModelling/NovelCovid19-ControlModellingGitHub"
    sub_path_3 <- "UNISON-UADY-VACCINATION-PRJ"
    sub_path_4 <- "r_sources/mcmc_parameter_estimation/runs"
prefix_time <- Sys.time()
prefix_time <- paste(date(prefix_time),
                    hour(prefix_time),
                    minute(prefix_time), sep="_")
file_name <- paste("run-", prefix_time,".RData", sep="")
runs_path <- 
    paste(sub_path_1, sub_path_2, sub_path_3, sub_path_4, file_name, sep = "/") 
save.image(file=runs_path)
#
#### Post analysis ####
#
mcmcm_post_analysis(nuts_sample = nuts_fit)
#
### Model Fit ####
# Model fitted values across the observed time period
#
model_fit()
#### Divergence analysis ####
# divergence_analysis()

