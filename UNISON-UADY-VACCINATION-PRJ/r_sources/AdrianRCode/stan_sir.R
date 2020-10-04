library(outbreaks)
library(tidyverse)
library(rstan, lib.loc = "/usr/lib/R/site-library")
library(StanHeaders, lib.loc = "/usr/lib/R/site-library")
library(gridExtra)

head(influenza_england_1978_school)
ggplot(data = influenza_england_1978_school) +
    geom_point(mapping = aes(x = date, y = in_bed)) +
    labs(y = "Number of students in bed")

rstan_options(auto_write = TRUE)
options (mc.cores = parallel::detectCores())
set.seed(3)

#### time series of cases ####
cases <- influenza_england_1978_school$in_bed  # Number of students in bed

# total count
N <- 763;

# times
n_days <- length(cases)
t <- seq(0, n_days, by = 1)
t0 = 0
t <- t[-1]

#### initial conditions ####
i0 <- 1
s0 <- N - i0
r0 <- 0
y0 = c(S = s0, I = i0, R = r0)

#### data for Stan ####
data_sir <- list(n_days = n_days,
                 y0 = y0,
                 t0 = t0,
                 ts = t,
                 N = N,
                 cases = cases)
####
niter <- 2000
model <- stan_model("sir.stan")
fit_sir_negbin <- sampling(model,
                           data = data_sir,
                           iter = niter,
                           chains = 4)
pars <- c('beta', 'gamma', "R0", "recovery_time")
print(fit_sir_negbin, pars = pars)
stan_dens(fit_sir_negbin, pars = pars, separate_chains = TRUE)
smr_pred <- cbind(as.data.frame(summary(
    fit_sir_negbin, pars = "pred_cases", probs = c(0.05, 0.5, 0.95))$summary), t, cases)
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names

ggplot(smr_pred, mapping = aes(x = t)) +
    geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "orange", alpha = 0.6) +
    geom_line(mapping = aes(x = t, y = X50.)) +
    geom_point(mapping = aes(y = cases)) +
    labs(x = "Day", y = "Number of students in bed")
