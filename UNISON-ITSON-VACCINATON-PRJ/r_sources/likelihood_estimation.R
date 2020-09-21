library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)
source("likelihood_function.R")
source("mckendrick_model.R")
# Load the data set
reported_data <- read.csv('./data/culiacan_prevalence_data.csv')
df <-data.frame(reported_data)
# INPUT
initial_state_values <- c(S = 762,
                          I = 1,
                          R = 0)
#
p_0 <- c(beta = 1.7, gamma = 0.45)
optimised <- optim(p_0,
                   fn = likelihood,
                   data = reported_data,
                   initial_state_values = initial_state_values,
                   rhs_ode_function = sir_model,
                   control = list(fnscale = -1)
)
par <- optimised$par
# ploting likelihood data fit and reported data
times <- seq(from = 0, to = 14, by = 0.1)
output <- as.data.frame(ode(y = initial_state_values,
                            times = times,
                            func = sir_model,
                            parms = par))
t_data <-output$time[output$time %in% reported_data$time]
I_t <- output$I[output$time %in% reported_data$time]
out_model <- data.frame(time = t_data, I = I_t)
# Plot of the data
ggplot() +
    geom_point(data = reported_data,
               aes(x = time, y = number_reported),
               colour = "red") +
    geom_point(data = out_model, aes(x = time, y = I)) +
    geom_line(data = output, aes(x=time, y=I)) +
    xlab("Time (days)") +
    ylab("Number of reported cases")

file_name <- "idm2_sir_reported_data.png"
ggsave(file_name)
