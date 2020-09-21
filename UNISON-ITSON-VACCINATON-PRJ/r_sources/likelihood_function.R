require(lubridate)
likelihood <- function (parameters,
                        data,
                        initial_state_values,
                        rhs_ode_function){
    last_index <- nrow(data)
    n_days <- time_length(
        difftime(as.Date(data[last_index, "FECHA_SINTOMAS"]),
                 as.Date(data[1, "FECHA_SINTOMAS"])),
                 "days"
        )
    times <- seq(from = 0, to = n_days, length.out = 1000)
    ### Calculate model output

    model_sample <- as.data.frame(ode(y = initial_state_values,
                                      times = times,
                                      func = rhs_ode_function,
                                      parms = parameters))
    model_sample_incidence <- model_sample$I[model_sample$time %in% data$time]
    # calculate likelihood function

    LL <- sum(
        dpois(
            x = data$number_reported,
            lambda = 0.6 * model_sample_incidence,
            log = TRUE
        )
    )
    return(LL)
}
