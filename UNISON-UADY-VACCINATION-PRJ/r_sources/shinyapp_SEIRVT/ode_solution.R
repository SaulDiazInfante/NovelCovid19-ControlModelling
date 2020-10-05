library("rjson")
require(shiny)
require(deSolve)
require(phaseR)

seirvt_mod = function(t, state, parameters){
    with(as.list(c(state, parameters)), {
    ####
    n_bar <- s + e + i_s + i_a + r + v + treat
    force_infection = (beta_s * i_s + beta_a * i_a) / n_bar
    rhs_s = mu * n_bar - force_infection * s - (mu + lambda_v) * s + delta_v * v
    rhs_e = force_infection * (epsilon * v + s) - (mu + delta_e) * e
    rhs_i_s = p * delta_e * e - (mu + mu_s + alpha_s + lambda_t) * i_s -
        (1.0 - q) * alpha_t * treat
    rhs_i_a = (1 - p) * delta_e * e - (mu + mu_a + alpha_a) * i_a
    rhs_r = alpha_s * i_s + alpha_a * i_a + q * alpha_t * treat - mu * r
    rhs_d = mu_s * i_s + mu_a * i_a
    rhs_v = lambda_v * s - epsilon * force_infection * v - (mu + delta_v) * v
    rhs_treat = lambda_t * i_s - (mu + alpha_t) * treat
    rhs = c(rhs_s, rhs_e, rhs_i_s, rhs_i_a, rhs_r, rhs_d, rhs_v, rhs_treat)
    return(list(rhs))
    })
}
#### load parameters from a json file####
loaded_paramters <- fromJSON(file = "model_vaccination_parameters.json")

#Plot1: renderPlot to be passed to UI tab 1
#output$plot1 = renderPlot({
parms = c(
    beta_s = loaded_paramters$beta_s,
    beta_a = loaded_paramters$beta_a,
    epsilon = loaded_paramters$epsilon,
    delta_e = loaded_paramters$delta_e,
    delta_v = loaded_paramters$delta_v,
    p = loaded_paramters$p,
    q = loaded_paramters$q,
    alpha_a = loaded_paramters$alpha_a,
    alpha_t = loaded_paramters$alpha_t,
    alpha_s = loaded_paramters$alpha_s,
    mu = loaded_paramters$mu,
    mu_s = loaded_paramters$mu_s,
    mu_a = loaded_paramters$mu_a,
    lambda_v = loaded_paramters$lambda_v,
    lambda_t = loaded_paramters$lambda_t,
    n_whole = loaded_paramters$n_whole,
    T = loaded_paramters$T
    )
####Initial Conditions####
n_whole = loaded_paramters$n_whole
initial_conditions = c(s = (n_whole - 2)/n_whole,
                       e = 0.0,
                       i_s = 1 / n_whole,
                       i_a = 1 / n_whole,
                       r = 0.0,
                       d = 0.0,
                       v = 0.0,
                       treat = 0.0
)
times = seq(0, 250, by = 1/1000)
#R0 = round(with(as.list(parms),
#                    beta_s + beta_a /
#                        (alpha_a + alpha_s + alpha_t + mu + mu_s + mu_a),
#    ),1
#    )
out = ode(y = initial_conditions,
              times = times,
              func = seirvt_mod,
              parms = parms)
out = as.data.frame(out)
#
#     #Plot1
#     #sel = out$time > input$T[1] & out$time < input$T[2]
plot(x = out$time,
          y = out$s,
          ylab = "fraction",
          xlab = "time",
          type = "l",
          #ylim = range(out[sel,-c(1,4)])
     )
     title(paste("R0=", R0))
     lines(x = out$time,  y = out$i_s, col = "red")
     lines(x = out$time, y = out$i_a, col = "green")
     legend("right",
            legend = c("I_s"),
            lty = c(1,1,1),
            col = c("red", "green")
)
#seirvt_mod(0, initial_conditions, parms)
