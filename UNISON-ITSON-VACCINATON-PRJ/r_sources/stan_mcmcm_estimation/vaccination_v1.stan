functions {
  real[] LSEIR(real t,
               real[] y,
               real[] theta,
               real[] x_r,
               int[] x_i) {
    real dy_dt[7];
    real mu_h = 0.000039139;
    real delta_e = 1 / 5.1;
    real alpha_a = 1 / 5.97;
    real alpha_s = 1 / 10.81;
    real delta_r = 1 / 365.0;
    real N_star = y[1] + y[2] + y[3] + y[4] + y[5];
    real FOI_1 = (theta[1] * y[4] + theta[2] * y[3]) * y[1] / N_star;
    //
    dy_dt[1] = mu_h * N_star - FOI_1 - mu_h * y[1] + delta_r * y[5];
    dy_dt[2] = FOI_1 - (delta_e + mu_h) * y[2];
    dy_dt[3] = theta[3] * delta_e * y[2] - (alpha_s + mu_h) * y[3];
    dy_dt[4] = (1.0 - theta[3]) * delta_e * y[2] - (alpha_a + mu_h) * y[4];
    dy_dt[5] = alpha_a * y[4] + alpha_s * y[3] - (delta_r + mu_h) * y[5];
    dy_dt[6] = 0;
    dy_dt[7] = theta[3] * delta_e * y[2];

    return dy_dt;
  }
  
}
data {
  int<lower = 1> n_obs;   // number of days observed
  int<lower = 1> n_theta; // number of model parameters
  int<lower = 1> n_difeq; // number of differential equations
  int<lower = 1> n_pop;   // population
  int y[n_obs];           // data, total number of infected individuals each day
  real t0;                // initial time point (zero)
  real ts[n_obs];         // time points observed
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real <lower = 0.05> theta[n_theta]; // model parameters {beta_s,beta_a,p}
  real <lower = 2, upper = 20> E0;    // initial value of exposed individuals
  real <lower = 2, upper = 10> Ia0;   // initial value of asymptomatic individuals
  real <lower = 0> phi_inv;
}

transformed parameters{
  real y_hat[n_obs, n_difeq]; // solution from the ODE solver
  real y_init[n_difeq];       // initial conditions for other variable
  real Is0 = 1;               // initial symptomatic infected by onset symptoms
  real phi = 1. / phi_inv;
  
  y_init[1] = n_pop - (Is0 + Ia0 + E0);
  y_init[2] = E0;
  y_init[3] = Is0;
  y_init[4] = Ia0;
  y_init[5] = 0;
  y_init[6] = 0;
  y_init[7] = Is0;
  
  y_hat = integrate_ode_rk45(LSEIR, y_init, t0, ts, theta, x_r, x_i);
}

model {
  real lambda[n_obs];      //negative binomial parameter
  //priors
  theta[1] ~ normal(1.0, 0.13);
  theta[2] ~ normal(1.0, 0.13);
  theta[3] ~ uniform(0, 0.25);
  E0  ~ uniform(2, 20);
  Ia0 ~ uniform(2, 10);
  phi_inv ~ exponential(4);
  
  //likelihood
  lambda[1] = y_hat[1, 7];
  
  for (i in 2:n_obs){
    lambda[i] = (y_hat[i, 7] - y_hat[(i - 1), 7]);
  }
  //y ~ poisson(lambda);
  y ~ neg_binomial_2(lambda, phi);
}

generated quantities {
  real R_0;      // Basic reproduction number
  real mu_h = 0.000039139;
  real delta_e = 1 / 5.1;
  real alpha_a = 1 / 5.97;
  real alpha_s = 1 / 10.81; 
  R_0 = (delta_e/(delta_e + mu_h)) * ((theta[2]*theta[3]/(alpha_s + mu_h)) 
                                      + ((1.0 - theta[3])*theta[1]/(alpha_a + mu_h)));
}
