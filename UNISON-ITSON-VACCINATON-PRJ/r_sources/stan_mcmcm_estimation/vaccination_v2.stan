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
    real beta_a  = 0.6435;
    real beta_s  = 0.9322;
    real rho     = 0.1227;
    real N_star = y[1] + y[2] + y[3] + y[4] + y[5];
    real FOI_1 = (beta_a * y[4] + beta_s * y[3]) * y[1] / N_star;
    //
    dy_dt[1] = mu_h * N_star - theta[1] * FOI_1 - mu_h * y[1] + delta_r * y[5];
    dy_dt[2] = theta[1] * FOI_1 - (delta_e + mu_h) * y[2];
    dy_dt[3] = rho * delta_e * y[2] - (alpha_s + mu_h) * y[3];
    dy_dt[4] = (1.0 - rho) * delta_e * y[2] - (alpha_a + mu_h) * y[4];
    dy_dt[5] = alpha_a * y[4] + 0.89 * alpha_s * y[3] - (delta_r + mu_h) * y[5];
    dy_dt[6] = 0.11 * alpha_s * y[3];
    dy_dt[7] = rho * delta_e * y[2];

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
  real <lower = 0.25, upper = 0.75> theta[n_theta]; // model parameters {epsilon}
  real <lower = 0> phi_inv;
}

transformed parameters{
  real y_hat[n_obs, n_difeq]; // solution from the ODE solver
  real y_init[n_difeq];       // initial conditions for other variable
  real E0   = 6587.5850; 
  real Is0  = 553.7035;
  real Ia0  = 3149.9240;
  real R0   = 3001.5470;
  real D0   = 3.0;
  real CIs0 = 77.0;
  real phi  = 1. / phi_inv;
  
  y_init[1] = n_pop - (Is0 + Ia0 + E0 + R0 + D0);
  y_init[2] = E0;
  y_init[3] = Is0;
  y_init[4] = Ia0;
  y_init[5] = R0;
  y_init[6] = D0;
  y_init[7] = CIs0;
  
  y_hat = integrate_ode_rk45(LSEIR, y_init, t0, ts, theta, x_r, x_i);
}

model {
  real lambda[n_obs];      //negative binomial parameter
  //priors
  theta[1] ~ uniform(0.25, 0.75);
  phi_inv  ~ exponential(4);
  
  //likelihood
  lambda[1] = y_hat[1, 7];
  
  for (i in 2:n_obs){
    lambda[i] = (y_hat[i, 7] - y_hat[(i - 1), 7]);
  }
  y ~ neg_binomial_2(lambda, phi);
}

generated quantities {
  real R_0;      // Basic reproduction number
  real mu_h = 0.000039139;
  real delta_e = 1 / 5.1;
  real alpha_a = 1 / 5.97;
  real alpha_s = 1 / 10.81;
  real beta_a  = 0.6435;
  real beta_s  = 0.9322;
  real rho     = 0.1227;
  
  R_0 = (delta_e/(delta_e + mu_h)) * ((rho*beta_s*theta[1]/(alpha_s + mu_h)) 
                                      + ((1.0 - rho)*beta_a*theta[1]/(alpha_a + mu_h)));
}
