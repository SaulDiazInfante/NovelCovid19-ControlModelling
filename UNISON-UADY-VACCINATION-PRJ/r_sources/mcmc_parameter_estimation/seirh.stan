functions {
    real[] seirvt(real t,  // time
             real[] y,
             real[] theta,
             real[] x_r,
             int[] x_i) {
//
    real dy_dt[8];
    real s = y[1];
    real e = y[2];
    real i_s = y[3];
    real i_a = y[4];
    real h = y[5];
    real r = y[6];
    real d = y[7];
    real n_star;
//
    real beta_s = theta[1];
    real beta_a = theta[2];
    real kappa = theta[3];
    real p = theta[4];
    real delta_h = theta[5];
    real mu_i_s = theta[6];
    real mu_h = theta[7];
    real gamma_s = theta[8];
    real gamma_a = theta[9];
    real gamma_h = theta[10];
//
    real force_infection;
    real mu = 3.913894e-05;
    n_star = s + e + i_s + i_a + h + r;

    force_infection = (beta_s * i_s + beta_a * i_a) / n_star;
    dy_dt[1] = mu * n_star - force_infection * s - mu * s;
    dy_dt[2] = force_infection * s - (mu + kappa) * e;
    dy_dt[3] = p * kappa * e - (gamma_s + mu_i_s + delta_h) * i_s;
    dy_dt[4] = (1.0 - p) * kappa * e - (gamma_a + mu) * i_a;
    dy_dt[5] = delta_h * i_s  - (gamma_h + mu_h + mu) * h;
    dy_dt[6] = gamma_s * i_s + gamma_a * i_a  + gamma_h * h - mu * r;
    dy_dt[7] = mu_i_s * i_s + mu_h * h;
    dy_dt[8] = p * kappa * i_s;
    return dy_dt;
  }
}
data {
    int<lower = 1> n_obs;       // number of days observed
    int<lower = 1> n_theta;     // number of model parameters
    int<lower = 1> n_difeq;     // number of differential equations
    int<lower = 1> n_pop;       // population
    int y[n_obs];               // data, total number of infected
    real t0;                    // initial time point (zero)
    real ts[n_obs];             // time points observed
}

transformed data {
  real x_r[0];
  int x_i[0];
}
parameters {
    real <lower = 5e-6> theta[n_theta];
    // real <lower = 0, upper = 1> E0;
    // real <lower = 0, upper=1> IA0;
    real <lower = 0, upper = n_pop> E0;
    real <lower = 0, upper = n_pop> IA0;
}
transformed parameters{
    real y_hat[n_obs, n_difeq]; // solution from the ODE solver
    real y_init[n_difeq];
    // initial conditions
    // real IS0 = 8.297217e-06;
    real IS0 = 22.0;
    real H0 = 0.0;
    real R0 = 0.0;
    real D0 = 0.0;
    real CIS0 = 74.0;
    y_init[1] = n_pop - (IS0 + IA0 + E0);
    y_init[2] = E0;
    y_init[3] = IS0;
    y_init[4] = IA0;
    y_init[5] = H0;
    y_init[6] = R0;
    y_init[7] = D0;
    y_init[8] = CIS0;
    y_hat = integrate_ode_rk45(seirvt, y_init, t0, ts, theta, x_r, x_i);
}
model {
    real lambda[n_obs];             // poisson parameter
                                    // priors
    theta[1] ~ normal(1.0, 0.1);    // beta_s
    theta[2] ~ normal(1.0, 0.1);    // beta_a
    theta[3] ~ gamma(10, 40);       // kappa
    theta[4] ~ uniform(0.1, 0.5);   // p
    theta[5] ~ gamma(10, 40);       // delta_h
    theta[6] ~ exponential(33.33333);   // mu_i_s
    theta[7] ~ exponential(25);   // mu_h
    theta[8] ~ gamma(10, 100);      // gamma_s
    theta[9] ~ gamma(10, 50);       // gamma_a
    theta[10] ~ gamma(10, 250);      // gamma_h
//
    // E0 ~ uniform(1.121246e-07, 3.363737e-05);
    // IA0 ~ uniform(7.848719e-06, 3.363737e-05);
    E0 ~ uniform(74, 2100);
    IA0 ~ uniform(74, 2100);
//    likelihood
    for (i in 1:n_obs){
       //lambda[i] = y_hat[i, 8] * n_pop;
        lambda[i] = y_hat[i, 8];
    }
    y ~ poisson(lambda);
}
generated quantities {
    real R_0;      // Basic reproduction number
    real mu = 3.913894e-05;
    real beta_s = theta[1];
    real beta_a = theta[2];
    real p = theta[4];
    real gamma_s = theta[4];
    real gamma_a = theta[5];
    real kappa = theta[6];
    R_0 = p * beta_s * kappa / ((gamma_s + mu) * (kappa + mu))
        + (1.0 - p) * beta_a * kappa / ((gamma_a + mu) * (kappa + mu));
}

