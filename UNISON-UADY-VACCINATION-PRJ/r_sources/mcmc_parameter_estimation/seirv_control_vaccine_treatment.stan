functions {
  real[] seirvt(real t,  // time
             real[] y,
             real[] theta,
             real[] x_r,
             int[] x_i) {

    real dy_dt[6];
    real s = y[1];
    real e = y[2];
    real i_s = y[3];
    real i_a = y[4];
    real r = y[5];
    real n_bar;
//
    real beta_s = theta[1];
    real beta_a = theta[2];
    real p = theta[3];
    real alpha_s = theta[4];
    real alpha_a = theta[5];
    real delta_e = theta[6];
//
    real force_infection;
    real mu = 3.913894e-05;
    n_bar = s + e + i_s + i_a + r;

    force_infection = (beta_s * i_s + beta_a * i_a) / n_bar;
    dy_dt[1] = mu * n_bar - force_infection *s - mu * s;
    dy_dt[2] = force_infection * s - (mu + delta_e) * e;
    dy_dt[3] = p * delta_e * e - (mu + alpha_s) * i_s;
    dy_dt[4] = (1.0 - p) * delta_e * e - (mu + alpha_a) * i_a;
    dy_dt[5] = alpha_a * i_s + alpha_a * i_a - mu * r;
    dy_dt[6] = p * delta_e * i_s;
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
    // real <lower = 0> theta[n_theta]; // model parameters {beta,gamma}
    real<lower = 0> theta[n_theta];
    real <lower = 0, upper = 1> E0;
    real <lower = 0, upper=1> IA0;
}
transformed parameters{
    real y_hat[n_obs, n_difeq]; // solution from the ODE solver
    real y_init[n_difeq];
    // initial conditions
    real IS0 = 4.15424288404158e-06;
    real R0 = 0.0;
    y_init[1] = 1.0 - (IS0 + IA0 + E0);
    y_init[2] = E0;
    y_init[3] = IS0;
    y_init[4] = IA0;
    y_init[5] = R0;
    y_init[6] = IS0;
    y_hat = integrate_ode_rk45(seirvt, y_init, t0, ts, theta, x_r, x_i);
}
model {
    real lambda[n_obs];      //poisson parameter
    //priors
    theta[1] ~ normal(1.0, 0.1);
    theta[2] ~ normal(1.0, 0.1);
    theta[3] ~ uniform(0.3, 0.8);
    theta[4] ~ gamma(10, 100);
    theta[5] ~ gamma(10, 50);
    theta[6] ~ gamma(10, 40);

    E0 ~ uniform(1.038561e-06, 1.038561e-05);
    IA0 ~ uniform(1.038561e-06, 1.038561e-05);
//    likelihood
    for (i in 1:n_obs){
       lambda[i] = y_hat[i, 6] * n_pop;
    }
    y ~ poisson(lambda);
}
generated quantities {
    real R_0;      // Basic reproduction number
    real mu = 3.913894e-05;
    real beta_s = theta[1];
    real beta_a = theta[2];
    real p = theta[3];
    real alpha_s = theta[4];
    real alpha_a = theta[5];
    real delta_e = theta[6];
    R_0 = p * beta_s * delta_e / ((alpha_s + mu) * (delta_e + mu))
        + (1.0 - p) * beta_a * delta_e / ((alpha_a + mu) * (delta_e + mu));
}

generated quantities {
    real R_0;      // Basic reproduction number
    real mu = 3.913894e-05;
    real beta_s = theta[1];
    real beta_a = theta[2];
    real p = theta[3];
    real alpha_s = theta[4];
    real alpha_a = theta[5];
    real delta_e = theta[6];
    R_0 = p * beta_s * delta_e / ((alpha_s + mu) * (delta_e + mu))
        + (1.0 - p) * beta_a * delta_e / ((alpha_a + mu) * (delta_e + mu));
}
