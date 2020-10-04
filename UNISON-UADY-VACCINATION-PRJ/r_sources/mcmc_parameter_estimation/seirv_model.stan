functions {
  real[] seirvt(real t,  // time
             real[] y,
             real[] theta,
             real[] x_r,
             int[] x_i) {

    real dy_dt[9];
    real n_bar;
    real force_infection;
    real mu = 3.913894e-05;
    real delta_e = 0.1960784;
    real mu_a = 0.125;

    n_bar = y[1] + y[2] +y[3] + y[4] + y[5] + y[6] + y[7] + y[8] + y[9];
    force_infection = (theta[2] * y[3] + theta[1] * y[4]) / n_bar;
    dy_dt[1] = mu * n_bar - force_infection * y[1] - (mu + theta[8]) * y[1] +
              theta[11] * y[7];
    dy_dt[2] = force_infection * (theta[4] * y[7] + y[1])
              - (mu + delta_e) * y[2];
    dy_dt[3] = theta[3] * delta_e * y[2]
              - (mu + theta[10] + theta[5] + theta[9]) * y[3]
              - (1.0 - theta[3]) * theta[7] * y[8];
    dy_dt[4] = (1.0 - theta[3]) * delta_e * y[2]
              - (mu + mu_a + theta[6]) * y[4];
    dy_dt[5] = theta[5] * y[3] + theta[6] * y[4]
                + theta[3] * theta[7] * y[8] - mu * y[5];
    dy_dt[6] = theta[10] * y[3] + mu_a * y[4];
    dy_dt[7] = theta[8] * y[1] - theta[4] * force_infection * y[7]
                - (mu + theta[11]) * y[7];
    dy_dt[8] = theta[9] * y[3] - (mu + theta[7]) * y[8];
    dy_dt[9] = theta[3] * delta_e * y[2];
    return dy_dt;
  }

}
data {
    int<lower = 1> n_obs;       // number of days observed
    int<lower = 1> n_theta;     // number of model parameters
    int<lower = 1> n_difeq;     // number of differential equations
    int<lower = 1> n_pop;       // population
    int y[n_obs];           // data, total number of infected individuals each day
    real t0;                // initial time point (zero)
    real ts[n_obs];         // time points observed
}

transformed data {
  real x_r[0];
  int x_i[0];
}
parameters {
    real <lower = 0> theta[n_theta]; // model parameters {beta,gamma}
    real <lower = 0, upper = 1> E0;  // initial fraction of exposed individuals
    real <lower=0, upper=1> IS0;
    real <lower=0, upper=1> IA0;
}
transformed parameters{
    real y_hat[n_obs, n_difeq]; // solution from the ODE solver
    real y_init[n_difeq];     // initial conditions

    y_init[1] = n_pop - (IS0 + IA0 + E0);
    y_init[2] = E0;
    y_init[3] = IS0;
    y_init[4] = IA0;
    y_init[5] = 0;
    y_init[6] = 0;
    y_init[7] = IS0;

    y_hat = integrate_ode_rk45(seirvt, y_init, t0, ts, theta, x_r, x_i);
}
model {
    real lambda[n_obs];      //poisson parameter
    //priors
    theta[1] ~ normal(1.0, 0.3);
    theta[2] ~ normal(1.0, 0.3);
    theta[3] ~ uniform(0, 0.5);
    theta[4] ~ gamma(10, 100);
    theta[5] ~ gamma(10, 50);
    theta[6] ~ gamma(10, 50);
    theta[7] ~ gamma(10, 50);
    theta[8] ~ gamma(10, 50);
    theta[9] ~ gamma(10, 50);
    theta[10] ~ gamma(10, 50);
    theta[11] ~ gamma(10, 50);

    E0 ~ uniform(0, 1.104167e-05);
    IS0 ~ uniform(0, 3.312501e-06);
    IA0 ~ uniform(0, 1.104167e-05);

  //likelihood
    for (i in 1:n_obs){
       lambda[i] = y_hat[i, 7] * n_pop;
    }
    y ~ poisson(lambda);
}
// generated quantities {
//  real R_0;      // Basic reproduction number
//  R_0 = theta[1]/theta[2];
//}
