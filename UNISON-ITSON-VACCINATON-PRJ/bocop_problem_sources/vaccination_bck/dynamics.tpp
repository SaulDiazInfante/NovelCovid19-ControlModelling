#include "header_dynamics"
{
    double beta_s = constants[0];
    double beta_a = constants[1];
    double epsilon = constants[2];
//
    double delta_e = constants[3];
    double delta_v = constants[4];
    double delta_r = constants[5];
//
    double p = constants[6];
//
    double alpha_a = constants[7];
    double alpha_t = constants[8];
    double alpha_s = constants[9];
//
    double mu = constants[10];
    double mu_s = constants[11];
    double mu_t = constants[12];
//
    double lambda_v = constants[13];
    double lambda_t = constants[14];
//
    double a_s = constants[15];
    double a_d = constants[16];
    double a_v = constants[17];
    double c_v = constants[18];
    double c_t = constants[19];

//
    Tdouble s = state[0];
    Tdouble e = state[1];
    Tdouble i_s = state[2];
    Tdouble i_a = state[3];
    Tdouble r = state[4];
    Tdouble d =  state[5];
    Tdouble v =  state[6];
    Tdouble treat =  state[7];
    Tdouble n_bar = state[0] + state[1] + state[2] +
                    state[3] + state[4] + state[6] +
                    state[7];
//
    Tdouble u_v = control[0];
    Tdouble u_t = control[1];
    Tdouble force_infection = (beta_s * i_s + beta_a * i_a) / n_bar;

    state_dynamics[0] = mu * n_bar - force_infection * s -
                        (mu + lambda_v * u_v) * s  + delta_v * v + delta_r * r;
    state_dynamics[1] = force_infection * (epsilon * v + s) - 
                        (mu + delta_e) * e;
    state_dynamics[2] = p * delta_e * e -
                        (mu + mu_s + alpha_s + lambda_t * u_t) * i_s;
 
    state_dynamics[3] = (1 - p) * delta_e * e - (mu + alpha_a) * i_a;
    state_dynamics[4] = alpha_s * i_s + alpha_a * i_a + alpha_t * treat -
                        (mu  + delta_r) * r;
    state_dynamics[5] = mu_s * i_s + mu_t * treat;
    state_dynamics[6] = lambda_v * u_v * s - 
                        epsilon * force_infection * v - (mu + delta_v) * v;
    state_dynamics[7] = lambda_t * u_t * i_s -
                        (alpha_t + mu_t + mu) * treat;

    state_dynamics[8] = lambda_v * u_v * (s + e + i_a + r);
    // Mayer form
    state_dynamics[9] = a_s * i_s + a_d * d +
                        0.5 * (c_v * u_v * u_v + c_t * u_t * u_t);                        
}
