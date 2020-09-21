// Function for the dynamics of the problem
// dy/dt = dynamics(y,u,z,p)
//
// The following are the input and output available variables 
// for the dynamics of your optimal control problem.
//
// Input :
// time : current time (t)
// normalized_time: t renormalized in [0,1]
// initial_time : time value on the first discretization point
// final_time : time value on the last discretization point
// dim_* is the dimension of next vector in the declaration
// state : vector of state variables
// control : vector of control variables
// algebraicvars : vector of algebraic variables
// optimvars : vector of optimization parameters
// constants : vector of constants
//
// Output :
// state_dynamics : vector giving the expression 
// of the dynamic of each state variable.
//
// The functions of your problem have to be written in C++ code
// Remember that the vectors numbering in C++ starts from 0
// (ex: the first component of the vector state is state[0])
//
// Tdouble variables correspond to 
// values that can change during optimization:
// states, controls, algebraic variables and optimization parameters.
// Values that remain constant during 
// optimization use standard types (double, int, ...).

#include "header_dynamics"
{
    
    // HERE : description of the function for the dynamics
    // Please give a function or a value for the 
    // dynamics of each state variable
    //
    double beta_s = constants[0];
    double beta_a = constants[1];
    double epsilon = constants[2];
    double delta_u = constants[3];
    double delta_e = constants[4];
    double delta_s = constants[5];
    double p = constants[6];
    double alpha_a = constants[7];
    double alpha_h = constants[8];
    double alpha_s = constants[9];
    double mu_s = constants[10];
    double mu_h = constants[11];
    double h_max = constants[12];
    double a_h = constants[13];
    double a_d = constants[14];
    double a_u = constants[15];
//   
    Tdouble l = state[0];
    Tdouble s = state[1];
    Tdouble e = state[2];
    Tdouble i_s = state[3];
    Tdouble i_a = state[4];
    Tdouble h = state[5];
    Tdouble r = state[6];
    Tdouble d =  state[7];
    Tdouble x_zero = state[8];
    Tdouble n_aster = 
        state[0] + state[1] + state[2] + state[3] + state[4] +
        state[5] + state[6];
    Tdouble u = control[0];
    //
    n_aster = l + s + e + i_s + i_a + h + r;
    Tdouble force_infection = (beta_s * i_s + beta_a * i_a) / n_aster;
    //   
    state_dynamics[0] = - l * (
                                epsilon * force_infection 
                                + delta_u * u);
    state_dynamics[1] = delta_u * l - force_infection * s;
    state_dynamics[2] = force_infection * (epsilon * l + s)
                        - delta_e * e;
    state_dynamics[3] = p * delta_e * e 
                        - (alpha_s + mu_s + delta_s) * i_s;
    //
    state_dynamics[4] = (1 - p) * delta_e * e - alpha_a * i_a ;
    state_dynamics[5] = delta_s * i_s - (alpha_h + mu_h) * h;
    state_dynamics[6] = alpha_s * i_s + alpha_a * i_a + alpha_h * h ;
    state_dynamics[7] = mu_s * i_s + mu_h * h;
    // Mayer form
    //state_dynamics[8] = a_d * d  + a_u * u * u + a_h * h  ; 
    //state_dynamics[8] =  a_d * d +  a_h * h + a_u * u ;
    state_dynamics[8] =  a_d * d +  a_h * h;
}

