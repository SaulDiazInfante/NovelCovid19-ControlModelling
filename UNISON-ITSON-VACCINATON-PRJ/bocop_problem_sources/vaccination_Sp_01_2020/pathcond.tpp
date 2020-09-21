// Function for the path constraints of the problem
// a <= g(t,y,u,z,p) <= b

// The following are the input and output available variables
// for the path constraints of your optimal control problem.

// Input :
// dim_path_constraints : number of path constraints
// time : current time (t)
// initial_time : time value on the first discretization point
// final_time : time value on the last discretization point
// dim_* is the dimension of next vector in the declaration
// state : vector of state variables
// control : vector of control variables
// algebraicvars : vector of algebraic variables
// optimvars : vector of optimization vector of optimization parameters
// constants : vector of constants

// Output :
// path_constraints :
// vector of path constraints expressions ("g" in the example above)

// The functions of your problem have to be written in C++ code
// Remember that the vectors numbering in C++ starts from 0
// (ex: the first component of the vector state is state[0])

// Tdouble variables correspond to values that can change during optimization:
// states, controls, algebraic variables and optimization parameters.
// Values that remain constant during optimization
// use standard types (double, int, ...).
//
#include "header_pathcond"
{
    Tdouble s = state[0];
    Tdouble e = state[1];
    Tdouble i_s = state[2];
    Tdouble i_a = state[3];
    Tdouble r = state[4];
    Tdouble d =  state[5];
    Tdouble v =  state[6];
    // Tdouble x =  state[7];
    
    double h_max = constants[20];
    double hospitalization_rate = constants[23];
    double d_max = constants[24];
    double n_pop = constants[25];
    double i_max = 0.65;
    // dynamic consistence
    
    path_constraints[0] = s + e + i_s + i_a + r + d + v - 1e0;
    path_constraints[1] = -s;
    path_constraints[2] = -e;
    path_constraints[3] = -i_s;
    path_constraints[4] = -i_a;
    path_constraints[5] = -r;
    path_constraints[6] = -d;
    path_constraints[7] = -v; 
    // control objectives
    path_constraints[8] = hospitalization_rate  * i_s - h_max;
    //path_constraints[9] = d - d_max;
    //path_constraints[8] = i_s + i_a - i_max;
    
}
