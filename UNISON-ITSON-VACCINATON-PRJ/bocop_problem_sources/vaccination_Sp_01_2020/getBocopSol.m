paht = '/home/saul/Desktop/WorikingProjects/Covid/COVID19-Sonora/bocop_problem_sources/vaccination_Aug_28_2020'...
verbpse = 1;
[dim_state,dim_control,dim_algebraic,dim_optimvars,dim_constant,...
    dim_boundarycond,dim_pathcond,time_steps,...
    time,stage,state,control,control_average,algebraic,constants,optimvars,...
    state_lb,state_ub,control_lb,control_ub,...
    boundarycond,boundarycond_lb,boundarycond_ub,...
    path_constraint, pathcond_lb,pathcond_ub,...
    dynamic_constraint,...
    adjoint,pathconstraint_mult,zl_state,zu_state,zl_control,zu_control,...
    objective,constraints_Infnorm,iterations] = readsolfile(path,verbose)