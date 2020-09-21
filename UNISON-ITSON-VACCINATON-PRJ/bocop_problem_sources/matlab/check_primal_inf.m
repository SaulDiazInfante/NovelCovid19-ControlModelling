% script for primal infeasibility diagnosis
% Pierre Martinon
% Inria and CMAP Ecole Polytechnique
% 2016



% read solution
verbose = 0;
path = input('Solution file path (default: problem.sol): ','s');
if isempty(path)
    path = 'problem.sol';
end
[dim_state,dim_control,dim_algebraic,dim_optimvars,dim_constant,...
    dim_boundarycond,dim_pathcond,time_steps,...
    time,stage,state,control,control_average,algebraic,constants,optimvars,...
    state_lb,state_ub,control_lb,control_ub,...
    boundarycond,boundarycond_lb,boundarycond_ub,...
    pathcond,pathcond_lb,pathcond_ub,...
    dynamic_constraint,...
    adjoint,pathconstraint_mult,zl_state,zu_state,zl_control,zu_control,...
    objective,constraints_Infnorm,iterations] = readsolfile2(path,verbose);

% check constraints norm
constraints_Infnorm

% state
for i=1:dim_state
state_err(:,i) = min(state(:,i) - state_lb(i),0) + max(state(:,i) - state_ub(i),0);
end
norm_state_err = norm(state_err,inf)

% control
for i=1:dim_control
control_err(:,i) = min(control(:,i) - control_lb(i),0) + max(control(:,i) - control_ub(i),0);
end
norm_control_err = norm(state_err,inf)

% boundary conditions
BC_err = min(boundarycond - boundarycond_lb,0) + max(boundarycond - boundarycond_ub,0);
norm_boundarycond_err = norm(BC_err,inf)

% path constraints
if (dim_pathcond > 0)
for i=1:dim_pathcond
PC_err = min(pathcond(:,i) - pathcond_lb(i),0) + max(pathcond(:,i) - pathcond_ub(i),0);
end
norm_pathcond = norm(PC_err,inf)
end

norm_dynamic_constraints = norm(dynamic_constraint,inf)

