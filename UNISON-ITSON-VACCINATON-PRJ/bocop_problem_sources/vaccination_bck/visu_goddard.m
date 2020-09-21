% visualization for Goddard problem
% Pierre Martinon
% Inria Saclay
% 2018

clear all;
close all;

% read solution
verbose = 0;
path = input('Solution file path (default: problem.sol): ','s');
if isempty(path)
    path = '/home/saul/Insync/saul.diazinfante@unison.mx/OneDrive Biz/UNISON/ARTICLES/Covid19/COVID19-Sonora/bocop_problem_sources/vaccination/unison_prj.sol';
end
[dim_state,dim_control,dim_algebraic,dim_optimvars,dim_constant,...
    dim_boundarycond,dim_pathcond,time_steps,...
    time,stage,state,control,control_average,algebraic,constants,optimvars,...
    state_lb,state_ub,control_lb,control_ub,...
    boundarycond,boundarycond_lb,boundarycond_ub,...
    path_constraint, pathcond_lb,pathcond_ub,...
    dynamic_constraint,...
    adjoint,pathconstraint_mult,zl_state,zu_state,zl_control,zu_control,...
    objective,constraints_Infnorm,iterations] = readsolfile2(path,verbose);
% plot trajectory
figure()
subplot(2,2,1); hold on
plot(time,state(:,1),'Linewidth',2)
title('POSITION')
subplot(2,2,2); hold on
plot(time,state(:,2),'Linewidth',2)
title('SPEED')
subplot(2,2,3); hold on
plot(time,state(:,3),'Linewidth',2)
title('MASS')
subplot(2,2,4); hold on
plot(time(1:end-1),control_average(:,1),'r','Linewidth',2)
title('ACCELERATION')
