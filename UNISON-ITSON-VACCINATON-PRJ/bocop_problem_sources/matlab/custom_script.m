% custom script for solution postprocessing (bocop nlp)
% Pierre Martinon, 2016


% read solution file
path = input('Solution name (default: problem.sol): ','s');
if isempty(path)
    path = 'problem.sol'
end
verbose = 0;
[time,stage,state,control,control_average,algebraic,constants,optimvars,...
    path_constraint,adjoint,pathconstraint_mult,zl_state,zu_state,zl_control,zu_control,...
    objective,constraints_Infnorm,iterations] = readsolfile2(path,verbose);

% graphs etc ...