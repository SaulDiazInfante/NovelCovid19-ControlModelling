% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% This code is published under the Eclipse Public License
% File: visualization_tree.m
% Authors: Daphne Giorgi, Pierre Martinon

% Tree view for visualization

function visualization_tree()

global free_time initial_time final_time 
global dim_state dim_control dim_algebraic dim_parameter dim_constant dim_boundarycond dim_constraint 
global time time_steps discretization_method time_stage 
global optimization_type batch_index batch_nrange batch_lowerbound batch_upperbound batch_directory
global initialization_type initialization_file
global paramid_type paramid_file paramid_separator
global state_names control_names algebraic_names parameters_names boundarycond_names constant_names constraint_names
global solution_file constants 
global ifbounds ifbounds_type bstates bstates_type bcontrols bcontrols_type balgebraics balgebraics_type
global bparameters bparameters_type bpathconstraints bpathconstraints_type
global initstate_type initstate initstatedisc_dim initcontrol_type initcontrol initcontroldisc_dim
global initalgebraic_type initalgebraic initalgebraicdisc_dim ninitparameter initparameter 
global objective constraints_L2norm constraints_Infnorm
global alltimes stage ntime nstage
global state control control_average algebraic parameter boundarycond  path_constraint dynamic_constraint 
global boundarycond_mult pathconstraint_mult adjoint_state 
global coef_a coef_b coef_c
global zl_states zl_control zl_algebraic zl_parameter zu_states zu_control zu_algebraic zu_parameter
global hpath filename

% First we read the solution file
filename = '/home/saul/Bocop-2.2.1-linux/examples/goddard';

readsolfile(filename);

fprintf('Available variables: time, state, control_average, algebraic, parameter, path_constraint\n');

% Variables node
Variables = uitreenode('v0', 'Variables', 'Variables', [], false);

StateVariables = uitreenode('v0', 'StateVariables', 'State Variables', [], false);
ControlVariables = uitreenode('v0', 'ControlVariables', 'Control Variables', [], false);
AlgebraicVariables = uitreenode('v0', 'AlgebraicVariables', 'Algebraic Variables', [], false);
OptimizationParameters = uitreenode('v0', 'OptimizationParameters', 'Optimization Parameters', [], false);

Variables.add(StateVariables);
Variables.add(ControlVariables);
Variables.add(AlgebraicVariables);
Variables.add(OptimizationParameters);

% Constraints node
Constraints = uitreenode('v0', 'Constraints', 'Constraints', [], false);

PathConstraints = uitreenode('v0','PathConstraints','Path Constraints',[],false);
DynamicConstraints = uitreenode('v0','DynamicConstraints','Dynamic Constraints',[],false);

Constraints.add(PathConstraints);
Constraints.add(DynamicConstraints);

% Boundary conditions node
BoundaryConditions = uitreenode('v0', 'BoundaryConditions', 'Boundary Conditions', [], false);

% Multipliers node
Multipliers = uitreenode('v0', 'Multipliers', 'Multipliers', [], false);

BoundaryConditionsMult = uitreenode('v0', 'BoundaryConditionsMult', 'Boundary Conditions', [], false);
PathConstraintsMult = uitreenode('v0', 'PathConstraintsMult', 'PathConstraints', [], false);
AdjointStates = uitreenode('v0', 'AdjointStates', 'Adjoint States', [], false);

Multipliers.add(BoundaryConditionsMult);
Multipliers.add(PathConstraintsMult);
Multipliers.add(AdjointStates);

% Root node
root = uitreenode('v0', 'PlottingTree', 'Variables for display', [], false);
root.add(Variables);
root.add(Constraints);
root.add(BoundaryConditions);
root.add(Multipliers);

% LOOPS ON THE VARIABLES TO GENERATE THE LEAVES

% State variables
for i = 1:dim_state
    StateVariables.add(uitreenode('v0',strcat('state',num2str(i)),state_names{i},[],false));
end

% Control variables
for i = 1:dim_control
    ControlVariables.add(uitreenode('v0',strcat('control',num2str(i)),control_names{i},[],false));
end

% Algebraic variables
for i = 1:dim_algebraic
    AlgebraicVariables.add(uitreenode('v0',strcat('algebraic',num2str(i)),algebraic_names{i},[],false));
end

% Optimization varaibles
for i = 1:dim_parameter
    OptimizationParameters.add(uitreenode('v0',strcat('parameter',num2str(i)),strcat(parameters_names{i},' = ',num2str(parameter(i))),[],false));
end
% LOOPS ON THE CONSTRAINTS TO GENERATE THE LEAVES

% Path constraints
for i = 1:dim_constraint
    PathConstraints.add(uitreenode('v0',strcat('path_constraint',num2str(i)),constraint_names{i},[],false));
end

% Dynamic constraints
for i = 1:dim_state
    DynamicConstraints.add(uitreenode('v0',strcat('dyn_constraint',num2str(i)),state_names{i},[],false));
end

% BOUNDARY CONDITIONS LEAVES
for i = 1:dim_boundarycond
    BoundaryConditions.add(uitreenode('v0',strcat('boundary_cond',num2str(i)),strcat(boundarycond_names{i},' = ',num2str(boundarycond(i))),[],false));
end

% LOOPS ON THE MULTIPLIERS TO GENERATE THE LEAVES

% Boundary conditions multipliers
for i = 1:dim_boundarycond
    BoundaryConditionsMult.add(uitreenode('v0',strcat('boundary_cond_mult',num2str(i)),strcat(boundarycond_names{i},' = ',num2str(boundarycond_mult(i))),[],false));
end

% Path constraits multipliers
for i = 1:dim_constraint
    PathConstraintsMult.add(uitreenode('v0',strcat('path_constraint_mult',num2str(i)),constraint_names{i},[],false));
end

% Adjoint states
for i = 1:dim_state
    AdjointStates.add(uitreenode('v0',strcat('adjoint_state',num2str(i)),state_names{i},[],false));
end

% TREE
f2=figure(2);
set(f2,'name','BOCOP - Display tree','numbertitle','off','pos',[800,300,400,400],'Menubar','none');
tree = uitree('v0', 'Root', root, 'SelectionChangeFcn',@callPlotFunctions);
tree.Position = [0,0,400,400];
tree.expand(root);

end

