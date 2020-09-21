% Simplified matlab interface for Bocop
% File: bocop.m
% Authors: Daphne Giorgi, Stephan Maindrault, Pierre Martinon
% Inria Saclay, 2013-2015

clear all
close all

%global variables
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

%define layout for user interface

%main frame
f = figure('Name','BOCOP - Solution loader','Position',[400,400,330,300],'Menubar','none');
set(f,'numbertitle','off')

%title
htext = uicontrol('Style','text','String','Load solution file','Fontsize',18,'Position',[40,230,250,40]);

%textfield for solution file name
lastopen = '../../examples/goddard/problem.sol';
hpath = uicontrol('Style','edit','String',lastopen,'Fontsize',16,'Position',[15,150,245,40]);
filename = get(hpath,'String');

%browser button
hbrowse = uicontrol('Style','pushbutton','String','...','Fontsize',17,'Position',[265,150,50,40],'Callback','launch');

%visualization tree button
hvisualisation = uicontrol('Style','pushbutton','String','Visualization','Fontsize',18,'Position',[15,80,300,50],'Callback','visualization_tree');

% save solution button
hsave = uicontrol('Style','pushbutton','String','Save variables','Fontsize',18,'Position',[15,20,300,50],'Callback','save_solution');


set(f,'Visible','on')


 