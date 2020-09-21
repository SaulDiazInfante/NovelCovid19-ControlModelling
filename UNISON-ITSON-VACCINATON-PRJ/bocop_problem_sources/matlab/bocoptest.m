% test for bocop function

%TODO: setters for X0
X0 = [{'state.0','constant',0.45};
    {'state.1','constant',0.47};
    {'state.2','constant',0.48};
    {'control.0','constant',0.51};
    {'optimvars','none',[0.1]}]

%TODO: setters for options
options_bocop = struct( 't0','0','tf','1','freetf','final','disc_steps','100','disc_method','gauss')
                 

options_ipopt = [{'max_iter','1000' };
                 {'tol','1e-10'};
                 {'mu_strategy','adaptive'};
                 {'print_level','5'};  
                 {'file_print_level','5'};
                 {'output_file','result.out'}] 

% try to merge these ones ?
dims = [4 3 1 0 1 1 6]
bounds = [1 1;0 0;1 1; 1.01 2e20; 0 2e20; -2e20 2e20; 0 2e20; 0 1; 0 2e20; -2e20 0]
constants = [3.5; 310; 500; 1; 7; 0.6]

[time,stage,state,control,optimvars,output] = bocoplauncher(X0,options_bocop,options_ipopt,dims,bounds,constants);

output