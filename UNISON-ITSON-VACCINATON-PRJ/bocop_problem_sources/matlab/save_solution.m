% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% This code is published under the Eclipse Public License
% File: save_solution.m
% Authors: Daphne Giorgi

% This function save the main solution vectors (time and variables) in a matlab format for a later load in matlab
function save_solution(times, state, average_control, algebraic, parameter, pathcond)

global time state control_average algebraic parameter path_constraint  filename

readsolfile(filename);

uisave({'time','state','control_average','algebraic','parameter','path_constraint'},'matlab_solution');

end