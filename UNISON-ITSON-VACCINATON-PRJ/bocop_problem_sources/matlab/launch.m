% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% File: launch.m
% Authors: Stephan Maindrault, Pierre Martinon

function launch()

%lauch file browser action
global hpath filename

%uigetfile2 remembers last open file
[file path]= uigetfile2('*.sol','Select any file');
filename=strcat(path,file);
set(hpath,'String',filename)

%load solution for visualization or use with custom scripts
visualization_tree()

end
