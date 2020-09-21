% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% File: variables_display.m
% Authors: Stephan Maindrault 

function variables_display()

global dim_state dim_control dim_algebraic
global state control_average algebraic parameter
global time stage initial_time final_time free_time
global  state_names control_names algebraic_names

dim_total=dim_state + dim_control +dim_algebraic;

if dim_total ~= 0
    row=0;
    col=0;
    
    f3 = figure(3);
    close(f3);
    f3 = figure(3);
    set(f3,'name','Variables','numbertitle','off')
    
    %display dimension
    
    sq=sqrt(dim_total);
    sq=round(sq);
    if (sq*sq < dim_total)
        col=sq+1;
        row=sq;
    else
        col=sq;
        row=sq;
    end
    
    %states display
    for i = 1 : dim_state
        subplot(row,col,i);
        plot(time,state(:,i),'linewidth',2);
        if strcmp(free_time,'none')
            xlim([initial_time final_time])
        else
            xlim([initial_time parameter(end)])
        end
        t=title(state_names(i)); set(t, 'FontSize', 15);
    end
    
    %controls display
    for i = 1 : dim_control
        subplot(row,col,dim_state + i);
        plot(time(1:end-1),control_average(:,i),'linewidth',2);
        if strcmp(free_time,'none')
            xlim([initial_time final_time])
         else
            xlim([initial_time parameter(end)])
         end
        t=title(control_names(i)); set(t, 'FontSize', 15);
    end
    
    %algebraics display
    for i = 1 : dim_algebraic
        subplot(row,col,dim_state + dim_control + i);
        plot(stage,algebraic(:,i),'linewidth',2);
        if strcmp(free_time,'none')
            xlim([initial_time final_time])
        else
            xlim([initial_time parameter(end)])
        end
        t=title(algebraic_names(i)); set(t, 'FontSize', 15);
    end
end
