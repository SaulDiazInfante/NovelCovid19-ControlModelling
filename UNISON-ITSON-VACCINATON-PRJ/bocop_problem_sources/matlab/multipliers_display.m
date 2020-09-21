% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% File: multipliers_display.m
% Authors: Stephan Maindrault 

function multipliers_display()

global dim_constraint dim_state
global pathconstraint_mult adjoint_state parameter
global stage time initial_time final_time free_time
global constraint_names state_names


dim_total=dim_constraint + dim_state;

if dim_total ~= 0
    
    row=0;
    col=0;
    
    f15 = figure(15);
    close(f15);
    f15 = figure(15);
    set(f15,'name','Multipliers','numbertitle','off')
    
    %display dimention
    sq=sqrt(dim_total);
    sq=round(sq);
    if (sq*sq < dim_total)
        col=sq+1;
        row=sq;
    else
        col=sq;
        row=sq;
    end
    
    %pathconstraint multiplier variables display
    for i = 1 : dim_constraint
        subplot(row,col,i);
        plot(stage,pathconstraint_mult(:,i),'linewidth',2);
        if strcmp(free_time,'none')
            xlim([initial_time final_time])
         else
            xlim([initial_time parameter(end)])
         end
        t=title(strcat('pathConstMult.',num2str(i-1))); set(t, 'FontSize', 16);
    end
    
    %adjoint state variables display
    for i = 1 : dim_state
        subplot(row,col,dim_constraint + i);
        plot(time(1:end-1),adjoint_state(:,i),'linewidth',2);
        if strcmp(free_time,'none')
            xlim([initial_time final_time])
         else
            xlim([initial_time parameter(end)])
         end
        t=title(strcat('adjointState.',num2str(i-1))); set(t, 'FontSize', 16);
    end
end
