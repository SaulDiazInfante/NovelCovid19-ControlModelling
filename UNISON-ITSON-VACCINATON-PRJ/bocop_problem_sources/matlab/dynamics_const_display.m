% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% File: dynamics_const_display.m
% Authors: Stephan Maindrault 

function dynamics_const_display(graphnum)

global dim_state
global dynamic_constraint parameter
global time initial_time final_time free_time
global state_names

if dim_state ~= 0
    
    col=0;
    row=0;

    if graphnum ~=0
        
    f14 = figure(14);
    close(f14)
    f14 = figure(14);
    set(f14,'name','Dynamic constraint','numbertitle','off')
    
    %dynamic constraint variable display
    for i = 1 : dim_state
        plot(time(1:end-1),dynamic_constraint(:,graphnum),'linewidth',2);
        if strcmp(free_time,'none')
            xlim([initial_time final_time])
         else
            xlim([initial_time parameter(end)])
        end
        t=title(strcat('dyncond.',num2str(graphnum-1))); set(t, 'FontSize', 16);
    end
    
    else
    f13 = figure(13);    
    close(f13)
    f13 = figure(13);
    set(f13,'name','Dynamic constraints','numbertitle','off')
    
    % dimension display
    sq=sqrt(dim_state);
    sq=round(sq);
    if (sq*sq < dim_state)
        col=sq+1;
        row=sq;
    else
        col=sq;
        row=sq;
    end
    
    %dynamic constraint variables display
     for i = 1 : dim_state
        subplot(row,col,i);
        plot(time(1:end-1),dynamic_constraint(:,i),'linewidth',2);
        if strcmp(free_time,'none')
            xlim([initial_time final_time])
        else
            xlim([initial_time parameter(end)])
        end
        t=title(strcat('dyncond.',num2str(i-1))); set(t, 'FontSize', 16);
     end
    end
end
