% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% File: states_display.m
% Authors: Stephan Maindrault 

function states_display(graphnum)

global dim_state
global state parameter
global time initial_time final_time free_time
global state_names

if dim_state ~= 0
    col=0;
    row=0;
    
    if graphnum ~= 0
 
    f5 = figure(5);
    close(f5);
    f5 = figure(5);
    set(f5,'name','State variable','numbertitle','off')
    
    %state variable display

    plot(time,state(:,graphnum),'linewidth',2);
     if strcmp(free_time,'none')
         xlim([initial_time final_time])
     else
         xlim([initial_time parameter(end)])
     end
    t=title(state_names(graphnum)); set(t, 'FontSize', 16);

    else  
    f4 = figure(4);
    close(f4);
    f4 = figure(4);
    set(f4,'name','State variables','numbertitle','off')
    
    %display dimension
    sq=sqrt(dim_state);
    sq=round(sq);
    if (sq*sq < dim_state)
        col=sq+1;
        row=sq;
    else
        col=sq;
        row=sq;
    end
    
    %state variables display
     for i = 1 : dim_state
        subplot(row,col,i);
        plot(time,state(:,i),'linewidth',2);
        if strcmp(free_time,'none')
            xlim([initial_time final_time])
        else
            xlim([initial_time parameter(end)])
        end
        t=title(state_names(i)); set(t, 'FontSize', 16);
     end
    end
end
