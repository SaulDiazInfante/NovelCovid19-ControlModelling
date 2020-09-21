% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% File: adjointstates_display.m
% Authors: Stephan Maindrault 

function adjointstates_display(graphnum)

global dim_state
global adjoint_state parameter
global time initial_time final_time free_time
global state_names


if dim_state ~=0
    
    row=0;
    col=0;
    
    if graphnum ~= 0
         
    f19 = figure(19);
    close(f19);
    f19 = figure(19);
    set(f19,'name','Adjoint state','numbertitle','off')
 
    %adjoint state variable display
    
    plot(time(1:end-1),adjoint_state(:,graphnum),'linewidth',2);
    if strcmp(free_time,'none')
           xlim([initial_time final_time])
    else
           xlim([initial_time parameter(end)])
    end
    
    t=title(strcat('adjointState.',num2str(graphnum-1))); set(t, 'FontSize', 16);
        
    else    
    f18 = figure(18);
    close(f18);
    f18 = figure(18);
    set(f18,'name','Adjoint states','numbertitle','off')
    
    %dimension display
    sq=sqrt(dim_state);
    sq=round(sq);
    if (sq*sq < dim_state)
        col=sq+1;
        row=sq;
    else
        col=sq;
        row=sq;
    end
    
    %adjoint state variables display
     for i = 1 : dim_state
        subplot(row,col,i);
        plot(time(1:end-1),adjoint_state(:,i),'linewidth',2);
         if strcmp(free_time,'none')
            xlim([initial_time final_time])
        else
            xlim([initial_time parameter(end)])
        end
        t=title(strcat('adjointState.',num2str(i-1))); set(t, 'FontSize', 16);
     end
    end
end
