% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% File: controls_display.m
% Authors: Stephan Maindrault 

function controls_display(graphnum)

global dim_control
global control_average parameter
global time initial_time final_time free_time
global control_names

if dim_control ~= 0    
    col=0;
    row=0;
    
    if graphnum ~=0
    
        f7 = figure(7);
        close(7);
        f7 = figure(7);
        set(f7,'name','Control variable','numbertitle','off')
    
        %control variable display
        plot(time(1:end-1),control_average(:,graphnum),'linewidth',2);
        if strcmp(free_time,'none')
            xlim([initial_time final_time])
        else
            xlim([initial_time parameter(end)])
        end
        t=title(control_names(graphnum)); set(t, 'FontSize', 16);
     
    else
        f6 = figure(6);
        close(f6);
        f6 = figure(6);
        set(f6,'name','Control variables','numbertitle','off')
        
        % display dimension
        sq=sqrt(dim_control);
        sq=round(sq);
        if (sq*sq < dim_control)
            col=sq+1;
            row=sq;
        else
            col=sq;
            row=sq;
        end
    
        %controls display
        for i = 1 : dim_control
            subplot(row,col,i);
            plot(time(1:end-1),control_average(:,i),'linewidth',2);
            if strcmp(free_time,'none')
                xlim([initial_time final_time])
            else
                xlim([initial_time parameter(end)])
            end
            t=title(control_names(i)); set(t, 'FontSize', 16)
        end
    end

end
