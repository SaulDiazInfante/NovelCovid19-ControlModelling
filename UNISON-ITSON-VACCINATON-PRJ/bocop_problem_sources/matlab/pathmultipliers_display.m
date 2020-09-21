% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% File: pathmultipliers_display.m
% Authors: Stephan Maindrault 

function pathmultipliers_display(graphnum)

global dim_constraint
global pathconstraint_mult parameter
global stage final_time initial_time free_time
global path_name

if dim_constraint ~= 0
    row=0;
    col=0;
    
    if graphnum ~= 0
        
        f17 = figure(17);
        close(f17);
        f17 = figure(17);
        set(f17,'name','Path constraint multiplier','numbertitle','off')
    
        %dynamic constraint variable display
   
        plot(stage,pathconstraint_mult(:,graphnum),'linewidth',2);
        if strcmp(free_time,'none')
            xlim([initial_time final_time])
        else
            xlim([initial_time parameter(end)])
        end
        t=title(strcat('pathConstMult.',num2str(graphnum-1)));set(t, 'FontSize', 16); 
       
        
    else    
        f16 = figure(16);
        close(f16);
        f16 = figure(16);
        set(f16,'name','Path constraint multipliers','numbertitle','off')
    
        %dimension display
        sq=sqrt(dim_constraint);
        sq=round(sq);
        if (sq*sq < dim_constraint)
            col=sq+1;
            row=sq;
        else
            col=sq;
            row=sq;
        end
    
        %adjoint state variables display
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
    end
end
