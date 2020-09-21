% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% File: constraints_display.m
% Authors: Stephan Maindrault 

function constraints_display()

global dim_constraint dim_state
global path_constraint dynamic_constraint bpathconstraints bpathconstraints_type parameter
global stage time initial_time final_time free_time
global constraint_names state_names

dim_total = dim_constraint + dim_state;
if dim_total ~= 0
    
    col=0;
    row=0;
    
    f10 = figure(10);
    close(f10);
    f10 = figure(10);
    set(f10,'name','Constraints','numbertitle','off')
    
    % dimension display
    sq=sqrt(dim_total);
    sq=round(sq);
    if (sq*sq < dim_total)
        col=sq+1;
        row=sq;
    else
        col=sq;
        row=sq;
    end
    
    % path constraint variables display
    for i = 1 : dim_constraint
        subplot(row,col,i);
        hold on
        plot(stage,path_constraint(:,i),'linewidth',2);
         if strcmp(bpathconstraints_type,'lower')
             plot(stage,ones(1,length(stage))*bpathconstraints(i,1),'r','linewidth',2)
         elseif strcmp(bpathconstraints_type,'upper')
             plot(stage,ones(1,length(stage))*bpathconstraints(i,2),'r','linewidth',2)
         elseif strcmp(bpathconstraints_type,'both')
             plot(stage,ones(1,length(stage))*bpathconstraints(i,1),'r','linewidth',2)
             plot(stage,ones(1,length(stage))*bpathconstraints(i,2),'r','linewidth',2)
         elseif strcmp(bpathconstraints_type,'equal')
             plot(stage,ones(1,length(stage))*bpathconstraints(i,1),'r','linewidth',2)
         end
         
         if strcmp(free_time,'none')
            xlim([initial_time final_time])
         else
            xlim([initial_time parameter(end)])
         end
         
         t=title(strcat('pathConstraint.',num2str(i-1))); set(t, 'FontSize', 16); 
         hold off
    end
     for i = 1 : dim_state
        subplot(row,col,dim_constraint + i);
        plot(time(1:end-1),dynamic_constraint(:,i),'linewidth',2);
        
        if strcmp(free_time,'none')
            xlim([initial_time final_time])
        else
            xlim([initial_time parameter(end)])
        end
         
        t=title(strcat('dyncond.',num2str(i-1))); set(t, 'FontSize', 16);
        
     end
    end


