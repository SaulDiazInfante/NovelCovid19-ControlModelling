% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% File: pathconstraints_display.m
% Authors: Stephan Maindrault 

function pathconstraints_display(graphnum)

global dim_constraint
global path_constraint bpathconstraints bpathconstraints_type parameter 
global stage initial_time final_time free_time
global constraint_names

if dim_constraint ~= 0
    
    col=0;
    row=0;
      
    if graphnum ~= 0
    
    f12 = figure(12);
    close(f12);
    f12 = figure(12);
    set(f12,'name','Path constraint variable','numbertitle','off')
    
    %Path constraint display
    hold on
    
    plot(stage,path_constraint(:,graphnum),'linewidth',2);
    if strcmp(bpathconstraints_type,'lower')
            plot(stage,ones(1,length(stage))*bpathconstraints(graphnum,1),'r','linewidth',2)
    elseif strcmp(bpathconstraints_type,'upper')
            plot(stage,ones(1,length(stage))*bpathconstraints(graphnum,2),'r','linewidth',2)
    elseif strcmp(bpathconstraints_type,'both')
            plot(stage,ones(1,length(stage))*bpathconstraints(graphnum,1),'r','linewidth',2)
            plot(stage,ones(1,length(stage))*bpathconstraints(graphnum,2),'r','linewidth',2)
    elseif strcmp(bpathconstraints_type,'equal')
            plot(stage,ones(1,length(stage))*bpathconstraints(graphnum,1),'r','linewidth',2)
    end
    
    if strcmp(free_time,'none')
        xlim([initial_time final_time])
    else
        xlim([initial_time parameter(end)])
    end
    
    t=title(strcat('pathConstraint.',num2str(graphnum-1))); set(t, 'FontSize', 16);
    hold off
    
    else
    f11 = figure(11);    
    close(f11);
    f11 = figure(11);
    set(f11,'name','Path constraint variables','numbertitle','off')
   
    % dimension display
    sq=sqrt(dim_constraint);
    sq=round(sq);
    if (sq*sq < dim_constraint)
        col=sq+1;
        row=sq;
    else
        col=sq;
        row=sq;
    end
    
    %path constraint variables display
     for i = 1 : dim_constraint
        subplot(row,col,i);
        plot(stage,path_constraint(:,i),'linewidth',2);
        hold on
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
    end
end
