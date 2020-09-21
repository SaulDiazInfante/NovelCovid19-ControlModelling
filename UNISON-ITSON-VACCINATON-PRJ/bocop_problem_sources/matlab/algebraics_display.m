% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% File: algebraics_display.m
% Authors: Stephan Maindrault 

function algebraics_display(graphnum)

global dim_algebraic
global algebraic
global stage
global algebraic_names


if dim_algebraic ~= 0
    col=0;
    row=0;
    
    if graphnum ~= 0
    
    f9 = figure(9);
    close(f9);    
    f9 = figure(9);
    set(f9,'name','Algebraic variable','numbertitle','off')
    
    %algebraic variable display
    
    plot(stage,algebraic(:,graphnum),'linewidth',2);
    t=title(algebraic_names(graphnum)); set(t, 'FontSize', 16);
       
    else
    f8 = figure(8);
    close(f8);
    f8 = figure(8);
    set(f8,'name','Algebraic variables','numbertitle','off')
    
    %dimension display
    sq = sqrt(dim_algebraic);
    sq = round(sq);
    if (sq*sq < dim_algebraic)
        col = sq+1;
        row = sq;
    else
        col = sq;
        row = sq;
    end
    
    %algebraics display
     for i = 1 : dim_algebraic
        subplot(row,col,i);
        plot(stage,algebraic(:,i),'linewidth',2);
        t=title(algebraic_names(i)); set(t, 'FontSize', 16);
     end
    end
end

