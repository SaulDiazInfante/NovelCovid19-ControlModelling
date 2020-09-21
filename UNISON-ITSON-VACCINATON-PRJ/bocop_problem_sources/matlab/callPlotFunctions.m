% Copyright (C) 2013 INRIA.
% All Rights Reserved.
% This code is published under the Eclipse Public License
% File: callPlotFunctions.m
% Authors: Daphne Giorgi

% This function calls the different plot functions, depending on the
% position in the tree. It makes a simple plot if we are on a leaf, and
% multiplots if we are on one of the nodes.

function callPlotFunctions(tree, value)

% We get the selected node
selectedNodes = tree.getSelectedNodes;
node = selectedNodes(1);
nodeName = node.getValue;

% If we are on the root we do nothing
if ~strcmp(nodeName,'PlottingTree')
    
    parentNode = node.getParent;
    parentName = parentNode.getValue;
    
    % If we are in a main category we plot all the variables under that
    % category
    if strcmp(parentName,'PlottingTree')
        plotAllVariables(nodeName);
    % If we are in a subcategory, we plot the specific variables of the
    % subcategory
    else
        subParentNode = parentNode.getParent;
        subParentName = subParentNode.getValue;
        
        if strcmp(subParentName,'PlottingTree')
            plotCategory(nodeName, parentName);
        % Else we are in a leaf and we plot one variable
        else    
            plotSimple(nodeName, parentName, subParentName);          
        end          
    end
end

end


% This function makes a multiplot with all the variables in a main category
function plotAllVariables(nodeName)

% If the main category is the variables category
% we display all the variables
if strcmp(nodeName, 'Variables')
    variables_display;
% else if the main category is the constraints category
% we display all the constraints
elseif strcmp(nodeName,'Constraints')
    constraints_display;
% else if the main category is the multipliers category
% we display all the multipliers
elseif strcmp(nodeName,'Multipliers')
    multipliers_display;
end

end


% This function makes a multiplot with the variables of a specific category
function plotCategory(nodeName, parentName)

% If the main category is the variables category
if strcmp(parentName, 'Variables')
    % If the subcategory is the states category
    % we display al the states
    if strcmp(nodeName,'StateVariables')
        states_display(0);
    % else if the subcategory is the controls category
    % we display al the controls
    elseif strcmp(nodeName,'ControlVariables')
        controls_display(0);
    % else if the subcategory is the algebraic variables category
    % we display al the algebraic variables
    elseif strcmp(nodeName,'AlgebraicVariables')
        algebraics_display(0);
    end   
% else if the main category is the constraints category  
elseif strcmp(parentName,'Constraints')
    % If the subcategory is the path constraints category
    % we display al the path constraints
    if strcmp(nodeName,'PathConstraints')
        pathconstraints_display(0);
    % else if the subcategory is the dynamic constraints category
    % we display all the dynamic constraints 
    elseif strcmp(nodeName,'DynamicConstraints')
         dynamics_const_display(0);
    end
% else if the main category is the multipliers category  
elseif strcmp(parentName,'Multipliers')
    % If the subcategory is the path constraints multipliers category
    % we display al the path constraints multipliers
    if strcmp(nodeName,'PathConstraintsMult')
        pathmultipliers_display(0);
    % else if  the subcategory is the adjoint states category
    % we display al the adjoint states
    elseif strcmp(nodeName,'AdjointStates')
        adjointstates_display(0);
    end
end

end

% This function make a simple plot of one variable (a leaf)
function plotSimple(nodeName, parentName, subParentName)

% If the main category is the variables category
if strcmp(subParentName, 'Variables')
    % If the subcategory is the states category
    % we take the index of the variable from its name
    % and we display it
    if strcmp(parentName,'StateVariables')
        index = str2num(strrep(nodeName,'state',''));
        states_display(index);
    % else if the subcategory is the controls category
    % we take the index of the variable from its name
    % and we display it
    elseif strcmp(parentName,'ControlVariables')
        index = str2num(strrep(nodeName,'control',''));
        controls_display(index);
    % else if the subcategory is the algebraics category
    % we take the index of the variable from its name
    % and we display it
    elseif strcmp(parentName,'AlgebraicVariables')
        index= str2num(strrep(nodeName,'algebraic',''));
        algebraics_display(index);
    end  
% else if the main category is the constraints category  
elseif strcmp(subParentName,'Constraints')
    % If the subcategory is the path constraints category
    % we take the index of the variable from its name
    % and we display it
    if strcmp(parentName,'PathConstraints')
        index = str2num(strrep(nodeName,'path_constraint',''));
        pathconstraints_display(index);
    % else if the subcategory is the dynamic constraints category
    % we take the index of the variable from its name
    % and we display it
    elseif strcmp(parentName,'DynamicConstraints')
        index  = str2num(strrep(nodeName,'dyn_constraint',''));
        dynamics_const_display(index);
    end
% else if the main category is the multipliers category
elseif strcmp(subParentName,'Multipliers')
    % If the subcategory is the path constraints multipliers category
    % we take the index of the variable from its name
    % and we display it
    if strcmp(parentName,'PathConstraintsMult')
        index = str2num(strrep(nodeName,'path_constraint_mult',''));
        pathmultipliers_display(index);
    % else if the subcategory is the adjoint states category
    % we take the index of the variable from its name
    % and we display it
    elseif strcmp(parentName,'AdjointStates')
        index = str2num(strrep(nodeName,'adjoint_state',''));
        adjointstates_display(index);
    end
end

end
