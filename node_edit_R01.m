function [MESHGRIDS_1,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2] = node_edit_R01(ID,action,overwrite,input,MESHGRIDS_1,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2)
% This function edits node properties and can also be used to change the
% node ID. This function does NOT allow for movement of the node. That can
% be performed with node_move. 
%
% Several different actions can be performed:
% action = 0: Delete the node. Overwrite must be set to 1.
% action = 1: Change node ID
% action = 2: Change thermal mass
% action = 3: Change initial temperature
%
% If there is a nonzero value already at the position, overwrite must be
% set to 1 to confirm overwrite.
%
% Version 1.0 completed 12/14/2024


%% Check if node exists
check = find(MESHGRIDS_1(:,1)==ID, 1);
if isempty(check)
    fprintf('ERROR (node_edit): Node does not exist.\n')
    return
end

%% Make edits
if action == 0 % Delete
    % Check if node is part of an element (and also a surface)
    indices = find(ELEMENTS_1(:,2:5)==ID, 1);
    if ~isempty(indices)
        fprintf('ERROR (node_edit): Node is part of an element and cannot be deleted.\n')
        return
    end

    % Check overwrite
    if overwrite ~= 1
        fprintf('ERROR (node_edit): Overwrite must be set to 1 to confirm deletion of node.\n')
        return
    end

    % MESHGRIDS_1
    MESHGRIDS_1(ID,:) = 0;
    
    % CONDUCTANCES_1
    CONDUCTANCES_1(ID,:) = 0;
    CONDUCTANCES_1(:,ID) = 0;

    % CONDUCTANCES_2
    if ~isempty(CONDUCTANCES_2)
        for ii = 1:length(CONDUCTANCES_2)
            if ~isempty(CONDUCTANCES_2{ii})
                mat = CONDUCTANCES_2{ii};
                [ind1,~] = find(mat(:,1:2)==ID);
                mat(ind1,:) = [];
                CONDUCTANCES_2{ii} = mat;
            end
        end
    end
    
    
elseif action == 1 % Node ID
    % Check if node exists
    check = find(MESHGRIDS_1(:,1)==input,1);
    if ~isempty(check) && overwrite ~= 1
        fprintf('ERROR (node_edit): New node ID already exists. Set overwrite=1 to overwrite the existing node.\n')
        return
    end
    
    % MESHGRIDS_1
    temp = MESHGRIDS_1(ID,:);
    MESHGRIDS_1(ID,:) = 0;
    MESHGRIDS_1(input,:) = temp;
    
    % ELEMENTS_1
    [ind1,ind2] = find(ELEMENTS_1(:,2:5)==ID);
    ind2 = ind2+1; % Fix the index since it skips one
    ELEMENTS_1(ind1,ind2) = input;
    
    % CONDUCTANCES_1
    temp1 = CONDUCTANCES_1(ID,:);
    temp2 = CONDUCTANCES_1(:,ID);
    CONDUCTANCES_1(ID,:) = 0;
    CONDUCTANCES_1(:,ID) = 0;
    CONDUCTANCES_1(input,:) = temp1;
    CONDUCTANCES_1(:,input) = temp2;
    
    % CONDUCTANCES_2
    if ~isempty(CONDUCTANCES_2)
        for ii = 1:length(CONDUCTANCES_2)
            if ~isempty(CONDUCTANCES_2{ii})
                mat = CONDUCTANCES_2{ii};
                [ind1,ind2] = find(mat(:,1:2)==ID);
                mat(ind1,ind2) = input;
                CONDUCTANCES_2{ii} = mat;
            end
        end
    end

elseif action == 2 % Thermal mass
    % Check if value is nonzero
    check = MESHGRIDS_1(ID,2);
    if check ~= 0 && overwrite ~= 1
        fprintf('ERROR (node_edit): Thermal mass already exists. Set overwrite=1 to overwrite the existing value.\n')
        return
    end
    MESHGRIDS_1(ID,2) = input;
elseif action == 3 % Initial temperature
    % Check if value is nonzero
    check = MESHGRIDS_1(ID,3);
    if check ~= 0 && overwrite ~= 1
        fprintf('ERROR (node_edit): Initial temperature already exists. Set overwrite=1 to overwrite the existing value.\n')
        return
    end
    MESHGRIDS_1(ID,3) = input;
else
    fprintf('ERROR (node_edit): Invalid action. Use 1 to change the node ID, 2 to change the thermal mass, 3 to change the initial temperature.\n')
    return
end


end