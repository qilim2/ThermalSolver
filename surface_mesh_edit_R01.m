function [MESHGRIDS_1,ELEMENTS_1,SURFACES_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = surface_mesh_edit_R01(ID,action,overwrite,input,MESHGRIDS_1,ELEMENTS_1,SURFACES_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2,THERMOPHYSICAL_1,THERMOOPTICAL_1)
% This function edits surface properties and can also be used to change the
% node IDs. This function does NOT allow for movement of the node. That can
% be performed with node_move. 
%
% Several different actions can be performed:
% action = 0: Delete the surface.
% action = 1: Change surface ID.
% action = 2: Change starting node ID (this will renumber all following
% nodes sequentially from this node). If any of the nodes already exists
% the function will abort regardless of the overwrite state.
% action = 3: Change the thermophysical property.
% action = 4: Change the topside thermooptical property.
% action = 5: Change the bottomside thermooptical property.
%
% If there is a nonzero value already at the position, overwrite must be
% set to 1 to confirm overwrite.
%
% Version 1.0 completed 12/14/2024


%% Check if surface exists
if ID > length(SURFACES_1)
    fprintf('ERROR (surface_mesh_edit): Surface does not exist.\n')
    return
end
if isempty(SURFACES_1{ID})
    fprintf('ERROR (surface_mesh_edit): Surface does not exist.\n')
    return
end

%% Get nodes
nodes = SURFACES_2{ID}.nodes;
nodecount = length(nodes);

%% Make edits
if action == 0 % Delete
    % Check overwrite
    if overwrite ~= 1
        fprintf('ERROR (surface_mesh_edit): Overwrite must be set to 1 to confirm deletion of surface.\n')
        return
    end

    % MESHGRIDS_1
    MESHGRIDS_1(nodes,:) = 0;

    % ELEMENTS_1
    [ind1,~] = find(ELEMENTS_1(:,1)==ID);
    ELEMENTS_1(ind1,:) = [];

    % SURFACES_1
    SURFACES_1{ID} = {};

    % SURFACES_2
    SURFACES_2{ID} = {};

    % CONDUCTANCES_1
    CONDUCTANCES_1(nodes,:) = 0;
    CONDUCTANCES_1(:,nodes) = 0;

    % CONDUCTANCES_2
    if ~isempty(CONDUCTANCES_2)
        for ii = 1:length(CONDUCTANCES_2)
            if ~isempty(CONDUCTANCES_2{ii})
                mat = CONDUCTANCES_2{ii};
                for jj = 1:length(nodes) % Iterate through list of nodes in the surface
                    [ind1,~] = find(mat(:,1:2)==nodes(jj)); % Find the rows the node exists in
                    mat(ind1,:) = []; % Delete the rows
                end
                CONDUCTANCES_2{ii} = mat;
            end
        end
    end

elseif action == 1 % Surface ID
    % Check if new surface ID exists
    if input <= length(SURFACES_1)
        if ~isempty(SURFACES_1{input}) && overwrite ~= 1
            fprintf('ERROR (surface_mesh_edit): New surface ID already exists. Set overwrite=1 to overwrite the existing node.\n')
            return
        end
    end

    % Update surface number in ELEMENTS_1
    [ind1,ind2] = find(ELEMENTS_1(:,1)==ID);
    ELEMENTS_1(ind1,ind2) = input;

    % Update SURFACES_1
    SURFACES_1{input} = SURFACES_1{ID};
    SURFACES_1{ID} = {};

    % Update SURFACES_2
    SURFACES_2{input} = SURFACES_2{ID};
    SURFACES_2{ID} = {};

elseif action == 2 % Start node ID
    % Get list of new nodes
    new_nodes = (input:input+nodecount-1)';

    % Check if any of new nodes already exists. Throw error and list
    % violating nodes if so.
    check = MESHGRIDS_1(new_nodes,1) ~= 0;
    if sum(check) > 0
        existing_nodes = new_nodes(check);
        fprintf('ERROR (surface_mesh_edit): The following nodes already exist. %s\n',mat2str(existing_nodes));
        return
    end

    % Assign nodes. Remove old nodes.
    MESHGRIDS_1(new_nodes,:) = MESHGRIDS_1(nodes,:);
    MESHGRIDS_1(nodes,:) = 0;
    SURFACES_2{ID}.nodes = new_nodes;

elseif action == 3 % Thermophysical property
    % Check if input is a string
    if ~isstring(input)
        fprintf('ERROR (surface_mesh_edit): Input must be a string.\n')
        return
    end

    % Check if input property exists
    indices = find(contains(THERMOPHYSICAL_1,input), 1); % Find indices in global storage matrix 
    if isempty(indices)
        fprintf('ERROR (surface_mesh_edit): Input property does not exist.\n')
        return
    end

    % Assign
    SURFACES_1{ID}.material = input;

elseif action == 4 % Topside thermooptical property
    % Check if input is a string
    if ~isstring(input)
        fprintf('ERROR (surface_mesh_edit): Input must be a string.\n')
        return
    end

    % Check if input property exists
    indices = find(contains(THERMOOPTICAL_1,input), 1); % Find indices in global storage matrix 
    if isempty(indices)
        fprintf('ERROR (surface_mesh_edit): Input property does not exist.\n')
        return
    end

    % Assign
    SURFACES_1{ID}.optical_top = input;

elseif action == 5 % Bottomside thermooptical property
    % Check if input is a string
    if ~isstring(input)
        fprintf('ERROR (surface_mesh_edit): Input must be a string.\n')
        return
    end

    % Check if input property exists
    indices = find(contains(THERMOOPTICAL_1,input), 1); % Find indices in global storage matrix 
    if isempty(indices)
        fprintf('ERROR (surface_mesh_edit): Input property does not exist.\n')
        return
    end

    % Assign
    SURFACES_1{ID}.optical_bot = input;
    
else
    fprintf('ERROR (node_edit): Invalid action. Use 1 to change the node ID, 2 to change the thermal mass, 3 to change the initial temperature.\n')
    return
end

end