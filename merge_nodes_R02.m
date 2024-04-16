function [MESHGRIDS_1,ELEMENTS_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2] = merge_nodes_R02(source_input,source_type,target_input,target_type,range_tolerance,keep_switch,MESHGRIDS_1,SURFACES_2,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2)
% This function takes in a two sets of inputs: sources and targets. Source
% nodes are checked against the target nodes using the range_tolerance, 
% merging the target with the source or vice-versa depending on the 
% keep_switch (0 = keep lower indexed node, 1 = keep higher indexed node).
% Source and target types are needed to identify what the sources and the
% targets are.
% 0 = all currently active nodes
% 1 = nodes
% 2 = elements
% 30 = surface
% 31 = surface south edge (input is the surface ID)
% 32 = surface east edge (input is the surface ID)
% 33 = surface north edge (input is the surface ID)
% 34 = surface west edge (input is the surface ID)

% Version 2.0 completed 12/7/2023
% Version 2.1 completed 12/11/2023

%% Pre-process
% Identify existing nodes and surfaces
nodes = find(MESHGRIDS_1(:,1)~=0);

%% Find source nodes
if source_type == 0 % All nodes currently active
    source_nodes = nodes;
elseif source_type == 1 % Nodes
    % Check if values exist in input
    for ii = 1:length(source_input)
        existcheck = MESHGRIDS_1(:,1) == source_input(ii);
        if sum(existcheck) < 1
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
    end
    source_nodes = source_input;
elseif source_type == 2 % Elements
    source_nodes = [];
    for ii = 1:length(source_input) % Grab nodes from each element and append to list
        existcheck = ELEMENTS_1(:,1) == source_input(ii);
        if sum(existcheck) < 1
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
        source_nodes = [source_nodes, ELEMENTS_1(source_input(ii),2:5)];
    end
    source_nodes = unique(source_nodes)'; % Remove duplicates and sort
elseif source_type == 30 % Entire surface
    source_nodes = [];
    for ii = 1:length(source_input)
        if isempty(SURFACES_2{source_input(ii)})
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
        source_nodes = [source_nodes; SURFACES_2{source_input(ii)}.nodes];
    end
    source_nodes = unique(source_nodes); % Remove duplicates and sort
elseif source_type == 31 % South edge of surface
    source_nodes = [];
    n_x1 = zeros(length(source_input),1);
    for ii = 1:length(source_input)
        if isempty(SURFACES_2{source_input(ii)})
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
        current_source = source_input(ii);
        surf_nodes = SURFACES_2{current_source}.nodes;
        n_x1(ii) = SURFACES_2{current_source}.n_x1;
        source_nodes = [source_nodes;surf_nodes(1:n_x1)];
    end
    source_nodes = unique(source_nodes); % Remove duplicates and sort
elseif source_type == 32 % East edge of surface
    source_nodes = [];
    n_x2 = zeros(length(source_input),1);
    for ii = 1:length(source_input)
        if isempty(SURFACES_2{source_input(ii)})
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
        current_source = source_input(ii);
        surf_nodes = SURFACES_2{current_source}.nodes;
        n_x1 = SURFACES_2{current_source}.n_x1;
        n_x2(ii) = SURFACES_2{current_source}.n_x2;
        for jj = 1:n_x2
            source_nodes = [source_nodes;surf_nodes(jj*n_x1)];
        end
    end
    source_nodes = unique(source_nodes); % Remove duplicates and sort
elseif source_type == 33 % North edge of surface
    source_nodes = [];
    n_x1 = zeros(length(source_input),1);
    for ii = 1:length(source_input)
        if isempty(SURFACES_2{source_input(ii)})
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
        current_source = source_input(ii);
        surf_nodes = SURFACES_2{current_source}.nodes;
        n_x1(ii) = SURFACES_2{current_source}.n_x1;
        n_x2 = SURFACES_2{current_source}.n_x2;
        source_nodes = [source_nodes;surf_nodes(n_x1*n_x2-n_x1+1:end)];
    end
    source_nodes = unique(source_nodes); % Remove duplicates and sort
elseif source_type == 34 % West edge of surface
    source_nodes = [];
    n_x2 = zeros(length(source_input),1);
    for ii = 1:length(source_input)
        if isempty(SURFACES_2{source_input(ii)})
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
        current_source = source_input(ii);
        surf_nodes = SURFACES_2{current_source}.nodes;
        n_x1 = SURFACES_2{current_source}.n_x1;
        n_x2(ii) = SURFACES_2{current_source}.n_x2;
        for jj = 1:n_x2
            source_nodes = [source_nodes;surf_nodes(jj*n_x1-n_x1+1)];
        end
    end
else
    fprintf('ERROR (merge_nodes): Input_type invalid. Valid input_types are 0, 1, 2, 30, 31, 32, 33, 34.\n')
    return
end
%disp(source_nodes)

%% Find target nodes
if target_type == 0 % All nodes currently active
    target_nodes = nodes;
elseif target_type == 1 % Nodes
    for ii = 1:length(source_input)
        existcheck = MESHGRIDS_1(:,1) == target_input(ii);
        if sum(existcheck) < 1
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
    end
    target_nodes = target_input;
elseif target_type == 2 % Elements
    target_nodes = [];
    for ii = 1:length(target_input) % Grab nodes from each element and append to list
        existcheck = ELEMENTS_1(:,1) == target_input(ii);
        if sum(existcheck) < 1
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
        target_nodes = [target_nodes, ELEMENTS_1(target_input(ii),2:5)];
    end
    target_nodes = unique(target_nodes)'; % Remove duplicates and sort
elseif target_type == 30 % Entire surface
    target_nodes = [];
    for ii = 1:length(target_input)
        if isempty(SURFACES_2{target_input(ii)})
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
        target_nodes = [target_nodes; SURFACES_2{target_input(ii)}.nodes];
    end
    target_nodes = unique(target_nodes); % Remove duplicates and sort
elseif target_type == 31 % South edge of surface
    target_nodes = [];
    n_x1 = zeros(length(target_input),1);
    for ii = 1:length(target_input)
        if isempty(SURFACES_2{target_input(ii)})
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
        current_target = target_input(ii);
        surf_nodes = SURFACES_2{current_target}.nodes;
        n_x1(ii) = SURFACES_2{current_target}.n_x1;
        target_nodes = [target_nodes;surf_nodes(1:n_x1)];
    end
    target_nodes = unique(target_nodes); % Remove duplicates and sort
elseif target_type == 32 % East edge of surface
    target_nodes = [];
    n_x2 = zeros(length(target_input),1);
    for ii = 1:length(target_input)
        if isempty(SURFACES_2{target_input(ii)})
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
        current_target = target_input(ii);
        surf_nodes = SURFACES_2{current_target}.nodes;
        n_x1 = SURFACES_2{current_target}.n_x1;
        n_x2(ii) = SURFACES_2{current_target}.n_x2;
        for jj = 1:n_x2
            target_nodes = [target_nodes;surf_nodes(jj*n_x1)];
        end
    end
    target_nodes = unique(target_nodes); % Remove duplicates and sort
elseif target_type == 33 % North edge of surface
    target_nodes = [];
    n_x1 = zeros(length(target_input),1);
    for ii = 1:length(target_input)
        if isempty(SURFACES_2{target_input(ii)})
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
        current_target = target_input(ii);
        surf_nodes = SURFACES_2{current_target}.nodes;
        n_x1(ii) = SURFACES_2{current_target}.n_x1;
        n_x2 = SURFACES_2{current_target}.n_x2;
        target_nodes = [target_nodes;surf_nodes(n_x1*n_x2-n_x1+1:end)];
    end
    target_nodes = unique(target_nodes); % Remove duplicates and sort
elseif target_type == 34 % West edge of surface
    target_nodes = [];
    n_x2 = zeros(length(target_input),1);
    for ii = 1:length(target_input)
        if isempty(SURFACES_2{target_input(ii)})
            fprintf('ERROR (merge_nodes): Value in input does not exist in the mesh.\n')
            return
        end
        current_target = target_input(ii);
        surf_nodes = SURFACES_2{current_target}.nodes;
        n_x1 = SURFACES_2{current_target}.n_x1;
        n_x2(ii) = SURFACES_2{current_target}.n_x2;
        for jj = 1:n_x2
            target_nodes = [target_nodes;surf_nodes(jj*n_x1-n_x1+1)];
        end
    end
else 
    fprintf('ERROR (merge_nodes): Input_type invalid. Valid input_types are 0, 1, 2, 30, 31, 32, 33, 34.\n')
    return
end
%disp(target_nodes)

%% Find nodes to keep and nodes to replace
source_count = length(source_nodes); % Get count of source nodes
source_dup = source_nodes; % Save original list of source nodes for later comparison
target_count = length(target_nodes); % Get count of target nodes
target_dup = target_nodes; % Save original list of target nodes for later comparison

% Keeping the higher index
if keep_switch == 1
    for ii = source_count:-1:1 % Iterate through source nodes
        current_source = source_nodes(ii);
        
        % Check if the range_tolerance is greater than 1/10th of the source node's parent element internodal lengths
        [ind1,ind2] = find(ELEMENTS_1(:,2:5)==current_source);
        ind2 = ind2 + 1; % Note: index is relative to ELEMENTS_1
        if ind2 < 5
            checknode = ind2+1; % Check against next node in element
        else
            checknode = ind2-1; % Check against previous node in element
        end
        checknode = ELEMENTS_1(ind1,checknode); % Get actual node number of the checknode
        temp = norm(MESHGRIDS_1(current_source,4:6)-MESHGRIDS_1(checknode,4:6));
        if range_tolerance > 0.1*temp
            fprintf('ERROR (merge_nodes): Range tolerance is too large and may find nodes within a source node''s surface. Please input a range tolerance less than %f.\n',0.1*temp);
            return
        end
        
        % Iterate through target nodes for the current source node
        for jj = target_count:-1:1
            current_target_node = target_nodes(jj);
            % Check if target node is within range
            nodal_length = norm(MESHGRIDS_1(current_source,4:6)-MESHGRIDS_1(current_target_node,4:6));
            if nodal_length <= range_tolerance
                if current_target_node > current_source
                    % Update thermal mass
                    MESHGRIDS_1(current_target_node,2) = MESHGRIDS_1(current_target_node,2) + MESHGRIDS_1(current_source,2);

                    % Keep the target node, replacing the source node in
                    % the list of source nodes and the list of target nodes
                    source_nodes(source_nodes == current_source) = current_target_node;
                    target_nodes(target_nodes == current_source) = current_target_node;
                elseif current_target_node < current_source
                    % Update thermal mass
                    MESHGRIDS_1(current_source,2) = MESHGRIDS_1(current_source,2) + MESHGRIDS_1(current_target_node,2);

                    % Keep the source node, replacing the target node in
                    % the list of target nodes and the list of source nodes
                    source_nodes(source_nodes == current_target_node) = current_source;
                    target_nodes(target_nodes == current_target_node) = current_source;
                end
            end
        end
    end
elseif keep_switch == 0
    for ii = 1:source_count % Iterate through source nodes
        current_source = source_nodes(ii);
        
        % Check if the range_tolerance is greater than 1/10th of the source node's parent element internodal lengths
        [ind1,ind2] = find(ELEMENTS_1(:,2:5)==current_source);
        ind2 = ind2 + 1; % Note: index is relative to ELEMENTS_1
        if ind2 < 5
            checknode = ind2+1; % Check against next node in element
        else
            checknode = ind2-1; % Check against previous node in element
        end
        checknode = ELEMENTS_1(ind1,checknode); % Get actual node number of the checknode
        temp = norm(MESHGRIDS_1(current_source,4:6)-MESHGRIDS_1(checknode,4:6));
        if range_tolerance > 0.1*temp
            fprintf('ERROR (merge_nodes): Range tolerance is too large and may find nodes within a source node''s surface. Please input a range tolerance less than %f.\n',0.1*temp);
            fprintf('If a connection between distant nodes/elements/surfaces are desired, consider using a Contactor or a Conductor.\n')
            return
        end
        
        % Iterate through target nodes for the current source node
        for jj = 1:target_count
            current_target_node = target_nodes(jj);
            % Check if target node is within range
            nodal_length = norm(MESHGRIDS_1(current_source,4:6)-MESHGRIDS_1(current_target_node,4:6));
            if nodal_length <= range_tolerance
                if current_target_node < current_source
                    % Update thermal mass
                    MESHGRIDS_1(current_target_node,2) = MESHGRIDS_1(current_target_node,2) + MESHGRIDS_1(current_source,2);

                    % Keep the target node, replacing the source node in
                    % the list of source nodes and the list of target nodes
                    source_nodes(source_nodes == current_source) = current_target_node;
                    target_nodes(target_nodes == current_source) = current_target_node;
                elseif current_target_node > current_source
                    % Update thermal mass
                    MESHGRIDS_1(current_source,2) = MESHGRIDS_1(current_source,2) + MESHGRIDS_1(current_target_node,2);

                    % Keep the source node, replacing the target node in
                    % the list of target nodes and the list of source nodes
                    source_nodes(source_nodes == current_target_node) = current_source;
                    target_nodes(target_nodes == current_target_node) = current_source;
                end
            end
        end
    end
else
    fprintf('ERROR (merge_nodes): Keep_switch invalid. Valid keep_switch values are 1 and 0.\n')
    return
end
%source_replace = [source_nodes,source_dup];
%target_replace = [target_nodes,target_dup];
%disp(source_replace)
%disp(target_replace)


%% Full replacement
% Replace nodes in MESHGRIDS_1, ELEMENTS_1, SURFACES_2, CONDUCTANCES_1
for ii = 1:source_count
    if source_nodes(ii) ~= source_dup(ii)
        % Find source_dup(ii) in MESHGRIDS_1, ELEMENTS_1, SURFACES_2,
        % CONDUCTANCES_1 and make the replacement or removal. Ignore if
        % already done.

        % Find source_dup(ii) in MESHGRIDS_1 and remove
        MESHGRIDS_1(MESHGRIDS_1(:,1)==source_dup(ii),:) = 0;
        
        % Find source_dup(ii) in ELEMENTS_1 and replace
        elem_surfs = ELEMENTS_1(:,1); % Keep first column since this may be overwritten in the next step
        change_indices = ELEMENTS_1==source_dup(ii); % This can be used to identify which elements have been changed, and therefore which surfaces have been changed
        ELEMENTS_1(change_indices) = source_nodes(ii);
        ELEMENTS_1(:,1) = elem_surfs; % Fix the surface IDs in the ELEMENTS_1 matrix
        
        % Find source_dup(ii) in SURFACES_1 and replace
        change_rows = sum(change_indices,2) > 0; % Find which elements have been changed
        change_surfs = elem_surfs(change_rows); % Extract the surfaces that have been changed
        %disp(change_surfs)
        for jj = 1:length(change_surfs) % Iterate through each changed surface and change nodes
            surf_nodes = SURFACES_2{change_surfs(jj)}.nodes;
            surf_nodes(surf_nodes==source_dup(ii))=source_nodes(ii);
            SURFACES_2{change_surfs(jj)}.nodes = surf_nodes;
        end

        % Find source_dup(ii) in CONDUCTANCES_1 and remove and replace.
        % Note that replacement only happens if the conductance between the
        % two nodes is zero.
        zero_locs_horz = CONDUCTANCES_1(source_nodes(ii),:) == 0;
        zero_locs_vert = CONDUCTANCES_1(:,source_nodes(ii)) == 0;
        CONDUCTANCES_1(source_nodes(ii),:) = CONDUCTANCES_1(source_nodes(ii),:) + CONDUCTANCES_1(source_dup(ii),:).*zero_locs_horz;
        CONDUCTANCES_1(:,source_nodes(ii)) = CONDUCTANCES_1(:,source_nodes(ii)) + CONDUCTANCES_1(:,source_dup(ii)).*zero_locs_vert;
        CONDUCTANCES_1(source_dup(ii),:) = 0;
        CONDUCTANCES_1(:,source_dup(ii)) = 0;

        % Find source_dup(ii) in CONDUCTANCES_2 and remove and replace.
        for jj = 1:length(CONDUCTANCES_2)
            if ~isempty(CONDUCTANCES_2{jj})
                mat = CONDUCTANCES_2{jj}; % Get data from cell
                [ind1,ind2,~] = find(mat(:,1:2) == source_dup(ii)); % Find where the to-be-replaced node IDs are
                mat(ind1,ind2) = source_nodes(ii); % Replace the node IDs
                CONDUCTANCES_2{jj} = mat; % Set values in holding matrix
            end
        end
    end
end

for ii = 1:target_count
    if target_nodes(ii) ~= target_dup(ii)
        % Find target_dup(ii) in MESHGRIDS_1, ELEMENTS_1, SURFACES_2,
        % CONDUCTANCES_1 and make the replacement or removal. Ignore if
        % already done.

        % Find target_dup(ii) in MESHGRIDS_1 and remove
        MESHGRIDS_1(MESHGRIDS_1(:,1)==target_dup(ii),:) = 0;
        
        % Find target_dup(ii) in ELEMENTS_1 and replace
        elem_surfs = ELEMENTS_1(:,1); % Keep first column since this may be overwritten in the next step
        change_indices = ELEMENTS_1==target_dup(ii); % This can be used to identify which elements have been changed, and therefore which surfaces have been changed
        ELEMENTS_1(change_indices) = target_nodes(ii);
        ELEMENTS_1(:,1) = elem_surfs; % Fix the surface IDs in the ELEMENTS_1 matrix
        
        % Find target_dup(ii) in SURFACES_1 and replace
        change_rows = sum(change_indices,2) > 0; % Find which elements have been changed
        change_surfs = elem_surfs(change_rows); % Extract the surfaces that have been changed
        %disp(change_surfs)
        for jj = 1:length(change_surfs) % Iterate through each changed surface and change nodes
            surf_nodes = SURFACES_2{change_surfs(jj)}.nodes;
            surf_nodes(surf_nodes==target_dup(ii))=target_nodes(ii);
            SURFACES_2{change_surfs(jj)}.nodes = surf_nodes;
        end

        % Find target_dup(ii) in CONDUCTANCES_1 and remove and replace
        CONDUCTANCES_1(target_nodes(ii),:) = CONDUCTANCES_1(target_nodes(ii),:) + CONDUCTANCES_1(target_dup(ii),:);
        CONDUCTANCES_1(:,target_nodes(ii)) = CONDUCTANCES_1(:,target_nodes(ii)) + CONDUCTANCES_1(:,target_dup(ii));
        CONDUCTANCES_1(target_dup(ii),:) = 0;
        CONDUCTANCES_1(:,target_dup(ii)) = 0;

        % Find source_dup(ii) in CONDUCTANCES_2 and remove and replace.
        for jj = 1:length(CONDUCTANCES_2)
            if ~isempty(CONDUCTANCES_2{jj})
                mat = CONDUCTANCES_2{jj}; % Get data from cell
                [ind1,ind2,~] = find(mat(:,1:2) == target_dup(ii)); % Find where the to-be-replaced node IDs are
                mat(ind1,ind2) = target_nodes(ii); % Replace the node IDs
                CONDUCTANCES_2{jj} = mat; % Set values in holding matrix
            end
        end
    end
end

end