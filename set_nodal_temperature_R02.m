function MESHGRIDS_1 = set_nodal_temperature_R02(T,nodes2edit,surfaces2edit,MESHGRIDS_1,SURFACES_2)
% Sets the initial temperature of the nodes or the surfaces input
% Set nodes or surfaces = -1 if not used
% Version 2.0 completed 9/25/2023

if nodes2edit~=-1
    for i=1:length(nodes2edit)
        indeces = find(MESHGRIDS_1(:,1)==nodes2edit(i), 1);
        if isempty(indeces) % Check if node exists
            % Error handling
            fprintf('ERROR (set_nodal_temperature): Node %i does not exist.\n',nodes2edit(i));
            return
        else
            MESHGRIDS_1(indeces(1,1),3) = T; % Set node at index with initial temperature
        end
    end
end
if surfaces2edit~=-1
    % Check if surface exists
    if isempty(SURFACES_2)
        fprintf('ERROR (set_nodal_temperature): No surfaces currently exist.\n');
        return
    end

    for i=1:length(surfaces2edit)
        surface_number = surfaces2edit(i);
        % Check if surface exists
        if SURFACES_2{surface_number}.ID~=surface_number
            % Error handling
            fprintf('ERROR (set_nodal_temperature): Surface %i does not exist.\n',surface_number);
            return
        else
            relevant_node_indeces = SURFACES_2{surface_number}.nodes; % Find indeces in the grids matrix of the surface nodes relevant to this surface
            for j=1:length(relevant_node_indeces) % Iterate through the list of relevant node indeces
                MESHGRIDS_1(relevant_node_indeces(j),3) = T; % Set the temperature of each node
            end
        end
    end
end


end