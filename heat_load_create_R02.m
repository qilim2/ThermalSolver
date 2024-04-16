function HL_vec = heat_load_create_R02(HL,HL_nodes,HL_surfaces,HL_surface_method,MESHGRIDS_1,SURFACES_2)
% Creates a heat load on either a node or a surface. Outputs as a vector 
% representing the heat on each node as produced by this heat load.
% HL - Heat load [W]
% Set HL_nodes, HL_surface = -1, HL_surface_method if not used
% HL_surface method: 1 to set all nodes in the surface as HL, 2 to divide HL among all nodes in the surface
% Version 2.0 completed 9/25/2023

nodes_total = height(MESHGRIDS_1); % Find total nodes in the mesh
HL_vec = [MESHGRIDS_1(:,1), zeros(nodes_total,1)]; % Initialize vector of zeros representing heat loads on each node matched with node number

if HL_nodes ~= -1
    for i = 1:length(HL_nodes)
        relevant_node = HL_nodes(i);
        if MESHGRIDS_1(relevant_node,1) == 0 || relevant_node > height(MESHGRIDS_1)
            % Error handling
            fprintf('ERROR (create_heat_load): Node %i does not exist.\n',relevant_node);
            return
        else
            HL_vec(relevant_node,2) = HL;
        end
    end
end
if HL_surfaces ~= -1
    for i=1:length(HL_surfaces)
        surface_number = HL_surfaces(i);
        % Check if surface exists
        if isempty(SURFACES_2{surface_number})
            % Error handling
            fprintf('ERROR (create_heat_load): Surface %i does not exist.\n',surface_number);
            return
        else
            relevant_node_indeces = SURFACES_2{surface_number}.nodes; % List relevant indeces as described in the surface
            for j=1:length(relevant_node_indeces) % Iterate through the list of relevant node indeces
                if HL_surface_method == 1
                    HL_vec(relevant_node_indeces(j),2) = HL; % Apply identical heat load to all nodes in the surface
                elseif HL_surface_method == 2
                    HL_vec(relevant_node_indeces(j),2) = HL/length(relevant_node_indeces); % Apply portion of the heat load to the nodes in the surface
                else
                    % Error handling
                    fprintf('ERROR (create_heat_load): Heat load application method not specified correctly.\n');
                    fprintf('Input -1 if a surface is not specified.\n'); 
                    fprintf('Input 1 to apply the same heat load across all nodes in the surface.\n');
                    fprintf('Input 2 to divide the heat load equivalently across the nodes in the surface\n')
                end
            end
        end
    end
end

end