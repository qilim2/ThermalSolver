function HL_vec = heat_flux_create_R01(HF,HF_surfaces,MESHGRIDS_1,SURFACES_2)
% Creates a heat load on either a node or a surface. Outputs as a vector 
% representing the heat on each node as produced by this heat load.
% HF - Heat flux [W/m2]
% Version 1.1 completed 9/25/2023


% Check if surfaces exist
if isempty(SURFACES_2)
    fprintf('ERROR (create_heat_flux): No surfaces currently exist.\n');
    return
end

mesh_grids_size = size(MESHGRIDS_1); % Find size of mesh
nodes_total = mesh_grids_size(1); % Find total nodes in the mesh
HL_vec = [MESHGRIDS_1(:,1), zeros(nodes_total,1)]; % Initialize vector of zeros representing heat loads on each node matched with node number


for i=1:length(HF_surfaces)
    surface_number = HF_surfaces(i);
    %fprintf('Adding heat flux to surface %i\n',surface_number);
    % Check if surface exists
    if isempty(SURFACES_2{surface_number})
        % Error handling
        fprintf('ERROR (create_heat_flux): Surface %i does not exist.\n',surface_number);
        return
    else
        surface_area = SURFACES_2{surface_number}.area; % Obtain the surface area of the chosen surface
        HL = HF*surface_area; % Create a heat load using the heat flux multiplied by the surface area
        relevant_node_indeces = SURFACES_2{surface_number}.nodes; % Find relevant indeces as described in the surface
        for j=1:length(relevant_node_indeces) % Iterate through the list of relevant node indeces
            HL_vec(relevant_node_indeces(j),2) = HL/length(relevant_node_indeces); % Apply portion of the heat load to the nodes in the surface
        end
    end
end


end