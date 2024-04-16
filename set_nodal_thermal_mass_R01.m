function MESHGRIDS_1 = set_nodal_thermal_mass_R01(surfaces2edit,thermophysical_name,THERMOPHYSICAL_1,MESHGRIDS_1,SURFACES_1,SURFACES_2)
% Assigns thermal mass to the nodes given the surface and the material.
% Version 1.0 completed 9/25/2023

% Check if surfaces exist
if isempty(SURFACES_2)
    fprintf('ERROR (set_nodal_thermal_mass): No surfaces currently exist.\n');
    return
end

for i=1:length(surfaces2edit)
    
    surface_number = surfaces2edit(i);

    % Check if surface exists
    if surface_number <= length(SURFACES_2)
        if isempty(SURFACES_2{surface_number})
            fprintf('ERROR (set_nodal_thermal_mass): Surface %i does not exist.\n',surface_number);
            return
        end
    else
        fprintf('ERROR (set_nodal_thermal_mass): Surface %i does not exist.\n',surface_number);
        return
    end
    
    % Extract data from relevant holding matrices
    surface_area = SURFACES_2{surface_number}.area;
    n_x1 = SURFACES_2{surface_number}.n_x1;
    n_x2 = SURFACES_2{surface_number}.n_x2;
    indeces = find(contains(THERMOPHYSICAL_1,thermophysical_name)); % Find indices of material in global storage matrix
    rho = str2double(THERMOPHYSICAL_1(indeces(1),2));
    cp = str2double(THERMOPHYSICAL_1(indeces(1),4));
    thickness = SURFACES_1{surface_number}.thickness;

    % Divide mass by nodes
    surf_node_total = n_x1*n_x2; % Total number of elements in the surface
    surface_mass = rho*surface_area*thickness; % Mass of surface [kg]
    int_node_total = (n_x1-2)*(n_x2-2); % Number of interior nodes
    corn_node_total = 4; % Number of corner nodes
    edge_node_total = surf_node_total-int_node_total-corn_node_total; % Number of edge nodes
    divisions_total = (n_x1-1)*(n_x2-1)*4; % Total number of divisions of the surface into quarter elements per node (num elems * 4)
    quarter_node_mass = surface_mass/divisions_total; % Mass per quarter node
    int_node_mass = 4*quarter_node_mass; % Mass of interior node
    corn_node_mass = quarter_node_mass; % Mass of corner node
    edge_node_mass = 2*quarter_node_mass; % Mass of edge node

    % Assign mass to nodes
    surface_nodes = SURFACES_2{surface_number}.nodes; % Find indeces in the grids matrix of the surface nodes relevant to this surface
    node_mass_vec = zeros(surf_node_total,1);
    node_mass_vec(1) = corn_node_mass; % First node in surface
    node_mass_vec(surf_node_total) = corn_node_mass; % Last node in surface
    node_mass_vec(n_x1) = corn_node_mass; % Node in bottom right corner
    node_mass_vec((n_x2-1)*n_x1+1) = corn_node_mass; % Node in top left corner
    if edge_node_total ~= 0
        node_mass_vec(2:n_x1-1) = edge_node_mass; % Edge nodes along bottom
        for j = 1:n_x2-2 % Edge nodes along left side
            node_mass_vec(j*n_x1+1) = edge_node_mass;
        end
        for j = 2:n_x2-1 % Edge nodes along right side
            node_mass_vec(j*n_x1) = edge_node_mass;
        end
        node_mass_vec((n_x2-1)*n_x1+2:surf_node_total-1) = edge_node_mass; % Edge nodes along top
    end
    if int_node_total ~= 0
        for k = 1:n_x2-2
            node_mass_vec(k*n_x1+2:k*n_x1+n_x1-1) = int_node_mass; % Interior nodes
        end
    end

    %disp(node_mass_vec)
    
    % Convert to thermal mass and assign to nodes in grid matrix
    for j=1:surf_node_total
        temp = surface_nodes(j);
        MESHGRIDS_1(temp,2) = node_mass_vec(j)*cp; % Set thermal mass of node
    end
end

end