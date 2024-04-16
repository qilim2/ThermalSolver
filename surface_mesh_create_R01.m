function [MESHGRIDS_1,SURFACES_2] = surface_mesh_create_R01(surface_number,origin,relative,vec_x1,vec_x2,n_x1,n_x2,start_node,MESHGRIDS_1,SURFACES_2)
% Generates a surface of nodes. The output matrix is of the following
% format:
% [Node number, Thermal mass, Current temperature, Parent surface, x1, x2, x3]
% Thermal mass and current temperature are defined separately.
% "Relative" is a switch that identifies whether the values given for
% vec_x1 and vec_x2 are nodal positions (not actually vectors) in the global 
% CSYS or are relative vectors from the origin point. (0 is global, 1 is relative)
% Version 1.3 completed 9/29/2023

% Check if nodes exist and set warning
if ~isempty(SURFACES_2) % If SURFACES_2 has data inside...
    if surface_number <= length(SURFACES_2) % If the input surface ID is not greater than the length of SURFACES_2...
        if ~isempty(SURFACES_2{surface_number}) % Check if the cell at the index provided by the surface ID is already occupied.
            fprintf('Warning: Surface %i already exists. Utilize the EDIT function to edit the surface. (NOTE to self: This is not yet implemented as of 9/29/2023.)\n',surface_number);
            return
        end
    end
end

if ~isempty(MESHGRIDS_1)
    temp = find(MESHGRIDS_1(:,1)==start_node, 1);
    if ~isempty(temp)
        % Error handling
        fprintf('ERROR (surface_mesh_create): Nodes overlap starting at %i. Please enter a different start node.\n',temp);
        %start_node = input("Please enter a different start node.");
        return
    end
end

if n_x1 < 2 || n_x2 < 2
    % Check if node count is reasonable
    fprintf('ERROR (surface_mesh_create): Too few nodes to create a surface! Node count along each edge must be 2 or greater.\n');
    return
else
    if relative ~= 1 % Convert nodal positions to relative vectors if global positions are given
        vec_x1 = vec_x1-origin;
        vec_x2 = vec_x2-origin;
    end
    vec_x3 = cross(vec_x1,vec_x2);
    l_x1 = norm(vec_x1); % Length in x1 direction [m]
    l_x2 = norm(vec_x2); % Length in x2 direction [m]
    l_x3 = norm(vec_x3); % Create magnitude of vec_x3 for normalization
    unit_x1 = vec_x1/l_x1; % Unit vector in x1 direction
    unit_x2 = vec_x2/l_x2; % Unit vector in x2 direction
    unit_x3 = vec_x3/l_x3; % Unit vector in x3 direction
    inner_angle = acosd(dot(vec_x1,vec_x2)/(norm(vec_x1)*norm(vec_x2))); % Inner angle of vectors
    if inner_angle ~= 90
        fprintf('ERROR (surface_mesh_create): Not a rectangular plate! Angle between given vectors must be 90 degrees.\n')
        return
    end
    
    % Defining nodes
    surface_node_total = n_x1*n_x2; % Total number of nodes in surface 01
    ones_vec = ones(surface_node_total,1); % Ones vector sized for node count
    
    % Local coordinates
    coords_x1 = linspace(0,l_x1,n_x1)'; % x1 coordinates of nodes [m]
    coords_x2 = linspace(0,l_x2,n_x2)'; % x2 coordinates of nodes [m]
    
    % Create nodes
    nodes = (start_node:start_node+surface_node_total-1)'; % Create list of nodes
    %fprintf('surface_mesh_create - Surface %i node count: %i\n',surface_number,surface_node_total);
    
    % Coords matrix
    current_local_node = 1; % Use for iterating through list of nodes
    coords_mat = ones(surface_node_total,3);
    for i=1:n_x2
        for j=1:n_x1
            coords_mat(current_local_node,:) = [coords_x1(j),coords_x2(i),0];
            current_local_node = current_local_node+1;
        end
    end
    
    % Convert local coords to global coords
    gcoords_mat = ones(3,surface_node_total); % Initialization
    for it = 1:surface_node_total
        gcoords_mat(:,it) = local2globalcoord(coords_mat(it,:)','rr',origin,[unit_x1,unit_x2,unit_x3]);
    end
    gcoords_mat = gcoords_mat';
    
    % Grid Matrix: [Node number, Thermal mass, Current temperature, x1, x2, x3]
    surface_grids = [nodes, NaN*ones_vec, NaN*ones_vec, gcoords_mat(:,1),gcoords_mat(:,2),gcoords_mat(:,3)];
    
    % Assign to MESHGRIDS_1
    for j = 1:height(surface_grids)
        MESHGRIDS_1(surface_grids(j,1),:) = surface_grids(j,:);
    end


    %MESHGRIDS_1 = [MESHGRIDS_1;surface_grids]; % Add surface grids to list of all grids
    %MESHGRIDS_1 = sortrows(MESHGRIDS_1); % Sort grids by node number
    
    surface.ID = surface_number; % Store surface ID
    surface.nodes = nodes; % Store nodes related
    surface.n_x1 = n_x1;
    surface.n_x2 = n_x2;
    surface.vec_x1 = vec_x1;
    surface.vec_x2 = vec_x2;
    surface.vec_x3 = unit_x3;
    surface.area = norm(cross(vec_x1,vec_x2));
    SURFACES_2{surface_number} = surface; % Store surface

end
end