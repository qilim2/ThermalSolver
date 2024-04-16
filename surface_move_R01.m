function MESHGRIDS_1 = surface_move_R01(ID,relative,x1,x2,x3,r1,r2,r3,MESHGRIDS_1,SURFACES_2)
% Moves and rotates a surface given the inputs.
% ID - Surface ID
% relative - Switch between relative movement from current origin (1) and absolute position of origin (0)
% x1 - movement in x1 [m]
% x2 - movement in x2 [m]
% x3 - movement in x3 [m]
% r1 - rotation about x1 [deg]
% r2 - rotation about x2 [deg]
% r3 - rotation about x3 [deg]
% Note that the x1, x2, and x3 axes are centered at the origin.

% Version 1.0 completed 12/12/2023

%% Check if surface exists
if ID > length(SURFACES_2)
    fprintf('ERROR (surface_move): Surface does not exist.\n')
    return
end
if isempty(SURFACES_2{ID})
    fprintf('ERROR (surface_move): Surface does not exist.\n')
    return
end

%% Translate
surf_nodes = SURFACES_2{ID}.nodes;
origin = surf_nodes(1);
origin_coords = MESHGRIDS_1(origin,4:6);
if relative == 0
    translation_vec = [x1,x2,x3] - origin_coords; % Find translation needed to get to the new origin point
elseif relative == 1
    translation_vec = [x1,x2,x3];
else
    fprintf('ERROR (surface_move): Invalid input for relative. Choose 1 for movement relative to initial position or 0 for absolute position.\n')
    return
end
for ii = 1:length(surf_nodes) % Iterate through each node in the surface and edit the coordinates
    current_node = surf_nodes(ii);
    MESHGRIDS_1(current_node,4:6) = MESHGRIDS_1(current_node,4:6) + translation_vec;
end
new_origin_coords = MESHGRIDS_1(origin,4:6);

%% Rotate
rotm = eul2rotm([deg2rad(r3),deg2rad(r2),deg2rad(r1)]); % Convert to rotation matrix (order is ZYX).
% Find vector between each point and the surface origin and apply rotation matrix
for ii = 2:length(surf_nodes)
    current_node = surf_nodes(ii);
    node_coords = MESHGRIDS_1(current_node,4:6);
    node_vec = node_coords - new_origin_coords;
    rotated_vec = rotm*node_vec';
    new_coords = new_origin_coords + rotated_vec';
    MESHGRIDS_1(current_node,4:6) = new_coords;
end

end