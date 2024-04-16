function INTERNODAL_AREA_1 = internodal_area_R03(surface_number, SURFACES_1, SURFACES_2)
% Calculates the internodal area between each node in a surface.
% Version 3 specifically creates a matrix of the same size as the surface input.
% Solves without for-loops or if-else statements.
% Version 3.0 completed 10/20/2023
 
% Extract data from relevant holding matrices
nodes = SURFACES_2{surface_number}.nodes;
n_x1 = SURFACES_2{surface_number}.n_x1;
n_x2 = SURFACES_2{surface_number}.n_x2;
v_x1 = SURFACES_2{surface_number}.vec_x1;
v_x2 = SURFACES_2{surface_number}.vec_x2;
thickness = SURFACES_1{surface_number}.thickness;

% Create widths
dx1 = norm(v_x1)/(n_x1-1);
dx2 = norm(v_x2)/(n_x2-1);

% Preallocate the matrix
node_count = length(nodes);
INTERNODAL_AREA_1 = zeros(node_count, node_count);

% Identify edge nodes for horizontal connections
nedge_horz = [1:n_x1, node_count-n_x1+1:node_count];

% Identify edge nodes for vertical connections
nedge_vert = [1:n_x1:node_count, n_x1:n_x1:node_count];

% Horizontal
n_horz = 1:node_count-1;
edge_horz = ismember(n_horz, nedge_horz);
n_horz_p1 = n_horz + 1;

% Input into mat
INTERNODAL_AREA_1(sub2ind([node_count, node_count], n_horz(edge_horz), n_horz_p1(edge_horz))) = thickness*dx2/2;
INTERNODAL_AREA_1(sub2ind([node_count, node_count], n_horz_p1(edge_horz), n_horz(edge_horz))) = thickness*dx2/2;

INTERNODAL_AREA_1(sub2ind([node_count, node_count], n_horz(~edge_horz), n_horz_p1(~edge_horz))) = thickness*dx2;
INTERNODAL_AREA_1(sub2ind([node_count, node_count], n_horz_p1(~edge_horz), n_horz(~edge_horz))) = thickness*dx2;

% For vertical connections
n_vert = 1:node_count-n_x1;
edge_vert = ismember(n_vert, nedge_vert);
n_vert_p1 = n_vert + n_x1;

INTERNODAL_AREA_1(sub2ind([node_count, node_count], n_vert(edge_vert), n_vert_p1(edge_vert))) = thickness*dx1/2;
INTERNODAL_AREA_1(sub2ind([node_count, node_count], n_vert_p1(edge_vert), n_vert(edge_vert))) = thickness*dx1/2;

INTERNODAL_AREA_1(sub2ind([node_count, node_count], n_vert(~edge_vert), n_vert_p1(~edge_vert))) = thickness*dx1;
INTERNODAL_AREA_1(sub2ind([node_count, node_count], n_vert_p1(~edge_vert), n_vert(~edge_vert))) = thickness*dx1;

end