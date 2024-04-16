function [elem_areas, elem_normals, nodal_areas] = nodal_area_elem_normals_R01(MESHGRIDS_1,ELEMENTS_1)
% Calculates the surface area partitioned to each node in the listed
% elements. Also calculates normals for each element.

% Version 1.0 completed 11/16/2023

%{
%% Method 1
% Calculate area of each element using area from surfaces 1
% Partition area of each element to each node, summing up with values from other elements

% Find non-empty surfaces
surface_list = find(~cellfun(@isempty,SURFACES_1));

% Extract areas of surfaces and input into holding matrix at the matching index
surface_areas = zeros(length(SURFACES_1));
for ii = 1:length(surface_list)
    current_surface = surface_list(ii);
    surface_areas(current_surface) = SURFACES_1{current_surface}.area;
end

% List surfaces the elements belong to and sort, then count number of elements share each surface
elems_temp = sortrows(ELEMENTS_1(:,1));
elems_per_surf = accumarray(elems_temp,1);

% Partition areas to each element
%}

%% Method 2
% Calculate area of each element using the element nodes

% Find nodes a, b, and d or each element
node_a = ELEMENTS_1(:,2);
node_b = ELEMENTS_1(:,3);
node_d = ELEMENTS_1(:,5);

% Find point info of nodes
a_points = MESHGRIDS_1(node_a,4:6);
b_points = MESHGRIDS_1(node_b,4:6);
d_points = MESHGRIDS_1(node_d,4:6);

% Get vectors
ab_vecs = b_points-a_points;
ad_vecs = d_points-a_points;

% Create unit vectors or element normal
elem_normals = cross(ab_vecs,ad_vecs);
elem_areas = vecnorm(elem_normals,3,2); % Calculate element areas (norm of normal vector)
elem_normals = rdivide(elem_normals,elem_areas); % Normalize to unit vectors

% Assign area to each node
quarter_areas = elem_areas/4; % Quarter of each element area
nodal_areas = zeros(height(MESHGRIDS_1),1);

for ii = 1:height(ELEMENTS_1)
    for jj = 2:5
        current_node = ELEMENTS_1(ii,jj);
        nodal_areas(current_node) = nodal_areas(current_node) + quarter_areas(ii);
    end
end

end