function CONDUCTANCES_1 = surface_conductance_R02(surface_number,CONDUCTANCES_1,SURFACES_1,SURFACES_2,THERMOPHYSICAL_1,SURFACELENGTHS_1)
% Creates a matrix listing all the conductances between each node in a
% surface.
% Version 1.0 completed 10/17/2023

material = SURFACES_1{surface_number}.material;
indices = find(contains(THERMOPHYSICAL_1,material)); % Find indices of material in global storage matrix
k = str2double(THERMOPHYSICAL_1(indices,3)); % Ignore matlab's suggestion here
relevant_nodes = SURFACES_2{surface_number}.nodes;

% Internodal Area
A = internodal_area_R03(surface_number,SURFACES_1,SURFACES_2);
%fprintf('(Within surface_conductance) Internodal area calculated.\n')

% Internodal Lengths
L = SURFACELENGTHS_1(relevant_nodes,relevant_nodes);

% Calculate local Conductances
G_local = k*A./L;

% Assign local G to global CONDUCTANCES_1 matrix
CONDUCTANCES_1(relevant_nodes,relevant_nodes) = G_local;

end