function surface_post_visualization_R03(timestep_number,MESHGRIDS_11,solver_temperature_output,elem_nodes_map)
% Plots elemental temperatures at the given timestep.

% MESHGRIDS_11 holds nodal coordinates and node numbers which do not match
% actual node indices.
% elem_nodes_map holds the indices of each node that maps to MESHGRIDS_11
% in each element.
% solver_temperature_output holds the temperatures at each instance in time

% Version 3.0 competed 4/8/2024

for i = 1:height(elem_nodes_map)
    relevant_node_indices = elem_nodes_map(i,:);
    x_coords = MESHGRIDS_11(relevant_node_indices,4);
    y_coords = MESHGRIDS_11(relevant_node_indices,5);
    z_coords = MESHGRIDS_11(relevant_node_indices,6);
    nodal_temperatures = solver_temperature_output(relevant_node_indices,timestep_number);
    fill3(x_coords,y_coords,z_coords,nodal_temperatures);
    hold on
end
pbaspect([1 1 1])
colormap(jet(256))
colorbar
%clim([0,20])
axis equal
xlabel('X'); ylabel('Y');zlabel('Z');
hold off
pause(0.000001)
end