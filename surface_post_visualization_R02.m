function surface_post_visualization_R02(MESHGRIDS,ELEMENTS)
% Plots the elements. Display normal vector with the back_switch.
% This version is specifically for finding nodes whose indices do not match
% their node number.

% Version 2 completed 10/17/2023

reshaped_indices = zeros(4,1);
for i = 1:height(ELEMENTS)
    relevant_nodes = ELEMENTS(i,2:5);
    for j = 1:length(relevant_nodes)
        reshaped_indices(j) = find(relevant_nodes(j)==MESHGRIDS(:,1));
    end
    x_coords = MESHGRIDS(reshaped_indices(:,1),4);
    y_coords = MESHGRIDS(reshaped_indices(:,1),5);
    z_coords = MESHGRIDS(reshaped_indices(:,1),6);
    nodal_temperatures = MESHGRIDS(reshaped_indices(:,1),3);
    fill3(x_coords,y_coords,z_coords,nodal_temperatures);
    hold on
end
pbaspect([1 1 1])
colormap(jet(256))
colorbar
%clim([0,20])
axis equal
hold off
pause(0.000001)
end