function surface_pre_visualization_R01(surfaces2plot,back_switch,marker_scale,MESHGRIDS_1,SURFACES_2)
% Plots the given surface. Display front/back with the back_switch.
% Only plots the surface. Does not display the elements.

% surface2plot: Desired surface to plot. Use a value <1 to select all surfaces.
% back_switch: Turns on the backface_color (0 for off, 1 for on)
% marker_scale: Size of the points in the plot

% Version 1.0 completed 10/2/2023
% Version 1.1 completed 12/4/2023 - Added input for multiple surfaces.
% Version 1.2 completed 12/8/2023 - Added input for selecting all surfaces.

if surfaces2plot < 1
    surfaces2plot = find(~cellfun(@isempty,SURFACES_2));
end

for ii = 1:length(surfaces2plot)
    surface2plot = surfaces2plot(ii);
    relevant_nodes = SURFACES_2{surface2plot}.nodes; % Make list of relevant nodes
    %disp(length(relevant_nodes))
    len = length(relevant_nodes); % Find length of list
    vec_x1 = SURFACES_2{surface2plot}.vec_x1;
    vec_x2 = SURFACES_2{surface2plot}.vec_x2;
    
    % Plot the nodes
    x_coords = zeros(length(relevant_nodes),1);
    y_coords = zeros(length(relevant_nodes),1);
    z_coords = zeros(length(relevant_nodes),1);
    for jj = 1:length(relevant_nodes)
        x_coords(jj) = MESHGRIDS_1(relevant_nodes(jj),4);
        y_coords(jj) = MESHGRIDS_1(relevant_nodes(jj),5);
        z_coords(jj) = MESHGRIDS_1(relevant_nodes(jj),6);
    end
    %disp(length(x_coords))
    scatter3(x_coords,y_coords,z_coords,marker_scale,'.');
    hold on
    offset_x = (x_coords(2)-x_coords(1))/3;
    offset_y = (y_coords(2)-y_coords(1))/3;
    offset_z = (z_coords(2)-z_coords(1))/3;
    textscatter3(x_coords+offset_x,y_coords+offset_y,z_coords+offset_z,string(relevant_nodes));
    %uistack(node_text,"top")
    
    % Plot frontface
    node_a = relevant_nodes(1);
    node_b = relevant_nodes(SURFACES_2{surface2plot}.n_x1);
    node_c = relevant_nodes(end);
    node_d = relevant_nodes(len-SURFACES_2{surface2plot}.n_x1+1);
    x_coords = [MESHGRIDS_1(node_a,4),MESHGRIDS_1(node_b,4),MESHGRIDS_1(node_c,4),MESHGRIDS_1(node_d,4)];
    y_coords = [MESHGRIDS_1(node_a,5),MESHGRIDS_1(node_b,5),MESHGRIDS_1(node_c,5),MESHGRIDS_1(node_d,5)];
    z_coords = [MESHGRIDS_1(node_a,6),MESHGRIDS_1(node_b,6),MESHGRIDS_1(node_c,6),MESHGRIDS_1(node_d,6)];
    c = [1 1 1 1];
    frontface = fill3(x_coords,y_coords,z_coords,c);
    frontface.FaceAlpha = 0.6;
    pbaspect([1 1 1]) % Set aspect ratio of plot to be 1:1:1
    
    % Plot backface
    if back_switch == 1
        
        vec_x3 = cross(vec_x1,vec_x2);
        unit_x3 = vec_x3/norm(vec_x3);
    
        back_off = 1e3; % Scale factor backface offset
        back_x = x_coords - unit_x3(1)/back_off;
        back_y = y_coords - unit_x3(2)/back_off;
        back_z = z_coords - unit_x3(3)/back_off;
    
        c = [0.5 0.5 0.5 0.5];
        backface = fill3(back_x,back_y,back_z,c);
        backface.FaceAlpha = 0.5;
    end
    
end
pbaspect([1 1 1]) % Set aspect ratio of plot to be 1:1:1
colorbar
axis equal

end