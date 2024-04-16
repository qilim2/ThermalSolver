function mesh_visualization_R01(back_switch,normal_switch,MESHGRIDS_1,ELEMENTS_1)
% Plots the entire mesh and displays normal vectors and backfaces upon request.

% Plot frontface
for i = 1:height(ELEMENTS_1)
    relevant_nodes = ELEMENTS_1(i,:);
    x_coords = MESHGRIDS_1(relevant_nodes(2:5),4);
    y_coords = MESHGRIDS_1(relevant_nodes(2:5),5);
    z_coords = MESHGRIDS_1(relevant_nodes(2:5),6);
    c = [1;1;1;1];
    fill3(x_coords,y_coords,z_coords,c);
    hold on

    % Plot backface
    if back_switch == 1
        vec_x1 = [x_coords(2)-x_coords(1);y_coords(2)-y_coords(1);z_coords(2)-z_coords(1)];
        vec_x2 = [x_coords(4)-x_coords(1);y_coords(4)-y_coords(1);z_coords(4)-z_coords(1)];
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

    % Plot normal vector
    if normal_switch == 1
        vec_x1 = [x_coords(2)-x_coords(1);y_coords(2)-y_coords(1);z_coords(2)-z_coords(1)];
        vec_x2 = [x_coords(4)-x_coords(1);y_coords(4)-y_coords(1);z_coords(4)-z_coords(1)];
        vec_x3 = cross(vec_x1,vec_x2);
        arrow_scale = 1;
        unit_x3 = vec_x3/norm(vec_x3)*arrow_scale;

        x_mid = mean(x_coords);
        y_mid = mean(y_coords);
        z_mid = mean(z_coords);
        quiver3(x_mid,y_mid,z_mid,unit_x3(1),unit_x3(2),unit_x3(3),'r'); % Plot normal vector
    end
end
pbaspect([1 1 1]) % Set aspect ratio of plot to be 1:1:1
colorbar
axis equal
end