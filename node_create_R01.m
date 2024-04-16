function MESHGRIDS_1 = node_create_R01(ID,thermal_mass,T_init,x_coord,y_coord,z_coord,MESHGRIDS_1)
% Creates an individual, unlinked node.
% Input the ID, the thermal mass, the initial temperature, and the global
% x, y, z coordinates.
% Version 1.0 completed 10/2/2023

if height(MESHGRIDS_1) >= ID
    if MESHGRIDS_1(ID,1) == ID
        fprintf('ERROR (create_node): Node ID already exists.\n')
        return
    end
end
MESHGRIDS_1(ID,:) = [ID, thermal_mass, T_init, x_coord, y_coord, z_coord];
end