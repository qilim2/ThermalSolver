function [MESHGRIDS_1,SURFACES_1] = surface_properties_assign_R02(surface_number,thermophysical_name,thermooptical_name_top,thermooptical_name_bottom,vec_x1,vec_x2,thickness,MESHGRIDS_1,SURFACES_1,SURFACES_2)
% Initiates the selected surface with defined properties

% INPUTS:
% surface_number - the number identifying the surface
% thermophysical_name - a string which is the identifier of the thermophysical property already created
% thermooptical_name_top - a string which is the identifier of the thermooptical property already created, identifying the thermooptical property to be applied to the top side of the surface
% thermooptical_name_bottom - a string which is the identifier of the thermooptical property already created, identifying the thermooptical property to be applied to the bottom side of the surface
% vec_x1 - the vector defining the x1 vector of the surface
% vec_x2 - the vector defining the x2 vector of the surface
% n_x1 - the number of nodes along the x1 direction
% n_x2 - the number of nodes along the x2 direction
% thickness - the thickness of the surface
% mesh_grids - global holding matrix of all nodes
% SURFACES_1 - global holding matrix of all surfaces
% THERMOPHYSICAL_1 - global holding matrix of all thermophysical properties
% THERMOOPTICAL_1 - global holding matrix of all thermo-optical properties

% OUTPUTS:
% mesh_grids - nodes are updated with thermal mass
% SURFACES_1 - surface properties are added to global holding matrix

% Version 2.0 completed 9/25/2023

if surface_number < 1
    fprintf('ERRROR (surface_properties_assign): Surface number must be 1 or greater.')
    return
end

surface.ID = surface_number; % Set number as the ID

% Check if surface exists
if isempty(SURFACES_2)
    fprintf('ERROR (surface_properties_assign): No surfaces currently exist.\n');
    return
else
    %temp = SURFACES_2{:}.ID==surface_number;
    if surface_number <= length(SURFACES_2)
        if isempty(SURFACES_2{surface_number})
            % Error handling
            fprintf('ERROR (surface_properties_assign): Surface does not exist.\n');
            return
        else
            % Create surface_props vector
            surface.area = norm(cross(vec_x1,vec_x2)); % Area of surface [m2]
            surface.material = thermophysical_name;
            surface.optical_top = thermooptical_name_top;
            surface.optical_bot = thermooptical_name_bottom;
            surface.thickness = thickness;
            surface.gtop = 1;
            surface.gbot = 1;
            surface.gtop_role = 0;
            surface.gbot_role = 0;
            surface.gtop_sn = 0;
            surface.gbot_sn = 0;
            SURFACES_1{surface.ID} = surface;
        end
    else
        fprintf('ERROR (surface_properties_assign): Surface does not exist.\n');
        return
    end
end
end