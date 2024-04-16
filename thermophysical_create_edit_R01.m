function THERMOPHYSICAL_1 = thermophysical_create_edit_R01(thermophysical_name,rho,k,cp,THERMOPHYSICAL_1)
% Creates a material thermophysical assignment specifying the density (rho [kg/m3]),
% conductivity (k [W/m/K]), and specific heat of a material (cp [J/kg/K]).
% Version 1.0 completed ~9/19/2023.

material_thermophysical = [thermophysical_name,rho,k,cp]; % Package the properties

indeces = find(contains(THERMOPHYSICAL_1,thermophysical_name)); % Find indeces in global storage matrix 
if isempty(indeces)
    THERMOPHYSICAL_1 = [THERMOPHYSICAL_1;material_thermophysical]; % If indeces cannot be found (i.e. the material is not already created), then add the material
else
    THERMOPHYSICAL_1(indeces(1),:) = material_thermophysical; % If indeces are found, replace the material properties
end
end