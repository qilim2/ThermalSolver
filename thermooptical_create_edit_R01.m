function THERMOOPTICAL_1 = thermooptical_create_edit_R01(thermooptical_name,alpha,epsilon,THERMOOPTICAL_1)
% Creates a material thermo-optical property specifying the solar absorptivity
% (alpha) and infrared emissivity (epsilon).
% Version 1.0 completed ~9/19/2023.

material_thermooptical = [thermooptical_name,alpha,epsilon]; % Package the properties

indeces = find(contains(THERMOOPTICAL_1,thermooptical_name)); % Find indeces in global storage matrix 
if isempty(indeces)
    THERMOOPTICAL_1 = [THERMOOPTICAL_1;material_thermooptical]; % If indeces cannot be found (i.e. the material is not already created), then add the material
else
    THERMOOPTICAL_1(indeces(1),:) = material_thermooptical; % If indeces are found, replace the material properties
end
end