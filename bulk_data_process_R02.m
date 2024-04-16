function bulk_data_process_R02(filename,MESHGRIDS_1,SURFACES_1,SURFACES_2,THERMOPHYSICAL_1,THERMOOPTICAL_1,ELEMENTS_1,CONDUCTANCES_1,CONDUCTANCES_2,HEATLOADS_1,~)
% Processes surfaces properties and sends information to nodes. This is
% done just before the actual solving. Packages the information into an
% input file.

% Version 2 runs with fprintf instead of writelines.
% Version 2.1 adds CONDUCTANCES_2

% Needs to also process heat loads and determine if they are valid (i.e.,
% are placed on existing nodes or surfaces)

% NECESSARY INFORMATION (for solving and plotting)
% Nodes and Temperatures (MESHGRIDS_1)
% Conductances (CONDUCTANCES_1)
% Elements (ELEMENTS_1)
% Surfaces and Thermo-optical Properties (SURFACES_1, SURFACES_2, THERMOOPTICAL_1)

% Version 2.1 completed 12/11/2023


%% Merge CONDUCTANCES_2 with CONDUCTANCES_1
for ii = 1:length(CONDUCTANCES_2) % Look through CONDUCTANCES_2 cell array
    if ~isempty(CONDUCTANCES_2{ii}) % Find nonempty cells
        mat = CONDUCTANCES_2{ii}; % Extract data from the current cell
        for jj = 1:height(mat) % Loop through each row of the data matrix
            CONDUCTANCES_1(mat(jj,1),mat(jj,2)) = mat(jj,3); % Set the conductance value between the two nodes
        end
        %disp(CONDUCTANCES_1)
    end
end

%% Reduce matrices and vectors
nodes = find(MESHGRIDS_1(:,1)~=0);
MESHGRIDS_11 = MESHGRIDS_1(nodes,:);
CONDUCTANCES_11 = CONDUCTANCES_1(nodes,nodes);
HEATLOADS_11 = HEATLOADS_1(nodes);

%% Check MESHGRIDS for open spots
temp = MESHGRIDS_11(:,1);
check1 = isnan(MESHGRIDS_11(:,2)); % Find remaining NaN values
check2 = isnan(MESHGRIDS_11(:,3)); % Find remaining NaN values
if sum(check1) > 0
    check11 = temp(check1); % List node IDs for NaN values
    fprintf('ERROR (bulk_data_process): Thermal mass at node ID(s) %s not defined.\n',mat2str(check11));
    return
end
if sum(check2) > 0
    check21 = temp(check2); % List node IDs for NaN values
    fprintf('ERROR (bulk_data_process): Initial temperature at node ID(s) %s not defined.\n',mat2str(check21));
    return
end

%% Store properties
for ii = 2:height(THERMOPHYSICAL_1)
    THERMOPHYSICAL_11(ii-1,:) = [ii-1,str2double(THERMOPHYSICAL_1(ii,2:4))];
end
for ii = 2:height(THERMOOPTICAL_1)
    THERMOOPTICAL_11(ii-1,:) = [ii-1,str2double(THERMOOPTICAL_1(ii,2:3))];
end

%% Reprocess surfaces and reassign properties
jj=1;
for ii = 1:length(SURFACES_1)
    if ~isempty(SURFACES_1{ii})
        SURFACES_11(jj,1) = SURFACES_1{ii}.ID;
        SURFACES_11(jj,2) = SURFACES_1{ii}.area;
        indices = find(contains(THERMOPHYSICAL_1,SURFACES_1{ii}.material));
        SURFACES_11(jj,3) = indices(1)-1;
        indices = find(contains(THERMOOPTICAL_1,SURFACES_1{ii}.optical_top));
        SURFACES_11(jj,4) = indices(1)-1;
        indices = find(contains(THERMOOPTICAL_1,SURFACES_1{ii}.optical_bot));
        SURFACES_11(jj,5) = indices(1)-1;
        SURFACES_11(jj,6) = SURFACES_1{ii}.thickness;
        SURFACES_11(jj,7:9) = SURFACES_2{ii}.vec_x3';
        SURFACES_11(jj,10) = SURFACES_1{ii}.gtop;
        SURFACES_11(jj,11) = SURFACES_1{ii}.gtop_role;
        SURFACES_11(jj,12) = SURFACES_1{ii}.gtop_sn;
        SURFACES_11(jj,13) = SURFACES_1{ii}.gbot;
        SURFACES_11(jj,14) = SURFACES_1{ii}.gbot_role;
        SURFACES_11(jj,15) = SURFACES_1{ii}.gbot_sn;
        jj=jj+1;
    end
end

%% Get other element info
[elem_areas, elem_normals, ~] = nodal_area_elem_normals_R01(MESHGRIDS_1,ELEMENTS_1);


%% Store in bulk data file with headers
% Create file
fid = fopen(filename,'w');

% Header
fprintf(fid,['THERMAL SOLVER BULK DATA FILE\n' ...
    'Qi Lim - University of Illinois Urbana-Champaign - M.S. Aerospace Engineering thesis\n' ...
    '%s\n'],char(datetime));

% Overview
fprintf(fid,'\nOVERVIEW\n');
fprintf(fid,'Node Count, Element Count, Surface Count, Thermophysical Properties Count, Thermo-optical Properties Count\n');
overview = [length(nodes),height(ELEMENTS_1),height(SURFACES_11),height(THERMOPHYSICAL_11),height(THERMOOPTICAL_11)];
fprintf(fid,'%s\n',mat2str(overview));
fprintf(fid,'END OF SECTION\n');

% Nodes
fprintf(fid,'\nNODES\n');
fprintf(fid,'Node ID, Thermal Mass, Initial Temperature, X Coordinate, Y Coordinate, Z Coordinate\n');
for ii = 1:height(MESHGRIDS_11)
    node_info = MESHGRIDS_11(ii,:);
    fprintf(fid,'%s\n',mat2str(node_info));
end
fprintf(fid,'END OF SECTION\n');

% Elements
fprintf(fid,'\nELEMENTS\n');
fprintf(fid,'Parent Surface, Node A ID, Node B ID, Node C ID, Node D ID, Elem Area, Elem Norm X, Elem Norm Y, Elem Norm Z\n');
for ii = 1:height(ELEMENTS_1)
    elem_info = [ELEMENTS_1(ii,:), elem_areas(ii), elem_normals(ii,:)];
    fprintf(fid,'%s\n',mat2str(elem_info));
end
fprintf(fid,'END OF SECTION\n');

% Thermophysical properties
fprintf(fid,'\nTHERMOPHYSICAL PROPERTIES\n');
fprintf(fid,'Thermophysical Property ID, Density, Conductivity, Specific Heat\n');
for ii = 1:height(THERMOPHYSICAL_11)
    thermo_info = THERMOPHYSICAL_11(ii,:);
    fprintf(fid,'%s\n',mat2str(thermo_info));
end
fprintf(fid,'END OF SECTION\n');

% Thermo-optical properties
fprintf(fid,'\nTHERMO-OPTICAL PROPERTIES\n');
fprintf(fid,'Thermo-optical Property ID, Solar Absorptivity, IR Emissivity\n');
for ii = 1:height(THERMOOPTICAL_11)
    optical_info = THERMOOPTICAL_11(ii,:);
    fprintf(fid,'%s\n',mat2str(optical_info));
end
fprintf(fid,'END OF SECTION\n');

% Surfaces
%   Surface properties
fprintf(fid,'\nSURFACES PROPERTIES\n');
fprintf(fid,['Surface ID, Surface Area, Thermophysical Material, Topside Thermo-optical Property, ' ...
    'Bottomside Thermo-optical Property, Surface Thickness, Surface Normal X, Surface Normal Y, ' ...
    'Surface Normal Z, Topside Radiation Group, Topside Rad Group Role, Topside Space Node, ' ...
    'Bottomside Radiation Group, Bottomside Rad Group Role, Bottomside Space Node\n']);
for ii = 1:height(SURFACES_11)
    surf_info = SURFACES_11(ii,:);
    fprintf(fid,'%s\n',mat2str(surf_info));
end
fprintf(fid,'END OF SECTION\n');
%   Surface definitions
fprintf(fid,'\nSURFACES DEFINITIONS\n');
fprintf(fid,'Surface ID, List of nodes in the surface\n');
for ii = 1:length(SURFACES_2)
    if ~isempty(SURFACES_2{ii})
        surf_info = [SURFACES_2{ii}.ID SURFACES_2{ii}.nodes'];
        fprintf(fid,'%s\n',mat2str(surf_info));
    end
end
fprintf(fid,'END OF SECTION\n');

% Conductances
fprintf(fid,'\nCONDUCTANCES\n');
fprintf(fid,'Conductances associated with every other node. Row defines which node is the source.\n');
for ii = 1:height(CONDUCTANCES_11)
    cond_info = CONDUCTANCES_11(ii,:);
    fprintf(fid,'%s\n',mat2str(cond_info));
end
fprintf(fid,'END OF SECTION\n');

% User heat loads
fprintf(fid,'\nUSER DEFINED HEAT LOADS\n');
fprintf(fid,'Heat loads applied to each node.\n');
cond_info = HEATLOADS_11;
fprintf(fid,'%s\n',mat2str(cond_info));
fprintf(fid,'END OF SECTION\n');

% Radiation Group Labels
% fprintf(fid,'\nRAD GROUP LABELS\n');
% for ii = 1:height(RADGROUPS_1)
%     cond_info = CONDUCTANCES_11(ii,:);
%     fprintf(fid,'%s\n',mat2str(cond_info));
% end

% Close file
fprintf(fid,'\nEND FILE');
fclose(fid);
end