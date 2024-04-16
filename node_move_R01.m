function MESHGRIDS_1 = node_move_R01(ID,relative,x1,x2,x3,MESHGRIDS_1,ELEMENTS_1)
% Moves a node as long as it is not part of a surface.
% ID - Node ID
% relative - Switch between relative movement from current origin (1) and absolute position of origin (0)
% x1 - movement in x1 [m]
% x2 - movement in x2 [m]
% x3 - movement in x3 [m]

% Version 1.0 completed 12/12/2023

%% Check if node exists
exist = MESHGRIDS_1(:,1) == ID;
exist = sum(exist,"all");
if exist == 0
    fprintf('ERROR(node_move): Node selected does not exist.\n')
    return
end

%% Check if node is part of a surface
in_surf = ELEMENTS_1(:,2:5) == ID;
in_surf = sum(in_surf,"all");
if in_surf ~= 0
    fprintf('ERROR (node_move): Node selected is part of a surface and cannot be moved.\n')
    return
end

%% Move
if relative == 0
    MESHGRIDS_1(ID,4:6) = [x1,x2,x3];
elseif relative == 1
    MESHGRIDS_1(ID,4:6) = MESHGRIDS(ID,4:6) + [x1,x2,x3];
else
    fprintf('ERROR (node_move): Invalid input for relative. Choose 1 for movement relative to initial position or 0 for absolute position.\n')
    return
end


end