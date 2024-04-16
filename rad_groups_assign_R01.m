function SURFACES_1 = rad_groups_assign_R01(ID,topbot,group_ID,role,sn,SURFACES_1)
% Assigns surfaces with radiation groups.
% ID - Surface ID
% topbot - Select both (1), topside (2), bottomside (3)
% group_ID - Group number
% role - Role of the surface (0=ignored, 1=active+blocker, 2=blocker)
% sn - Boolean for turning on/off spacenode vision (0=off, 1=on)

% Version 1.0 completed 1/18/2024

if role ~= 0 && role ~= 1 && role ~= 2
    fprintf('ERROR (rad_groups_assign): role must be 9 (ignored), 1 (active++blocker) or 2 (blocker).\n')
    return
end

if ID > length(SURFACES_1) || isempty(SURFACES_1{ID}) || ID < 1
    fprintf('ERROR (rad_groups_assign): Invalid surface ID.\n')
    return
end

if sn ~= 1 && sn ~= 0
    fprintf('ERROR (rad_groups_assign): Space node boolean should be 0 (off) or 1 (on). This describes the surface side(s) capabilitiy to see the space node.\n')
    % Note: Does not necessarily need to be active. Can be a blocker and
    % therefore prevents light from reaching other objects.
    return
end

% Assign surface to group
if topbot == 1
    SURFACES_1{ID}.gtop = group_ID;
    SURFACES_1{ID}.gtop_role = role;
    SURFACES_1{ID}.gtop_sn = sn;
    SURFACES_1{ID}.gbot = group_ID;
    SURFACES_1{ID}.gbot_role = role;
    SURFACES_1{ID}.gbot_sn = sn;
elseif topbot == 2
    SURFACES_1{ID}.gtop = group_ID;
    SURFACES_1{ID}.gtop_role = role;
    SURFACES_1{ID}.gtop_sn = sn;
elseif topbot == 3
    SURFACES_1{ID}.gbot = group_ID;
    SURFACES_1{ID}.gbot_role = role;
    SURFACES_1{ID}.gbot_sn = sn;
else
    fprintf('ERROR (rad_groups_assign): topbot should be 1 (both sides), 2 (topside), 3 (bottomside).\n')
    return
end
end