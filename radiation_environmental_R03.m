function [Qdot_total_S,Qdot_total_A,Qdot_total_I,S] = radiation_environmental_R03(T_body,AF,vS,v_body,g,enormals,ecoords,rotms,alpha,epsilon,sab,elem_solar_proj_area,illum_t,illum_b,elem_body_t_VF,elem_body_b_VF,eclipse,elem_areas,R_body)
% Calculates the heating rate of elements based on incoming radiation
% (solar, albedo, Earth IR/planetshine). Solves one radiation group at a
% time.

% vS = vector from s/c center to Sun
% v_body = vector from s/c center to orbital body
% enormals = matrix of normal vectors for ALL elements
% rotms = cell array of rotation matrices needed to flatten each element to 2D

% g = struct with the following information (group information)
% .elems = list of element IDs in this group
% .tsn = boolean vector of length(.elems) indicating which elems have sn-active topsides
% .bsn = boolean vector of length(.elems) indicating which elems have sn-active bottomsides
% .trole = vector of length(.elems) indicating the role of each element topside
% .brole = vector of length(.elems) indicating the role of each element bottomside
% .centers = Nx3 matrix of coordinates of each element center

% sab = solar_angle_bool = boolean vector indicating which side of an element is illuminated by the sun

% ecoords = cell array of the coordinates of each element in 3x4 matrix format

% Version 3.0 completed 2/26/2024

%% Initialization
%TE = 255; % Earth average temperature [K]
SBC = 5.670374419e-8; % Stefan-Boltzmann Constant [W/m2/K4]
S = solar_flux_R01(norm(vS)); % Solar flux [W/m2]
vRayS = repmat(vS,length(g.elems),1) - g.centers; % From element to sun
vRay_body = repmat(v_body,length(g.elems),1) - g.centers; % From element to earth
dS = vRayS./vecnorm(vRayS,2,2); % Get unit vectors along rays
dbody = vRay_body./vecnorm(vRay_body,2,2); % Get unit vectors along rays

%% Element Prep
% Normal unit vectors
n_t = enormals(g.elems,:); % Get normal unit vectors of each element
% n_b = -n_t; % Reverse normal unit vectors for bottomside
sab = boolean(sab(g.elems)); % Reduce to relevant elements (solar angle bool)
illum_t = illum_t(g.elems); % Reduce to relevant elements (need top and bottom for orbiting body)
illum_b = illum_b(g.elems); % Reduce to relevant elements

%%% Solar
% Pare down number of elements needed to do shadowing checks by removing those that are not sn-active and illuminated
elems_t_bool_1 = boolean(g.tsn) & sab & ~boolean(eclipse); % Create boolean vector indicating if each element topside is both sn-active and illuminated
elems_b_bool_1 = boolean(g.bsn) & ~sab & ~boolean(eclipse); % Create boolean vector indicating if each element bottomside is both sn-active and illuminated
elems_t_solar_shadecheck = g.elems(elems_t_bool_1); % Reduce list of elements to only the elements that are both sn-active and illuminated
elems_b_solar_shadecheck = g.elems(elems_b_bool_1); % Reduce list of elements to only the elements that are both sn-active and illuminated

%%% Orbital body
% Pare down number of elements needed to do shadowing checks by removing those that are not sn-active and illuminated
elems_t_bool_2 = boolean(g.tsn) & illum_t; % Create boolean vector indicating if each element topside is both sn-active and illuminated
elems_b_bool_2 = boolean(g.bsn) & illum_b; % Create boolean vector indicating if each element bottomside is both sn-active and illuminated
elems_t_body_shadecheck = g.elems(elems_t_bool_2); % Reduce list of elements to only the elements that are both sn-active and illuminated
elems_b_body_shadecheck = g.elems(elems_b_bool_2); % Reduce list of elements to only the elements that are both sn-active and illuminated

%%% Blockers
blockers = g.elems(g.trole ~= 0 | g.brole ~= 0); % Get list of blockers by finding elements that are not ignored, including both topside or bottomside
blocker_centers = g.centers(blockers,:);

%%% Reduce vectors for shadowing calculations (also avoids broadcasting)
dS_t_1 = dS(elems_t_bool_1,:); % Reduce list of unit vectors
dS_b_1 = dS(elems_b_bool_1,:); % Reduce list of unit vectors
dbody_t_1 = dbody(elems_t_bool_2,:); % Reduce list of unit vectors
dbody_b_1 = dbody(elems_b_bool_2,:); % Reduce list of unit vectors

%% Solar
%%% Shadowing
parfor ii = 1:length(blockers) % Check through each source elem
    blocker = blockers(ii);
    tc = ecoords{blocker};
    tr = rotms{blocker};
    p = blocker_centers(ii,:)';
    n = n_t(ii,:)';
    
    o = g.centers(elems_t_solar_shadecheck,:)'; % Origin of the rays
    d = dS_t_1'; % Ray unit vector out to the sun
    
    [~,~,hit_check(ii,:)] = ray_elem_intersect_R04(o,d,p,n,tc,tr); % Perform element hit check
end
unshadowed_t_solar = ~boolean(sum(hit_check,1)); % 0 = in shadow, 1 = not in shadow
hit_check = []; % Reset hit check

parfor ii = 1:length(blockers)
    blocker = blockers(ii);
    tc = ecoords{blocker};
    tr = rotms{blocker};
    p = blocker_centers(ii,:)';
    n = n_t(ii,:)';
    
    o = g.centers(elems_b_solar_shadecheck,:)'; % Origin of the rays
    d = dS_b_1'; % Ray unit vector out to the sun
    
    [~,~,hit_check(ii,:)] = ray_elem_intersect_R04(o,d,p,n,tc,tr); % Perform element hit check
end
unshadowed_b_solar = ~boolean(sum(hit_check,1)); % 0 = in shadow, 1 = not in shadow
hit_check = []; % Reset hit check

% Elements that are unshadowed, have active-sn, have active-role (i.e. elements that should have Qdot):
elems2calc_t_solar = elems_t_solar_shadecheck(boolean(unshadowed_t_solar));
elems2calc_b_solar = elems_b_solar_shadecheck(boolean(unshadowed_b_solar));

%%% Heating calculation
Qdot_t_solar = S*alpha(elems2calc_t_solar,1).*elem_solar_proj_area(elems2calc_t_solar);
Qdot_b_solar = S*alpha(elems2calc_b_solar,2).*elem_solar_proj_area(elems2calc_b_solar);

%% Orbital body
%%% Shadowing
parfor ii = 1:length(blockers) % Check through each source elem
    blocker = blockers(ii);
    tc = ecoords{blocker};
    tr = rotms{blocker};
    p = blocker_centers(ii,:)';
    n = n_t(ii,:)';
    
    % o = g.centers(elems_t_body_shadecheck,:)'; % Origin of the rays
    o = g.centers(elems_t_bool_2,:)';
    n_sources = enormals(elems_t_body_shadecheck,:)';
    body_edges = n_sources*R_body;
    d = (body_edges - o);
    d = d./(vecnorm(d')');
    % d = dbody_t_1'; % Ray unit vector out to the body
    
    [~,~,hit_check(ii,:)] = ray_elem_intersect_R04(o,d,p,n,tc,tr); % Perform element hit check
end
unshadowed_t_body = ~boolean(sum(hit_check,1)); % 0 = in shadow, 1 = not in shadow
hit_check = []; % Reset hit check

parfor ii = 1:length(blockers)
    blocker = blockers(ii);
    tc = ecoords{blocker};
    tr = rotms{blocker};
    p = blocker_centers(ii,:)';
    n = n_t(ii,:)';
    
    % o = g.centers(elems_b_body_shadecheck,:)'; % Origin of the rays
    o = g.centers(elems_b_bool_2,:)';
    n_sources = -enormals(elems_b_body_shadecheck,:)';
    body_edges = n_sources*R_body;
    d = (body_edges - o);
    d = d./(vecnorm(d')');
    % d = dbody_b_1'; % Ray unit vector out to the sun
    
    [~,~,hit_check(ii,:)] = ray_elem_intersect_R04(o,d,p,n,tc,tr); % Perform element hit check
end
unshadowed_b_body = ~boolean(sum(hit_check,1)); % 0 = in shadow, 1 = not in shadow


% Elements that are unshadowed, have active-sn, have active-role (i.e. elements that should have Qdot)
elems2calc_t_body = elems_t_body_shadecheck(boolean(unshadowed_t_body));
elems2calc_b_body = elems_b_body_shadecheck(boolean(unshadowed_b_body));

% Bug fixing
% elems2calc_t_body = elems_t_body_shadecheck;
% elems2calc_b_body = elems_b_body_shadecheck;

%%% Heating calculation - Albedo
vSBody = v_body-vS;
reflection_angle = vector_angles_R03(v_body,vSBody);
Qdot_t_albedo = S*AF*alpha(elems2calc_t_body,1).*elem_body_t_VF(elems2calc_t_body).*elem_areas(elems2calc_t_body)*abs(cos(reflection_angle))*(eclipse~=1);
Qdot_b_albedo = S*AF*alpha(elems2calc_b_body,2).*elem_body_b_VF(elems2calc_b_body).*elem_areas(elems2calc_b_body)*abs(cos(reflection_angle))*(eclipse~=1);

%%% Heating calculation - IR
Qdot_t_IR = SBC*(T_body^4)*epsilon(elems2calc_t_body,1).*elem_body_t_VF(elems2calc_t_body).*elem_areas(elems2calc_t_body);
Qdot_b_IR = SBC*(T_body^4)*epsilon(elems2calc_b_body,2).*elem_body_b_VF(elems2calc_b_body).*elem_areas(elems2calc_b_body);


%% Qdot Sum
Qdot_total_S = zeros(height(enormals),1);
Qdot_total_A = zeros(height(enormals),1);
Qdot_total_I = zeros(height(enormals),1);
for ii = 1:length(elems2calc_t_solar)
    elem = elems2calc_t_solar(ii);
    Qdot_total_S(elem) = Qdot_total_S(elem) + Qdot_t_solar(ii); % Add to current value
end
for ii = 1:length(elems2calc_b_solar)
    elem = elems2calc_b_solar(ii);
    Qdot_total_S(elem) = Qdot_total_S(elem) + Qdot_b_solar(ii); % Add to current value
end
for ii = 1:length(elems2calc_t_body)
    elem = elems2calc_t_body(ii);
    Qdot_total_A(elem) = Qdot_total_A(elem) + Qdot_t_albedo(ii); % Add to current value
end
for ii = 1:length(elems2calc_t_body)
    elem = elems2calc_t_body(ii);
    Qdot_total_I(elem) = Qdot_total_I(elem) + Qdot_t_IR(ii); % Add to current value
end
for ii = 1:length(elems2calc_b_body)
    elem = elems2calc_b_body(ii);
    Qdot_total_A(elem) = Qdot_total_A(elem) + Qdot_b_albedo(ii); % Add to current value
end
for ii = 1:length(elems2calc_b_body)
    elem = elems2calc_b_body(ii);
    Qdot_total_I(elem) = Qdot_total_I(elem) + Qdot_b_IR(ii); % Add to current value
end

% Bug checking
% Qdot_total_S = Qdot_total_S/4;
% Qdot_total_B = Qdot_total_B/4;

end
