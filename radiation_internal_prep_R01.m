function out = radiation_internal_prep_R01(g,ecoords,ray_count,rotms,epsilon,cutoff)
% Prepares element connections for internal radiation of each radiation
% group.

% g = struct giving information of the current rad group
% ecoords = coordinates of each node in the element
% ray_count = number of rays to be shot for the monte-carlo ray tracing
% (used for view factor calculations)
% epsilon = Nx2 matrix giving IR emissivity/absorptivity values for the top
% and bottom sides of each element N.
% cutoff = tolerance that kills reflection function once an emitted ray has
% negligible energy

% Get elements
elems = g.elems; % ID of elements in the radiation group
% N = length(elems); % Number of elements in the radiation group

% Get epsilon for only the group elements
eps = epsilon(elems,:);
eps_tb = [eps(:,1);eps(:,2)]'; % Concatenate top and bottom into one vector
% rho_tb = 1-eps_tb; % Get reflectivities

% Get view factors
[VF_t,VF_b,FB_t,FB_b] = internal_view_factors_R02(g,ecoords,ray_count,rotms);
% VF_t = view factor that the topside of each element in the group has of each element in the group, including itself.
% VF_b = view factor that the bottomside of each element in the group has of each element in the group, including itself.
% FB_t = (frontback) indicates if the view factor in VF_f is to the topside of an element (1) or the bottomside (0) or itself (0)
% FB_b = (frontback) indicates if the view factor in VF_b is to the topside of an element (1) or the bottomside (0) or itself (0)

assignin('base','VF_t',VF_t);
assignin('base','VF_b',VF_b);
assignin('base','FB_t',FB_t);
assignin('base','FB_b',FB_b);
% fprintf('View factors topsides of each element (uncaring of seeing top or bottom)\n')
% checktop = [elems,sum(VF_t,2)]
% checktop_sum = sum(checktop)
% fprintf('View factors bottomsides of each element (uncaring of seeing top or bottom)\n')
% checkbot = [elems,sum(VF_b,2)]
% checkbot_sum = sum(checkbot)

% Find only active elements that give and take radiation
tactive = g.trole == 1; bactive = g.brole == 1;
active = [tactive,bactive];

% Get matrices for topside viewing and bottomside viewing
% Four matrices of the same size each
F_tt = VF_t.*FB_t; % NxN matrix (same size as the VF matrices but now with zeros for target bottomsides). Shows the view of each element topside to the topside of other elements.
F_tb = VF_t.*(~FB_t); % NxN matrix (same size as the VF matrices but now with zeros for target topsides). Shows the view of each element topside to the bottomside of other elements.
F_bt = VF_b.*FB_b; % NxN matrix (same size as the VF matrices but now with zeros for target bottomsides). Shows the view of each element bottomside to the topside of other elements.
F_bb = VF_b.*(~FB_b); % NxN matrix (same size as the VF matrices but now with zeros for target topsides). Shows the view of each element bottomside tot he bottomside of other elements.

% Handle emissivity and reflectivity values
% F_tt_eps = F_tt.*repmat(eps(:,1)',N,1); % Get topside epsilons (Nx1) and rotate to match the columns (1xN), then repeat to create NxN matrix, then multiply by element.
% F_tb_eps = F_tb.*repmat(eps(:,2)',N,1);
% F_bt_eps = F_tt.*repmat(eps(:,1)',N,1);
% F_bb_eps = F_tb.*repmat(eps(:,2)',N,1);
% F_tt_rho = F_tt.*(1-repmat(eps(:,1)',N,1)); % Get topside epsilons (Nx1) and rotate to match the columns (1xN), then repeat to create NxN matrix, then subtract from 1 to get reflectivities, then multiply by element.
% F_tb_rho = F_tb.*(1-repmat(eps(:,2)',N,1));
% F_bt_rho = F_tt.*(1-repmat(eps(:,1)',N,1));
% F_bb_rho = F_tb.*(1-repmat(eps(:,2)',N,1));

% F_eps = [F_tt_eps, F_tb_eps; F_bt_eps, F_bb_eps]; % Concatenate to make emissivity matrix
% F_rho = [F_tt_rho, F_tb_rho; F_bt_rho, F_bb_rho]; % Concatenate to make reflectivity matrix
F = [F_tt,F_tb;F_bt,F_bb]; % Concatenate to get 2Nx2N matrix

% Run through reflections to get proportionality matrix
RPM = rad_proportions_R01(eps_tb,F,cutoff);
out = RPM(active,active); % Prune out inactive elements

end