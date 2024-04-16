function [vf, o_final, r_final, frontback] = mcrt_R03(scoords,tcoords,ray_count,srotm,trotms,topbot)
% Performs Monte Carlo Ray Tracing (MCRT) to find the view factor between
% the source element and the target elements. This includes the shadowing
% of each target element upon the others.

%%% Inputs
% scoords = Source element coordinate matrix given as a 3x4 matrix with row
%           1 giving the x coordinates of each point, row 2 giving the y coordinates
%           of each poinot, and row 3 giving the z coordinates of each point. Points
%           are defined in order of lower left then ccw around the rectangular
%           element.
% tcoords = Target element coordinate cell array containing the coordinates
%           of each node in 3x4 matrix form/ 
%           Example: tcoords{elem_ID}.a = [xcoord;ycoord;zcoord]
% ray_count = Number of rays to be shot from the source element.
% srotms = Rotation matrix to rotate position vectors into the source
%          element's frame of reference.
% trotms = Cell array giving the rotation matrix to rotate position vectors
%          from s/c frame of reference into each target element's frame of 
%          reference. Same length as tcoords.
% topbot = Choose to shoot rays from topside or bottomside of the element
%          (1 or -1)

%%% Outputs
% vf = View factors of each target element.
% o_final = Coordinates of each ray origin that survives.
% r_final = Location of intersect of each ray with each blocker.
% frontback = Indicates which side is seen (1=front, 0=back). Must be used 
%             in conjunction with vf since it does not care if the view factor is zero.

% Version 1.1 completed 2/19/2024

%% Prep
A = norm(scoords(:,1) - scoords(:,2)); % Length - Node a to node b
B = norm(scoords(:,1) - scoords(:,4)); % Length - Node a to node d
tcount = length(tcoords);

%% Generate rays
% Ray origin point
u = A*rand(1,ray_count); % 1xray_count random u pos
v = B*rand(1,ray_count); % 1xray_count random v pos
s = [u;v;zeros(1,ray_count)]; % Position vectors of each source point in the reference frame of the source element.
o = srotm'*s + repmat(scoords(:,1),1,ray_count);
% o = bsxfun(@plus,srotm'*s,scoords(:,1)); % Rotate position vectors to s/c body frame. Add values to origin of source element to get position.

% Ray vector
angle1 = asin(sqrt(rand(1,ray_count))); % Get phi angle
angle2 = rand(1,ray_count)*2*pi; % Get theta angle
u_x = cos(angle2).*sin(angle1); % Unit vector x
u_y = sin(angle2).*sin(angle1); % Unit vector y
u_z = cos(angle1); % Unit vector z
ray_unit_vec = topbot*[u_x;u_y;u_z]./vecnorm([u_x;u_y;u_z]); % Ray unit vector in the reference frame of the source element.
d = srotm'*ray_unit_vec; % Rotate ray unit vector into the s/c body frame.
%d = reshape(d,3,ray_count); % Reshape to get 3xN matrix of unit vectors where each column is a different, random unit vector.


%% Ray hit/intersection checks
%%% Find t for each element
t_keep = zeros(tcount,ray_count);
r_keep = zeros(3,ray_count,tcount);
tn = zeros(tcount,3);
parfor ii = 1:tcount
    tc = tcoords{ii}; % Get coordinates of target element
    p = tc(:,1); % Get a point on the target element
    n = cross(tc(:,2)-tc(:,1),tc(:,4)-tc(:,1)); % Get normal vector of the target element
    n = n/norm(n); % Make into unit vector
    tn(ii,:) = n'; % Save normal vector into Nx3 matrix
    tr = trotms{ii}; % Get rotation matrix of target element
    [t_out,r_out,~] = ray_elem_intersect_R04(o,d,p,n,tc,tr); % Perform element hit check
    %   t output in 1 x ray_count vector format with Infs at indices that fail or do not hit.
    %   r output in 3 x ray_count matrix format

    t_keep(ii,:) = t_out;
    r_keep(:,:,ii) = r_out;
end

%%% Check which has lower t
%   t_keep is a tcount x ray_count matrix
[val_check,min_target] = min(t_keep,[],1); % Identify which was the first target element hit for each ray
%[target_elems,~,~] = unique(min_target)

% Split out ray hits and count
hit_count = zeros(1,tcount);
parfor ii = 1:tcount
    look_for = ii; % Element to look for
    indices_bool = min_target == look_for; % Get indices of min_target that match current element
    vals_at_indices = val_check(indices_bool); % Get values of t for the rays hitting the current element
    vals_at_indices(isinf(vals_at_indices) | isnan(vals_at_indices)) = []; % Prune out Infs and NaNs
    hit_count(ii) = numel(vals_at_indices); % Get count of elements that are not Infs that hit the element
    indices_keep(ii,:) = indices_bool;
end
vf = hit_count./ray_count;

%% Visualization prep
o_final = cell(1,ray_count);
r_final = cell(1,ray_count);
parfor ii = 1:tcount
    indices_bool = indices_keep(ii,:); % Indices of relevant rays
    
    o_final{ii} = o(:,indices_bool); % Get origins of rays that hit each element, filtered to each element
    r_temp = r_keep(:,:,ii);
    r_final{ii} = r_temp(:,indices_bool);
end

%% See topside or bottomside?
% Check if the source element sees the topside or the bottomside of each
% target element

% Get source element normal vector
sx1 = scoords(:,2) - scoords(:,1);
sx2 = scoords(:,4) - scoords(:,1);
sn = topbot*cross(sx1,sx2);
sn = repmat(sn,1,tcount)';
inc_angles = vector_angles_R03(sn,tn);
frontback = (inc_angles > pi/2)';

end
