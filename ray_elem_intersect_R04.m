function [t_out,r_out,hit_check] = ray_elem_intersect_R04(o,d,p,n,tc,tr)
% Checks if the launched rays hit the target element.

% o = origins of rays
% d = unit vectors of rays
% p = point on target element
% n = normal vector of target element
% tc = coordinates of target element
% tr = rotation matrix of target element (s/c frame to target frame)

[s1,s2] = size(d); % Get size of unit vector matrix (need to know how many rays are shot)
n = repmat(n,1,s2); % Repeat normal vector to form matrix of appropriate size
dotdn = dot(d,n); % Get dot product between each ray unit vector and the normal vector
parallelism_check = abs(dotdn) < 1e-6; % Check if too small to create threshold for numerical stability

t = dot(p-o,n)./dotdn; % Create t values for each ray (where each ray intersects plane of the target element)
tbool = t>0; % Check if greater t is in front of source element
%ptbool = ~parallelism_check & tbool; 
trep = repmat(t,s1,1); % Repeat t vector for easier elementwise multiplication
r = o + trep.*d; % Get r for each value t
r_2D = tr*r; % Transform position vector from s/c body frame to target element frame
r_2D(3,:) = []; % Remove third row
tc_2D = tr*tc; % Transform coordinates of target element from s/c body frame to target element frame
[in,on] = inpolygon(r_2D(1,:),r_2D(2,:),tc_2D(1,:),tc_2D(2,:)); % Check if points r are within the bounds of the target element
in_on = in | on; % Combine inside points with points on the edge

hit_check = (~parallelism_check & tbool) & in_on; % Create boolean to indicate which rays pass checks and hit the element
t_out = t;
t_out(~hit_check) = Inf; % Set failed rays
r_out = r;
r_out(:,~hit_check) = NaN; % Give bad values for bad rays
end