function [VF_t,VF_b,FB_t,FB_b] = internal_view_factors_R02(g,ecoords,ray_count,rotms)
% Calculates the view factor between elements in a radiation group. Output
% is a matrix giving the view factor of the source element (row) of the
% target element (column).

% g = struct with the following information (group information)
% .elems = list of element IDs in this group
% .tsn = boolean vector of length(.elems) indicating which elems have sn-active topsides
% .bsn = boolean vector of length(.elems) indicating which elems have sn-active bottomsides
% .trole = vector of length(.elems) indicating the role of each element topside
% .brole = vector of length(.elems) indicating the role of each element bottomside
% .centers = Nx3 matrix of coordinates of each element center

% ecoords = coordinates of all elements
% ray_count = number of rays to be shot
% rotms = cell array of rotation matrices that rotate from s/c body frame to local element reference frame

% This version creates a VF matrix with front and back sides of each
% element handled separately. This allows for the radiation proportionality
% matrix to handle elements with sides of different emissivities.

% Version 2.0 completed 2/21/2024

elems1 = g.elems;
N1 = length(elems1); % Get number of elements in the rad group

% Pare down list of elements to only elements that are active or blocking
role_check = g.trole == 1 | g.trole == 2 | g.brole == 1 | g.brole == 2;
% elems2 = elems1(role_check);
% N2 = length(elems2); % Get number of elements to check
% gcoords = ecoords(elems2); % Get coordinates of only the elements in the rad group that are active or blocking
% grotms = rotms{elems2};


% Iterate through list of elements in the group
VF_t = zeros(N1); VF_b = zeros(N1);
FB_t = zeros(N1); FB_b = zeros(N1);
for ii = 1:N1
    selem = elems1(ii); % Get current element ID
    if g.trole(ii) == 1 || g.brole(ii) == 1 % Only perform calcs for elements that radiate out
        
        scoords = ecoords{selem}; % Get coordinates of the source element
        srotm = rotms{selem}; % Get the rotation matrix of the source element

        telems = elems1; telems(ii) = NaN; 
        telems = telems(role_check); % Keep only elements that can block or absorb
        telems(isnan(telems)) = []; % Exclude source element
        telems_indices = 1:N1; telems_indices(ii) = NaN;
        telems_indices = telems_indices(role_check);
        telems_indices(isnan(telems_indices)) = []; % Get list of indices relating telems to list of elems in rad group

        tcoords = ecoords(telems); % Get coordinates of the target elements (elements that can possibly block or absorb radiation)
        trotms = rotms(telems); % Get the rotation matrices of the target elements
        % tcoords(ii) = []; % Remove coordinates of the source element from the list of target elements coords
        % trotms = grotms; % Get the rotation matrices of the target elements
        % trotms(ii) = []; % Remove the rotation matrix of the source element from the list of target element rotms

        [vf_t, ~, ~, frontback_of_topside] = mcrt_R03(scoords,tcoords,ray_count,srotm,trotms,1); % Get source-topside view factors
        [vf_b, ~, ~, frontback_of_backside] = mcrt_R03(scoords,tcoords,ray_count,srotm,trotms,-1); % Get source-backside view factors
        % frontback indicates if the source element sees the front or the back
        % side of the target elements. This must be used with the condition of
        % the view factor (vf1/2) since it does not care if the view factor is
        % zero.
        vf_t(isnan(vf_t)) = 0; vf_b(isnan(vf_b)) = 0;
        vf_t(isinf(vf_t)) = 0; vf_b(isinf(VF_b)) = 0;
        
        for jj = 1:length(telems) % Iterate through list of outputs (same length as list of target elements)
            telem_index = telems_indices(jj); % Get index of current target element, relative to list of elements in rad group
            VF_t(ii,telem_index) = vf_t(jj);
            VF_b(ii,telem_index) = vf_b(jj);
            FB_t(ii,telem_index) = frontback_of_topside(jj);
            FB_b(ii,telem_index) = frontback_of_backside(jj);
        end

        % This sets VF as a matrix specifically with elements of the group and
        % handles self viewing
        % if ii > 1 && ii < N1
        %     tempvf1t = vf_t(1:ii); tempvf2t = vf_t(ii+1:end);
        %     tempvf1b = vf_b(1:ii); tempvf2b = vf_b(ii+1:end);
        %     tempfb1t = frontback_of_topside(1:ii); tempfb2t = frontback_of_topside(ii+1:end);
        %     tempfb1b = frontback_of_backside(1:ii); tempfb2b = frontback_of_backside(ii+1:end);
        %     VF_t(ii,:) = [tempvf1t, false, tempvf2t];
        %     VF_b(ii,:) = [tempvf1b, false, tempvf2b];
        %     FB_t(ii,:) = [tempfb1t, false, tempfb2t];
        %     FB_b(ii,:) = [tempfb1b, false, tempfb2b];
        % elseif ii == 1
        %     VF_t(ii,:) = [false, vf_t];
        %     VF_b(ii,:) = [false, vf_b];
        %     FB_t(ii,:) = [false, frontback_of_topside];
        %     FB_b(ii,:) = [false, frontback_of_backside];
        % elseif ii == N1
        %     VF_t(ii,:) = [vf_t, false];
        %     VF_b(ii,:) = [vf_b, false];
        %     FB_t(ii,:) = [frontback_of_topside, false];
        %     FB_b(ii,:) = [frontback_of_backside, false];
        % end
    end
end


end