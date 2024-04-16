function CONDUCTANCES_2 = contactor_create_R01(ID,source_input,source_type,target_input,target_type,cond,cond_type,area,range_tolerance,priority_switch,MESHGRIDS_1,SURFACES_1,SURFACES_2,CONDUCTANCES_2)
% Creates a contact conductance between two surfaces or edges. This value
% is partitioned out to the participating nodes. Turn on priority_switch to
% force each source node to seek out only the closest target node within
% tolerance range. Area value uses the following logic:
% area = -2: use area of target surface
% area = -1: use area of source surface
% area = 0: find and use smaller of the two areas
% area > 0: use this specified area in any calculations requiring it

% Result is CONDUCTANCES 2 which holds: [node 1, node 2, conductance,
% length between nodes].

% Version 1.0 completed 12/11/2023


%% Pre-check
if source_input == target_input
    fprintf('ERROR (contactor_create): Source and target cannot be the same.\n')
    return
end

%% Surface existence check
if source_type == 30 || source_type == 31 || source_type == 32 || source_type == 33 || source_type == 34
    for ii = 1:length(source_input)
        if isempty(SURFACES_2{source_input(ii)})
            fprintf('ERROR (contactor_create): Surface %i does not exist.\n',source_input(ii))
            return
        end
    end
end

if target_type == 30 || target_type == 31 || target_type == 32 || target_type == 33 || target_type == 34
    for ii = 1:length(target_input)
        if isempty(SURFACES_2{target_input(ii)})
            fprintf('ERROR (contactor_create): Surface %i does not exist.\n',target_input(ii))
            return
        end
    end
end

%% Source acquisition
n_x1 = SURFACES_2{source_input}.n_x1;
n_x2 = SURFACES_2{source_input}.n_x2;
surf_nodes = SURFACES_2{source_input}.nodes;
if source_type == 30 % Surface
    source_nodes = surf_nodes;
    source_area = SURFACES_2{source_input}.area;
elseif source_type == 31 % South edge
    source_nodes = surf_nodes(1:n_x1);
    source_area = SURFACES_1{source_input}.thickness*n_x1;
elseif source_type == 32 % East edge
    source_nodes = zeros(n_x2,1);
    for ii = 1:n_x2
        source_nodes(ii) = surf_nodes(jj*n_x1);
    end
    source_area = SURFACES_1{source_input}.thickness*n_x2;
elseif source_type == 33 % North edge
    source_nodes = surf_nodes(n_x1*n_x2-n_x1+1:end);
    source_area = SURFACES_1{source_input}.thickness*n_x1;
elseif source_type == 34 % West edge
    source_nodes = zeros(n_x2,1);
    for ii = 1:n_x2
        source_nodes(ii) = surf_nodes(ii*n_x1-n_x1+1);
    end
    source_area = SURFACES_1{source_input}.thickness*n_x2;
else % Error handling
    fprintf(['ERROR (contactor_create): Source_type does not exist. ' ...
        'Viable types are 30 (surf), 31 (south edge), ' ...
        '32 (east edge), 33 (north edge), 34 (west edge).\n'])
    return
end

%% Target acquisition
n_x1 = SURFACES_2{target_input}.n_x1;
n_x2 = SURFACES_2{target_input}.n_x2;
surf_nodes = SURFACES_2{target_input}.nodes;
if target_type == 30 % Surface
    target_nodes = surf_nodes;
    target_area = SURFACES_2{target_input}.area;
elseif target_type == 31 % South edge
    target_nodes = surf_nodes(1:n_x1);
    target_area = SURFACES_1{target_input}.thickness*n_x1;
elseif target_type == 32 % East edge
    target_nodes = zeros(n_x2,1);
    for ii = 1:n_x2
        target_nodes(ii) = surf_nodes(jj*n_x1);
    end
    target_area = SURFACES_1{target_input}.thickness*n_x2;
elseif target_type == 33 % North edge
    target_nodes = surf_nodes(n_x1*n_x2-n_x1+1:end);
    target_area = SURFACES_1{target_input}.thickness*n_x1;
elseif target_type == 34 % West edge
    target_nodes = zeros(n_x2,1);
    for ii = 1:n_x2
        target_nodes(ii) = surf_nodes(ii*n_x1-n_x1+1);
    end
    target_area = SURFACES_1{target_input}.thickness*n_x2;
else % Error handling
    fprintf(['ERROR (contactor_create): target_type does not exist. ' ...
        'Viable types are 30 (surf), 31 (south edge), ' ...
        '32 (east edge), 33 (north edge), 34 (west edge).\n'])
    return
end

%% Area handling
if area > 0
    faying_area = area;
elseif area == 0
    if source_area < target_area
        faying_area = source_area;
    else
        faying_area = target_area;
    end
elseif area == -1
    faying_area = source_area;
elseif area == -2
    faying_area = target_area;
else
    fprintf('ERROR (contactor_create): Area input issue. Please input a different value.\n')
    return
end

%% Conduction type handling
if cond_type == 2 % Area conductance, h
    G = cond*faying_area;
elseif cond_type == 3 % Total conductance, G
    G = cond;
else % Error handling
    fprintf('ERROR (contactor_create): Cond_type not accepted. Use 2 for area conductance [W/m2] or 3 for total conductance [W/K].')
end

%% Connection finding
% Find nodal coords
source_coords = MESHGRIDS_1(source_nodes,4:6);
target_coords = MESHGRIDS_1(target_nodes,4:6);

% Check each source node against the target nodes
dist = custom_pdist2(source_coords, target_coords);

% Check against range tolerance
[u,v] = find(dist <= range_tolerance);
w = zeros(length(u),1);
for ii = 1:length(u)
    w(ii) = dist(u(ii),v(ii));
end
u = source_nodes(u); % Source nodes numbers in range
v = target_nodes(v); % Target nodes numbers in range
dup_check = u ~= v;
u = u(dup_check); v = v(dup_check); w = w(dup_check); % Remove self-self connections if they exist
G_split = G/length(u); % Conductance value to be assigned to each connection
result = [u,v,G_split*ones(length(u),1),w];

if priority_switch == 1 % Link only to shortest node pair that can be found
    % Find unique values of first column
    [unique_vals, ~, ~] = unique(result(:,1), 'stable');
    
    % Run through each unique value and compare to the duplicates
    rows2keep = [];
    for ii = 1:length(unique_vals)
        % Find indices of duplicates
        dup_indices = find(result(:,1)==unique_vals(ii));
    
        % Compare fourth value of unique val and the duplicates
        dist_vals = result(dup_indices,4);
        [~,I] = min(dist_vals); % Find the subindex of the minimum of the duplicates distances amongst the duplicates
        min_index = dup_indices(I); % Get the actual index of the row
        rows2keep = [rows2keep;min_index];
    end
    result = result(rows2keep,:);
end

% Add the diagonal flip
result = [result; result(:,2),result(:,1),result(:,3:4)];

CONDUCTANCES_2{ID} = result; % Assign to contactor cell array


end


function distances = custom_pdist2(X, Y)
    % Check input dimensions
    [m, n] = size(X);
    [p, q] = size(Y);

    if n ~= q
        error('Input matrices must have the same number of columns.');
    end

    % Compute pairwise distances using Euclidean distance formula
    distances = sqrt(sum((reshape(X, m, 1, n) - reshape(Y, 1, p, n)).^2, 3));
end
