function CONDUCTANCES_2 = conductor_create_R02(ID,source,target,target_type,cond,cond_type,area,split,MESHGRIDS_1,ELEMENTS_1,SURFACES_1,SURFACES_2,CONDUCTANCES_1,CONDUCTANCES_2)
% Creates a conductance between a node and another entity (node, element,
% edge, surface).

% Source must be a node
% Target can be a node, an element, an edge, or a surface

% ID - ID of conductor
% source - Source node number
% target - Target ID
% target_type - Specify if the target is a node (1), an element (2), a surface (30), or a surface edge (31, 32, 33, 34)
% cond - conduction value
% cond_type - Specify what type of conduction is used (1-conductivity (k), 2-thermal/area conductance (h), 3-total conductance (G))
% area - Specify area to calculate total conductance. Required for cond_type 1 and 2. Input 0 or less for targets of type 2+ to have the function calculate this based on element or surface area/thickness.
% split - Specify if conduction value is split between targets (1) or is per node in the target (0)

% Version 2.0 completed 11/16/2023

if split == 1 && length(target) == 1
    fprintf('Warning: Split is set to 1 while there is only a single target.\n')
end

if split == 1 % Split divides cond value by number of targets, not necessarily number of nodes.
    cond = cond/length(target);
elseif split ~= 0 && split ~=1
    fprintf('Invalid split input. Please choose 0 (no split) or 1 (split cond value between the target nodes.\n')
    return
end

if cond_type ~= 1 && cond_type ~= 2 && cond_type ~= 3
    fprintf('Invalid conduction type. Please input 1, 2, or 3 for cond_type with appropriate ancillary information.\n');
    return
end

if target_type == 1 && cond_type ~= 3 && area == 0
    fprintf('Warning: Area set to zero for a node-node connection without given total conductance. The conductance value will be zero.\n')
end

if cond_type == 3 && area ~= 0
    fprintf('Note: Area given but not used. cond_type is set to 3, Total Conductance.\n')
end

if (target_type == 2 || target_type == 31 || target_type == 32 || target_type == 33 || target_type == 33) && area > 0 && cond_type ~= 3
    fprintf('Note: Area given and used in priority over existing element or surface area. Input area=0 to use the existing area.\n')
end

if target_type == 1 % Target is node(s)
    mat1 = zeros(length(target),3);
    mat2 = zeros(length(target),3);
    for ii = 1:length(target)
        % Check target
        if source == target(ii)
            fprintf('Target node cannot be the same as the source node.\n')
            return
        end
        if height(MESHGRIDS_1) < target(ii) || target(ii) < 0
        %if isempty(isempty(find(MESHGRIDS_1(:,1)==target(ii), 1)))
            fprintf('Target node %i does not exist. Select a different target node.\n',target(ii));
            return
        end

        if cond_type == 1 % Thermal conductivity
            L = norm(MESHGRIDS_1(source,4:6)-MESHGRIDS_1(target(ii),4:6)); % Calculate length
            G = cond*area/L; % Calculate total conductance
        elseif cond_type == 2 % Thermal area conductance
            G = cond*area; % Calculate total conductance
        elseif cond_type == 3 % Total thermal conductance
            G = cond;
        end
        % CONDUCTANCES_1(source,target(ii)) = G;
        % CONDUCTANCES_1(target(ii),source) = G;
        mat1(ii,:) = [source,target(ii),G];
        mat2(ii,:) = [target(ii),source,G];
    end
    mat = [mat1;mat2];
elseif target_type == 2 % Target is element(s)
    % Check target
    for ii = 1:length(target)
        if target(ii) > height(ELEMENTS_1)
            fprintf('Target element %i does not exist. Select a different target element.\n',target(ii));
            return
        end
    end

    areas = zeros(length(target));
    for ii = 1:length(target) % Iterate through target elements
        if area <= 0 % Calculate area of element if needed. If area is given, then this overrides area of all selected elements.
            nodes = ELEMENTS_1(target(ii),2:5);
            if sum(find(nodes==source)) > 0 % Check if source is within element
                fprintf('Source node must not be within the target set.\n')
                return
            end
            vec_1 = MESHGRIDS_1(nodes(2),4:6)-MESHGRIDS_1(nodes(1),4:6);
            vec_2 = MESHGRIDS_1(nodes(4),4:6)-MESHGRIDS_1(nodes(1),4:6);
            areas(ii) = norm(cross(vec_1,vec_2));
        else
            areas(ii) = area;
        end
    end
    
    COND_temp = [];
    if cond_type == 1 % Thermal conductivity
        for ii = 1:length(target) % Iterate through target elements
            nodes = ELEMENTS_1(target(ii),2:5);
            for jj = 1:length(nodes) % Iterate through nodes
                L = norm(MESHGRIDS_1(source,4:6)-MESHGRIDS_1(nodes(jj),4:6)); % Calculate length
                G = cond*areas(ii)/L; % Calculate total conductance
                COND_temp = [COND_temp; nodes(jj), G]; % Assign to temporary array
            end
        end
    elseif cond_type == 2 % Thermal area conductance
        for ii = 1:length(target) % Iterate through target elements
            nodes = ELEMENTS_1(target(ii),2:5);
            for jj = 1:length(nodes) % Iterate through nodes
                G = cond*areas(ii);
                COND_temp = [COND_temp; nodes(jj), G/length(nodes)]; % Assign to temporary array
            end
        end
    elseif cond_type == 3 % Total thermal conductance
        for ii = 1:length(target) % Iterate through target elements
            nodes = ELEMENTS_1(target(ii),2:5);
            for jj = 1:length(nodes) % Iterate through nodes
                G = cond;
                COND_temp = [COND_temp; nodes(jj), G/length(nodes)]; % Assign to temporary array
            end
        end
    end
    
    % Find duplicates and select greater value
    sorted_matrix = sortrows(COND_temp, -2); % Sort the input matrix by the second column in descending order.
    [~, unique_indices, ~] = unique(sorted_matrix(:,1), 'stable'); % Find the indices of the duplicate values in the first column of the sorted matrix.
    COND_temp = sorted_matrix(unique_indices, :); % Keep only the duplicate rows with the highest corresponding value in the second column.
    %fprintf('Cond_temp = %s\n',mat2str(COND_temp));

    % Assign to matrix
    mat1 = zeros(height(COND_temp),3);
    mat2 = zeros(height(COND_temp),3);
    for ii = 1:height(COND_temp)
        % CONDUCTANCES_1(source,COND_temp(ii,1)) = COND_temp(ii,2);
        % CONDUCTANCES_1(COND_temp(ii,1),source) = COND_temp(ii,2);
        mat1(ii,:) = [source,COND_temp(ii,1),G];
        mat2(ii,:) = [COND_temp(ii,1),source,G];
    end
    mat = [mat1;mat2];

elseif target_type == 30 % Target is surface(s)

    areas = zeros(length(target),1);
    for ii = 1:length(target)
        current_target = target(ii);
        if area <= 0 % Calculate area of surface if needed. If area is given, then this overrides area of all selected elements.
            areas(ii) = SURFACES_2{current_target}.area; 
        else
            areas(ii) = area;
        end
    end

    % Iterate through targets
    for ii = 1:length(target)
        % Find nodes in target surface
        current_target = target(ii);
        nodes = SURFACES_2{current_target}.nodes;
        mat1 = [];
        mat2 = [];
        if cond_type == 1 % Thermal conductivity
            for jj = 1:length(nodes)
                L = norm(MESHGRIDS_1(source,4:6)-MESHGRIDS_1(nodes(jj),4:6)); % Calculate length
                G = cond*areas(ii)/L/length(nodes); % Calculate total conductance of each node in target
                % CONDUCTANCES_1(source,nodes(jj)) = G;
                % CONDUCTANCES_1(nodes(jj),source) = G;
                mat1 = [mat1; source,nodes(jj),G];
                mat2 = [mat2; nodes(jj),source,G];
            end
        elseif cond_type == 2 % Thermal area conductance
            cond = cond*areas(ii); % Convert h to G based on area of target
            for jj = 1:length(nodes)
                G = cond/length(nodes);
                % CONDUCTANCES_1(source,nodes(jj)) = cond/length(nodes);
                % CONDUCTANCES_1(nodes(jj),source) = cond/length(nodes);
                mat1 = [mat1; source,nodes(jj),G];
                mat2 = [mat2; nodes(jj),source,G];
            end
        elseif cond_type == 3 % Total thermal conductance
            for jj = 1:length(nodes)
                G = cond/length(nodes);
                % CONDUCTANCES_1(source,nodes(jj)) = cond/length(nodes);
                % CONDUCTANCES_1(nodes(jj),source) = cond/length(nodes);
                mat1 = [mat1; source,nodes(jj),G];
                mat2 = [mat2; nodes(jj),source,G];
            end
        end
    end
    mat = [mat1;mat2];
elseif target_type == 31 % Target is surface(s) south edge

    areas = zeros(length(target));
    for ii = 1:length(target) % Iterate through target edges
        current_target = target(ii);
        if area <= 0 % Calculate area of element if needed. If area is given, then this overrides area of all selected elements.
            areas(ii) = norm(SURFACES_2{current_target}.vec_x1)*SURFACES_1{current_target}.thickness;
        else
            areas(ii) = area;
        end
    end
    
    % Find nodes
    edge_nodes = [];
    n_x1 = zeros(length(target),1);
    for ii = 1:length(target)
        current_target = target(ii);
        surf_nodes = SURFACES_2{current_target}.nodes;
        n_x1(ii) = SURFACES_2{current_target}.n_x1;
        edge_nodes = [edge_nodes;surf_nodes(1:n_x1), areas(ii)/n_x1(ii)*ones(n_x1(ii),1), n_x1(ii)*ones(n_x1(ii),1)]; % Includes area for each node in second column
    end

    % Calculate conductance for each node
    if cond_type == 1
        G = zeros(height(edge_nodes),1);
        for ii = 1:length(edge_nodes)
            L = norm(MESHGRIDS_1(source,4:6)-MESHGRIDS_1(edge_nodes(ii,1),4:6)); % Calculate length
            G(ii) = cond*edge_nodes(ii,2)/L;
        end
    elseif cond_type == 2
        G = zeros(height(edge_nodes),1);
        for ii = 1:length(edge_nodes)
            G(ii) = cond*edge_nodes(ii,2);
        end
    elseif cond_type == 3
        G = zeros(height(edge_nodes),1);
        for ii = 1:length(edge_nodes)
            G(ii) = cond/edge_nodes(ii,3);
        end
    end
    %disp(G)

    COND_temp = [edge_nodes(:,1), G]; % Concatenate into one matrix
    
    % Find duplicates and select greater value
    sorted_matrix = sortrows(COND_temp, -2); % Sort the input matrix by the second column in descending order.
    [~, unique_indices, ~] = unique(sorted_matrix(:,1), 'stable'); % Find the indices of the duplicate values in the first column of the sorted matrix.
    COND_temp = sorted_matrix(unique_indices, :); % Keep only the duplicate rows with the highest corresponding value in the second column.
    
    % Assign to matrix
    mat1 = zeros(length(nodes),3);
    mat2 = zeros(length(nodes),3);
    for ii = 1:height(COND_temp)
        % CONDUCTANCES_1(source,COND_temp(ii,1)) = COND_temp(ii,2);
        % CONDUCTANCES_1(COND_temp(ii,1),source) = COND_temp(ii,2);
        mat1(ii) = [source,nodes(jj),G];
        mat2(ii) = [nodes(jj),source,G];
    end
    mat = [mat1;mat2];

elseif target_type == 32 % Target is surface(s) east edge
    areas = zeros(length(target));
    for ii = 1:length(target) % Iterate through target edges
        current_target = target(ii);
        if area <= 0 % Calculate area of surface if needed. If area is given, then this overrides area of all selected elements.
            areas(ii) = norm(SURFACES_2{current_target}.vec_x2)*SURFACES_1{current_target}.thickness;
        else
            areas(ii) = area;
        end
    end
    
    % Find nodes
    edge_nodes = [];
    n_x2 = zeros(length(target),1);
    for ii = 1:length(target)
        current_target = target(ii);
        surf_nodes = SURFACES_2{current_target}.nodes;
        n_x1 = SURFACES_2{current_target}.n_x1;
        n_x2(ii) = SURFACES_2{current_target}.n_x2;
        for jj = 1:n_x2
            edge_nodes = [edge_nodes;surf_nodes(jj*n_x1), areas(ii)/n_x2(ii), n_x2(ii)];
        end
    end

    % Calculate conductance for each node
    if cond_type == 1
        G = zeros(height(edge_nodes),1);
        for ii = 1:length(edge_nodes)
            L = norm(MESHGRIDS_1(source,4:6)-MESHGRIDS_1(edge_nodes(ii,1),4:6)); % Calculate length
            G(ii) = cond*edge_nodes(ii,2)/L;
        end
    elseif cond_type == 2
        G = zeros(height(edge_nodes),1);
        for ii = 1:length(edge_nodes)
            G(ii) = cond*edge_nodes(ii,2);
        end
    elseif cond_type == 3
        G = zeros(height(edge_nodes),1);
        for ii = 1:length(edge_nodes)
            G(ii) = cond/edge_nodes(ii,3);
        end
    end
        
    COND_temp = [edge_nodes(:,1), G]; % Concatenate into one matrix
    
    % Find duplicates and select greater value
    sorted_matrix = sortrows(COND_temp, -2); % Sort the input matrix by the second column in descending order.
    [~, unique_indices, ~] = unique(sorted_matrix(:,1), 'stable'); % Find the indices of the duplicate values in the first column of the sorted matrix.
    COND_temp = sorted_matrix(unique_indices, :); % Keep only the duplicate rows with the highest corresponding value in the second column.

    % Assign to matrix
    mat1 = zeros(length(nodes),3);
    mat2 = zeros(length(nodes),3);
    for ii = 1:height(COND_temp)
        % CONDUCTANCES_1(source,COND_temp(ii,1)) = COND_temp(ii,2);
        % CONDUCTANCES_1(COND_temp(ii,1),source) = COND_temp(ii,2);
        mat1(ii) = [source,nodes(jj),G];
        mat2(ii) = [nodes(jj),source,G];
    end
    mat = [mat1;mat2];

elseif target_type == 33 % Target is surface(s) north edge
    areas = zeros(length(target));
    for ii = 1:length(target) % Iterate through target edges
        current_target = target(ii);
        if area <= 0 % Calculate area of element if needed. If area is given, then this overrides area of all selected elements.
            areas(ii) = norm(SURFACES_2{current_target}.vec_x1)*SURFACES_1{current_target}.thickness;
        else
            areas(ii) = area;
        end
    end
    
    % Find nodes
    edge_nodes = [];
    n_x1 = zeros(length(target),1);
    for ii = 1:length(target)
        current_target = target(ii);
        surf_nodes = SURFACES_2{current_target}.nodes;
        n_x1(ii) = SURFACES_2{current_target}.n_x1;
        n_x2 = SURFACES_2{current_target}.n_x2;
        edge_nodes = [edge_nodes;surf_nodes(n_x1*n_x2-n_x1+1:end), areas(ii)/n_x1(ii)*ones(n_x1(ii),1), n_x1(ii)*ones(n_x1(ii),1)];
    end

    % Calculate conductance for each node
    if cond_type == 1
        G = zeros(height(edge_nodes),1);
        for ii = 1:length(edge_nodes)
            L = norm(MESHGRIDS_1(source,4:6)-MESHGRIDS_1(edge_nodes(ii,1),4:6)); % Calculate length
            G(ii) = cond*edge_nodes(ii,2)/L;
        end
    elseif cond_type == 2
        G = zeros(height(edge_nodes),1);
        for ii = 1:length(edge_nodes)
            G(ii) = cond*edge_nodes(ii,2);
        end
    elseif cond_type == 3
        G = zeros(height(edge_nodes),1);
        for ii = 1:length(edge_nodes)
            G(ii) = cond/edge_nodes(ii,3);
        end
    end
    
    COND_temp = [edge_nodes(:,1), G]; % Concatenate into one matrix
    
    % Find duplicates and select greater value
    sorted_matrix = sortrows(COND_temp, -2); % Sort the input matrix by the second column in descending order.
    [~, unique_indices, ~] = unique(sorted_matrix(:,1), 'stable'); % Find the indices of the duplicate values in the first column of the sorted matrix.
    COND_temp = sorted_matrix(unique_indices, :); % Keep only the duplicate rows with the highest corresponding value in the second column.

    % Assign to matrix
    mat1 = zeros(length(nodes),3);
    mat2 = zeros(length(nodes),3);
    for ii = 1:height(COND_temp)
        % CONDUCTANCES_1(source,COND_temp(ii,1)) = COND_temp(ii,2);
        % CONDUCTANCES_1(COND_temp(ii,1),source) = COND_temp(ii,2);
        mat1(ii) = [source,nodes(jj),G];
        mat2(ii) = [nodes(jj),source,G];
    end
    mat = [mat1;mat2];


elseif target_type == 34 % Target is surface(s) west edge
    areas = zeros(length(target));
    for ii = 1:length(target) % Iterate through target surfaces
        current_target = target(ii);
        if ~isempty(SURFACES_1{current_target})
            if area <= 0 % Calculate area of element if needed. If area is given, then this overrides area of all selected elements.
                areas(ii) = norm(SURFACES_2{current_target}.vec_x2)*SURFACES_1{current_target}.thickness;
            else
                areas(ii) = area;
            end
        end
    end
    
    % Find nodes
    edge_nodes = [];
    n_x2 = zeros(length(target),1);
    for ii = 1:length(target)
        current_target = target(ii);
        surf_nodes = SURFACES_2{current_target}.nodes;
        n_x1 = SURFACES_2{current_target}.n_x1;
        n_x2(ii) = SURFACES_2{current_target}.n_x2;
        for jj = 1:n_x2
            edge_nodes = [edge_nodes;surf_nodes(jj*n_x1-n_x1+1), areas(ii)/n_x2(ii),n_x2(ii)];
        end
    end
    %disp(edge_nodes)

    % Calculate conductance for each node
    if cond_type == 1
        G = zeros(height(edge_nodes),1);
        for ii = 1:length(edge_nodes)
            L = norm(MESHGRIDS_1(source,4:6)-MESHGRIDS_1(edge_nodes(ii,1),4:6)); % Calculate length
            G(ii) = cond*edge_nodes(ii,2)/L;
        end
    elseif cond_type == 2
        G = zeros(height(edge_nodes),1);
        for ii = 1:length(edge_nodes)
            G(ii) = cond*edge_nodes(ii,2);
        end
    elseif cond_type == 3
        G = zeros(height(edge_nodes),1);
        for ii = 1:length(edge_nodes)
            G(ii) = cond/edge_nodes(ii,3);
        end
    end
        
    COND_temp = [edge_nodes(:,1), G]; % Concatenate into one matrix
    
    % Find duplicates and select greater value
    sorted_matrix = sortrows(COND_temp, -2); % Sort the input matrix by the second column in descending order.
    [~, unique_indices, ~] = unique(sorted_matrix(:,1), 'stable'); % Find the indices of the duplicate values in the first column of the sorted matrix.
    COND_temp = sorted_matrix(unique_indices, :); % Keep only the duplicate rows with the highest corresponding value in the second column.

    % Assign to matrix
    mat1 = zeros(length(nodes),3);
    mat2 = zeros(length(nodes),3);
    for ii = 1:height(COND_temp)
        % CONDUCTANCES_1(source,COND_temp(ii,1)) = COND_temp(ii,2);
        % CONDUCTANCES_1(COND_temp(ii,1),source) = COND_temp(ii,2);
        mat1(ii) = [source,nodes(jj),G];
        mat2(ii) = [nodes(jj),source,G];
    end
    mat = [mat1;mat2];

else
    fprintf('Invalid target type. Please input 1 for node(s), 2 for element(s), 30 for surface(s),');
    fprintf('31 for surface(s) south edge, 32 for surface(s) east edge, 33 for surface(s) north edge, 34 for surface(s) west edge.\n');
    return
end

CONDUCTANCES_2{ID} = mat;

end