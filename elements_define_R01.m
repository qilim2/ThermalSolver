function ELEMENTS_1 = elements_define_R01(MESHGRIDS_1,SURFACES_2)
% Creates a matrix of elements agnostic of the surfaces.
ELEMENTS_1 = NaN*ones(height(MESHGRIDS_1),5);
current_elem = 1;

% Version 1.0 completed 10/2/2023

for current_surf = 1:length(SURFACES_2) % Iterate through all the surfaces to find the relevant nodes
    if ~isempty(SURFACES_2{current_surf}) % Check if the container is empty
        
        % Extract information from SURFACES_2 about the current surface
        relevant_nodes = SURFACES_2{current_surf}.nodes;
        n_x1 = SURFACES_2{current_surf}.n_x1;
        n_x2 = SURFACES_2{current_surf}.n_x2;
        
        % Identify node a of each element
        indices2remove = ones(n_x2-1,1); % Initialize without index in top row
        for i = 1:n_x2-1
            indices2remove(i) = i*n_x1;
        end
        relevant_indices = 1:length(relevant_nodes);
        relevant_indices(end-n_x1+1:end) = []; % Remove top row
        relevant_indices(indices2remove) = []; % Remove last column
        
        % Create elements
        for i=1:length(relevant_indices)
            current_index = relevant_indices(i);
            node_a = relevant_nodes(current_index);
            node_b = relevant_nodes(current_index+1);
            node_c = relevant_nodes(current_index+1+n_x1);
            node_d = relevant_nodes(current_index+n_x1);
            ELEMENTS_1(current_elem,:) = [current_surf,node_a,node_b,node_c,node_d];
            current_elem = current_elem + 1; % Iterate up
        end
    end
end

for i = height(ELEMENTS_1):-1:1
    if isnan(ELEMENTS_1(i,1))
        ELEMENTS_1(i,:) = [];
    end
end


end