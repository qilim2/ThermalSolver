function [Qdot_out_nodes,Qdot_out_elems] = radiation_out_R03(Tn1,pre_rad_out_coefficient,elem_nodes_map)
% Calculates Qdot for input elements depending on the given coefficient and
% current temperature of the nodes. This function MUST work with nodes
% since the temperature of each node may be different.

% Tn1 = Vector of current temperatures for each node.
% pre_rad_out_coefficient = Truncated vector of precalculated values for each element. Repeated four times to match list of nodes.
% elem_nodes_map = Nx4 matrix of nodes listing the four nodes of each input element.

% Version 3.0 completeed 4/3/2024

% Reshaped list of nodes that is in order based on elements
reshaped_nodes_vec = reshape(elem_nodes_map',[numel(elem_nodes_map),1]);

% Preallocate
Qdot_out_nodes = zeros(max(reshaped_nodes_vec),1);
Qdot_out_elems = zeros(height(elem_nodes_map),1);

% Solve
Qdot_tempval = pre_rad_out_coefficient.*(Tn1(reshaped_nodes_vec).^4);

% Sum values at nodes
for ii = 1:length(Qdot_tempval)
    % Assign to vector whose shape matches nodes vector
    current_node = reshaped_nodes_vec(ii); % Get node index
    Qdot_out_nodes(current_node,1) = Qdot_out_nodes(current_node,1) + Qdot_tempval(ii); % Assign/add to overall vector

    % Assign to parent element
    current_elem = ceil(ii/4);
    Qdot_out_elems(current_elem,1) = Qdot_out_elems(current_elem,1) + Qdot_tempval(ii);
end

end