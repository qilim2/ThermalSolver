function Qdot_nodes = elem2nodes_R01(Qdot_elems,elem_nodes_map,node_count)
% Spreads heat for each element out to the nodes

Qdot_nodes = zeros(node_count,1);
for ii = 1:length(Qdot_elems)
    Qdot_add = Qdot_elems(ii)/4;
    element_nodes = elem_nodes_map(ii,:);
    Qdot_nodes(element_nodes) = Qdot_add + Qdot_nodes(element_nodes);
end
end