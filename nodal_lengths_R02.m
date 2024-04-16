function SURFACELENGTHS_1 = nodal_lengths_R02(surfaces2edit,MESHGRIDS_1,ELEMENTS_1,SURFACELENGTHS_1)
% Generates a holding matrices of the lengths between each node. Works only
% within a surface.
% Only calculates for nodes that connect.

% Version 2.0 completed 10/19/2023

for kk = 1:length(surfaces2edit)
    surface_number = surfaces2edit(kk);
    relevant_elems = find(ELEMENTS_1(:,1)==surface_number);
    for ii = 1:length(relevant_elems)
        elem_nodes = [ELEMENTS_1(relevant_elems(ii),2:5), ELEMENTS_1(relevant_elems(ii),2)]; % Element nodes with node_a wrapped to the back too
        for jj = 1:4
            L = norm(MESHGRIDS_1(elem_nodes(jj),4:6)-MESHGRIDS_1(elem_nodes(jj+1),4:6));
            SURFACELENGTHS_1(elem_nodes(jj),elem_nodes(jj+1)) = L;
            SURFACELENGTHS_1(elem_nodes(jj+1),elem_nodes(jj)) = L;
        end
    end
end
end