function out = rect_elem_center_R01(elem_index,ELEMENTS_11,MESHGRIDS_11)
% Calculates center of a rectangular element given its index.
nodeA = ELEMENTS_11(elem_index,2); % Get node A ID
nodeC = ELEMENTS_11(elem_index,4); % Get node C ID
nodeAind = MESHGRIDS_11(:,1)==nodeA; % Get boolean index of node A
nodeCind = MESHGRIDS_11(:,1)==nodeC; % Get boolean index of node C
nodeAcoords = MESHGRIDS_11(nodeAind,4:6); % Get coords of node A
nodeCcoords = MESHGRIDS_11(nodeCind,4:6); % Get coords of node C
out = mean([nodeAcoords;nodeCcoords]); % Average
end