function Qdot_e2e = radiation_internal_R01(g,RPM,elem_count,Qdot_out_elems_top,Qdot_out_elems_bot)
% Calculates the surface-to-surface radiation heat transfer for a radiation
% group. Outputs to full matrix showing element to element radiation
% exchange rates.
% https://www.afs.enea.it/project/neptunius/docs/fluent/html/th/node116.htm

% RPM gives the radiation proportion matrix for each radiation group,
% listed only for the active elements. Isolate active elements using the
% trole and brole boolean vectors and the list of elements in each group.

% Version 1.0 completed 2/26/2024

elems = g.elems; % Get list of elements in rad group
tactive = g.trole == 1; bactive = g.brole == 1;
elems_t_active = elems(tactive); % Get elements that have topside active
elems_b_active = elems(bactive); % Get list of elements that have bottomside active
elems_active = [elems_t_active;elems_b_active];
Qdot_out = [Qdot_out_elems_top(elems_t_active);Qdot_out_elems_bot(elems_b_active)];


Qdot_temp = sum(RPM.*Qdot_out,1)';
Qdot_e2e = zeros(elem_count,1);
for ii = 1:length(elems_active)
    elem = elems_active(ii);
    Qdot_e2e(elem) = Qdot_e2e(elem) + Qdot_temp(ii);
end

end