function v = cp_R(p)
% CP_R - Species non-dimensional heat capacities.
%
%        This method returns an array containing the pure species
%        standard-state heat capacities at constant pressure.
%
v = thermo_get(p.tp_id,38);
