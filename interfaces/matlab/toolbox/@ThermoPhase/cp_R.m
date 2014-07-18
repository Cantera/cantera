function v = cp_R(tp)
% CP_R - Species non-dimensional heat capacities.
%
%        This method returns an array containing the pure species
%        standard-state heat capacities at constant pressure.
%

v = thermo_get(tp.tp_id, 38);
