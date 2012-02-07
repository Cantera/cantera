function v = electricPotential(a)
% ELECTRICPOTENTIAL - the electric potential of the phase
%
v = thermo_get(a.tp_id,25);
