function v = electricPotential(tp)
% ELECTRICPOTENTIAL - the electric potential of the phase
%

v = thermo_get(tp.tp_id, 25);
