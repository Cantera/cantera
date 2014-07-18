function setElectricPotential(tp, phi)
% SETELECTRICPOTENTIAL  Set the electric potential [V].
%

thermo_set(tp.tp_id, 2, phi);
