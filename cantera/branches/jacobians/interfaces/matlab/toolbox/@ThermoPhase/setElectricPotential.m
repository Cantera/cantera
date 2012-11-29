function a = setElectricPotential(a,phi)
% SETELECTRICPOTENTIAL  Set the electric potential [V].
%
thermo_set(a.tp_id,2,phi);
