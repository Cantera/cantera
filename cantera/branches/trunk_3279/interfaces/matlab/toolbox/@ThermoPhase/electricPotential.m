function v = electricPotential(tp)
% ELECTRICPOTENTIAL  Get the electric potential.
% v = electricPotential(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     The electric potential of the phase. Units: V
%

v = thermo_get(tp.tp_id, 25);
