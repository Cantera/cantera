function v = critDensity(tp)
% CRITDENSITY  Get the critical density.
% v = critDensity(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Critical density. Units: kg/m**3
%

v = thermo_get(tp.tp_id, 21);
