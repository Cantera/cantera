function v = cp_mole(tp)
% CP_MOLE  Get the molar-basis specific heat at constant pressure.
% v = cp_mole(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Molar basis specific heat of the mixture at
%     constant pressure. Units: J/kmol-K
%

v = thermo_get(tp.tp_id, 6);
