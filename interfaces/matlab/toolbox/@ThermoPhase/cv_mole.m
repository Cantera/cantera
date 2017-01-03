function v = cv_mole(tp)
% CV_MOLE  Get the molar-basis specific heat at constant volume.
% v = cv_mole(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Molar basis specific heat of the mixture at
%     constant volume. Units: J/kmol-K
%

v = thermo_get(tp.tp_id, 7);
