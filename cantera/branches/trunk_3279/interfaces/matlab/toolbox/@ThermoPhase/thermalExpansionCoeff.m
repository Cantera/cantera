function a = thermalExpansionCoeff(tp)
% THERMALEXPANSIONCOEFF  Get the thermal expansion coefficient.
% a = thermalExpansionCoeff(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :return:
%     Thermal Expansion Coefficient. Units: 1/K
%

a = thermo_get(tp.tp_id, 27);
