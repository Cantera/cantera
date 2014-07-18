function n = molarDensity(tp)
% MOLARDENSITY  Get the molar density.
% n = molarDensity(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Molar density. Units: kmol/m^3
%

n = phase_get(tp.tp_id, 3);
