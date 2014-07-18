function b = isothermalCompressibility(tp)
% ISOTHERMALCOMPRESSIBILITY  Get the isothermal compressibility.
% b = isothermalCompressibility(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Isothermal Compressibility. Units: 1/Pa
%

b = thermo_get(tp.tp_id, 26);

