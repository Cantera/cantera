function setState_Psat(tp, px)
% SETSTATE_PSAT  Set the fluid state using the given pressure and quality.
% setState_Psat(tp,px)
% The fluid state will be set to a saturated liquid-vapor state using the
% input pressure and vapor fraction (quality) as the independent,
% intensive variables.
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param px:
%     Vector of length 2 containing the desired values for the pressure (Pa)
%     and the vapor fraction
%

thermo_set(tp.tp_id, 24, px);
