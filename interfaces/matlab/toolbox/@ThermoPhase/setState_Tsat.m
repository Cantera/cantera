function setState_Tsat(tp, tx)
% SETSTATE_TSAT  Set the fluid state using the given temperature and quality.
% setState_Tsat(tp,tx)
% The fluid state will be set to a saturated liquid-vapor state using the
% input temperature and vapor fraction (quality) as the independent,
% intensive variables.
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param tx:
%     Vector of length 2 containing the desired values for the temperature (K)
%     and the vapor fraction (quality)
%

thermo_set(tp.tp_id, 25, tx);
