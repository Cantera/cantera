function setState_Tsat(tp, tx)
% SETSTATE_TSAT  Set the fluid to a saturated state at a given temperature.
% setState_Tsat(tp,tx)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param tx:
%     Temperature. Units: K
%

thermo_set(tp.tp_id, 25, tx);
