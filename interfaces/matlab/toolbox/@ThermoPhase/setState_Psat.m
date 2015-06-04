function setState_Psat(tp, px)
% SETSTATE_PSAT  Set the fluid to a saturated state at a given pressure.
% setState_Psat(tp,px)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param px:
%     Pressure. Units: Pa
%
thermo_set(tp.tp_id, 24, px);
