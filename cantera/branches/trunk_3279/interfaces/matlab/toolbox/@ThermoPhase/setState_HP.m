function setState_HP(tp, hp)
% SETSTATE_HP  Set the specific enthalpy and pressure.
% setState_HP(tp,hp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param hp:
%     Vector of length 2 containing the desired values for the specific
%     enthalpy (J/kg) and pressure (Pa).
%

if hp(2) <= 0.0
    error('The pressure must be positive.');
end

thermo_set(tp.tp_id, 20, hp);
