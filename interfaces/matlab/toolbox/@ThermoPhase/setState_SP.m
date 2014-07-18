function setState_SP(tp, sp)
% SETSTATE_SP  Set the specific entropy and pressure.
% setState_SP(tp,sp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param sp:
%     Vector of length 2 containing the desired values for the specific
%     entropy (J/kg-K) and pressure (Pa).
%

if sp(1) <= 0.0
    error('The specific entropy must be positive.');
end
if sp(2) <= 0.0
    error('The pressure must be positive.');
end
thermo_set(tp.tp_id, 23, sp);
