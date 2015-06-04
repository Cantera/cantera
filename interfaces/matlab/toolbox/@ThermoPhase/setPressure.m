function setPressure(tp, p)
% SETPRESSURE  Set the pressure.
% setPressure(tp,p)
% The pressure is set by changing the density holding the
% temperature and chemical composition fixed.
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param p:
%     Pressure. Units: Pa
%

if p <= 0.0
    error('The pressure must be positive.')
end

thermo_set(tp.tp_id, 1, p);
