function setState_UV(tp, uv)
% SETSTATE_UV  Set the specific internal energy and specific volume.
% setState_UV(tp,uv)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param uv:
%     Vector of length 2 containing
%     the desired values for the specific internal energy (J/kg) and
%     specific volume (m**3/kg).
%

if uv(2) <= 0.0
    error('The specific volume must be positive.');
end
thermo_set(tp.tp_id, 21, uv);
