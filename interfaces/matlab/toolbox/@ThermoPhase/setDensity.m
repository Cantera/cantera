function setDensity(tp, rho)
% SETDENSITY  Set the density.
% setDensity(tp,rho)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param rho:
%     Density. Units: kg/m**3
%

if rho <= 0.0
    error('The density must be positive.');
end

phase_set(tp.tp_id, 2, rho);
