function a = setDensity(a,rho)
% SETDENSITY - Set the density [kg/m^3].
%
%   setDensity(phase, 0.01);
%
if rho <= 0.0
    error('the density must be positive');
end

phase_set(a.tp_id,2,rho);
