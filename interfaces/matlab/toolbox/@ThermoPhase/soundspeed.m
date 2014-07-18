function c = soundspeed(tp)
% SOUNDSPEED  Get the speed of sound.
% c = soundspeed(tp)
% If the phase is an ideal gas, the speed of sound
% is calculated by:
%
% .. math:: c = \sqrt{\gamma * R * T}
%
% where :math:`\gamma` is the ratio of specific heats, :math:`R` is
% the specific gas constant, and :math:`T` is the temperature. If the
% phase is not an ideal gas, the speed of sound is calculated by
%
% .. math:: c = \sqrt{\left(\frac{\partial p}{\partial \rho}\right)_s}
%
% where :math:`p` is the pressure and :math:`\rho` is the density,
% and the subscript :math:`s` indicates constant entropy. This is
% approximated by slightly increasing the density at constant entropy
% and computing the change in pressure.
%
% .. math:: c = \sqrt{\frac{p_1 - p_0}{\rho_1-\rho_0}}
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :return:
%     The speed of sound. Units: m/s
%

if isIdealGas(tp)
    gamma = cp_mass(tp)/cv_mass(tp);
    wtm = meanMolecularWeight(tp);
    r = gasconstant/wtm;
    c = sqrt(gamma*r*temperature(tp));
else
    rho0 = density(tp);
    p0 = pressure(tp);
    s0 = entropy_mass(tp);
    rho1 = 1.001*rho0;
    set(tp, 'Density', rho1, 'Entropy', s0);
    p1 = pressure(tp);
    dpdrho_s = (p1 - p0)/(rho1 - rho0);
    c = sqrt(dpdrho_s);
end
