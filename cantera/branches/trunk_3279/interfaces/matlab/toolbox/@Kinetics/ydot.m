function v = ydot(a)
% YDOT  Get the mass production rates of the species.
% v = ydot(a)
% Evaluates the source term :math:`\dot{\omega}_k M_k /\rho`
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which the ydots are desired.
% :return:
%     Returns a vector of length nSpecies. Units: kg/s
%

v = kinetics_get(a.id, 24, 0);
