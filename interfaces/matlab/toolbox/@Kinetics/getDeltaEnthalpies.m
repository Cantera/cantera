function dH = getDeltaEnthalpies(a)
% GETDELTAENTHALPIES  Get the enthalpy of reaction for each reaction.
% dH = getDeltaEnthalpies(a)
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics) for which the enthalpies of
%     reaction are desired.
% :return:
%     Returns a vector of the enthalpy of reaction for each
%     reaction. Units: J/kmol
%

dH = kinetics_get(a.id, 17, 0);
