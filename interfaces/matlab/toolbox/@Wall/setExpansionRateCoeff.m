function setExpansionRateCoeff(w, k)
% SETEXPANSIONRATECOEFF  Set the expansion rate coefficient.
% setExpansionRateCoeff(w, k)
% The expansion rate coefficient
% determines the velocity of the wall for a given pressure
% differential between the left and right reactors, according to the
% formula
%
% .. math:: v = K(P_{left}-P_{right})
%
% where :math:`v` is velocity, :math:`K` is the expansion rate
% coefficient, and :math:`P` is the pressure.
%
% :param w:
%     Instance of class :mat:func:`Wall`
% :param k:
%     Double, wall expansion rate coefficient. Units: m/(s-Pa)
%

wallmethods(9, wall_hndl(w), k);
