function n = Nitrogen()
% NITROGEN  Return an object representing nitrogen.
% n = Nitrogen()
% The object returned by this method implements an accurate equation of
% state for nitrogen that can be used in the liquid, vapor, saturated
% liquid/vapor, and supercritical regions of the phase diagram. The
% equation of state is taken from
%
% Reynolds, W. C. *Thermodynamic Properties in SI: graphs, tables, and
% computational equations for forty substances* Stanford: Stanford
% University, 1979. Print.
%
% :return:
%     Instance of class :mat:func:`Solution`
%

n = Solution('liquidvapor.cti', 'nitrogen');
