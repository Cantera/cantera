function o = Oxygen()
% OXYGEN  Return an object representing oxygen.
% o = Oxygen()
% The object returned by this method implements an accurate equation of
% state for oxygen that can be used in the liquid, vapor, saturated
% liquid/vapor, and supercritical regions of the phase diagram.  The
% equation of state is taken from
%
% Reynolds, W. C. *Thermodynamic Properties in SI: graphs, tables, and
% computational equations for forty substances* Stanford: Stanford
% University, 1979. Print.
%
% :return:
%     Instance of class :mat:func:`Solution`
%

o = Solution('liquidvapor.yaml', 'oxygen');
