function h = Hydrogen()
% HYDROGEN  Return an object representing hydrogen.
% h = Hydrogen()
% The object returned by this method implements an accurate equation of
% state for hydrogen that can be used in the liquid, vapor, saturated
% liquid/vapor, and supercritical regions of the phase diagram. The
% equation of state is taken from
%
% Reynolds, W. C. *Thermodynamic Properties in SI: graphs, tables, and
% computational equations for forty substances* Stanford: Stanford
% University, 1979. Print.
%
% For more details, see classes Cantera::PureFluid and tpx::hydrogen in the
% Cantera C++ source code documentation.
%
% :return:
%     Instance of class :mat:func:`Solution`
%

h = Solution('liquidvapor.yaml', 'hydrogen');
