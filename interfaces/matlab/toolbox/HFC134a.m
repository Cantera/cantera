function h = HFC134a()
% HFC134A  Return an object representing refrigerant HFC134a.
% h = HFC134a()
% The object returned by this method implements an accurate equation of
% state for refrigerant HFC134a (R134a) that can be used in the liquid,
% vapor, saturated liquid/vapor, and supercritical regions of the phase
% diagram. Implements the equation of state given in:
% R. Tillner-Roth and H. D. Baehr. "An International Standard Formulation for
% The Thermodynamic Properties of 1,1,1,2-Tetrafluoroethane (HFC-134a) for
% Temperatures From 170 K to 455 K and Pressures up to 70 MPa". J. Phys.
% Chem. Ref. Data, Vol. 23, No. 5, 1994. pp. 657--729.
% http://dx.doi.org/10.1063/1.555958
%
% For more details, see classes Cantera::PureFluid and tpx::HFC134a in the
% Cantera C++ source code documentation.
%
% :return:
%     Instance of class :mat:func:`Solution`
%

h = Solution('liquidvapor.yaml', 'HFC-134a');
