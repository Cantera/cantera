function w = Water()
% WATER - Return an object representing water.
%
%   The object returned by this method implements an accurate equation of
%   state for water that can be used in the liquid, vapor, saturated
%   liquid/vapor, and supercritical regions of the phase diagram.  The
%   equation of state is taken from W. C. Reynolds, "Thermodynamic
%   Properties in SI."
%
%   For more details, see classes Cantera::PureFluid and tpx::water in the
%   Cantera C++ source code documentation.
%
w = importPhase('liquidvapor.cti','water');
