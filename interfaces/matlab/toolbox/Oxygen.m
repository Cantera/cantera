function n = Oxygen()
% OXYGEN - Return an object representing oxygen.
%
%   The object returned by this method implements an accurate equation of
%   state for oxygen that can be used in the liquid, vapor, saturated
%   liquid/vapor, and supercritical regions of the phase diagram.  The
%   equation of state is taken from W. C. Reynolds, "Thermodynamic
%   Properties in SI."
%
n = importPhase('liquidvapor.cti','oxygen');
