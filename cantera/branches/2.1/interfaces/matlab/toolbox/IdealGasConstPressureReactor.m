function r = IdealGasConstPressureReactor(contents)
% IDEALGASCONSTPRESSUREREACTOR - Create an IdealGasConstPressureReactor object.
%
%    An IdealGasConstPressureReactor is an instance of class Reactor where the
%    pressure is held constant. The volume is not a state variable, but
%    instead takes on whatever value is consistent with holding the pressure
%    constant. Additionally, its governing equations are specialized for the
%    ideal gas equation of state (and do not work correctly with other
%    thermodynamic models).
%
%        r1 = IdealGasConstPressureReactor      % an empty reactor
%        r2 = IdealGasConstPressureReactor(gas) % a reactor containing a gas
%
%   See also: Reactor
%
if nargin == 0
    contents = 0;
end
r = Reactor(contents, 6);
