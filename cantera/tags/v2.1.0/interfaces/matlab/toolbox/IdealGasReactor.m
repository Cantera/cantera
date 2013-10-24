function r = IdealGasReactor(contents)
% IDEALGASREACTOR - Create an IdealGasReactor object.
%
%    An IdealGasReactor is an instance of class Reactor where the governing
%    equations are specialized for the ideal gas equation of state (and do not
%    work correctly with other thermodynamic models).
%
%        r1 = IdealGasReactor         % an empty reactor
%        r2 = IdealGasReactor(gas)    % a reactor containing a gas
%
%   See also: Reactor
%
if nargin == 0
    contents = 0;
end
r = Reactor(contents, 5);
