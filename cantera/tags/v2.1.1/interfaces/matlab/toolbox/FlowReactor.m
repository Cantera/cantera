function r = FlowReactor(contents)
% FLOWREACTOR - Create a FlowReactor object.
%
%    A reactor representing adiabatic plug flow in a constant-area duct.
%
%        r1 = FlowReactor         % an empty reactor
%        r2 = FlowReactor(gas)    % a reactor containing a gas
%
%   See also: Reactor
%
if nargin == 0
    contents = 0;
end
r = Reactor(contents, 3);
