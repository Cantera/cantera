function x = ReactorNet(reactors)
% REACTOR - Create a ReactorNet object.
%
%    A ReactorNet object is a container that holds one or more
%    Reactor objects. 
%
%   See also: Reservoir
%
if nargin == 1
else
  error('wrong number of arguments');
end

x.index = reactornetmethods(0);
if x.index < 0
  error(geterr);
end
x = class(x,'ReactorNet');

% add reactors
unfinished


