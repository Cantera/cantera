function x = ReactorNet(reactors)
% REACTOR - Create a ReactorNet object.
%
%    A ReactorNet object is a container that holds one or more
%    Reactor objects.
%
if nargin == 1
else
    error('wrong number of arguments to ReactorNet constructor');
end

if isa(reactors, 'Reactor')
    % Allow simpler syntax for creating a network with one reactor
    reactors = {reactors};
end

x.index = reactornetmethods(0,0);
if x.index < 0
    error(geterr);
end
x = class(x,'ReactorNet');

% add reactors
nr = length(reactors);
for i = 1:nr
    addReactor(x,reactors{i});
end
