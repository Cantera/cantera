function x = ReactorNet(reactors)
% REACTORNET  ReactorNet class constructor.
% x = ReactorNet(reactors)
% A :mat:func:`ReactorNet` object is a container that holds one
% or more :mat:func:`Reactor` objects in a network. :mat:func:`ReactorNet`
% objects are used to simultaneously advance the state of one or
% more coupled reactors.
%
% Example::
%
%     >> r1 = Reactor(gas1)
%     >> r2 = Reactor(gas2)
%     >> <... install walls, inlets, outlets, etc...>
%
%     >> reactor_network = ReactorNet({r1, r2})
%     >> advance(reactor_network, time)
%
% See also: :mat:func:`Reactor`
%
% :param reactors:
%     A single instance of :mat:func:`Reactor` or a cell array
%     of instances of :mat:func:`Reactor`
% :return:
%     Instance of class :mat:func:`ReactorNet`
%

if nargin ~= 1
    error('Wrong number of arguments to ReactorNet constructor.');
end

if isa(reactors, 'Reactor')
    % Allow simpler syntax for creating a network with one reactor
    reactors = {reactors};
end

x.index = reactornetmethods(0, 0);
if x.index < 0
    error(geterr);
end
x = class(x, 'ReactorNet');

% add reactors
nr = length(reactors);
for i = 1:nr
    addReactor(x, reactors{i});
end
