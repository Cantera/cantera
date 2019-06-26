function x = Reactor(contents, typ)
% REACTOR  Reactor class constructor.
% x = Reactor(contents, typ)
% A :mat:func:`Reactor` object simulates a perfectly-stirred reactor.
% It has a time-dependent state, and may be coupled to other reactors
% through flow lines or through walls that may expand or
% contract and/or conduct heat.
%
% .. code-block:: matlab
%
%     >> r1 = Reactor      % an empty reactor
%     >> r2 = Reactor(gas) % a reactor containing a phase
%
% See also: :mat:func:`Reservoir`, :mat:func:`IdealGasReactor`,
% :mat:func:`IdealGasConstPressureReactor`, :mat:func:`ConstPressureReactor`
%
% :param contents:
%     Instance of class :mat:func:`Solution` representing the contents of the
%     reactor
% :param typ:
%     Character array, reactor type. Options are:
%
%     'Reservoir'
%     'Reactor'
%     'FlowReactor'
%     'ConstPressureReactor'
%     'IdealGasReactor'
%     'IdealGasConstPressureReactor'
%
% :return:
%     Instance of class :mat:func:`Reactor`
%

if nargin == 0
    contents = 0;
    typ = 'Reactor';
elseif nargin == 1
    typ = 'Reactor';
elseif nargin > 2
    error('too many arguments');
end

if isa(typ, 'double')
    warning('Definition via integer type to be deprecated after Cantera 2.5')
    reactor_types = {'Reservoir' 'Reactor' 'FlowReactor' ...
                     'ConstPressureReactor' 'IdealGasReactor' ...
                     'IdealGasConstPressureReactor'};
    typ = reactor_types(typ);
end

x.type = char(typ);
x.index = reactormethods(0, x.type);
if x.index < 0
    error(geterr);
end
x.contents = contents;
x = class(x, 'Reactor');

if isa(contents, 'Solution')
    insert(x, contents);
elseif ~(isa(contents, 'double') && contents == 0)
    error('Reactor contents must be an object of type "Solution"')
end
