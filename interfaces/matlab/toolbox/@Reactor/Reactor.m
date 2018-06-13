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
%     Integer, reactor type. Options are:
%
%     1. Reservoir
%     2. Reactor
%     3. Flow Reactor
%     4. Constant Pressure Reactor
%     5. Ideal Gas Reactor
%     6. Ideal Gas Constant Pressure Reactor
%
% :return:
%     Instance of class :mat:func:`Reactor`
%

if nargin == 0
    contents = 0;
    typ = 2;
elseif nargin == 1
    typ = 2;
elseif nargin > 2
    error('too many arguments');
end

x.index = reactormethods(0, typ);
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
