function r = Reservoir(contents)
% RESERVOIR  Create a Reservoir object.
% r = Reservoir(contents)
% A :mat:func:`Reservoir` is an instance of class :mat:func:`Reactor`
% configured so that its intensive state is constant in time. A
% reservoir may be thought of as infinite in extent, perfectly mixed,
% and non-reacting, so that fluid may be extracted or added without
% changing the composition or thermodynamic state. Note that even
% if the reaction mechanism associated with the fluid in the
% reactor defines reactions, they are disabled within
% reservoirs. Examples:
%
% .. code-block:: matlab
%
%     r1 = Reservoir         % an empty reservoir
%     r2 = Reservoir(gas)    % a reservoir containing a gas
%
% See also: :mat:func:`Reactor`
%
% :param contents:
%     Cantera :mat:func:`Solution` to be set as the contents of the
%     reactor
% :return:
%     Instance of class :mat:func:`Reactor`
%

if nargin == 0
    contents = 0;
end
r = Reactor(contents, 'Reservoir');
