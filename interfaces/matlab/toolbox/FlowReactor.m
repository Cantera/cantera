function r = FlowReactor(contents)
% FLOWREACTOR  Create a flow reactor object.
% r = FlowReactor(contents)
% A reactor representing adiabatic plug flow in a constant-area
% duct. Examples:
%
% .. code-block:: matlab
%
%     r1 = FlowReactor         % an empty reactor
%     r2 = FlowReactor(gas)    % a reactor containing a gas
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
r = Reactor(contents, 'FlowReactor');
