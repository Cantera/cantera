function r = ConstPressureReactor(contents)
% CONSTPRESSUREREACTOR  Create a constant pressure reactor object.
% r = ConstPressureReactor(contents)
% A :mat:func:`ConstPressureReactor` is an instance of class
% :mat:func:`Reactor` where the pressure is held constant. The volume
% is not a state variable, but instead takes on whatever value is
% consistent with holding the pressure constant. Examples:
%
% .. code-block:: matlab
%
%     r1 = ConstPressureReactor         % an empty reactor
%     r2 = ConstPressureReactor(contents)    % a reactor containing contents
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
r = Reactor(contents, 'ConstPressureReactor');
