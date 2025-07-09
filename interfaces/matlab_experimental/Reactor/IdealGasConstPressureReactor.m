classdef IdealGasConstPressureReactor < Reactor
    % Create a constant pressure reactor with an ideal gas. ::
    %
    %     >> r = IdealGasConstPressureReactor(contents, name)
    %
    % An :mat:class:`IdealGasConstPressureReactor` is an instance of
    % :mat:class:`Reactor` where the pressure is held constant.
    % The volume is not a state variable, but instead takes on
    % whatever value is consistent with holding the pressure constant.
    % Additionally, its governing equations are specialized for the
    % ideal gas equation of state (and do not work correctly with other
    % thermodynamic models). Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = IdealGasConstPressureReactor(gas) % a reactor containing a gas
    %
    % See also: :mat:class:`Reactor`
    %
    % :param contents:
    %     Cantera :mat:class:`Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :return:
    %     Instance of class :mat:class:`IdealGasConstPressureReactor`.

    methods

        function r = IdealGasConstPressureReactor(contents, name)
            % Constructor

            if nargin < 2
                name = '(none)';
            end

            r@Reactor(contents, 'IdealGasConstPressureReactor', name);
        end

    end
end
