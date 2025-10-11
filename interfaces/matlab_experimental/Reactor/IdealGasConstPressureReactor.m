classdef IdealGasConstPressureReactor < Reactor
    % Create a constant pressure reactor with an ideal gas. ::
    %
    %     >> r = IdealGasConstPressureReactor(phase, name)
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
    % :param phase:
    %     Cantera :mat:class:`Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :return:
    %     Instance of class :mat:class:`IdealGasConstPressureReactor`.

    methods

        function r = IdealGasConstPressureReactor(phase, name)
            % Constructor

            if nargin < 2
                name = '(none)';
            end

            r@Reactor(phase, 'IdealGasConstPressureReactor', name);
        end

    end
end
