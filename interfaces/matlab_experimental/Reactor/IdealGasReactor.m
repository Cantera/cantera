classdef IdealGasReactor < Reactor
    % Create a reactor with an ideal gas. ::
    %
    %     >> r = IdealGasReactor(phase, name)
    %
    % An :mat:class:`IdealGasReactor` is an instance of :mat:class:`Reactor` where
    % the governing equations are specialized for the ideal gas equation of state
    % (and do not work correctly with other thermodynamic models). Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = IdealGasReactor(gas)    % a reactor containing a gas
    %
    % See also: :mat:class:`Reactor`
    %
    % :param phase:
    %     Cantera :mat:class:`Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :return:
    %     Instance of class :mat:class:`IdealGasReactor`.

    methods

        function r = IdealGasReactor(phase, name)
            % Constructor

            if nargin < 2
                name = '(none)';
            end

            r@Reactor(phase, 'IdealGasReactor', name);
        end

    end

end
