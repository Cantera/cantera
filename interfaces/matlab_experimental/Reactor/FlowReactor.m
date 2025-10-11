classdef FlowReactor < Reactor
    % Create a flow reactor object. ::
    %
    %     >> r = FlowReactor(phase, name)
    %
    % A reactor representing adiabatic plug flow in a constant-area
    % duct. Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = FlowReactor(gas)    % a reactor containing a gas
    %
    % See also: :mat:class:`Reactor`
    %
    % :param phase:
    %     Cantera :mat:class:`Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :return:
    %     Instance of class :mat:class:`FlowReactor`.

    methods

        function r = FlowReactor(phase, name)
            % Constructor

            if nargin < 2
                name = '(none)';
            end

            r@Reactor(phase, 'FlowReactor', name);
        end

    end
end
