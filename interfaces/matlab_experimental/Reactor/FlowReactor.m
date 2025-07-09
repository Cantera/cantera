classdef FlowReactor < Reactor
    % Create a flow reactor object. ::
    %
    %     >> r = FlowReactor(contents, name)
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
    % :param contents:
    %     Cantera :mat:class:`Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :return:
    %     Instance of class :mat:class:`FlowReactor`.

    methods

        function r = FlowReactor(contents, name)
            % Constructor

            if nargin < 2
                name = '(none)'
            end

            r@Reactor(contents, 'FlowReactor', name);
        end

    end
end
