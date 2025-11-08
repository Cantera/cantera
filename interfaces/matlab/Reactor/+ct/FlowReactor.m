classdef FlowReactor < ct.ReactorBase
    % Create a flow reactor object. ::
    %
    %     >> r = ct.FlowReactor(phase, name, clone)
    %
    % A reactor representing adiabatic plug flow in a constant-area
    % duct. Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = ct.FlowReactor(gas)    % a reactor containing a gas
    %
    % See also: :mat:class:`ct.ReactorBase`
    %
    % :param phase:
    %     Cantera :mat:class:`ct.Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :param clone:
    %    Determines whether to clone `phase` so that the internal state of
    %    this reactor is independent of the original Solution object and
    %    any Solution objects used by other reactors in the network.
    %    (optional; default is true).

    properties (SetAccess = public)

        massFlowRate % Mass flow rate [kg/s].

    end

    methods

        function r = FlowReactor(phase, name, clone)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "(none)"
                clone (1,1) logical = true
            end

            id = ct.impl.call('mReactor_new', 'FlowReactor', phase.solnID, clone, name);
            r@ct.ReactorBase(id);
        end

        %% FlowReactor Get Methods

        function flag = get.massFlowRate(obj)
            rate = ct.impl.call('mReactor_massFlowRate', obj.id);
        end

        %% FlowReactor Set Methods

        function set.massFlowRate(obj, MFR)
            ct.impl.call('mReactor_setMassFlowRate', obj.id, MFR);
            obj.massFlowRate = MFR;
        end

    end
end
